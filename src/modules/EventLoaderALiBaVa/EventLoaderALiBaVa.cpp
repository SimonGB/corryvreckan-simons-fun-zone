/**
 * @file
 * @brief Implementation of module EventLoaderALiBaVa
 *
 * @copyright Copyright (c) 2021-2022 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <dirent.h>

//#include "objects/Pixel.hpp"
//#include "ALiBaVa/auxfunctions.h"
#include "EventLoaderALiBaVa.h"
//#include "ALiBaVa/HDFRoot.h"

using namespace corryvreckan;

EventLoaderALiBaVa::EventLoaderALiBaVa(Configuration& config, std::shared_ptr<Detector> detector)
    : Module(config, detector), detector_(detector) {}

void EventLoaderALiBaVa::initialize() {

    // Take input directory, run number,
    // lower timecut, upper timecut,        --- in nanoseconds
    // first X ALiBaVa events to ignore, channel range,
    // from global parameters in the config.
    // Default values are set first

    config_.setDefault<int>("run", 0);
    config_.setDefault<double>("LowerTimecut", 0);
    config_.setDefault<double>("UpperTimecut", std::numeric_limits<double>::max());
    config_.setDefault<int>("IgnoreEvents", 1);
    config_.setDefault<double>("Chargecut", std::numeric_limits<double>::max());
    config_.setDefault<double>("CalibrationConstant", 1.0);
    config_.setDefault<bool>("CrosstalkCorrection", false);
    config_.setDefaultArray<unsigned int>("ROI", {0, 255});
    config_.setDefault<int>("Polarity", -1);

    m_inputDirectory = config_.getPath("input_directory");
    m_run = config_.get<int>("run");
    m_timecut_lower = config_.get<double>("LowerTimecut");
    m_timecut_upper = config_.get<double>("UpperTimecut");
    m_ignore_events = config_.get<int>("IgnoreEvents");
    m_chargecut = config_.get<double>("Chargecut");
    m_calibration_constant = config_.get<double>("CalibrationConstant");
    m_correct_crosstalk = config_.get<bool>("CrosstalkCorrection");
    m_roi = config_.getArray<unsigned int>("ROI");
    m_polarity = config_.get<int>("Polarity");
    if(!m_correct_crosstalk) {
        config_.setDefault<double>("b_one", 0.);
        config_.setDefault<double>("b_two", 0.);
    }
    m_b_one = config_.get<double>("b_one");
    m_b_two = config_.get<double>("b_two");

    // Open the input directory
    DIR* directory = opendir(m_inputDirectory.c_str());
    if(directory == nullptr) {
        LOG(ERROR) << "Directory \'" << m_inputDirectory << "\' does not exist or was not supplied.";
        return;
    }

    // Read the run-files (data, pedestal and calibration) in the folder
    dirent* entry;
    while(entry = readdir(directory)) {
        if(entry->d_type == DT_REG) {
            std::string entryName = entry->d_name;
            if(entryName.find(std::to_string(m_run) + ".dat") != std::string::npos ||
               entryName.find("dat_run" + std::to_string(m_run) + ".hdf") != std::string::npos) {
                m_datafilename = m_inputDirectory + "/" + entryName;
            }
            if(entryName.find(std::to_string(m_run) + ".ped") != std::string::npos ||
               entryName.find("ped_run" + std::to_string(m_run) + ".hdf") != std::string::npos) {
                m_pedestalfilename = m_inputDirectory + "/" + entryName;
            }
            if(entryName.find(std::to_string(m_run) + ".cal") != std::string::npos ||
               entryName.find("cal_run" + std::to_string(m_run) + ".hdf") != std::string::npos) {
                m_calibrationfilename = m_inputDirectory + "/" + entryName;
            }
        }
    }

    // Log errors in case the files aren't found in the folder.
    // The datafile can also be supplied directly in the config.
    if(m_datafilename.length() == 0) {
        LOG(ERROR) << "No data file was found for ALiBaVa in " << m_inputDirectory;
        return;
    }
    if(m_pedestalfilename.length() == 0) {
        LOG(WARNING) << "No pedestal file was found."
                     << "\n"
                     << "Datafile will be used for pedestal";
    }

    if(m_calibrationfilename.length() == 0) {
        LOG(WARNING) << "No calibration file was found."
                     << "\n"
                     << "Results will be uncalibrated: ADC = charge.";
    }

    // Create histograms
    std::string title = "Charge of signal; Charge [e]; # entries";
    hChargeSignal = new TH1F("chargeSignal", title.c_str(), 100, -96000, 32000);

    title = "ADC of signal; Signal [ADC]; # entries";
    hADCSignal = new TH1F("ADCSignal", title.c_str(), 100, -200, 600);

    title = "Signal to noise ratio; SNRatio; # entries";
    hSNR = new TH1F("SNRatio", title.c_str(), 100, -50, 200);

    title = "Uncorrected pedestal; # channel; Pedestal[ADC]";
    hPedestal = new TH1F("pedestal", title.c_str(), 256, 0, 256);

    title = "Corrected pedestal; # channel; Pedestal[ADC]";
    hPedestalCorrect = new TH1F("pedestalCorrect", title.c_str(), 256, 0, 256);

    title = "Uncorrected Noise; # channel; Noise[ADC]";
    hNoise = new TH1F("noise", title.c_str(), 256, 0, 256);

    title = "Corrected Noise; # channel; Noise[ADC]";
    hNoiseCorrect = new TH1F("noiseCorrect", title.c_str(), 256, 0, 256);

    title = "Time profile; Time [ns], Ave. signal highest channel [ADC]";
    hTimeProfile = new TProfile("timeProfile", title.c_str(), 35, 0, 35, 0, 200);

    // Create a pointer with the data file.
    ALiBaVaPointer = DataFileRoot::OpenFile(m_datafilename.c_str());
    // DEBUG stuff
    // print the vector
    // for (auto i: m_roi)
    // std::cout << i << std::endl;
    // print the size of the vector
    // std::cout << m_roi.size() << std::endl;
    // Sort vector to avoid errors later on
    std::sort(m_roi.begin(), m_roi.end());
    // Set the region of interest
    ALiBaVaPointer->set_ROI(m_roi);
    // Get the vector with all ROI channels listed to use here
    m_roi_ch = ALiBaVaPointer->get_ROI();

    ALiBaVaPointer->set_polarity(m_polarity);

    const char* ped_f = "alibava_ped.ped";
    const char* cal_f = "alibava_cal.cal";
    // Create a pointer with the pedestal file
    DataFileRoot* PedestalPointer = DataFileRoot::OpenFile(m_pedestalfilename.c_str());
    PedestalPointer->set_ROI(m_roi);

    // Calculate the pedestals, and compute and apply the common mode noise correction
    PedestalPointer->compute_pedestals_alternative();

    for(int chan : m_roi_ch) {
        double ped_val, noise_val;
        ped_val = PedestalPointer->ped(chan);
        noise_val = PedestalPointer->noise(chan);
        hPedestal->SetBinContent(chan, ped_val);
        hNoise->SetBinContent(chan, noise_val);
    }

    PedestalPointer->compute_cmmd_alternative();
    for(int chan : m_roi_ch) {
        double ped_val, noise_val;
        ped_val = PedestalPointer->ped(chan);
        noise_val = PedestalPointer->noise(chan);
        hPedestalCorrect->SetBinContent(chan, ped_val);
        hNoiseCorrect->SetBinContent(chan, noise_val);
    }

    // Save the calculated pedestal information in a temporary file
    PedestalPointer->save_pedestals(ped_f);
    PedestalPointer->close();
    delete PedestalPointer;
    // Load the calculated pedestal info into the original datafile
    ALiBaVaPointer->load_pedestals(ped_f, kTRUE);

    int columns, rows;

    columns = detector_->nPixels().X();
    rows = detector_->nPixels().Y();

    if(columns > rows) {
        m_horizontal = false;
    } else {
        m_horizontal = true;
    }

    std::vector<unsigned int> all_ch(256);
    std::iota(all_ch.begin(), all_ch.end(), 0);

    std::vector<unsigned int> mask_ch;
    std::set_difference(
        all_ch.begin(), all_ch.end(), m_roi_ch.begin(), m_roi_ch.end(), std::inserter(mask_ch, mask_ch.begin()));

    for(int i : mask_ch) {
        if(m_horizontal) {
            detector_->maskChannel(0, i);
        } else {
            detector_->maskChannel(i, 0);
        }
    }

    // Set the timecuts
    ALiBaVaPointer->set_timecut(m_timecut_lower, m_timecut_upper);

    // Ignore the first X events to ensure synchronisation, default is X = 1 which ignores the first event.
    for(int ievt = 0; ievt < m_ignore_events; ievt++) {
        ALiBaVaPointer->read_event();
    }
}

StatusCode EventLoaderALiBaVa::run(const std::shared_ptr<Clipboard>& clipboard) {
    // During running, every Corryvreckan event will get one ALiBaVa event

    // Create the pixelvector
    PixelVector pixels;
    // Get the event that already exists on the Corryvreckan clipboard
    auto event = clipboard->getEvent();

    // Read a data event from the ALiBaVa data file
    // Give feedback according to return code
    int return_code = ALiBaVaPointer->read_event();
    if(return_code == -1) {
        // LOG(WARNING) << "End of data file reached.";
        return StatusCode::EndRun;
        // Not sure if this is end of run or something else. Need to see difference between HDF5 and binary
    } else if(return_code == 0) {
        LOG(ERROR) << "There\'s something wrong (0) with the datafile";
        LOG(ERROR) << "Terminating run";
        return StatusCode::EndRun;
    } else if(return_code == 1) {
        // This means the event was read properly.
    } else if(return_code == 2) {
        LOG(ERROR) << "There\'s something wrong (2) with the datafile";
        LOG(ERROR) << "Terminating run";
        return StatusCode::EndRun;
    } else if(return_code == 3) {
        // LOG(ERROR) << "There\'s something wrong (3) with the datafile at event number " << iEvent-1;
        // LOG(ERROR) << "Terminating run";
        // return StatusCode::EndRun;
        // Still need to figure this out... I think it's just end of run. Need to see difference between HDF5 and binary
        // LOG(WARNING) << "End of data file reached.";
        return StatusCode::EndRun;
    } else if(return_code == 4) {
        LOG(ERROR) << "There\'s something wrong (4) with the HDF5 datafile";
        LOG(ERROR) << "Terminating run";
        return StatusCode::EndRun;
    } else {
        LOG(ERROR) << "This shouldn\'t happen.";
        LOG(ERROR) << "Terminating run";
        return StatusCode::EndRun;
    }
    // Calculate the common mode for the signal in this event
    ALiBaVaPointer->calc_common_mode_signal();
    // Process the opened data event, i.e. pedestal correction, common mode noise correction
    ALiBaVaPointer->process_event();
    // This gets the TDC time from the event, allowing timecuts around the event peak
    // The timecut is set in the ALiBaVa_loader() function.
    double TDCTime = ALiBaVaPointer->time();
    if(!ALiBaVaPointer->valid_time(TDCTime)) {
        clipboard->putData(pixels, detector_->getName());
        return StatusCode::NoData;
    }

    double trigger_ts;

    if(!clipboard->getEvent()->triggerList().empty()) {
        trigger_ts = clipboard->getEvent()->triggerList().begin()->second;
        LOG(DEBUG) << "Using trigger timestamp " << Units::display(trigger_ts, "us") << " as cluster timestamp.";

    } else {
        trigger_ts = 0;
        LOG(DEBUG) << "No trigger timestamp setting 0 ns as cluster timestamp.";
    }

    double channels_Sig_corrected[m_roi_ch.size()];
    /*
    if(m_correct_crosstalk){
      double channels_Sig[m_upper_channel-m_lower_channel+1];
      double b_one = m_b_one;
      double b_two = m_b_two;

      for(int chan = m_lower_channel; chan <= m_upper_channel; chan++){
        channels_Sig[chan-m_lower_channel] = ALiBaVaPointer->ADC_signal(chan);
      }

      channels_Sig_corrected[0] = (1+b_one+b_two)*channels_Sig[0];
      channels_Sig_corrected[1] = (1+b_one+b_two)*channels_Sig[1]-b_one*channels_Sig[0];
      for(int chan = m_lower_channel+2; chan <= m_upper_channel-2; chan++){
        channels_Sig_corrected[chan-m_lower_channel] =
    (1+b_one+b_two)*channels_Sig[chan-m_lower_channel]-b_one*channels_Sig[chan-m_lower_channel-1]-b_two*channels_Sig[chan-m_lower_channel-2];
      }
      channels_Sig_corrected[m_upper_channel-1] =
    (1+b_one)*channels_Sig[m_upper_channel-1]-b_one*channels_Sig[m_upper_channel-2]-b_two*channels_Sig[m_upper_channel-3];
      channels_Sig_corrected[m_upper_channel] =
    (1)*channels_Sig[m_upper_channel]-b_one*channels_Sig[m_upper_channel-1]-b_two*channels_Sig[m_upper_channel-2];
    }
    */
    double max_signal = 0;
    // This loops over the channels in the current ALiBaVa event
    for(int chan : m_roi_ch) {
        // In order, these are the calibration factor, the calibrated signal, and the ADC signal
        // If a pedestal file is supplied, pedestal and common mode error correction is applied to all
        // double calibration = ALiBaVaPointer->get_gain(chan);
        // double CalSignal = ALiBaVaPointer->signal(chan);
        // double ADCSignal = ALiBaVaPointer->signal(chan)*calibration;
        // double SNRatio = ALiBaVaPointer->sn(chan);
        double ADCSignal = 0;
        double SNRatio = 0;

        if(m_correct_crosstalk) {
            ADCSignal = channels_Sig_corrected[chan];
            SNRatio = channels_Sig_corrected[chan] / ALiBaVaPointer->noise(chan);
        } else {
            ADCSignal = ALiBaVaPointer->ADC_signal(chan);
            SNRatio = ALiBaVaPointer->sn(chan);
        }
        double CalSignal = ADCSignal * m_calibration_constant;

        if(ADCSignal > max_signal) {
            max_signal = ADCSignal;
            // std::cout << chan << std::endl;
        }

        // The chargecut is applied here
        // std::cout << CalSignal << " greater than " << m_chargecut << "?\n";
        // Not needed anymore but left in for now, other wise stuff explodes in Correlation
        if(CalSignal > m_chargecut) {
            // std::cout << "Here should be an entry" << CalSignal << std::endl;
            // if(true){
            // Create a pixel for every channel in this event with all the information and put it in the vector.
            // The value in the pixel reserved for the ADC value is used for the S/N ratio multiplied by 100000.

            std::shared_ptr<Pixel> pixel;

            if(m_horizontal) {
                pixel = std::make_shared<Pixel>(detector_->getName(), 0, chan, SNRatio * 100000, CalSignal, trigger_ts);
            } else {
                pixel = std::make_shared<Pixel>(detector_->getName(), chan, 0, SNRatio * 100000, CalSignal, trigger_ts);
            }

            pixels.push_back(pixel);

            // Fill the histograms
            hChargeSignal->Fill(CalSignal);
            hADCSignal->Fill(ADCSignal);
            hSNR->Fill(SNRatio);
        }
    }
    // std::cout << "Maximum signal is: " << max_signal << std::endl;
    hTimeProfile->Fill(TDCTime, max_signal, 1);
    // Put the created vector of pixels on the clipboard.
    clipboard->putData(pixels, detector_->getName());

    // If the pixels vector is empty, report this to Corryvreckan
    if(pixels.empty()) {
        return StatusCode::NoData;
    }

    // std::vector<double> event_time;
    // std::vector<double> highest_signal;

    // event_time.push_back(TDCTime);
    // highest_signal.push_back(max_signal);

    // Report the end of this event to Corryvreckan
    return StatusCode::Success;
}

void EventLoaderALiBaVa::finalize(const std::shared_ptr<ReadonlyClipboard>&) {
    ALiBaVaPointer->close();
    delete ALiBaVaPointer;
}
