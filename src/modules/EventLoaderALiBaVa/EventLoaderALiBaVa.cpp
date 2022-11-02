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
#include "EventLoaderALiBaVa.h"

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
    config_.setDefault<double>("timecut_low", 0);
    config_.setDefault<double>("timecut_up", std::numeric_limits<double>::max());
    config_.setDefault<int>("ignore_events", 1);
    config_.setDefault<double>("chargecut", std::numeric_limits<double>::max());
    config_.setDefault<double>("calibration_constant", 1.0);
    config_.setDefaultArray<unsigned int>("ROI", {0, 255});
    config_.setDefault<int>("polarity", -1);

    std::string input_directory = config_.getPath("input_directory");
    int run = config_.get<int>("run");
    double timecut_low = config_.get<double>("timecut_low");
    double timecut_up = config_.get<double>("timecut_up");
    int ignore_events = config_.get<int>("ignore_events");
    m_chargecut = config_.get<double>("chargecut");
    m_calibration_constant = config_.get<double>("calibration_constant");
    std::vector<unsigned int> roi = config_.getArray<unsigned int>("ROI");
    int polarity = config_.get<int>("polarity");

    // Open the input directory
    DIR* directory = opendir(input_directory.c_str());
    if(directory == nullptr) {
        LOG(ERROR) << "Directory \'" << input_directory << "\' does not exist or was not supplied.";
        return;
    }

    std::string datafilename;
    std::string pedestalfilename;
    
    // Read the run-files (data, pedestal and calibration) in the folder
    dirent* entry;
    while(entry = readdir(directory)) {
        if(entry->d_type == DT_REG) {
            std::string entryName = entry->d_name;
            if(entryName.find(std::to_string(run) + ".dat") != std::string::npos ||
               entryName.find("dat_run" + std::to_string(run) + ".hdf") != std::string::npos) {
                datafilename = input_directory + "/" + entryName;
            }
            if(entryName.find(std::to_string(run) + ".ped") != std::string::npos ||
               entryName.find("ped_run" + std::to_string(run) + ".hdf") != std::string::npos) {
                pedestalfilename = input_directory + "/" + entryName;
            }
        }
    }

    // Log errors in case the files aren't found in the folder.
    // The datafile can also be supplied directly in the config.
    if(datafilename.length() == 0) {
        LOG(ERROR) << "No data file was found for ALiBaVa in " << input_directory;
        return;
    }
    if(pedestalfilename.length() == 0) {
        LOG(WARNING) << "No pedestal file was found."
                     << "\n"
                     << "Datafile will be used for pedestal";
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
    ALiBaVaPointer = DataFileRoot::OpenFile(datafilename.c_str());
    // Sort vector to avoid errors later on
    std::sort(roi.begin(), roi.end());
    // Set the region of interest
    ALiBaVaPointer->set_ROI(roi);
    // Get the vector with all ROI channels listed to use here
    m_roi_ch = ALiBaVaPointer->get_ROI();

    ALiBaVaPointer->set_polarity(polarity);

    const char* ped_f = "alibava_ped.ped";
    const char* cal_f = "alibava_cal.cal";
    // Create a pointer with the pedestal file
    DataFileRoot* PedestalPointer = DataFileRoot::OpenFile(pedestalfilename.c_str());
    PedestalPointer->set_ROI(roi);

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
    ALiBaVaPointer->set_timecut(timecut_low, timecut_up);

    // Ignore the first X events to ensure synchronisation, default is X = 1 which ignores the first event.
    for(int ievt = 0; ievt < ignore_events; ievt++) {
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
        // Still need to figure this out... I think it's just end of run. Need to see difference between HDF5 and binary
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

    double max_signal = 0;
    // This loops over the channels in the current ALiBaVa event
    for(int chan : m_roi_ch) {
        double ADCSignal = ALiBaVaPointer->ADC_signal(chan);
        double SNRatio = ALiBaVaPointer->sn(chan);
        double CalSignal = ADCSignal * m_calibration_constant;

        if(ADCSignal > max_signal) {
            max_signal = ADCSignal;
        }

        // The chargecut is applied here
        // Not needed anymore but left in for now, other wise stuff explodes in Correlation
        if(CalSignal > m_chargecut) {
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
    hTimeProfile->Fill(TDCTime, max_signal, 1);
    // Put the created vector of pixels on the clipboard.
    clipboard->putData(pixels, detector_->getName());

    // If the pixels vector is empty, report this to Corryvreckan
    if(pixels.empty()) {
        return StatusCode::NoData;
    }

    // Report the end of this event to Corryvreckan
    return StatusCode::Success;
}

void EventLoaderALiBaVa::finalize(const std::shared_ptr<ReadonlyClipboard>&) {
    ALiBaVaPointer->close();
    delete ALiBaVaPointer;
}
