/**
 * @file
 * @brief Implementation of module EventLoaderALiBaVa
 *
 * @copyright Copyright (c) 2021-2022 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "EventLoaderALiBaVa.h"
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <dirent.h>

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
    config_.setDefault<double>("calibration_constant", 1.0);
    config_.setDefault<double>("chargecut", 0);
    config_.setDefault<int>("polarity", -1);

    std::string input_directory = config_.getPath("input_directory");
    int run = config_.get<int>("run");
    double timecut_low = config_.get<double>("timecut_low");
    double timecut_up = config_.get<double>("timecut_up");
    int ignore_events = config_.get<int>("ignore_events");
    m_calibration_constant = config_.get<double>("calibration_constant");
    m_chargecut = config_.get<double>("chargecut");
    int polarity = config_.get<int>("polarity");

    // Check if input directory exists
    std::filesystem::path directory = input_directory;
    if(!std::filesystem::exists(directory)) {
        throw InvalidValueError(config_, "input_directory", "The directory could not be found");
    }

    std::string datafilename;
    std::string pedestalfilename;

    // Read the run-files (data, pedestal and calibration) in the folder
    for(auto const& dir_entry : std::filesystem::directory_iterator{directory}) {
        std::string entryName = dir_entry.path();

        if(entryName.find(std::to_string(run) + ".dat") != std::string::npos ||
           entryName.find("dat_run" + std::to_string(run) + ".hdf") != std::string::npos) {
            datafilename = entryName;
        }
        if(entryName.find(std::to_string(run) + ".ped") != std::string::npos ||
           entryName.find("ped_run" + std::to_string(run) + ".hdf") != std::string::npos) {
            pedestalfilename = entryName;
        }
    }

    // Log errors in case the files aren't found in the folder.
    // The datafile can also be supplied directly in the config.
    if(datafilename.length() == 0) {
        throw InvalidValueError(config_, "run", "No valid ALiBaVa data file with this run number was found in directory");
        return;
    }
    if(pedestalfilename.length() == 0) {
        LOG(WARNING) << "No pedestal file was found." << std::endl << "Datafile will be used for pedestal";
        pedestalfilename = datafilename;
    }

    // Create histograms
    hChargeSignal = new TH1F("chargeSignal", "Charge of signal; Charge [e]; # entries", 100, -96000, 32000);

    hADCSignal = new TH1F("ADCSignal", "ADC of signal; Signal [ADC]; # entries", 100, -200, 600);

    hSNR = new TH1F("SNRatio", "Signal to noise ratio; SNRatio; # entries", 100, -50, 200);

    hPedestal = new TH1F("pedestal", "Uncorrected pedestal; # channel; Pedestal[ADC]", 256, -0.5, 255.5);

    hPedestalCorrect = new TH1F("pedestalCorrect", "Corrected pedestal; # channel; Pedestal[ADC]", 256, -0.5, 255.5);

    hNoise = new TH1F("noise", "Uncorrected Noise; # channel; Noise[ADC]", 256, -0.5, 255.5);

    hNoiseCorrect = new TH1F("noiseCorrect", "Corrected Noise; # channel; Noise[ADC]", 256, -0.5, 255.5);

    hTimeProfile =
        new TProfile("timeProfile", "Time profile; Time [ns]; Ave. signal highest channel [ADC]", 35, 0, 35, 0, 200);

    hPedestalCorrect2D = new TH2F(
        "pedestalCorrect2D", "Corrected pedestal in 2D; # columns; # rows; Pedestal[ADC]", 256, -0.5, 255.5, 1, -0.5, 0.5);

    hNoiseCorrect2D = new TH2F(
        "noiseCorrect2D", "Corrected Noise in 2D; # columns; # rows; Pedestal[ADC]", 256, -0.5, 255.5, 1, -0.5, 0.5);

    // Create a shared pointer with the data file.
    m_alibava.reset(DataFileRoot::OpenFile(datafilename.c_str()));

    // Find all non masked channels from the detector config and put them into a vector
    for(unsigned int col = 0; col < 256; col++) {
        if(!detector_->masked(static_cast<int>(col), 0)) {
            m_roi_ch.push_back(col);
        }
    }

    // Set the region of interest
    m_alibava->set_ROI(m_roi_ch);

    // Set the polarity of the signal
    m_alibava->set_polarity(polarity);

    // Create a pointer with the pedestal file
    DataFileRoot* PedestalPointer = DataFileRoot::OpenFile(pedestalfilename.c_str());
    PedestalPointer->set_ROI(m_roi_ch);

    // Calculate the pedestals, and compute and apply the common mode noise correction
    PedestalPointer->compute_pedestals_alternative();

    for(auto chan : m_roi_ch) {
        double ped_val, noise_val;
        ped_val = PedestalPointer->ped(chan);
        noise_val = PedestalPointer->noise(chan);
        hPedestal->SetBinContent(static_cast<int>(chan + 1), ped_val);
        hNoise->SetBinContent(static_cast<int>(chan + 1), noise_val);
    }

    PedestalPointer->compute_cmmd_alternative();
    for(auto chan : m_roi_ch) {
        double ped_val, noise_val;
        ped_val = PedestalPointer->ped(chan);
        noise_val = PedestalPointer->noise(chan);
        hPedestalCorrect->SetBinContent(static_cast<int>(chan + 1), ped_val);
        hNoiseCorrect->SetBinContent(static_cast<int>(chan + 1), noise_val);
        hPedestalCorrect2D->SetBinContent(static_cast<int>(chan + 1), 1, ped_val);
        hNoiseCorrect2D->SetBinContent(static_cast<int>(chan + 1), 1, noise_val);
    }

    // Save the calculated pedestal information in a temporary file
    const std::string ped_f = "alibava_ped.ped";
    PedestalPointer->save_pedestals(ped_f.c_str());
    PedestalPointer->close();
    delete PedestalPointer;
    // Load the calculated pedestal info into the original datafile
    m_alibava->load_pedestals(ped_f.c_str(), kTRUE);

    // Set the timecuts
    m_alibava->set_timecut(timecut_low, timecut_up);

    // Ignore the first X events to ensure synchronisation, default is X = 1 which ignores the first event.
    for(int ievt = 0; ievt < ignore_events; ievt++) {
        m_alibava->read_event();
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
    int return_code = m_alibava->read_event();

    if(return_code == 1) {
        LOG(DEBUG) << "Successfully read event from ALiBaVa file";
    } else if(return_code == -1) {
        LOG(DEBUG) << "Reached end of the ALiBaVa file, requesting end of run";
        return StatusCode::EndRun;
    } else {
        throw ModuleError("Issue with ALiBaVa file, return code " + std::to_string(return_code));
    }

    // Calculate the common mode for the signal in this event
    m_alibava->calc_common_mode_signal();
    // Process the opened data event, i.e. pedestal correction, common mode noise correction
    m_alibava->process_event();
    // This gets the TDC time from the event, allowing timecuts around the event peak
    // The timecut is set in the ALiBaVa_loader() function.
    double TDCTime = m_alibava->time();
    if(!m_alibava->valid_time(TDCTime)) {
        LOG(DEBUG) << "Event time of " << TDCTime << " ns outside of timecut limits; ignoring event";
        return StatusCode::DeadTime;
    }

    double trigger_ts;

    if(!clipboard->getEvent()->triggerList().empty()) {
        trigger_ts = clipboard->getEvent()->triggerList().begin()->second;
        LOG(DEBUG) << "Using trigger timestamp " << Units::display(trigger_ts, "us") << " as event timestamp.";

    } else {
        trigger_ts = 0;
        LOG(DEBUG) << "No trigger timestamp setting 0 ns as cluster timestamp.";
    }

    double max_signal = 0;
    // This loops over the channels in the current ALiBaVa event
    for(auto chan : m_roi_ch) {
        double ADCSignal = m_alibava->ADC_signal(chan);
        double SNRatio = m_alibava->sn(chan);
        double CalSignal = ADCSignal * m_calibration_constant;

        if(ADCSignal > max_signal) {
            max_signal = ADCSignal;
        }

        // The chargecut is applied here
        // Not needed anymore but left in for now, other wise stuff explodes in Correlation
        if(CalSignal > m_chargecut) {
            // Create a pixel for every channel in this event with all the information and put it in the vector.
            // The value in the pixel reserved for the ADC value is used for the S/N ratio multiplied by 100000.

            std::shared_ptr<Pixel> pixel =
                std::make_shared<Pixel>(detector_->getName(), chan, 0, SNRatio * 100000, CalSignal, trigger_ts);

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
    m_alibava->close();
    // delete m_alibava;
    m_alibava.reset();
}
