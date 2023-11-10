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
    config_.setDefault<std::string>("charge_calibration_formula", "1.0");
    config_.setDefault<double>("chargecut", 0);
    config_.setDefault<int>("polarity", -1);

    std::string input_directory = config_.getPath("input_directory");
    int run = config_.get<int>("run");
    double timecut_low = config_.get<double>("timecut_low");
    double timecut_up = config_.get<double>("timecut_up");
    int ignore_events = config_.get<int>("ignore_events");
    chargecut_ = config_.get<double>("chargecut");
    int polarity = config_.get<int>("polarity");

    charge_calibration_formula_ = parse_formula();

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
    hChargeSignal = new TH1F("chargeSignal", "Charge of signal; Charge [e]; # entries", 200, -2000, 50000);

    hADCSignal = new TH1F("ADCSignal", "ADC of signal; Signal [ADC]; # entries", 140, -20.5, 119.5);

    hSNR = new TH1F("SNRatio", "Signal to noise ratio; SNRatio; # entries", 100, -50, 200);

    hPedestal = new TH1F("pedestal", "Pedestal; # channel; Pedestal [ADC]", 256, -0.5, 255.5);

    hPedestalCorrect = new TH1F("pedestalCorrect", "Corrected pedestal; # channel; Pedestal [ADC]", 256, -0.5, 255.5);

    hNoise = new TH1F("noise", "Uncorrected Noise; # channel; Noise [ADC]", 256, -0.5, 255.5);

    hNoiseCorrect = new TH1F("noiseCorrect", "Corrected Noise; # channel; Noise [ADC]", 256, -0.5, 255.5);

    hNoiseElectrons = new TH1F("noiseElectrons", "Corrected Noise in electrons; # channel; Noise [e]", 256, -0.5, 255.5);

    hTimeProfile =
        new TProfile("timeProfile", "Time profile; Time [ns]; Ave. signal highest channel [ADC]", 35, 0, 35, 0, 200);

    hPedestalCorrect2D = new TH2F(
        "pedestalCorrect2D", "Corrected pedestal in 2D; # columns; # rows; Pedestal [ADC]", 256, -0.5, 255.5, 1, -0.5, 0.5);

    hNoiseCorrect2D = new TH2F(
        "noiseCorrect2D", "Corrected Noise in 2D; # columns; # rows; Pedestal [ADC]", 256, -0.5, 255.5, 1, -0.5, 0.5);

    hNoiseElectrons2D = new TH2F("noiseElectrons2D",
                                 "Corrected Noise in 2D in electrons; # columns; # rows; Noise [e]",
                                 256,
                                 -0.5,
                                 255.5,
                                 1,
                                 -0.5,
                                 0.5);

    hTemperature =
        new TProfile("tempProfile", "Temperature vs Time; Time [s]; Temperature [K]", 120, 0, 7200, 223.15, 323.15);

    // Create a shared pointer with the data file.
    alibava_.reset(DataFileRoot::OpenFile(datafilename.c_str()));

    // Find all non masked channels from the detector config and put them into a vector
    for(unsigned int col = 0; col < 256; col++) {
        if(!detector_->masked(static_cast<int>(col), 0)) {
            roi_ch_.push_back(col);
        }
    }

    // Set the region of interest
    alibava_->set_ROI(roi_ch_);

    // Set the polarity of the signal
    alibava_->set_polarity(polarity);

    // Create a pointer with the pedestal file
    DataFileRoot* PedestalPointer = DataFileRoot::OpenFile(pedestalfilename.c_str());
    PedestalPointer->set_ROI(roi_ch_);

    // Calculate the pedestals, and compute and apply the common mode noise correction
    PedestalPointer->compute_pedestals_alternative();

    for(auto chan : roi_ch_) {
        double ped_val, noise_val;
        ped_val = PedestalPointer->ped(chan);
        noise_val = PedestalPointer->noise(chan);
        hPedestal->SetBinContent(static_cast<int>(chan + 1), ped_val);
        hNoise->SetBinContent(static_cast<int>(chan + 1), noise_val);
    }

    PedestalPointer->compute_cmmd_alternative();
    double mean_ped_temp = PedestalPointer->get_mean_pedestal_temp();
    for(auto chan : roi_ch_) {
        double ped_val, noise_val;
        ped_val = PedestalPointer->ped(chan);
        noise_val = PedestalPointer->noise(chan);
        hPedestalCorrect->SetBinContent(static_cast<int>(chan + 1), ped_val);
        hNoiseCorrect->SetBinContent(static_cast<int>(chan + 1), noise_val);
        hNoiseElectrons->SetBinContent(static_cast<int>(chan + 1),
                                       noise_val * charge_calibration_formula_->Eval(mean_ped_temp));
        hPedestalCorrect2D->SetBinContent(static_cast<int>(chan + 1), 1, ped_val);
        hNoiseCorrect2D->SetBinContent(static_cast<int>(chan + 1), 1, noise_val);
        hNoiseElectrons2D->SetBinContent(
            static_cast<int>(chan + 1), 1, noise_val * charge_calibration_formula_->Eval(mean_ped_temp));
    }

    // Save the calculated pedestal information in a temporary file
    const std::string ped_f = "alibava_ped.ped";
    PedestalPointer->save_pedestals(ped_f.c_str());
    PedestalPointer->close();
    delete PedestalPointer;
    // Load the calculated pedestal info into the original datafile
    alibava_->load_pedestals(ped_f.c_str(), kTRUE);

    // Set the timecuts
    alibava_->set_timecut(timecut_low, timecut_up);

    // Ignore the first X events to ensure synchronisation, default is X = 1 which ignores the first event.
    for(int ievt = 0; ievt < ignore_events; ievt++) {
        alibava_->read_event();
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
    int return_code = alibava_->read_event();

    if(return_code == 1) {
        LOG(DEBUG) << "Successfully read event from ALiBaVa file";
    } else if(return_code == -1) {
        LOG(DEBUG) << "Reached end of the ALiBaVa file, requesting end of run";
        return StatusCode::EndRun;
    } else {
        throw ModuleError("Issue with ALiBaVa file, return code " + std::to_string(return_code));
    }

    // Calculate the common mode for the signal in this event
    alibava_->calc_common_mode_signal();
    // Process the opened data event, i.e. pedestal correction, common mode noise correction
    alibava_->process_event();
    // This gets the TDC time from the event, allowing timecuts around the event peak
    // The timecut is set in the ALiBaVa_loader() function.
    double TDCTime = alibava_->time();
    if(!alibava_->valid_time(TDCTime)) {
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

    // Get temperature and convert it into Kelvin
    double temp = alibava_->temp() + 273.15;
    hTemperature->Fill(static_cast<double>(Units::convert(event->start(), "s")), temp, 1);

    double max_signal = 0;
    // This loops over the channels in the current ALiBaVa event
    for(auto chan : roi_ch_) {
        double ADCSignal = alibava_->ADC_signal(chan);
        double SNRatio = alibava_->sn(chan);
        double CalSignal = ADCSignal * charge_calibration_formula_->Eval(temp);

        if(ADCSignal > max_signal) {
            max_signal = ADCSignal;
        }

        // The chargecut is applied here
        // Not needed anymore but left in for now, other wise stuff explodes in Correlation
        if(CalSignal > chargecut_) {
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
    alibava_->close();
    // delete alibava_;
    alibava_.reset();
}

std::shared_ptr<TFormula> EventLoaderALiBaVa::parse_formula() {

    const std::string formula_title = "charge_calibration_formula";
    auto function = config_.get<std::string>(formula_title);

    // Create the formula for the charge calibration
    std::shared_ptr<TFormula> formula;
    formula = std::make_shared<TFormula>(formula_title.c_str(), function.c_str(), false);

    // Check if it is a valid formula
    if(!formula->IsValid()) {
        throw InvalidValueError(config_, formula_title, "Invalid formula");
    } else {
        LOG(DEBUG) << "Found charge calibration formula";
    }

    // Check if we have the correct number of dimensions
    if(formula->GetNdim() > 1) {
        throw InvalidValueError(config_, formula_title, "Invalid number of dimensions, only 1d is supported");
    }

    auto n_params = static_cast<size_t>(formula->GetNpar());
    if(n_params != 0) {
        LOG(DEBUG) << "Found charge calibration formula with " << n_params << " parameters";

        auto parameters = config_.getArray<double>("charge_calibration_parameters");
        if(n_params != parameters.size()) {
            throw InvalidValueError(
                config_,
                "charge_calibration_parameters",
                "The number of provided parameters does not line up with the sum of parameters in the function.");
        }

        for(auto i = 0; i < formula->GetNpar(); i++) {
            formula->SetParameter(i, parameters.at(static_cast<size_t>(i)));
        }
    }
    return formula;
}