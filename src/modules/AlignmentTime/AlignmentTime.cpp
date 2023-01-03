/**
 * @file
 * @brief Implementation of module AlignmentTime
 *
 * @copyright Copyright (c) 2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "AlignmentTime.h"

using namespace corryvreckan;

AlignmentTime::AlignmentTime(Configuration& config, std::shared_ptr<Detector> detector)
    : Module(config, std::move(detector)) {

    // get the name of the detector used as time reference
    config_.setDefault<std::string>("reference_name", "");
    reference_name_ = config_.get<std::string>("reference_name");
    if(reference_name_.empty()) {
        LOG(WARNING) << "Module called without reference_name " << reference_name_.c_str();
    }
}

void AlignmentTime::initialize() {

    timestamps_[reference_name_] = {};

    // Reference time stamp histograms
    std::string title = reference_name_ + ";pixel timestamps [ms]; # entries";
    hTimeStampsRef = new TH1D("timeStampsRef", title.c_str(), 3e6, 0, 3e3);
    title = reference_name_ + ";pixel timestamps [s]; # entries";
    hTimeStampsRef_long = new TH1D("timeStampsRef_long", title.c_str(), 3e6, 0, 3e3);

    // Iterate through detectors we want to align
    for(auto& detector : get_detectors()) {
        std::string detectorName = detector->getName();
        LOG(DEBUG) << "Initialise for detector " + detectorName;

        timestamps_[detectorName] = {};

        // Detector time stamp histograms
        title = detectorName + ";pixel timestamps [ms]; # entries";
        hTimeStamps[detectorName] = new TH1D("timeStamps", title.c_str(), 3e6, 0, 3e3);
        title = detectorName + ";pixel timestamps [s]; # entries";
        hTimeStamps_long[detectorName] = new TH1D("timeStamps_long", title.c_str(), 3e6, 0, 3e3);
    }

    // Initialise member variables
    m_eventNumber = 0;
}

StatusCode AlignmentTime::run(const std::shared_ptr<Clipboard>& clipboard) {

    reference_filled_ = 0;

    // Loop over all detectors
    for(auto& detector : get_detectors()) {
        // Get the detector name
        std::string detectorName = detector->getName();
        LOG(DEBUG) << "Detector with name " << detectorName;

        // Get all pixels for this detector
        auto pixels = clipboard->getData<Pixel>(detectorName);
        if(pixels.empty()) {
            LOG(DEBUG) << "Detector " << detectorName << " does not have any pixels on the clipboard";
        }
        LOG(DEBUG) << "Picked up " << pixels.size() << " pixels for device " << detectorName;

        // Iterate pixels
        for(size_t iP = 0; iP < pixels.size(); iP++) {
            Pixel* pixel = pixels[iP].get();
            // Fill the timestamps of this pixel into container
            timestamps_[detectorName].emplace_back(pixel->timestamp());
        }
    }

    // Fill also reference timestamps, once
    if(reference_filled_) {
        LOG(DEBUG) << "Reference detector already filled";
    } else {
        // Get all pixels for this detector
        auto pixels = clipboard->getData<Pixel>(reference_name_);
        if(pixels.empty()) {
            LOG(DEBUG) << "Detector " << reference_name_ << " does not have any pixels on the clipboard";
        }
        LOG(DEBUG) << "Picked up " << pixels.size() << " pixels for device " << reference_name_;

        // Iterate pixels
        for(size_t iP = 0; iP < pixels.size(); iP++) {
            Pixel* pixel = pixels[iP].get();
            // Fill the timestamps of this pixel into container
            timestamps_[reference_name_].emplace_back(pixel->timestamp());
        }

        reference_filled_ = 1;
    }

    // Increment event counter
    m_eventNumber++;

    // Return value telling analysis to keep running
    return StatusCode::Success;
}

void AlignmentTime::finalize(const std::shared_ptr<ReadonlyClipboard>&) {

    // Loop over time stamps in reference detector
    for(auto ts : timestamps_[reference_name_]) {
        hTimeStampsRef->Fill(static_cast<double>(Units::convert(ts, "ms")));
        hTimeStampsRef_long->Fill(static_cast<double>(Units::convert(ts, "s")));
    }

    // Loop over all detectors other
    for(auto& detector : get_detectors()) {
        // Get the detector name
        std::string detectorName = detector->getName();
        LOG(DEBUG) << "Detector with name " << detectorName;

        for(auto ts : timestamps_[detectorName]) {
            hTimeStamps[detectorName]->Fill(static_cast<double>(Units::convert(ts, "ms")));
            hTimeStamps_long[detectorName]->Fill(static_cast<double>(Units::convert(ts, "s")));
        }
    }

    LOG(DEBUG) << "Analysed " << m_eventNumber << " events";
}
