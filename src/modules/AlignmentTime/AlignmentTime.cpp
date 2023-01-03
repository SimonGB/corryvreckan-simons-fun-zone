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
    config_.setDefault<std::string>("reference_name", "noname");
    reference_name_ = config_.get<std::string>("reference_name");
    if(strcmp(reference_name_.c_str(), "noname")){
        LOG(WARNING) << "Module called without reference_name";
    }

}

void AlignmentTime::initialize() {

    for(auto& detector : get_detectors()) {
        std::string detectorName = detector->getName();
        LOG(DEBUG) << "Initialise for detector " + detectorName;

        timestamps_[detectorName] = {};

        // time stamp histograms
        std::string title = detectorName + ";pixel timestamps [ms]; # entries";
        hTimeStamps[detectorName] = new TH1D("timeStamps", title.c_str(), 3e6, 0, 3e3);
        title = detectorName + ";pixel timestamps [s]; # entries";
        hTimeStamps_long[detectorName] = new TH1D("timeStampsRef_long", title.c_str(), 3e6, 0, 3e3);
        title = reference_name_ + ";pixel timestamps [ms]; # entries";
        hTimeStampsRef[detectorName] = new TH1D("timeStampsRef", title.c_str(), 3e6, 0, 3e3);
        title = reference_name_ + ";pixel timestamps [s]; # entries";
        hTimeStampsRef_long[detectorName] = new TH1D("timeStampsRef_long", title.c_str(), 3e6, 0, 3e3);

    }

    // Initialise member variables
    m_eventNumber = 0;
}

StatusCode AlignmentTime::run(const std::shared_ptr<Clipboard>&) {

    // Loop over all detectors
    for(auto& detector : get_detectors()) {
        // Get the detector name
        std::string detectorName = detector->getName();
        LOG(DEBUG) << "Detector with name " << detectorName;
    }

    // Increment event counter
    m_eventNumber++;

    // Return value telling analysis to keep running
    return StatusCode::Success;
}

void AlignmentTime::finalize(const std::shared_ptr<ReadonlyClipboard>&) {

    LOG(DEBUG) << "Analysed " << m_eventNumber << " events";
}
