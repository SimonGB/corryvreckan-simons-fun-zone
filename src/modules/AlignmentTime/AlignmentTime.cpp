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

    // Get the name of the detector used as time reference
    config_.setDefault<std::string>("reference_name", "");
    reference_name_ = config_.get<std::string>("reference_name");
    if(reference_name_.empty()) {
        LOG(WARNING) << "Module called without reference_name. Output plots will be empty";
    }

    // Get the scan parameters
    // TODO: Should I make better use of the configuration class?
    config_.setDefault<double>("shift_start", 0);
    config_.setDefault<double>("shift_end", 0);
    config_.setDefault<int>("shift_n", 0);
    config_.setDefault<double>("time_scale", 0);
    config_.setDefault<int>("time_nbins", 0);
    shift_start_ = config_.get<double>("shift_start");
    shift_end_ = config_.get<double>("shift_end");
    shift_n_ = config_.get<int>("shift_n");
    time_scale_ = config_.get<double>("time_scale");
    time_nbins_ = config_.get<int>("time_nbins");

    // Checking user input
    shift_user_ = true;
    LOG(INFO) << "Configured to scan from shift_start = " << shift_start_;
    LOG(INFO) << "to shift_end = " << shift_end_;
    LOG(INFO) << "in shift_n = " << shift_n_ << " steps.";
    if(shift_start_ == shift_end_ || shift_n_ == 0) {
        shift_user_ = false;
        LOG(INFO) << "Attempting to guess reasonable scan parameters.";
    }
    time_user_ = true;
    LOG(INFO) << "Using time scale of time_scale = " << time_scale_;
    LOG(INFO) << "and time_nbins = " << time_nbins_;
    if(time_scale_ == 0 || time_nbins_ == 0) {
        time_user_ = false;
        LOG(INFO) << "Attempting to guess reasonable range and binning.";
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

    // check if reference timestamps are filled
    if(timestamps_[reference_name_].size() == 0) {
        LOG(ERROR) << "No timestamps found for time reference";
        return;
    }
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

        if(timestamps_[detectorName].size() == 0) {
            LOG(ERROR) << "No timestamps found for " << detectorName;
            // TODO: There are some histograms which will not be initialized. Will this cause crashes?
            continue;
        }
        for(auto ts : timestamps_[detectorName]) {
            hTimeStamps[detectorName]->Fill(static_cast<double>(Units::convert(ts, "ms")));
            hTimeStamps_long[detectorName]->Fill(static_cast<double>(Units::convert(ts, "s")));
        }

        // calculate final scan parameters and perform scan
        calculateParameters(detectorName);
        scanDelay(detectorName);
    }

    LOG(DEBUG) << "Analysed " << m_eventNumber << " events";
}

// Calculating parameters from user input, or guess.
void AlignmentTime::calculateParameters(std::string detectorName) {

    // Steps need to be smaller than the trigger period.
    // Take difference between first and last and divide by size.
    double start = timestamps_[detectorName].at(0);
    double end = timestamps_[detectorName].at(timestamps_[detectorName].size() - 1);
    double period = (end - start) / static_cast<double>(timestamps_[detectorName].size());

    if(shift_user_) {
        shift_step_ = (shift_end_ - shift_start_) / static_cast<double>(shift_n_);
    } else {
        // And make it 20 times smaller (my guess is as good as yours).
        shift_step_ = period / 20.;
        // Lets make 200 steps around 0.
        shift_n_ = 200;
        shift_start_ = -shift_step_ * static_cast<double>(shift_n_) / 2.;
        shift_end_ = shift_step_ * static_cast<double>(shift_n_) / 2.;
        // And tell the world.
        LOG(INFO) << "Calculated to scan from shift_start = " << shift_start_;
        LOG(INFO) << "to shift_end = " << shift_end_;
        LOG(INFO) << "in shift_n = " << shift_n_ << " steps.";
    }

    if(!time_user_) {
        // Same as range above
        time_scale_ = period * 5.;
        time_nbins_ = 200;
        // Tell the world.
        LOG(INFO) << "Using calculated time scale of time_scale = " << time_scale_;
        LOG(INFO) << "and time_nbins = " << time_nbins_;
    }

    // TODO: Should probably catch cases where user try to use overly large histograms.

    return;
}

// Scan delay
void AlignmentTime::scanDelay(std::string detectorName) {

    // Create histogram
    std::string title = detectorName + ";time shift [ms]; #Deltat [ms]; # entries";
    hResidualVsShift[detectorName] = new TH2D(
        "hResidualVsShift", title.c_str(), shift_n_, shift_start_, shift_end_, time_nbins_, -time_scale_, time_scale_);

    // Scanning the shift
    uint64_t counter = 0;
    for(auto shift = shift_start_; shift < shift_end_; shift += shift_step_) {
        // Satisfy my impatiance
        if(0 == counter % 10) {
            LOG(DEBUG) << "  testing shift " << shift;
        }

        // Iterate hits in the detector
        for(auto detector_ts : timestamps_[detectorName]) {
            // Apply shift
            auto detector_ts_shifted = detector_ts - shift;
            // Fill difference between best match w.r.t. reference time stamp
            hResidualVsShift[detectorName]->Fill(
                shift, detector_ts_shifted - findClosest(timestamps_[reference_name_], detector_ts_shifted));
        }
        counter++;
    }
    return;
}

// Modiefied Binary search algorithm adapted from
// https://www.geeksforgeeks.org/cpp-program-to-find-closest-number-in-array/
// Time Complexity: O(log(n))
// Auxiliary Space: O(log(n)) (implicit stack is created due to recursion)
// Returns the array element closest to the target value.
double AlignmentTime::findClosest(std::vector<double> const& arr, double target) {
    uint64_t n = arr.size();

    // If we hit the edge, by chance
    if(target <= arr[0])
        return arr[0];
    if(target >= arr[n - 1])
        return arr[n - 1];

    // Doing binary search
    uint64_t i = 0, j = n, mid = 0;
    while(i < j) {
        mid = (i + j) / 2;

        // Eeturn if we hit
        if(arr[mid] == target)
            return arr[mid];

        // If target is less than array element, then search in left
        if(target < arr[mid]) {
            // If target is greater than previous to mid, return closest of two
            if(mid > 0 && target > arr[mid - 1])
                return whichCloser(arr[mid - 1], arr[mid], target);
            // Iteratively repeat for left half
            j = mid;
        }
        // If target is greater than mid, then search right in the same way
        else {
            if(mid < n - 1 && target < arr[mid + 1])
                return whichCloser(arr[mid], arr[mid + 1], target);
            i = mid + 1;
        }
    } // While

    // Only single element left after search
    return arr[mid];
}
// Helper for findClosest.
// Compares two values to target and returns the closer one.
double AlignmentTime::whichCloser(double val1, double val2, double target) {
    if(target - val1 >= val2 - target)
        return val2;
    else
        return val1;
}
