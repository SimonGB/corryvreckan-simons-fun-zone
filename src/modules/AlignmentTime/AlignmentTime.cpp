/**
 * @file
 * @brief Implementation of module AlignmentTime
 *
 * @copyright Copyright (c) 2017-2022 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "AlignmentTime.h"

using namespace corryvreckan;

AlignmentTime::AlignmentTime(Configuration& config, std::shared_ptr<Detector> detector)
    : Module(config, std::move(detector)) {

    // Check if the user wants to directly apply the determined correction
    config_.setDefault<bool>("update_time_offset", false);
    update_time_offset = config_.get<bool>("update_time_offset");

    // Get the scan parameters
    // TODO: Should I make better use of the configuration class?
    config_.setDefault<double>("shift_start", 1);
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
    LOG(INFO) << "Configured to scan from shift_start = " << Units::display(shift_start_, {"s", "ms", "us"});
    LOG(INFO) << "to shift_end = " << Units::display(shift_end_, {"s", "ms", "us"});
    LOG(INFO) << "in shift_n = " << shift_n_ << " steps.";
    if(shift_start_ > shift_end_ || shift_n_ == 0) {
        shift_user_ = false;
        LOG(INFO) << "Attempting to guess reasonable scan parameters.";
    }
    time_user_ = true;
    LOG(INFO) << "Using time_scale = " << Units::display(shift_end_, {"s", "ms", "us"});
    LOG(INFO) << "and time_nbins = " << time_nbins_;
    if(time_scale_ == 0 || time_nbins_ == 0) {
        time_user_ = false;
        LOG(INFO) << "Attempting to guess reasonable range and binning.";
    }
}

void AlignmentTime::initialize() {

    // Get the name of the detector used as time reference
    std::shared_ptr<Detector> reference = get_reference();
    config_.setDefault<std::string>("time_reference_name", reference->getName());
    time_reference_name_ = config_.get<std::string>("time_reference_name");
    LOG(INFO) << "Using " << time_reference_name_ << " as reference.";

    timestamps_[time_reference_name_] = {};

    // Reference time stamp histograms
    std::string title = time_reference_name_ + ";pixel timestamps [ms]; # entries";
    hTimeStampsRef = new TH1D("timeStampsRef", title.c_str(), 3e6, 0, 3e3);
    title = time_reference_name_ + ";pixel timestamps [s]; # entries";
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
        auto pixels = clipboard->getData<Pixel>(time_reference_name_);
        if(pixels.empty()) {
            LOG(DEBUG) << "Detector " << time_reference_name_ << " does not have any pixels on the clipboard";
        }
        LOG(DEBUG) << "Picked up " << pixels.size() << " pixels for device " << time_reference_name_;

        // Iterate pixels
        for(size_t iP = 0; iP < pixels.size(); iP++) {
            Pixel* pixel = pixels[iP].get();
            // Fill the timestamps of this pixel into container
            timestamps_[time_reference_name_].emplace_back(pixel->timestamp());
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
    if(timestamps_[time_reference_name_].size() == 0) {
        LOG(ERROR) << "No timestamps found for time reference";
        return;
    }
    // Loop over time stamps in reference detector
    for(auto ts : timestamps_[time_reference_name_]) {
        hTimeStampsRef->Fill(static_cast<double>(Units::convert(ts, "ms")));
        hTimeStampsRef_long->Fill(static_cast<double>(Units::convert(ts, "s")));
    }

    // Loop over all other detectors
    for(auto& detector : get_detectors()) {
        // Get the detector name
        std::string detectorName = detector->getName();
        LOG(DEBUG) << "Detector with name " << detectorName;

        if(timestamps_[detectorName].size() == 0) {
            LOG(ERROR) << "No timestamps found for " << detectorName;
            continue;
        }

        LOG(DEBUG) << "Filling time stamp histograms";
        for(auto ts : timestamps_[detectorName]) {
            hTimeStamps[detectorName]->Fill(static_cast<double>(Units::convert(ts, "ms")));
            hTimeStamps_long[detectorName]->Fill(static_cast<double>(Units::convert(ts, "s")));
        }

        // calculate final scan parameters and perform scan
        calculate_parameters(detectorName);
        scan_delay(detectorName);
        if(update_time_offset) {
            find_delay(detectorName);
        }
    }

    LOG(DEBUG) << "Analysed " << m_eventNumber << " events";
}

// Calculating parameters from user input, or guess.
void AlignmentTime::calculate_parameters(std::string detectorName) {

    // Steps need to be smaller than the trigger period.
    // Take difference between first and last and divide by size.
    double start = timestamps_[detectorName].at(0);
    double end = timestamps_[detectorName].at(timestamps_[detectorName].size() - 1);
    double period = (end - start) / static_cast<double>(timestamps_[detectorName].size());

    if(shift_user_) {
        shift_step_ = (shift_end_ - shift_start_) / static_cast<double>(shift_n_);
    } else {
        // Calculate period and make it 20 times smaller (my guess is as good as yours).
        shift_step_ = period / 20.;
        // Lets make 200 steps around 0.
        shift_n_ = 200;
        shift_start_ = -shift_step_ * static_cast<double>(shift_n_) / 2.;
        shift_end_ = shift_step_ * static_cast<double>(shift_n_) / 2.;
        // And tell the world.
        LOG(INFO) << "Calculated to scan from shift_start = " << Units::display(shift_start_, {"s", "ms", "us"});
        LOG(INFO) << "to shift_end = " << Units::display(shift_end_, {"s", "ms", "us"});
        LOG(INFO) << "in shift_n = " << shift_n_ << " steps.";
    }

    if(!time_user_) {
        // Same as range above
        time_scale_ = period * 5.;
        time_nbins_ = 200;
        // Tell the world.
        LOG(INFO) << "Using calculated time scale of time_scale = " << Units::display(time_scale_, {"s", "ms", "us"});
        LOG(INFO) << "and time_nbins = " << time_nbins_;
    }

    if(shift_n_ * time_nbins_ > 1e6) {
        LOG(WARNING) << "Using large number of bins in 2D histogram: shift_n = " << shift_n_
                     << " times time_nbins = " << time_nbins_;
        LOG(WARNING) << "This might cause crashes if there is not enough memory. Consider adjustment!";
    }

    return;
}

// Scan delay
void AlignmentTime::scan_delay(std::string detectorName) {

    LOG(INFO) << "Starting delay scan";

    // Create histogram
    std::string title = detectorName + ";time shift [ms]; #Deltat [ms]; # entries";
    hResidualVsShift[detectorName] = new TH2D("hResidualVsShift",
                                              title.c_str(),
                                              shift_n_,
                                              shift_start_ / 1e6 - shift_step_ / 2e6,
                                              shift_end_ / 1e6 - shift_step_ / 2e6,
                                              time_nbins_,
                                              -time_scale_ / 1e6,
                                              time_scale_ / 1e6);

    // Scanning the shift
    uint64_t counter = 0;
    for(auto shift = shift_start_; shift < shift_end_; shift += shift_step_) {
        // Satisfy my impatiance
        if(0 == counter % 10) {
            LOG(DEBUG) << "  testing shift " << Units::display(shift, {"s", "ms", "us"});
        }

        // Iterate hits in the detector
        for(auto detector_ts : timestamps_[detectorName]) {
            // Apply shift
            auto detector_ts_shifted = detector_ts + shift;
            // Calculate difference between shifted ts and best matching reference ts.
            auto residual = detector_ts_shifted - find_closest(timestamps_[time_reference_name_], detector_ts_shifted);
            hResidualVsShift[detectorName]->Fill(static_cast<double>(Units::convert(shift, "ms")),
                                                 static_cast<double>(Units::convert(residual, "ms")));
        }
        counter++;
    }

    return;
}

// Find delay
void AlignmentTime::find_delay(std::string detectorName) {

    LOG(INFO) << "Trying to estimate best delay, i.e. the maximum in hResidualVsShift";

    // If the scan parameters are good,
    // the maximum of the histogram indicates the right shift
    int max, tmp;
    hResidualVsShift[detectorName]->GetBinXYZ(hResidualVsShift[detectorName]->GetMaximumBin(), max, tmp, tmp);
    double best_shift = hResidualVsShift[detectorName]->GetXaxis()->GetBinCenter(max);

    // Histogram is in ms, convert back do base unit
    best_shift /= static_cast<double>(Units::convert(1., "ms"));

    // Now update detector, to adjust geometry file
    auto detector = get_detector(detectorName);
    LOG(INFO) << "Updating time offset for detector " << detectorName;
    LOG(INFO) << "Old: " << Units::display(detector->timeOffset(), {"s", "ms", "us"});
    detector->setTimeOffset(detector->timeOffset() + best_shift);
    LOG(INFO) << "New: " << Units::display(detector->timeOffset(), {"s", "ms", "us"}) << " with best shift of "
              << Units::display(best_shift, {"s", "ms", "us"});

    return;
}

// Modified Binary search algorithm adapted from
// https://www.geeksforgeeks.org/cpp-program-to-find-closest-number-in-array/
// Time Complexity: O(log(n))
// Auxiliary Space: O(log(n)) (implicit stack is created due to recursion)
// Returns the array element closest to the target value.
double AlignmentTime::find_closest(std::vector<double> const& arr, double target) {
    uint64_t n = arr.size();

    // If we hit the edge, by chance
    if(target <= arr[0]) {
        return arr[0];
    }
    if(target >= arr[n - 1]) {
        return arr[n - 1];
    }

    // Helper for findClosest as lambda expression.
    // Compares two values to target and returns the closer one.
    auto which_closer = [](double val1, double val2, double targ) {
        if(targ - val1 >= val2 - targ) {
            return val2;
        } else {
            return val1;
        }
    };

    // Doing binary search
    uint64_t i = 0, j = n, mid = 0;
    while(i < j) {
        mid = (i + j) / 2;

        // Return if we hit
        if(arr[mid] == target) {
            return arr[mid];
        }

        // If target is less than array element, then search in left
        if(target < arr[mid]) {
            // If target is greater than previous to mid, return closest of two
            if(mid > 0 && target > arr[mid - 1]) {
                return which_closer(arr[mid - 1], arr[mid], target);
            }
            // Iteratively repeat for left half
            j = mid;
        }
        // If target is greater than mid, then search right in the same way
        else {
            if(mid < n - 1 && target < arr[mid + 1]) {
                return which_closer(arr[mid], arr[mid + 1], target);
            }
            i = mid + 1;
        }
    } // While

    // Only single element left after search
    return arr[mid];
}
