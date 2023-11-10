/**
 * @file
 * @brief Implementation of module AhcalTextWriter
 *
 * @copyright Copyright (c) 2019-2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "AhcalTextWriter.h"

using namespace corryvreckan;

AhcalTextWriter::AhcalTextWriter(Configuration& config, std::vector<std::shared_ptr<Detector>> detectors)
    : Module(config, std::move(detectors)) {}

void AhcalTextWriter::initialize() {

    for(auto& detector : get_detectors()) {
        LOG(DEBUG) << "Initialise for detector " + detector->getName();
    }

    // Create output file
    config_.setDefault<std::string>("file_name", "data");
    config_.setDefault<uint32_t>("runnr", 0);
    config_.setDefault<double>("incident_z", Units::get<double>(2000, "mm"));

    m_incident_z = config_.get<double>("incident_z");
    m_runnr = config_.get<uint32_t>("runnr");

    LOG(DEBUG) << "parameter incident_z = " << m_incident_z;

    output_file_name_ = createOutputFile(config_.get<std::string>("file_name"), "txt", true);
    output_file_ = std::make_unique<std::ofstream>(output_file_name_);

    *output_file_ << "# AHCAL incident position from the track" << std::endl;
    *output_file_ << "#Runnr\tevt\ttrigID\tx\ty\tchi2ndof" << std::endl;

    // Read include and exclude list
    if(config_.has("include") && config_.has("exclude")) {
        throw InvalidValueError(config_, "exclude", "include and exclude parameter are mutually exclusive");
    } else if(config_.has("include")) {
        auto inc_arr = config_.getArray<std::string>("include");
        include_.insert(inc_arr.begin(), inc_arr.end());
    } else if(config_.has("exclude")) {
        auto exc_arr = config_.getArray<std::string>("exclude");
        exclude_.insert(exc_arr.begin(), exc_arr.end());
    }

    // Initialise member variables
    m_eventNumber = 0;
}

StatusCode AhcalTextWriter::run(const std::shared_ptr<Clipboard>& clipboard) {

    if(!clipboard->isEventDefined()) {
        ModuleError("No Clipboard event defined, cannot continue");
    }

    // Print the current event:
    //    *output_file_ << "=== " << m_eventNumber << " ===" << std::endl;
    //	*output_file_ << *clipboard->getEvent() <<  std::endl;
    auto trigger_multiplicity = clipboard->getEvent()->triggerList().size();
    uint32_t trigger_number;
    if(trigger_multiplicity == 1) {
        trigger_number = clipboard->getEvent()->triggerList().cbegin()->first;
        LOG(TRACE) << "Trigger number is " << trigger_number;
    }
    auto data = clipboard->getAll();
    LOG(DEBUG) << "Clipboard has " << data.size() << " different object types.";

    if(trigger_multiplicity != 1) {
        LOG(TRACE) << "skipping event " << m_eventNumber++ << ", because it has " << trigger_multiplicity << " triggers.";
        *output_file_ << "#skipping event " << m_eventNumber << " because of incorrect number of triggers ("
                      << trigger_multiplicity << "), tracks (" << clipboard->getData<Track>().size() << ")" << std::endl;
        m_eventNumber++;
        return StatusCode::Success;
    }

    auto tracks = clipboard->getData<Track>();
    if(tracks.size() != 1) {
        *output_file_ << "#skipping event " << m_eventNumber << " with trigger " << trigger_number
                      << " because of incorrect number of tracks (" << tracks.size() << ")" << std::endl;
        m_eventNumber++;
        return StatusCode::Success;
    }
    //	const auto track=tracks.cbegin();
    *output_file_ << m_runnr << "\t" << m_eventNumber << "\t" << trigger_number << "\t";
    for(auto& track : tracks) {
        auto globalIntercept = track->getIntercept(m_incident_z);
        *output_file_ << globalIntercept.X() << "\t" << globalIntercept.Y() << "\t";
        auto Chi2ndof = track->getChi2ndof();
        *output_file_ << Chi2ndof << "\t";
        *output_file_ << std::endl;
    }
    //	for (auto &block : data) {
    //		try {
    //			auto type_idx = block.first;
    //			auto class_name = corryvreckan::demangle(type_idx.name());
    //			auto class_name_full = corryvreckan::demangle(type_idx.name(), true);
    //			LOG(TRACE) << "Received objects of type \"" << class_name << "\", in full \"" << class_name_full << "\"";
    //
    //			//AHCAL print track incident position
    //			if (class_name.find("Track") != std::string::npos) {
    //				LOG(TRACE) << "Found track!";
    //				auto &tracks = block.second;
    //				LOG(TRACE) << "DEBUG size of the 'track' vector = " << tracks.size();
    //				const auto objvec = std::static_pointer_cast<ObjectVector>(block.second.cbegin()->second);
    //				if (objvec->size() != 1) {
    //					*output_file_ << "#skipping event " << m_eventNumber << " with trigger " << trigger_number << " because of
    //incorrect number of tracks ("
    //							<< objvec->size() << std::endl;
    //					LOG(DEBUG) << "skipping event with " << objvec->size() << " tracks. TrigN=" << trigger_number;
    //					continue;
    //				}
    //				LOG(TRACE) << "Trigger " << trigger_number << "has size of objvec = " << objvec->size();
    //
    //				const auto &detector_block = block.second.cbegin();
    //				const auto &obj = std::static_pointer_cast<ObjectVector>(block.second.cbegin()->second)->cbegin();
    //				LOG(TRACE) << "DEBUG obj print = " << *obj;
    ////				StraightLineTrack
    //			} else
    //				continue;
    //
    //			// Check if these objects should be read
    //			if ((!include_.empty() && include_.find(class_name) == include_.end()) || (!exclude_.empty() &&
    //exclude_.find(class_name) != exclude_.end())) { 				LOG(TRACE) << "Ignoring object " <<
    //corryvreckan::demangle(type_idx.name()) << " because it has been excluded or not explicitly included"; 				continue;
    //			}
    //			continue;
    //			for (auto &detector_block : block.second) {
    //				// Get the detector name
    //				std::string detector_name;
    //				if (!detector_block.first.empty()) {
    //					detector_name = detector_block.first;
    //				} else {
    //					detector_name = "<global>";
    //				}
    //				LOG(TRACE) << "Processing name " << detector_name;
    //
    //				*output_file_ << "--- " << detector_name << " ---" << std::endl;
    //
    //				auto objects = std::static_pointer_cast<ObjectVector>(detector_block.second);
    //				for (auto &object : *objects) {
    //					*output_file_ << *object << std::endl;
    //					LOG(TRACE) << "DEBUG object print " << *object;
    //				}
    //			}
    //		} catch (...) {
    //			LOG(WARNING) << "Cannot process object of type" << corryvreckan::demangle(block.first.name());
    //			return StatusCode::NoData;
    //		}
    //	}

    // Increment event counter
    m_eventNumber++;

    // Return value telling analysis to keep running
    return StatusCode::Success;
}

void AhcalTextWriter::finalize(const std::shared_ptr<ReadonlyClipboard>&) {

    LOG(DEBUG) << "Analysed " << m_eventNumber << " events";
}
