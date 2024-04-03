/**
 * @file
 * @brief Implementation of module TextWriter
 *
 * @copyright Copyright (c) 2019-2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "TextWriter.h"

using namespace corryvreckan;

TextWriter::TextWriter(Configuration& config, std::vector<std::shared_ptr<Detector>> detectors)
    : Module(config, std::move(detectors)) {}

void TextWriter::initialize() {

    // Create output file
    config_.setDefault<std::string>("file_name", "data");
    output_file_name_ = createOutputFile(config_.get<std::string>("file_name"), "txt", true);
    output_file_ = std::make_unique<std::ofstream>(output_file_name_);

    *output_file_ << "# Corryvreckan ASCII data" << std::endl << std::endl;

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

StatusCode TextWriter::run(const std::shared_ptr<Clipboard>& clipboard) {

    if(!clipboard->isEventDefined()) {
        ModuleError("No Clipboard event defined, cannot continue");
    }

    // Print the current event:
    *output_file_ << "=== " << m_eventNumber << " ===" << std::endl;
    *output_file_ << *clipboard->getEvent() << std::endl;

    auto data = clipboard->getAll();
    LOG(DEBUG) << "Clipboard has " << data.size() << " different object types.";

    for(auto& block : data) {
        try {
            auto type_idx = block.first;
            auto class_name = corryvreckan::demangle(type_idx.name());
            auto class_name_full = corryvreckan::demangle(type_idx.name(), true);
            LOG(TRACE) << "Received objects of type \"" << class_name << "\"";

            // Check if these objects should be read
            if((!include_.empty() && include_.find(class_name) == include_.end()) ||
               (!exclude_.empty() && exclude_.find(class_name) != exclude_.end())) {
                LOG(TRACE) << "Ignoring object " << corryvreckan::demangle(type_idx.name())
                           << " because it has been excluded or not explicitly included";
                continue;
            }

            for(auto& detector_block : block.second) {
                // Get the detector name
                std::string detector_name;
                if(!detector_block.first.empty()) {
                    detector_name = detector_block.first;
                } else {
                    detector_name = "<global>";
                }

                *output_file_ << "--- " << detector_name << " ---" << std::endl;

                auto objects = std::static_pointer_cast<ObjectVector>(detector_block.second);
                for(auto& object : *objects) {
                    *output_file_ << *object << std::endl;
                }
            }
        } catch(...) {
            LOG(WARNING) << "Cannot process object of type" << corryvreckan::demangle(block.first.name());
            return StatusCode::NoData;
        }
    }

    // Increment event counter
    m_eventNumber++;

    // Return value telling analysis to keep running
    return StatusCode::Success;
}

void TextWriter::finalize(const std::shared_ptr<ReadonlyClipboard>&) {

    LOG(DEBUG) << "Analysed " << m_eventNumber << " events";
}
