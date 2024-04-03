/**
 * @file
 * @brief Implementation of module EventLoaderMuPixTelescope
 *
 * @copyright Copyright (c) 2019-2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "EventLoaderMuPixTelescope.h"
#include <string>
#include "dirent.h"
#include "objects/Cluster.hpp"
#include "objects/Pixel.hpp"
#include "objects/Track.hpp"

#include "TDirectory.h"
using namespace corryvreckan;

EventLoaderMuPixTelescope::EventLoaderMuPixTelescope(Configuration& config, std::vector<std::shared_ptr<Detector>> detectors)
    : Module(config, std::move(detectors)), blockFile_(nullptr) {

    config_.setDefault<bool>("is_sorted", false);
    config_.setDefault<bool>("ts2_is_gray", false);
    config_.setDefault<unsigned>("buffer_depth", 1000);
    config_.setDefault<double>("time_offset", 0.0);
    config_.setDefault<double>("reference_frequency", 125.);
    config_.setDefault<uint>("nbits_timestamp", 10);
    config_.setDefault<uint>("nbits_tot", 6);
    config_.setDefault<uint>("ckdivend", 0);
    config_.setDefault<uint>("ckdivend2", 7);
    config_.setDefault<bool>("use_both_timestamps", false);

    use_both_timestamps_ = config_.get<bool>("use_both_timestamps");
    nbitsTS_ = config_.get<uint>("nbits_timestamp");
    nbitsToT_ = config_.get<uint>("nbits_tot");
    ckdivend_ = config_.get<uint>("ckdivend");
    ckdivend2_ = config_.get<uint>("ckdivend2");
    refFrequency_ = config_.get<double>("reference_frequency");
    inputDirectory_ = config_.getPath("input_directory");
    buffer_depth_ = config.get<unsigned>("buffer_depth");
    isSorted_ = config_.get<bool>("is_sorted");
    if(config.count({"run", "input_file"}) > 1) {
        throw InvalidCombinationError(config, {"run", "input_file"}, "run and input_file are mutually exclusive.");
    } else if(config_.has("input_file")) {
        input_file_ = config_.get<std::string>("input_file");
    } else {
        runNumber_ = config_.get<int>("run");
    }

    // simplifying calculations:
    multiplierToT_ = (1. + static_cast<double>(ckdivend2_)) / (1. + static_cast<double>(ckdivend_)) * 2.;
    // timestamp is calculated with respect to a 4ns base, tot wrt 8ns
    timestampMask_ = ((0x1) << nbitsTS_) - 1; // raw timestamp from data
    // timestamp after shift and 4ns base change
    timestampMaskExtended_ = ((0x1) << ((ckdivend_ + 1) * (nbitsTS_ + 1))) - 1;
    totMask_ = ((0x1) << nbitsToT_) - 1;
    clockToTime_ = 4. / refFrequency_ * 125.;
    maxToT_ =
        static_cast<double>(((totMask_ + 1) * static_cast<uint>(multiplierToT_)) & timestampMaskExtended_) * clockToTime_;
}

void EventLoaderMuPixTelescope::initialize() {
    // extract the tag from each detetcor name
    // find the corresponding detetctors:
    for(auto d : get_detectors()) {
        if(typeString_to_typeID.find(d->getType()) != typeString_to_typeID.end()) {
            LOG(INFO) << "----- Added detetctor:" << d->getName();
            detectors_.push_back(d);
        } else
            LOG(INFO) << "------ NOT Added detetctor:" << d->getName();
    }
    for(auto detector : detectors_) {
        std::string tag = detector->getName();
        if(tag.find("_") < tag.length())
            tag = tag.substr(tag.find("_") + 1);
        uint tag_ = uint(stoi(tag, nullptr, 16));
        tags_.push_back(tag_);
        LOG(INFO) << detector->getName() << " is using the fpga link tag " << std::hex << tag_;
        if(typeString_to_typeID.find(detector->getType()) == typeString_to_typeID.end()) {
            throw KeyValueParseError("tag " + std::to_string(tag_), "Sensor tag not supported");
        }

        removed_[tag_] = 0;
        stored_[tag_] = 0;
        counterHits_[tag_] = 0;
        // take the offset from the geometry file
        timeOffset_[tag_] = detector->timeOffset();
        types_[tag_] = typeString_to_typeID.at(detector->getType());
        LOG(INFO) << "Detector " << detector->getType() << "is assigned to type id " << types_.at(tag_);
        names_[tag_] = detector->getName();
        pixelbuffers_[tag_]; // create an empty map entry
        pixels_[tag_];
    }

    std::stringstream ss;
    if(input_file_.size() == 0) {
        ss << std::setw(6) << std::setfill('0') << runNumber_;
        std::string s = ss.str();
        input_file_ = "telescope_run_" + s + ".blck";
    }

    // check the if folder and file do exist
    dirent* entry;
    bool foundFile = false;
    DIR* directory = opendir(inputDirectory_.c_str());
    if(directory == nullptr) {
        throw MissingDataError("Cannot open directory: " + inputDirectory_);
    }
    while((entry = readdir(directory))) {
        if(entry->d_name == input_file_) {
            foundFile = true;
            break;
        }
    }
    if(!foundFile) {

        throw MissingDataError("Cannot open data file: " + input_file_);
    } else
        LOG(INFO) << "File " << input_file_ << " found";
    std::string file = (inputDirectory_ + "/" + entry->d_name);
    LOG(INFO) << "reading " << file;
    blockFile_ = new BlockFile(file);
    if(!blockFile_->open_read()) {
        throw MissingDataError("Cannot read data file: " + input_file_);
    }
    TDirectory* dir = getROOTDirectory();

    // create the histograms for all sensor
    for(auto& detector : detectors_) {
        auto name = detector->getName();
        TDirectory* local_directory = dir->mkdir(name.c_str());

        if(local_directory == nullptr) {
            throw RuntimeError("Cannot create or access local ROOT directory for module " + this->getUniqueName());
        }
        local_directory->cd();
        std::string title = name + "_hitMap; column; row";
        hHitMap[name] = new TH2F("hitMap",
                                 title.c_str(),
                                 detector->nPixels().x(),
                                 -.05,
                                 detector->nPixels().x() - .5,
                                 detector->nPixels().y(),
                                 -.05,
                                 detector->nPixels().y() - .5);
        title = name + "hitMap of out of event hits; column; row";
        hdiscardedHitmap[name] = new TH2F("discardedhitMap",
                                          title.c_str(),
                                          detector->nPixels().x(),
                                          -.05,
                                          detector->nPixels().x() - .5,
                                          detector->nPixels().y(),
                                          -.05,
                                          detector->nPixels().y() - .5);
        title = name + "pixelToT; ToT / ns;";
        hPixelToT[name] = new TH1F("pixelToT", title.c_str(), 2 * 2048, 10 * (-1024.5), 10 * 1023.5);
        title = name + "pixelTS; TS in clock cycles; ";
        hTimeStamp[name] = new TH1F("pixelTS", title.c_str(), 1024, -0.5, 1023.5);

        title = name + "ts_ToT; TS ToT in clock cycles; ";
        hts_ToT[name] = new TH1F("ts_ToT", title.c_str(), 1024, -0.5, 1023.5);

        title = name + "hHitsEvent; # hits per event; ";
        hHitsEvent[name] = new TH1F("hHitsEvent", title.c_str(), 300, -.5, 299.5);
        title = name + "hitsper1kevents; corry events /1k; hits per 1k events";
        hitsPerkEvent[name] = new TH1F("hHitsPerkEvent", title.c_str(), 1000, -.5, 999.5);
        title = name + "fpga vs chip clock;chip clock;fpga clock";
        raw_fpga_vs_chip[name] = new TH2F("raw_fpga_vs_chip", title.c_str(), 1024, 0, 1023, 2048, 0, 2047);
        title = name + "fpga vs chip clock;chip clock;fpga clock";
        raw_fpga_vs_chip_corrected[name] =
            new TH2F("raw_fpga_vs_chip_corrected", title.c_str(), 1024, 0, 1023, 2048, 0, 2047);

        title = name + "ts tst tot;ts1;ts tot";
        ts_TS1_ToT[name] = new TH2F("ts_ts1_tot", title.c_str(), 1024, 0, 1023, 2048, 0, 2047);
        title = "ts1_ts2";
        ts1_ts2[name] = new TH2F("ts1_ts2", title.c_str(), 1024, 0, 1023, 1024, 0, 1023);

        title = name + "Delay of chip events wrt. telescope frame;fpga clock@ chip clock 0;#events";
        chip_delay[name] = new TH1F("chip_delay", title.c_str(), 2048, -1023, 1023);
    }
}

void EventLoaderMuPixTelescope::finalize(const std::shared_ptr<ReadonlyClipboard>&) {

    for(auto d : tags_)
        LOG(INFO) << names_.at(d) << ": Number of hits put to clipboard: " << stored_.at(d)
                  << " and number of removed (not fitting in an event) hits: " << removed_.at(d);
    if(!isSorted_)
        LOG(INFO) << "Increasing the buffer depth might reduce this number.";
}

StatusCode EventLoaderMuPixTelescope::read_sorted(const std::shared_ptr<Clipboard>& clipboard) {
    PixelVector hits;
    if(!blockFile_->read_next(tf_)) {
        return StatusCode::EndRun;
    }
    for(uint i = 0; i < tf_.num_hits(); ++i) {
        RawHit h = tf_.get_hit(i);
        auto tag = h.tag() & uint(~0x3);
        if(names_.count(tag) != 1) {
            throw RuntimeError("Unknown pixel tag read in data: " + to_string(tag));
        }
        pixels_[tag].push_back(read_hit(h, tag, tf_.timestamp()));
    }
    // If no event is defined create one
    if(!clipboard->isEventDefined()) {
        // The readout FPGA creates frames with a length of 128 time stamps, each 8ns. The int division cuts of lowest bits
        int begin = int(pixels_.begin()->second.front()->timestamp()) / 1024;
        clipboard->putEvent(std::make_shared<Event>(double(begin * 1024), double((begin + 1) * 1024)));
    }
    return StatusCode::Success;
}

StatusCode EventLoaderMuPixTelescope::read_unsorted(const std::shared_ptr<Clipboard>& clipboard) {

    // Fill buffer
    fillBuffer();

    // check if an event is defined - if not we make one
    if(!clipboard->isEventDefined()) {
        auto minTS = 1e40;
        for(auto t : tags_) {
            if((pixelbuffers_.at(t).size() > 0) && (minTS > pixelbuffers_.at(t).top()->timestamp()))
                minTS = pixelbuffers_.at(t).top()->timestamp();
        }
        LOG(DEBUG) << "Defining Event around " << minTS;

        minTS = minTS - (static_cast<int>(minTS) % 8192);
        clipboard->putEvent(std::make_shared<Event>(minTS, minTS + 8192));
        LOG(DEBUG) << clipboard->getEvent()->start() << "\t" << clipboard->getEvent()->end();
    } else {
        LOG(TRACE) << "Event is defined";
    }

    for(auto t : tags_) {
        while(!pixelbuffers_.at(t).empty()) {
            // Always fill the buffer:
            fillBuffer();

            auto pixel = pixelbuffers_.at(t).top();
            auto event = clipboard->getEvent();
            auto position = event->getTimestampPosition(pixel->timestamp());

            if(position == Event::Position::BEFORE) {
                LOG(DEBUG) << " Old hit found: " << Units::display(pixel->timestamp(), "us") << " vs prev end ("
                           << eventNo_ - 1 << ")\t" << Units::display(prev_event_end_, "us") << " and current start \t"
                           << Units::display(clipboard->getEvent()->start(), "us")
                           << " and duration: " << Units::display(clipboard->getEvent()->duration(), {"us", "ns"})
                           << " and number of triggers: " << clipboard->getEvent()->triggerList().size();
                LOG(TRACE) << pixel->timestamp() << " < " << clipboard->getEvent()->start();
                removed_.at(t)++;
                hdiscardedHitmap.at(names_.at(t))->Fill(pixel->column(), pixel->row());
                pixelbuffers_.at(t).pop(); // remove top element
                continue;
            } else if(position == Event::Position::DURING) {
                LOG(DEBUG) << " Adding pixel hit: " << Units::display(pixel->timestamp(), "us") << " vs prev end ("
                           << eventNo_ - 1 << ")\t" << Units::display(prev_event_end_, "us") << " and current start \t"
                           << Units::display(clipboard->getEvent()->start(), "us")
                           << " and duration: " << Units::display(clipboard->getEvent()->duration(), {"us", "ns"});
                int col = pixel->column();
                int row = pixel->row();

                if(detectors_.back()->masked(col, row)) {
                    LOG(DEBUG) << "Masking pixel (col, row) = (" << col << ", " << row << ")" << std::endl;
                    removed_.at(t)++;
                    pixelbuffers_.at(t).pop();
                    continue;
                } else {
                    LOG(DEBUG) << "Storing pixel (col, row) = (" << col << ", " << row << ")" << std::endl;
                }
                pixels_.at(t).push_back(pixel);
                hHitMap.at(names_.at(t))->Fill(pixel.get()->column(), pixel.get()->row());
                hPixelToT.at(names_.at(t))->Fill(pixel.get()->raw());
                // hPixelToT.at(names_.at(t))->Fill(pixel.get()->raw());
                // display the 10 bit timestamp distribution
                hTimeStamp.at(names_.at(t))->Fill(fmod((pixel.get()->timestamp() / 8.), pow(2, 10)));
                pixelbuffers_.at(t).pop(); // remove top element
            } else {
                LOG(DEBUG) << "Pixel position with respect to event: " << corryvreckan::to_string(position);
                break;
            }
        }
    }
    prev_event_end_ = clipboard->getEvent()->end();

    // Return value telling analysis to keep running
    for(auto t : tags_) {
        if(pixels_.at(t).size() > 0)
            return StatusCode::Success;
    }
    bool data_in_buffer = false;
    for(auto t : tags_) {
        if(pixelbuffers_.at(t).size() > 0) {
            data_in_buffer = true;
        }
    }

    if(!data_in_buffer) {
        return StatusCode::EndRun;
    }

    return StatusCode::NoData;
}

std::shared_ptr<Pixel> EventLoaderMuPixTelescope::read_hit(const RawHit& h, uint tag, long unsigned int corrected_fpgaTime) {

    uint16_t time = 0x0;
    // TS can be sampled on both edges - keep this optional
    auto detector = detectors_.back();
    if((h.get_ts2() == uint16_t(-1)) || (!use_both_timestamps_)) {
        (detector->getType() == "telepix2" ? time = (timestampMask_ & h.timestamp_raw())
                                           : time = ((timestampMask_ & h.timestamp_raw()) << 1));

    } else if(h.timestamp_raw() > h.get_ts2()) {
        time = ((timestampMask_ & h.timestamp_raw()) << 1);
    } else {
        time = ((timestampMask_ & h.timestamp_raw()) << 1) + 0x1;
    }

    // in case the ts clock is divided
    time *= (ckdivend_ + 1);

    double time_shifted = static_cast<double>(time) * static_cast<double>(ckdivend_ + 1);

    auto name = detector->getName();

    ts1_ts2[name]->Fill(h.get_ts2(), h.timestamp_raw());
    double px_timestamp =
        clockToTime_ * (static_cast<double>((corrected_fpgaTime >> 1) & (0xFFFFFFFFFFFFF - timestampMask_)) + time_shifted) -
        static_cast<double>(timeOffset_.at(tag));

    hts_ToT[name]->Fill(static_cast<double>(h.tot_decoded()));

    // store the ToT information if reasonable
    double tot_timestamp = clockToTime_ * (static_cast<double>(h.tot_decoded()) * multiplierToT_);

    ts_TS1_ToT[name]->Fill(static_cast<double>((static_cast<uint>(px_timestamp / 8)) & timestampMaskExtended_),
                           (static_cast<double>(static_cast<uint>(tot_timestamp / 8) & timestampMaskExtended_)));

    double tot = tot_timestamp - (time_shifted * clockToTime_);

    // catch lapse of ToT time stamp
    while(tot < 0)
        tot += maxToT_;

    return std::make_shared<Pixel>(names_.at(tag), h.column(), h.row(), tot, tot, px_timestamp);
}

void EventLoaderMuPixTelescope::fillBuffer() {
    long unsigned int corrected_fpgaTime = 0;
    unsigned int raw_time = 0;
    // here we need to check quite a number of cases
    bool buffers_full = false;
    // do not fill the buffer if it is already full anyways
    for(auto t : tags_) {
        if(pixelbuffers_.at(t).size() < buffer_depth_)
            buffers_full = false;
    }
    if(buffers_full) {
        LOG(DEBUG) << "Buffer still fully filled";
        return;
    }
    while(!buffers_full) {
        if(blockFile_->read_next(tf_)) {
            //      std::cout << "Reading hit: "<< tf_.timestamp() <<std::endl;
            // no hits in data - can only happen if the zero suppression is switched off, skip the event
            if(tf_.num_hits() == 0) {
                continue;
            }
            // need to determine the sensor layer that is identified by the tag
            RawHit h = tf_.get_hit(0);
            // tag does not match - continue reading if data is not sorted
            uint tag = h.tag() & (~0x3);
            if(names_.count(tag) != 1) {
                throw RuntimeError("Unknown pixel tag read in data: " + to_string(tag));
            }
            // all hits in one frame are from the same sensor. Copy them
            for(uint i = 0; i < tf_.num_hits(); ++i) {
                h = tf_.get_hit(i, types_.at(tag));
                LOG(TRACE) << "Filling buffer with " << h;

                // this assumes a few things:
                // time from fpga is using a 500MHz clock (4 times the clock used for the hit timestamp
                corrected_fpgaTime = (tf_.timestamp() >> 2);
                // just take 10 bits from the hit timestamp
                raw_time = h.timestamp_raw() & 0x3FF;
                // get the fpga time +1bit just for plots
                raw_fpga_vs_chip.at(names_.at(tag))->Fill(raw_time, static_cast<double>(corrected_fpgaTime & 0x7FF));
                chip_delay.at(names_.at(tag))->Fill(static_cast<double>((corrected_fpgaTime & 0x3FF) - raw_time));
                // if the chip timestamp is smaller than the fpga we have a bit flip on the 11th bit
                if(((corrected_fpgaTime & 0x3FF) < raw_time)) { // && (corrected_fpgaTime>1024)) {
                    corrected_fpgaTime -= 1024;
                }
                raw_fpga_vs_chip_corrected.at(names_.at(tag))
                    ->Fill(raw_time, static_cast<double>(corrected_fpgaTime & 0x7FF));

                pixelbuffers_.at(tag).push(read_hit(h, tag, (corrected_fpgaTime * 4)));
            }
            buffers_full = true;
            for(auto t : tags_) {
                if(pixelbuffers_.at(t).size() < buffer_depth_)
                    buffers_full = false;
            }
        } else {
            LOG(INFO) << "Reached eof";
            for(auto t : tags_)
                LOG(INFO) << "Hits buffered in " << names_.at(t) << ": " << pixelbuffers_.at(t).size();
            eof_ = true;
            break;
        }
    }
}
std::map<std::string, int> EventLoaderMuPixTelescope::typeString_to_typeID = {{"mupix8", MP8_SORTED_TS2},
                                                                              {"mupix9", MP10_SORTED_TS2},
                                                                              {"mupix10", MP10_UNSORTED_GS1_GS2},
                                                                              {"mupix11", MP11_UNSORTED_GS1_GS2},
                                                                              {"atlaspix3", AP3_UNSORTED_GS1_GS2},
                                                                              {"run2020v1", R20V1_UNSORTED_GS1_GS2_GS3},
                                                                              {"run2020v2", R20V2_UNSORTED_GS1_GS2_GS3},
                                                                              {"run2020v3", R20V3_UNSORTED_GS1_GS2_GS3},
                                                                              {"run2020v4", R20V4_UNSORTED_GS1_GS2_GS3},
                                                                              {"run2020v5", R20V5_UNSORTED_GS1_GS2_GS3},
                                                                              {"run2020v6", R20V6_UNSORTED_GS1_GS2_GS3},
                                                                              {"run2020v7", R20V7_UNSORTED_GS1_GS2_GS3},
                                                                              {"run2020v8", R20V8_UNSORTED_GS1_GS2_GS3},
                                                                              {"telepix2", TELEPIX2_UNSORTED_GS1_GS2_GS3}};

StatusCode EventLoaderMuPixTelescope::run(const std::shared_ptr<Clipboard>& clipboard) {
    eventNo_++;

    for(auto t : tags_)
        pixels_.at(t).clear();
    // get the hits
    StatusCode result = (isSorted_ ? read_sorted(clipboard) : read_unsorted(clipboard));
    LOG(TRACE) << "StatusCode: " << static_cast<int>(result);
    for(auto t : tags_) {

        hHitsEvent.at(names_.at(t))->Fill(double(pixels_.at(t).size()));
        counterHits_.at(t) += pixels_.at(t).size();
        if(eventNo_ % 1000 == 0) {
            int point = eventNo_ / 1000;
            hitsPerkEvent.at(names_.at(t))->Fill(point, double(counterHits_.at(t)));
            counterHits_.at(t) = 0;
        }

        if(pixels_.at(t).size() > 0)
            clipboard->putData(pixels_.at(t), names_.at(t));
        stored_.at(t) += pixels_.at(t).size();
    }
    return result;
}
