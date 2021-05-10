/**
 * @file
 * @brief Implementation of module EventLoaderMuPixTelescope
 *
 * @copyright Copyright (c) 2019-2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "EventLoaderMuPixTelescope.h"
#include <string>
#include "dirent.h"
#include "objects/Cluster.hpp"
#include "objects/Pixel.hpp"
#include "objects/Track.hpp"
using namespace corryvreckan;

EventLoaderMuPixTelescope::EventLoaderMuPixTelescope(Configuration& config, std::vector<std::shared_ptr<Detector> > detectors)
    : Module(config, detectors), blockFile_(nullptr) {

    config_.setDefault<bool>("is_sorted", false);
    config_.setDefault<bool>("ts2_is_gray", false);
    config_.setDefault<unsigned>("buffer_depth", 1000);
    config_.setDefault<double>("time_offset", 0.0);
    inputDirectory_ = config_.getPath("input_directory");
    buffer_depth_ = config.get<unsigned>("buffer_depth");
    isSorted_ = config_.get<bool>("is_sorted");
    if(config.count({"run", "input_file"}) > 1) {
        throw InvalidCombinationError(config, {"run", "input_file"}, "run and input_file are mutually exclusive.");
    } else if(config_.has("input_file")) {
        input_file_ = config_.get<string>("input_file");
    } else {
        runNumber_ = config_.get<int>("run");
    }
    // find the corresponding detetctors:
    for(auto d : detectors){
        if(typeString_to_typeID.find(d->getType()) != typeString_to_typeID.end()) {
            detectors_.push_back(d);
        }
    }
}

void EventLoaderMuPixTelescope::initialize() {
    // extract the tag from each detetcor name
    for(auto detector : detectors_)
    {
        string tag = detector->getName();
        if(tag.find("_") < tag.length())
            tag = tag.substr(tag.find("_") + 1);
        uint tag_ = uint(stoi(tag, nullptr, 16));
        LOG(DEBUG) << detector->getName() << " is using the fpga link tag " << hex << tag_;
        if(typeString_to_typeID.find(detector->getType()) == typeString_to_typeID.end()) {
            throw KeyValueParseError("tag " + std::to_string(tag_), "Sensor tag not supported");
        }
        removed_[tag] = 0;
        stored_[tag] = 0;
        // take the offset from the geometry file
        timeOffset_[tag_] = detector->timeOffset();
        types_[tag_] = typeString_to_typeID.at(detector->getType());
        LOG(INFO) << "Detector " << detector->getType() << "is assigned to type id " << types_.at(tag_);
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
        LOG(INFO) << "File found" << endl;
    string file = (inputDirectory_ + "/" + entry->d_name);
    LOG(INFO) << "reading " << file;
    blockFile_ = new BlockFile(file);
    if(!blockFile_->open_read()) {
        throw MissingDataError("Cannot read data file: " + input_file_);
    }

    // create the histograms for all sensor
    for(auto & detector : detectors_){
        auto name = detector->getName();
        TDirectory* directory = getROOTDirectory();
        TDirectory* local_directory = directory->mkdir(name.c_str());

        if(local_directory == nullptr) {
            throw RuntimeError("Cannot create or access local ROOT directory for module " + this->getUniqueName());
        }
        local_directory->cd();
        std::string title = name+"_hitMap; column; row";
        hHitMap[name] = new TH2F("hitMap",
                           title.c_str(),
                           detector->nPixels().x(),
                           -.05,
                           detector->nPixels().x() - .5,
                           detector->nPixels().y(),
                           -.05,
                           detector->nPixels().y() - .5);
        title = name+"hitMap of out of event hits; column; row";
        hdiscardedHitmap[name] = new TH2F("discardedhitMap",
                                    title.c_str(),
                                    detector->nPixels().x(),
                                    -.05,
                                    detector->nPixels().x() - .5,
                                    detector->nPixels().y(),
                                    -.05,
                                    detector->nPixels().y() - .5);
        title = name+"pixelToT; ToT in TS2 clock cycles.;";
        hPixelToT[name] = new TH1F("pixelToT",title.c_str(), 64, -0.5, 63.5);
        title = name+"pixelTS; TS in clock cycles; ";
        hTimeStamp[name] = new TH1F("pixelTS", title.c_str(), 1024, -0.5, 1023.5);
        title = name+"hHitsEvent; # hits per event; ";
        hHitsEvent[name] = new TH1F("hHitsEvent", title.c_str(), 300, -.5, 299.5);
        title = name+"hitsper1kevents; corry events /1k; hits per 1k events";
        hitsPerkEvent[name] = new TH1F("hHitsPerkEvent", title.c_str(), 1000, -.5, 999.5);
        title = name+ "fpga vs chip clock;chip clock;fpga clock";
        raw_fpga_vs_chip[name] =
                new TH2F("raw_fpga_vs_chip", title.c_str(), 1024, 0, 1023, 2048, 0, 2047);
        title = name+"fpga vs chip clock;chip clock;fpga clock";
        raw_fpga_vs_chip_corrected[name] =
                new TH2F("raw_fpga_vs_chip_corrected", title.c_str(), 1024, 0, 1023, 2048, 0, 2047);
        title = name+"Delay of chip events wrt. telescope frame;fpga clock@ chip clock 0;#events";
        chip_delay[name] = new TH1F(
                    "chip_delay", title.c_str(), 2048, -1023, 1023);
    }
}

void EventLoaderMuPixTelescope::finalize(const std::shared_ptr<ReadonlyClipboard>&) {

    for(auto d : detectors_) LOG(INFO) << "Number of hits put to clipboard: " << stored_.at(d->getName())
              << " and number of removed (not fitting in an event) hits: " << removed_.at(d->getName());
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
        auto tag = h.tag() &uint(~0x3);
        h = tf_.get_hit(i,tag);
        // Convert time stamp to ns - i'd like to do this already on the mupix8_DAQ side, but have not found the time yet,
        // assuming 10bit ts
        double px_timestamp = 8 * static_cast<double>(((tf_.timestamp() >> 2) & 0xFFFFFFFFFFC00) + h.timestamp_raw());
        // setting tot and charge to zero here - needs to be improved
        pixels_[tag].push_back(std::make_shared<Pixel>(names_.at(tag), h.column(), h.row(), 0, 0, px_timestamp));
    }
    // If no event is defined create one
    if(clipboard->getEvent() == nullptr) {
        // The readout FPGA creates frames with a length of 128 time stamps, each 8ns. The int division cuts of lowest bits
        int begin = int(pixels_.front()->timestamp()) / 1024;
        clipboard->putEvent(std::make_shared<Event>(double(begin * 1024), double((begin + 1) * 1024)));
    }
    return StatusCode::Success;
}

StatusCode EventLoaderMuPixTelescope::read_unsorted(const std::shared_ptr<Clipboard>& clipboard) {
    if(!eof_)
        fillBuffer();
    else
        return StatusCode::EndRun;
    while(true) {
        if(pixelbuffer_.size() == 0)
            break;
        auto pixel = pixelbuffer_.top();
        if((pixel->timestamp() < clipboard->getEvent()->start())) {
            LOG(DEBUG) << " Old hit found: " << Units::display(pixel->timestamp(), "us") << " vs prev end (" << eventNo_ - 1
                       << ")\t" << Units::display(prev_event_end_, "us") << " and current start \t"
                       << Units::display(clipboard->getEvent()->start(), "us")
                       << " and duration: " << clipboard->getEvent()->duration()
                       << "and number of triggers: " << clipboard->getEvent()->triggerList().size();
            removed_++;
            hdiscardedHitmap->Fill(pixel->column(), pixel->row());
            pixelbuffer_.pop(); // remove top element
            continue;
        }
        if(pixelbuffer_.size() && (pixel->timestamp() < clipboard->getEvent()->end()) &&
                (pixel->timestamp() > clipboard->getEvent()->start())) {
            LOG(DEBUG) << " Adding pixel hit: " << Units::display(pixel->timestamp(), "us") << " vs prev end ("
                       << eventNo_ - 1 << ")\t" << Units::display(prev_event_end_, "us") << " and current start \t"
                       << Units::display(clipboard->getEvent()->start(), "us")
                       << " and duration: " << Units::display(clipboard->getEvent()->duration(), "us");
            pixels_.push_back(pixel);
            hHitMap->Fill(pixel.get()->column(), pixel.get()->row());
            hPixelToT->Fill(pixel.get()->raw());
            // display the 10 bit timestamp distribution
            hTimeStamp->Fill(fmod((pixel.get()->timestamp() / 8.), pow(2, 10)));
            pixelbuffer_.pop();
        } else {
            break;
        }
    }
    if(pixelbuffer_.size() < buffer_depth_)
        fillBuffer();
    // Return value telling analysis to keep running
    if(pixelbuffer_.size() == 0)
        return StatusCode::NoData;
    prev_event_end_ = clipboard->getEvent()->end();
    return StatusCode::Success;
}

void EventLoaderMuPixTelescope::fillBuffer() {
    long unsigned int temp_fpga_time = 0;
    unsigned int raw_time = 0;
    unsigned int overlap_fpga = 0;
    // here we need to check quite a number of cases
    while(pixelbuffer_.size() < buffer_depth_) {
        if(blockFile_->read_next(tf_)) {
            // no hits in data - can only happen if the zero suppression is switched off, skip the event
            if(tf_.num_hits() == 0) {
                continue;
            }
            // need to determine the sensor layer that is identified by the tag
            RawHit h = tf_.get_hit(0);
            // tag does not match - continue reading if data is not sorted
            if(((h.tag() & uint(~0x3)) != tag_)) {
                continue;
            }
            // all hits in one frame are from the same sensor. Copy them
            for(uint i = 0; i < tf_.num_hits(); ++i) {
                h = tf_.get_hit(i, type_);

                // this assumes a few things:
                // time from fpga is using a 500MHz clock (4 times the clock used for the hit timestamp
                temp_fpga_time = (tf_.timestamp() >> 2);
                // just take 10 bits from the hit timestamp
                raw_time = h.timestamp_raw() & 0x3FF;
                // get the fpga time +1bit just for plots
                raw_fpga_vs_chip->Fill(raw_time, static_cast<double>(temp_fpga_time & 0x7FF));
                chip_delay->Fill(static_cast<double>((temp_fpga_time & 0x3FF) - raw_time));
                // if the chip timestamp is smaller than the fpga we have a bit flip on the 11th bit
                if((temp_fpga_time & 0x3FF) < raw_time) {
                    temp_fpga_time -= 1024;
                }
                raw_fpga_vs_chip_corrected->Fill(raw_time, static_cast<double>(temp_fpga_time & 0x7FF));

                // convert timestamp to ns - i'd like to do this already on the mupix8_DAQ side, but have not found the time
                // yet, assuming 10bit ts
                double px_timestamp = 8 * static_cast<double>((temp_fpga_time & 0xFFFFFFFFFFC00) + raw_time) - timeOffset_;
                LOG(TRACE) << "Pixel timestamp " << px_timestamp;
                // setting tot and charge to zero here - needs to be improved
                pixelbuffer_.push(std::make_shared<Pixel>(detector_->getName(), h.column(), h.row(), 0, 0, px_timestamp));
            }
        } else {
            eof_ = true;
            break;
        }
    }
}
std::map<std::string, int> EventLoaderMuPixTelescope::typeString_to_typeID = {{"mupix8", MP8_SORTED_TS2},
                                                                              {"mupix9", MP10_SORTED_TS2},
                                                                              {"mupix10", MP10_UNSORTED_GS1_GS2},
                                                                              {"run2020v1", R20V1_UNSORTED_GS1_GS2_GS3},
                                                                              {"run2020v2", R20V2_UNSORTED_GS1_GS2_GS3},
                                                                              {"run2020v3", R20V3_UNSORTED_GS1_GS2_GS3},
                                                                              {"run2020v4", R20V4_UNSORTED_GS1_GS2_GS3},
                                                                              {"run2020v5", R20V5_UNSORTED_GS1_GS2_GS3},
                                                                              {"run2020v6", R20V6_UNSORTED_GS1_GS2_GS3},
                                                                              {"run2020v7", R20V7_UNSORTED_GS1_GS2_GS3},
                                                                              {"run2020v8", R20V8_UNSORTED_GS1_GS2_GS3},
                                                                              {"run2020v9", R20V9_UNSORTED_GS1_GS2_GS3}};

StatusCode EventLoaderMuPixTelescope::run(const std::shared_ptr<Clipboard>& clipboard) {
    eventNo_++;
    pixels_.clear();
    // get the hits
    StatusCode result = (isSorted_ ? read_sorted(clipboard) : read_unsorted(clipboard));
    hHitsEvent->Fill(double(pixels_.size()));
    counterHits_ += pixels_.size();
    if(eventNo_ % 1000 == 0) {
        int point = eventNo_ / 1000;
        hitsPerkEvent->Fill(point, double(counterHits_));
        counterHits_ = 0;
    }
    if(pixels_.size() > 0)
        clipboard->putData(pixels_, detector_->getName());
    stored_ += pixels_.size();
    return result;
}
