/**
 * @file
 * @brief Implementation of module EventLoaderALiBaVa
 *
 * @copyright Copyright (c) 2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <dirent.h>

#include "objects/Pixel.hpp"
#include "ALiBaVa/auxfunctions.h"
#include "EventLoaderALiBaVa.h"
#include "ALiBaVa/HDFRoot.h"

using namespace corryvreckan;

EventLoaderALiBaVa::EventLoaderALiBaVa(Configuration& config, std::shared_ptr<Detector> detector)
    : Module(config, detector) {
      m_detector = detector;
    }

void EventLoaderALiBaVa::initialize() {

  // Take input directory from global parameters
  std::string inputDirectory = config_.getPath("input_directory");

  // LOG(DEBUG) << "input directory read";

  // Open the root directory
  DIR* directory = opendir(inputDirectory.c_str());
  if(directory == nullptr) {
      LOG(ERROR) << "Directory " << inputDirectory << " does not exist";
      return;
  }
  dirent* entry;

  // LOG(DEBUG) << "input directory opened";

  // Read the files (data, pedestal and calibration) in the folder
  while(entry = readdir(directory)) {
      if(entry->d_type == DT_REG) {
          std::string entryName = entry->d_name;
          if(entryName.find("dat") != std::string::npos){
            m_datafilename = inputDirectory + "/" + entryName;
          }
          if(entryName.find("ped") != std::string::npos){
            m_pedestalfilename = inputDirectory + "/" + entryName;
          }
          if(entryName.find("cal") != std::string::npos){
            m_calibrationfilename = inputDirectory + "/" + entryName;
          }
      }
  }

  // LOG(DEBUG) << "input files read";

  if(m_datafilename.length() == 0) {
      LOG(ERROR) << "No data file was found for ALiBaVa in " << inputDirectory;
      return;
  }

  if(m_pedestalfilename.length() == 0) {
      LOG(ERROR) << "No pedestal file was found." << "\n" << "Cannot continue without pedestal information.";
      return;
  }

  if(m_calibrationfilename.length() == 0) {
      LOG(WARNING) << "No calibration file was found." << "\n" << "Results will be uncalibrated.";
  }

  // LOG(DEBUG) << "passed file checks";

  // Open the ALiBaVa data
  ALiBaVaPointer = DataFileRoot::OpenFile(m_datafilename.c_str(), m_pedestalfilename.c_str(), m_calibrationfilename.c_str());
  // LOG(DEBUG) << "pointer created to alibava file";
  ALiBaVa_loader(ALiBaVaPointer, m_datafilename.c_str(), m_pedestalfilename.c_str(), m_calibrationfilename.c_str());
  // LOG(DEBUG) << "Initializer finished";
  nEvents = ALiBaVaPointer->nevents();
  ALiBaVaPointer->rewind();

  while(!ALiBaVaPointer->read_event()){
    ALiBaVaPointer->process_event();
    double currentAlibavaTimeStamp = ALiBaVaPointer->time();
    if(currentAlibavaTimeStamp > lastAlibavaTimeStamp){
      lastAlibavaTimeStamp = currentAlibavaTimeStamp;
    }
  }
  // LOG(STATUS) << lastAlibavaTimeStamp;

  ALiBaVaPointer->rewind();
  LOG(STATUS) << "It's not stuck, it's just really slow...";
}

StatusCode EventLoaderALiBaVa::run(const std::shared_ptr<Clipboard>& clipboard) {
  PixelVector pixels;
  auto event = clipboard->getEvent();
  // if(event->start() == 6) return StatusCode::EndRun;
  // LOG(STATUS) << "Corry EVENT START = " << event->start()/billion << "s";
  // LOG(STATUS) << "Corry EVENT END = "<< event->end()/billion<< "s";

  if(event->getTimestampPosition(lastAlibavaTimeStamp*billion) == Event::Position::BEFORE)
  {
    return StatusCode::EndRun;
  }

  while(!ALiBaVaPointer->read_event()){
      ALiBaVaPointer->process_event();
      // LOG(STATUS) << ALiBaVaPointer->time();
      // LOG(STATUS) << "----------";
      // LOG(STATUS) << ALiBaVaPointer->time();
      // LOG(STATUS) << ALiBaVaPointer->time()*billion;
      // LOG(STATUS) << ALiBaVaPointer->time()*billion/billion;
      // LOG(STATUS) << "----------";
      double AlibavaEventTime = ALiBaVaPointer->time()*billion;
      auto position = event->getTimestampPosition(AlibavaEventTime);
      // debug
      // LOG(STATUS) << "EVENT START " << event->start();
      // LOG(STATUS) << "EVENT END "<< event->end();
      // if(position == Event::Position::BEFORE) LOG(STATUS) << "BEFORE";
      // if(position == Event::Position::DURING) LOG(STATUS) << "DURING";
      // if(position == Event::Position::AFTER) LOG(STATUS) << "AFTER";
      // if(position == Event::Position::UNKNOWN) LOG(STATUS) << "UNKNOWN";
      // if(position != Event::Position::UNKNOWN && position != Event::Position::AFTER && position != Event::Position::DURING && position != Event::Position::BEFORE) LOG(STATUS) << "Wtf";

      if(position == Event::Position::DURING){
        iEvent++;
        // LOG(STATUS) << "Alibava event in Corry event: ALIBAVATIME = " << AlibavaEventTime/billion << "s";
        // LOG(STATUS) << "Alibava event in corry event";
        for(int chan = 0; chan < ALiBaVaPointer->nchan(); chan++){
          double calibration = ALiBaVaPointer->get_gain(chan);
          double CalSignal = ALiBaVaPointer->signal(chan);
          double ADCSignal = ALiBaVaPointer->signal(chan)*calibration;
          auto pixel = std::make_shared<Pixel>(m_detector->getName(), chan, 0, ADCSignal, CalSignal, AlibavaEventTime);
          pixels.push_back(pixel);
        }
      }
      // if(position == Event::Position::AFTER){
      //   LOG(STATUS) << "Alibava event after Corry event: ALIBAVATIME = " << AlibavaEventTime/1000000000 << "s";
      //   LOG(STATUS) << "Reverting to last Alibava state.";
      //   // LOG(STATUS) << "Alibava event after corry event";
      //   ALiBaVaPointer->restore();
      //   clipboard->putData(pixels, m_detector->getName());
      //   if(pixels.empty()){
      //     return StatusCode::NoData;
      //   }
      //   return StatusCode::Success;
      // }

      // ALiBaVaPointer->save();

  }
  ALiBaVaPointer->rewind();
  clipboard->putData(pixels, m_detector->getName());
  if(pixels.empty()){
    return StatusCode::NoData;
  }
  return StatusCode::Success;
}


void EventLoaderALiBaVa::finalize(const std::shared_ptr<ReadonlyClipboard>& clipboard) {
    ALiBaVaPointer->close();
    delete ALiBaVaPointer;
    LOG(STATUS) << "Analysed " << iEvent << " of total " << nEvents << " ALiBaVa events";
}
