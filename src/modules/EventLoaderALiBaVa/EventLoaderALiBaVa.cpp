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

  LOG(DEBUG) << "input directory read";

  // Open the root directory
  DIR* directory = opendir(inputDirectory.c_str());
  if(directory == nullptr) {
      LOG(ERROR) << "Directory " << inputDirectory << " does not exist";
      return;
  }
  dirent* entry;

  LOG(DEBUG) << "input directory opened";

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
      LOG(WARNING) << "No pedestal file was found." << "\n" << "Datafile will be used for pedestal";
  }

  if(m_calibrationfilename.length() == 0) {
      LOG(WARNING) << "No calibration file was found." << "\n" << "Results will be uncalibrated.";
  }

  LOG(DEBUG) << "passed file checks";

  // Open the ALiBaVa data
  // ALiBaVaPointer = DataFileRoot::OpenFile(m_datafilename.c_str(), m_pedestalfilename.c_str(), m_calibrationfilename.c_str());
  ALiBaVaPointer = DataFileRoot::OpenFile(m_datafilename.c_str());
  // ALiBaVaPointer = DataFileRoot::OpenFile(m_pedestalfilename.c_str());
  LOG(DEBUG) << "pointer created to alibava file";
  ALiBaVa_loader(ALiBaVaPointer, m_datafilename.c_str(), m_pedestalfilename.c_str(), m_calibrationfilename.c_str());
  LOG(DEBUG) << "Initializer finished";
  nEvents = ALiBaVaPointer->nevents();
  ALiBaVaPointer->rewind();

  // while(!ALiBaVaPointer->read_event()){
  //   ALiBaVaPointer->process_event(kFALSE);
  //   double currentAlibavaTimeStamp = ALiBaVaPointer->clock_counter();
  //   if(currentAlibavaTimeStamp > lastAlibavaTimeStamp){
  //     lastAlibavaTimeStamp = currentAlibavaTimeStamp;
  //   }
  // }
  // LOG(STATUS) << lastAlibavaTimeStamp;

  // ALiBaVaPointer->rewind();
  // LOG(STATUS) << "It's not stuck, it's just really slow...";
  // LOG(DEBUG) << "The reading out of ALiBaVa events doesn't happen in-order, so I'm forced to iterate over ALL Alibava events for every single Corryvreckan event...";
}

StatusCode EventLoaderALiBaVa::run(const std::shared_ptr<Clipboard>& clipboard) {
  PixelVector pixels;
  auto event = clipboard->getEvent();
  // if(event->start() == 6) return StatusCode::EndRun;
  // LOG(STATUS) << "Corry EVENT START = " << event->start()/billion << "s";
  // LOG(STATUS) << "Corry EVENT END = "<< event->end()/billion<< "s";

  if(event->getTimestampPosition(lastAlibavaTimeStamp) == Event::Position::BEFORE)
  {
    return StatusCode::EndRun;
  }

  // LOG(DEBUG) << "----------DEBUG----------";
  bool test=false;
  if(test){
    int i = 0;
    while(i < 4){
      ALiBaVaPointer->read_event();
      ALiBaVaPointer->process_event(kFALSE);
      LOG(DEBUG) << "Event: " << i;
      LOG(DEBUG) << "nr. of channels: " << 128;
      int activechannels = 0;
      for(int chan = 0; chan < 128; chan++){
        double calibration = ALiBaVaPointer->get_gain(chan);
        if(calibration != 0){
          double CalSignal = ALiBaVaPointer->signal(chan);
          double ADCSignal = ALiBaVaPointer->signal(chan)*calibration;
          double noise = ALiBaVaPointer->noise(chan);
          double pedestal = ALiBaVaPointer->ped(chan);
          // double AlibavaEventTime = ALiBaVaPointer->time()*billion;
          double AlibavaEventTime = ALiBaVaPointer->clock_counter();
          // double AlibavaEventTime = ALiBaVaPointer->time()*billion;
          double temp = ALiBaVaPointer->temp();
          double data = ALiBaVaPointer->data(chan);
          LOG(DEBUG) << "\t Channel: " << chan << " \t Calibration: " << calibration << " \t Pedestal: " << pedestal << " \t Signal: " << CalSignal << " \t ADC: " << ADCSignal << " \t Timestamp: " << AlibavaEventTime << "\t Temp: " << temp << "\t Raw data: " << data;
          activechannels++;
        }
      }
      LOG(DEBUG) << "Nr. of active channels: " << activechannels;
      i++;
    }
    return StatusCode::EndRun;
  }
  // LOG(DEBUG) << "----------DEBUG----------";

  while(!ALiBaVaPointer->read_event()){
      ALiBaVaPointer->process_event(kFALSE);
      // LOG(STATUS) << ALiBaVaPointer->time();
      double AlibavaEventTime = ALiBaVaPointer->clock_counter();
      double TDCTime = ALiBaVaPointer->time();
      if(!ALiBaVaPointer->valid_time(TDCTime)) continue;
      // double AlibavaEventTime = ALiBaVaPointer->time()*billion;
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
    LOG(DEBUG) << "Analysed " << iEvent << " of total " << nEvents << " ALiBaVa events";
}
