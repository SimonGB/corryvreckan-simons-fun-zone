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
}

StatusCode EventLoaderALiBaVa::run(const std::shared_ptr<Clipboard>& clipboard) {
  // LOG(STATUS) << "Starting run";
  // std::string detectorName = "TEST";
  PixelVector pixels;
  while(!ALiBaVaPointer->read_event()){
      // LOG(DEBUG) << "reading event";
      iEvent++;
      // LOG(DEBUG) << "iEvent incremented";
      // LOG(DEBUG) << "Pixelvector created";
      ALiBaVaPointer->process_event();
      // LOG(DEBUG) << "event processed";
      double eventTime = ALiBaVaPointer->time();
      // LOG(DEBUG) << "time retrieved";

      for(int chan = 0; chan < ALiBaVaPointer->nchan(); chan++){
      // for(int chan = 0; chan < 5; chan++){

        // LOG(DEBUG) << "Starting channel " << chan;
        double calibration = ALiBaVaPointer->get_gain(chan);
        // LOG(DEBUG) << calibration;
        // LOG(DEBUG) << "calibration defined";
        double CalSignal = ALiBaVaPointer->signal(chan);
        // LOG(DEBUG) << "signal defined";
        double ADCSignal = ALiBaVaPointer->signal(chan)*calibration;
        // LOG(DEBUG) << "raw signal defined";
        // LOG(DEBUG) << chan;
        // LOG(DEBUG) << ADCSignal;
        // LOG(DEBUG) << CalSignal;
        // LOG(DEBUG) << eventTime;
        // LOG(DEBUG) << detectorName;
        auto pixel = std::make_shared<Pixel>(m_detector->getName(), chan, 0, ADCSignal, CalSignal, eventTime);
        // LOG(DEBUG) << "pixel created";
        pixels.push_back(pixel);
        // LOG(DEBUG) << "pixel pushed back into vector";
        // LOG(STATUS) << "Channel number " << chan;
      }
      // LOG(STATUS) << "Event number " << iEvent;
      break;
      // if(iEvent > 10)
      // {
      //   // LOG(STATUS) << "###########################################################################################################";
      //   break;
      // }
  }


  clipboard->putData(pixels, m_detector->getName());
  if(iEvent >= nEvents){
    return StatusCode::EndRun;
  }
  if(pixels.empty()){
    return StatusCode::NoData;
  }
  // LOG(STATUS) << "Ending run";
  return StatusCode::Success;
}


void EventLoaderALiBaVa::finalize(const std::shared_ptr<ReadonlyClipboard>& clipboard) {
    ALiBaVaPointer->close();
    delete ALiBaVaPointer;
    LOG(STATUS) << "Analysed " << iEvent << " of total " << nEvents << " events";
}
