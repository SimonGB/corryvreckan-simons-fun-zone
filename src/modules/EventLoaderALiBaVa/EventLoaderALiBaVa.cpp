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

using namespace corryvreckan;

EventLoaderALiBaVa::EventLoaderALiBaVa(Configuration& config, std::shared_ptr<Detector> detector)
    : Module(config, move(detector)) {
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

  LOG(DEBUG) << "input files read";

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

  LOG(DEBUG) << "passed file checks";

  // Open the ALiBaVa data
  //########### TEST AREA ############
  LOG(DEBUG) << "####################################################";
  // DataFileRoot * test = DataFileRoot::OpenFile(m_datafilename.c_str(),m_pedestalfilename.c_str(),m_calibrationfilename.c_str());
  LOG(DEBUG) << m_datafilename.c_str();
  LOG(DEBUG) << m_calibrationfilename.c_str();
  LOG(DEBUG) << m_pedestalfilename.c_str();

  DataFileRoot * test = DataFileRoot::OpenFile(m_datafilename.c_str());
  // test->open(m_pedestalfilename.c_str());
  LOG(DEBUG) << test->valid();
  LOG(DEBUG) << test->type();
  LOG(DEBUG) << "####################################################";
  //##################################
  // ALiBaVaPointer = DataFileRoot::OpenFile(m_datafilename.c_str(), m_pedestalfilename.c_str(), m_calibrationfilename.c_str());
  // LOG(DEBUG) << "pointer created to alibava file";
  // ALiBaVa_loader(ALiBaVaPointer, m_datafilename.c_str(), m_pedestalfilename.c_str(), m_calibrationfilename.c_str());

  LOG(DEBUG) << "Initializer finished";
}

StatusCode EventLoaderALiBaVa::run(const std::shared_ptr<Clipboard>& clipboard) {
  LOG(DEBUG) << "Starting run";
  while(!ALiBaVaPointer->read_event()){
      LOG(DEBUG) << "reading event";
      nEvents++;
      LOG(DEBUG) << "nEvents incremented";
      PixelVector pixels;
      LOG(DEBUG) << "Pixelvector created";
      ALiBaVaPointer->process_event();
      LOG(DEBUG) << "event processed";
      double eventTime = ALiBaVaPointer->time();
      LOG(DEBUG) << "time retrieved";

      for(int chan = 0; chan < ALiBaVaPointer->nchan(); chan++){
        LOG(DEBUG) << "Starting channel " << chan;
        double calibration = ALiBaVaPointer->get_gain(chan);
        LOG(DEBUG) << "calibration defined";
        double CalSignal = ALiBaVaPointer->signal(chan);
        LOG(DEBUG) << "signal defined";
        double ADCSignal = ALiBaVaPointer->signal(chan)*calibration;
        LOG(DEBUG) << "raw signal defined";
        LOG(DEBUG) << chan << 0 << ADCSignal << CalSignal << eventTime;
        LOG(DEBUG) << m_detector->getName();
        auto pixel = std::make_shared<Pixel>(m_detector->getName(), chan, 0, ADCSignal, CalSignal, eventTime);
        LOG(DEBUG) << "pixel created";
        pixels.push_back(pixel);
        LOG(DEBUG) << "pixel pushed back into vector";
      }
      clipboard->putData(pixels, m_detector->getName());
  }
  return StatusCode::Success;
}


void EventLoaderALiBaVa::finalize(const std::shared_ptr<ReadonlyClipboard>& clipboard) {
    ALiBaVaPointer->close();
    delete ALiBaVaPointer;
    LOG(DEBUG) << "Analysed " << nEvents << " events";
}
