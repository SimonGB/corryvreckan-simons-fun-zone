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

#include "objects/Pixel.hpp"
#include "ALiBaVa/auxfunctions.h"
#include "EventLoaderALiBaVa.h"

using namespace corryvreckan;
using namespace DataFileRoot;
using namespace std;

EventLoaderALiBaVa::EventLoaderALiBaVa(Configuration& config, vector<shared_ptr<Detector>> detectors)
    : Module(config, move(detector)) {
      m_detector = detector;
    }

void EventLoaderALiBaVa::initialize() {

  string m_datafilename;
  string m_pedestalfilename;
  string m_calibrationfilename;

  // Take input directory from global parameters
  string inputDirectory = config_.getPath("input_directory");

  // Open the root directory
  DIR* directory = opendir(inputDirectory.c_str());
  if(directory == nullptr) {
      LOG(ERROR) << "Directory " << inputDirectory << " does not exist";
      return;
  }
  dirent* entry;

  // Read the files (data, pedestal and calibration) in the folder
  while(entry = readdir(directory)) {
      if(entry->d_type == isFile) {
          string filename = inputDirectory + "/" + entry->d_name;
          if(filename.find("dat") != string::npos){
            m_datafilename = filename;
          }
          if(filename.find("ped") != string::npos){
            m_pedestalfilename = filename;
          }
          if(filename.find("cal") != string::npos){
            m_calibrationfilename = filename;
          }
      }
  }

  if(m_datafilename.length() == 0) {
      LOG(ERROR)("No data file was found for ALiBaVa in " + inputDirectory);
      return;
  }

  if(m_pedestalfilename.length() == 0) {
      LOG(ERROR) << "No pedestal file was found." << "\n" << "Cannot continue without pedestal information.";
      return;
  }

  if(m_calibrationfilename.length() == 0) {
      LOG(WARNING) << "No calibration file was found." << "\n" << "Results will be uncalibrated.";
  }

  // Open the ALiBaVa data
  ALiBaVaPointer = DataFileRoot::OpenFile();
  ALiBaVa_loader(ALiBaVaPointer, m_datafilename, m_pedestalfilename, m_calibrationfilename);

}

StatusCode EventLoaderALiBaVa::run(const shared_ptr<Clipboard>&) {

  while(!ALiBaVaPointer->read_event()){
      nEvents++;
      PixelVector pixels;
      ALiBaVaPointer->process_event();
      double eventTime = ALiBaVaPointer->time();

      for(int chan = 0; chan < ALiBaVaPointer->nchan(); chan++){
        calibration = ALiBaVaPointer->get_gain();
        CalSignal = ALiBaVaPointer->signal(chan);
        ADCSignal = ALiBaVaPointer->signal(chan)*calibration;
        auto pixel = make_shared<Pixel>(m_detector->getName(), chan, 0, ADCSignal, CalSignal, eventTime);
        pixels.push_back(pixel);
      }

  }

  clipboard->putData(pixels, m_detector->getName());
}

void EventLoaderALiBaVa::finalize(const shared_ptr<ReadonlyClipboard>&) {
    ALiBaVaPointer.close();
    delete ALiBaVaPointer;
    LOG(DEBUG) << "Analysed " << nEvents << " events";
}
