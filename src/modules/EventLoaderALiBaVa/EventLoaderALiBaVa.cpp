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

  // Take input directory, run number,
  // lower timecut, upper timecut,        --- in nanoseconds
  // first X ALiBaVa events to ignore, channel range,
  // from global parameters in the config.
  // Default values are set first

  config_.setDefault<bool>("get_time_residuals", false);
  config_.setDefault<int>("run", 0);
  config_.setDefault<int>("LowerTimecut", 0);
  config_.setDefault<int>("UpperTimecut", std::numeric_limits<int>::max());
  config_.setDefault<int>("IgnoreEvents", 1);
  config_.setDefault<int>("LowerChannel", 0);
  config_.setDefault<int>("UpperChannel", 255);
  config_.setDefault<int>("Chargecut", std::numeric_limits<int>::max());

  std::string m_inputDirectory = config_.getPath("input_directory");
  int m_run = config_.get<int>("run");
  int m_timecut_lower = config_.get<int>("LowerTimecut");
  int m_timecut_upper = config_.get<int>("UpperTimecut");
  int m_ignore_events = config_.get<int>("IgnoreEvents");
  int m_lower_channel = config_.get<int>("LowerChannel");
  int m_upper_channel = config_.get<int>("UpperChannel");
  int m_chargecut = config_.get<int>("Chargecut");

  // Open the input directory
  DIR* directory = opendir(m_inputDirectory.c_str());
  if(directory == nullptr) {
      LOG(ERROR) << "Directory \'" << m_inputDirectory << "\' does not exist or was not supplied";
      return;
  }

  if(m_run == EMPTY){

  }

  // Read the run-files (data, pedestal and calibration) in the folder
  dirent* entry;
  while(entry = readdir(directory)) {
      if(entry->d_type == DT_REG) {
          std::string entryName = entry->d_name;
          if(entryName.find((std::to_string(m_run)+"*.dat").c_str()) != std::string::npos){
            m_datafilename = m_inputDirectory + "/" + entryName;
          }
          if(entryName.find(std::to_string(m_run)+"*.ped") != std::string::npos){
            m_pedestalfilename = m_inputDirectory + "/" + entryName;
          }
          if(entryName.find(std::to_string(m_run)+"*.cal") != std::string::npos){
            m_calibrationfilename = m_inputDirectory + "/" + entryName;
          }
      }
  }

  // Log errors in case the files aren't found in the folder.
  // The datafile can also be supplied directly in the config.
  if(m_datafilename.length() == 0) {
      LOG(WARNING) << "No data file was found for ALiBaVa in " << m_inputDirectory;
      m_datafilename = config_.get<std::string>("filename");
      if(m_datafilename.length() == 0) {
        LOG(ERROR) << << "No data filepath was supplied in the config. There is no data file.";
        return;
      }
      else{
        LOG(WARNING) << "Opening" << m_datafilename << " directly from filepath in config.";
      }
  }

  if(m_pedestalfilename.length() == 0) {
      LOG(WARNING) << "No pedestal file was found." << "\n" << "Datafile will be used for pedestal";
  }

  if(m_calibrationfilename.length() == 0) {
      LOG(WARNING) << "No calibration file was found." << "\n" << "Results will be uncalibrated: ADC = charge.";
  }

  // Create a pointer with the data file.
  ALiBaVaPointer = DataFileRoot::OpenFile(m_datafilename.c_str());
  // Call the loader function, which sets everything up.
  // For some reason in ALiBaVa's code, it needs the pointer AND the data file again: AliBaVa and the 40 redundancies?
  ALiBaVa_loader(ALiBaVaPointer, m_datafilename.c_str(), m_pedestalfilename.c_str(), m_calibrationfilename.c_str());
  // Set the timecuts
  ALiBaVaPointer->set_timecut(m_timecut_lower, m_timecut_upper)
  // Check how many events the data file contains.
  nEvents = ALiBaVaPointer->nevents();

// Ignore the first X events to ensure synchronisation, default is X = 1 which ignores the first event.
for(int ievt = 0; ievt < m_ignore_events; ievt++){
    ALiBaVaPointer->read_event()
  }

  // Create histograms
  chargeHist = new TH1F("charge","charge",5000,-2500,2500);
}

StatusCode EventLoaderALiBaVa::run(const std::shared_ptr<Clipboard>& clipboard) {
  // During running, every Corryvreckan event will get one ALiBaVa event
  // Increment the event counter
  iEvent++;
  // Create the pixelvector
  PixelVector pixels;
  // Get the event that already exists on the Corryvreckan clipboard
  auto event = clipboard->getEvent();
  // Read a data event from the ALiBaVa data file
  // give feedback according to return code
  int return_code = ALiBaVaPointer->read_event();
  if(return_code == -1){
    LOG(WARNING) << "End of data file reached.";
    return StatusCode::EndRun;
  }
  elif(return_code == 0){
    LOG(ERROR) << "There's something wrong (0) with the datafile at event number " << iEvent-1;
    LOG(ERROR) << "Terminating run";
    return StatusCode::EndRun;
  }
  elif(return_code == 1){
    // This means the event was read properly.
  }
  elif(return_code == 2 ){
    LOG(ERROR) << "There's something wrong (2) with the datafile at event number " << iEvent-1;
    LOG(ERROR) << "Terminating run";
    return StatusCode::EndRun;
  }
  elif(return_code == 3){
    LOG(ERROR) << "There's something wrong (3) with the datafile at event number " << iEvent-1;
    LOG(ERROR) << "Terminating run";
    return StatusCode::EndRun;
  }
  elif(return_code == 4){
    LOG(ERROR) << "There's something wrong (4) with the HDF5 datafile at event number " << iEvent-1;
    LOG(ERROR) << "Terminating run";
    return StatusCode::EndRun;
  }
  else(){
    LOG(ERROR) << "This shouldn't happen. Current event number: " << iEvent-1;";
    LOG(ERROR) << "Terminating run";
    return StatusCode::EndRun;
  }

  // Process the opened data event, i.e. pedestal correction, common mode noise correction
  ALiBaVaPointer->process_event();

  // This gets the TDC time from the event, allowing timecuts around the event peak
  // The timecut is set in the ALiBaVa_loader() function.
  double TDCTime = ALiBaVaPointer->time();
  if(!ALiBaVaPointer->valid_time(TDCTime)){
    clipboard->putData(pixels, m_detector->getName());
    return StatusCode::NoData;
  }

  // This loops over the channels in the current ALiBaVa event
  for(int chan = m_lower_channel; chan <= m_upper_channel; chan++){
    // In order, these are the calibration factor, the calibrated signal, and the ADC signal
    // If a pedestal file is supplied, pedestal and common mode error correction is applied to all
    double calibration = ALiBaVaPointer->get_gain(chan);
    double CalSignal = ALiBaVaPointer->signal(chan);
    double ADCSignal = ALiBaVaPointer->signal(chan)*calibration;

    // The chargecut is applied here
    if(CalSignal < m_chargecut){
      // Create a pixel for every channel in this event with all the information and put it in the vector.
      auto pixel = std::make_shared<Pixel>(m_detector->getName(), chan, 0, ADCSignal, CalSignal, 0);
      pixels.push_back(pixel);

      // Fill the histograms
      chargeHist->Fill(CalSignal);
    }
  }

  // Put the created vector of pixels on the clipboard.
  clipboard->putData(pixels, m_detector->getName());

  // If the pixels vector is empty, report this to Corryvreckan
  if(pixels.empty()){
    return StatusCode::NoData;
  }

  // Report the end of this event to Corryvreckan
  return StatusCode::Success;
}


void EventLoaderALiBaVa::finalize(const std::shared_ptr<ReadonlyClipboard>& clipboard) {
    ALiBaVaPointer->close();
    delete ALiBaVaPointer;
    LOG(DEBUG) << "Analysed " << iEvent << " of total " << nEvents << " ALiBaVa events";
}
