/**
 * @file
 * @brief Implementation of module ClusteringSeed
 *
 * @copyright Copyright (c) 2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */
//test
#include "ClusteringSeed.h"
#include "objects/Pixel.hpp"

using namespace corryvreckan;
using namespace std;

ClusteringSeed::ClusteringSeed(Configuration& config, std::shared_ptr<Detector> detector)
    : Module(config, detector), m_detector(detector) {}

void ClusteringSeed::initialize() {

    for(auto& detector : get_detectors()) {
        LOG(DEBUG) << "Initialise for detector " + detector->getName();
    }

    config_.setDefault<double>("SeedThreshold", 5.);
    config_.setDefault<double>("NeighbourThreshold", 2.5);
    config_.setDefault<int>("LowerChannel", 0);
    config_.setDefault<int>("UpperChannel", 127);
    config_.setDefault<bool>("CalculateCrosstalk", true);

    m_seedThreshold = config_.get<double>("SeedThreshold");
    m_neighbourThreshold = config_.get<double>("NeighbourThreshold");
    m_lower_channel = config_.get<int>("LowerChannel");
    m_upper_channel = config_.get<int>("UpperChannel");
    m_calculate_crosstalk = config_.get<bool>("CalculateCrosstalk");

    std::string title = m_detector->getName() + " Cluster size;cluster size [#];events";
    clusterSize = new TH1F("clusterSize", title.c_str(), 15, 0, 15);

    title = m_detector->getName() + " Cluster seed charge;cluster seed charge [e];events";
    clusterSeedCharge = new TH1F("clusterSeedCharge", title.c_str(), 40, 2500, 8000);

    title = m_detector->getName() + " Cluster Charge;cluster charge [e];events";
    clusterCharge = new TH1F("clusterCharge", title.c_str(), 80, 2000, 35000);

    title = m_detector->getName() + " Cluster Position;x [px];events";
    clusterPosition = new TH1F("clusterPosition", title.c_str(), 128, 0, 128);

    title = m_detector->getName() + " Eta Distribution;eta distribution [eta];events";
    etaDistribution = new TH1F("etaDistribution", title.c_str(), 100, 0, 1);

    title = m_detector->getName() + " Signal left 1;relative charge to seed ;events";
    leftSignal = new TH1F("leftSignal", title.c_str(), 100, -2, 2);
    title = m_detector->getName() + " Signal right 1;relative charge to seed ;events";
    rightSignal = new TH1F("rightSignal", title.c_str(), 100, -2, 2);
    title = m_detector->getName() + " Signal left 2;relative charge to seed ;events";
    leftBOSignal = new TH1F("leftBOSignal", title.c_str(), 100, -2, 2);
    title = m_detector->getName() + " Signal right 2;relative charge to seed ;events";
    rightBOSignal = new TH1F("rightBOSignal", title.c_str(), 100, -2, 2);

    title = m_detector->getName() + " Debug event SNR;channel ;SNR";
    DEBUG_event_SNR = new TH1F("DEBUG_event_SNR", title.c_str(), 128, 0, 128);
    title = m_detector->getName() + " Debug event charge;channel ;charge";
    DEBUG_event_charge = new TH1F("DEBUG_event_charge", title.c_str(), 128, 0, 128);
}

StatusCode ClusteringSeed::run(const std::shared_ptr<Clipboard>& clipboard) {

  // Get the pixels in the strip
  auto pixels = clipboard->getData<Pixel>(m_detector->getName());
  if(pixels.empty()) {
      LOG(DEBUG) << "Detector " << m_detector->getName() << " does not have any pixels on the clipboard";
      return StatusCode::Success;
  }

  // Prepare cluster vector and clusters
  // For strip detectors we only define 1 cluster at most per event
  // Strictly speaking the clustervector is not necessary, but it's easier to work like this
  ClusterVector deviceClusters;
  auto cluster = std::make_shared<Cluster>();
  PixelVector neighbourCandidates; // for the neighbours
  std::shared_ptr<Pixel> seed;

  double greatestSeedSNRatio = 0.;
  int seedCoordinate = -1;
  int nChannels = 0;
  double clusterChargeVal = 0.;
  double seedCharge = 0.;
  int clusterSizeVal = 0;

  // The seed and neighbours have a signal-noise ratio requirement
  // The SNratio*100000 should be saved as the raw pixel value in the eventLoader
  // Which is the case for the ALiBaVa eventLoader
  for(auto pixel : pixels){
    neighbourCandidates.push_back(pixel);
    nChannels++;
    if(nChannels-1 == 0 || nChannels-1 == m_upper_channel-m_lower_channel) continue;
    double pixelSNRatio = pixel->raw()/100000.0;
    if(pixelSNRatio > m_seedThreshold && pixelSNRatio > greatestSeedSNRatio){
      greatestSeedSNRatio = pixelSNRatio;
      seedCoordinate = nChannels-1;
      seed = pixel;
      seedCharge = pixel->charge();
      clusterChargeVal = seedCharge;
      clusterSizeVal = 1;
    }
  }

  if(seedCoordinate != -1){ // i.e. if at least one seed is found

    cluster->addPixel(&*seed);
    bool right_neighbourFound = false;
    bool left_neighbourFound = false;
    bool rightBO_neighbourFound = false;
    bool leftBO_neighbourFound = false;

    bool clusterContact = true;
    for(int rightChan = seedCoordinate+1; rightChan < nChannels; rightChan++){
      if(!clusterContact) break;
      clusterContact = false;
      std::shared_ptr<Pixel> rNeighbourCandidate = neighbourCandidates[rightChan];
      double rNeighCandSNRatio = rNeighbourCandidate->raw()/100000.0;
      if(rNeighCandSNRatio > m_neighbourThreshold){
        if(right_neighbourFound) rightBO_neighbourFound = true;
        right_neighbourFound = true;
        clusterContact = true;
        cluster->addPixel(&*rNeighbourCandidate);
        clusterChargeVal += rNeighbourCandidate->charge();
        clusterSizeVal++;
      }
    }

    clusterContact = true;
    for(int leftChan = seedCoordinate-1; leftChan >= 0; leftChan--){
      if(!clusterContact) break;
      clusterContact = false;
      std::shared_ptr<Pixel> lNeighbourCandidate = neighbourCandidates[leftChan];
      double lNeighCandSNRatio = lNeighbourCandidate->raw()/100000.0;
      if(lNeighCandSNRatio > m_neighbourThreshold){
        if(left_neighbourFound) leftBO_neighbourFound = true;
        left_neighbourFound = true;
        clusterContact = true;
        cluster->addPixel(&*lNeighbourCandidate);
        clusterChargeVal += lNeighbourCandidate->charge();
        clusterSizeVal++;
      }
    }

    // This is for the calculation of the eta distribution
    // and to add the neighbours to the cluster
    double seedCharge=seed->charge();
    double leftNeighbourCharge;
    double rightNeighbourCharge;
    double eta;
    bool equalCharges = false;

     // if(left_neighbourFound && right_neighbourFound){
     //     QLeftOverQSeed += leftNeighbourCharge/seedCharge;
     //     QRightOverQSeed += rightNeighbourCharge/seedCharge;
     //   }

     // This could be shortened considerably and made more efficient
     // but I changed it to this because it's more clear for debugging
     if(leftBO_neighbourFound && rightBO_neighbourFound && left_neighbourFound && right_neighbourFound){
        QLeftOverQSeed += neighbourCandidates[seedCoordinate-1]->charge()/seedCharge;
        QRightOverQSeed += neighbourCandidates[seedCoordinate+1]->charge()/seedCharge;
        QLeftBOOverQSeed += neighbourCandidates[seedCoordinate-2]->charge()/seedCharge;
        QRightBOOverQSeed += neighbourCandidates[seedCoordinate+2]->charge()/seedCharge;
        leftSignal->Fill(neighbourCandidates[seedCoordinate-1]->charge()/seedCharge);
        rightSignal->Fill(neighbourCandidates[seedCoordinate+1]->charge()/seedCharge);
        leftBOSignal->Fill(neighbourCandidates[seedCoordinate-2]->charge()/seedCharge);
        rightBOSignal->Fill(neighbourCandidates[seedCoordinate+2]->charge()/seedCharge);


        // This is just debug stuff
        debugVar++;
        if(debugVar == 50){
          for(int chan = 0; chan <= 128; chan++){
            if(chan>m_lower_channel && chan<=m_upper_channel+1){
              DEBUG_event_SNR->SetBinContent(chan, neighbourCandidates[chan-m_lower_channel-1]->raw()/100000.0);
              DEBUG_event_charge->SetBinContent(chan, neighbourCandidates[chan-m_lower_channel-1]->charge());
            }
            else{
              DEBUG_event_SNR->SetBinContent(chan, 0);
              DEBUG_event_charge->SetBinContent(chan, 0);
            }
          }
        }
        // End debug stuff
      }

      if(left_neighbourFound && right_neighbourFound){
        double leftNeighbourCharge = neighbourCandidates[seedCoordinate-1]->charge();
        double rightNeighbourCharge = neighbourCandidates[seedCoordinate+1]->charge();
        if(leftNeighbourCharge > rightNeighbourCharge){
          eta = leftNeighbourCharge/(leftNeighbourCharge+seedCharge);
        }
        else if(leftNeighbourCharge < rightNeighbourCharge){
          eta = seedCharge/(seedCharge+rightNeighbourCharge);
        }
        else equalCharges = true;
      }
      else if(left_neighbourFound && !right_neighbourFound){
        eta = neighbourCandidates[seedCoordinate-1]->charge()/(neighbourCandidates[seedCoordinate-1]->charge()+seedCharge);
      }
      else if(!left_neighbourFound && right_neighbourFound){
        eta = seedCharge/(seedCharge+neighbourCandidates[seedCoordinate+1]->charge());
      }
      else{}



     //  if(left_neighbourFound){
     //      leftNeighbourCharge = neighbourCandidates[seedCoordinate-1]->charge();
     //  }
     //  if(right_neighbourFound){
     //      rightNeighbourCharge = neighbourCandidates[seedCoordinate+1]->charge();
     //   }
     //
     // if(leftNeighbourCharge > rightNeighbourCharge){
     //   eta = leftNeighbourCharge/(leftNeighbourCharge+seedCharge);
     // }
     // else if(leftNeighbourCharge < rightNeighbourCharge){
     //   eta = seedCharge/(seedCharge+rightNeighbourCharge);
     // }
     // else{
     //   eta = -1.0;
     // }

    // Fill the histograms
    clusterSize->Fill(clusterSizeVal);
    clusterSeedCharge->Fill(seedCharge);
    clusterCharge->Fill(clusterChargeVal);
    clusterPosition->Fill(m_lower_channel+seedCoordinate);
    if((right_neighbourFound || left_neighbourFound) & !equalCharges) etaDistribution->Fill(eta);
    // If there are values of -1 in the eta distribution, this means that
    // there are neighbours on both sides with the exact same charge...
    // Put the cluster on the clipboard
    deviceClusters.push_back(cluster);
    clipboard->putData(deviceClusters, m_detector->getName());
  }

  return StatusCode::Success;
}


StatusCode ClusteringSeed::runTEMP(const std::shared_ptr<Clipboard>& clipboard){
if(false){
  // Get the pixels in the strip
  auto pixels = clipboard->getData<Pixel>(m_detector->getName());
  if(pixels.empty()) {
      LOG(DEBUG) << "Detector " << m_detector->getName() << " does not have any pixels on the clipboard";
      return StatusCode::Success;
  }

  std::shared_ptr<Pixel> orderedStrip[128];
  // The seed and neighbours have a signal-noise ratio requirement
  // The SNratio*100000 should be saved as the raw pixel value in the eventLoader
  // Which is the case for the ALiBaVa eventLoader
  for(auto pixel : pixels){
    orderedStrip[pixel->column()] = pixel;
  }

  std::vector<PixelVector> pixelChainCollection;
  bool lastPixelNeighbour = false;
  bool seedFound = false;
  PixelVector pixelChain;
  for(int ichan = 0; ichan < 128; ichan++){
    if(!lastPixelNeighbour){
      if(seedFound){
        pixelChainCollection.push_back(pixelChain);
      }
      seedFound = false;
      pixelChain.clear();
    }
    lastPixelNeighbour = false;
    double pixelSNRatio = orderedStrip[ichan]->raw()/100000.0;
    if(pixelSNRatio > m_neighbourThreshold){
      pixelChain.push_back(orderedStrip[ichan]);
      if(pixelSNRatio > m_seedThreshold){
        seedFound = true;
      }
      lastPixelNeighbour = true;
    }
  }
}
}


void ClusteringSeed::finalize(const std::shared_ptr<ReadonlyClipboard>&) {
    if(m_calculate_crosstalk) LOG(WARNING) << "Crosstalk correction calculation results: b1 = " << QLeftOverQSeed-QRightOverQSeed << " and b2 = " << QLeftBOOverQSeed-QRightBOOverQSeed;
    LOG(DEBUG) << "Analysed " << m_eventNumber << " events";
}
