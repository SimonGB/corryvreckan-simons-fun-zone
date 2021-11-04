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

    config_.setDefault<double>("SeedThreshold", 3.0);
    config_.setDefault<double>("NeighbourThreshold", 1.8);

    double m_seedThreshold = config_.get<double>("SeedThreshold");
    double m_neighbourThreshold = config_.get<double>("NeighbourThreshold");


    std::string title = m_detector->getName() + " Cluster size;cluster size [#];events";
    clusterSize = new TH1F("clusterSize", title.c_str(), 3, 1, 4);

    title = m_detector->getName() + " Cluster seed charge;cluster seed charge [e];events";
    clusterSeedCharge = new TH1F("clusterSeedCharge", title.c_str(), 35, 0, 3500);

    title = m_detector->getName() + " Cluster Charge;cluster charge [e];events";
    clusterCharge = new TH1F("clusterCharge", title.c_str(), 80, 0, 10500);

    title = m_detector->getName() + " Cluster Position;x [px];events";
    clusterPosition = new TH1F("clusterPosition", title.c_str(), m_detector->nPixels().X(), -0.5, m_detector->nPixels().X()-0.5);

    title = m_detector->getName() + " Eta Distribution;eta distribution [eta];events";
    etaDistribution = new TH1F("etaDistribution", title.c_str(), 100, -0.5, 1.5);
}

StatusCode ClusteringSeed::run(const std::shared_ptr<Clipboard>& clipboard) {

  // Get the pixels in the strip
  auto pixels = clipboard->getData<Pixel>(m_detector->getName());
  if(pixels.empty()) {
      LOG(DEBUG) << "Detector " << m_detector->getName() << " does not have any pixels on the clipboard";
      return StatusCode::Success;
  }

  // Prepare cluster vector and clusters
  // For strip detectors we only define 1 cluster at most per events
  // Strictly speaking the vector is not necessary, but it's easier to work like this
  ClusterVector deviceClusters;
  auto cluster = std::make_shared<Cluster>();
  PixelVector neighbourCandidates; // for the neighbours
  std::shared_ptr<Pixel> seed;

  double greatestSeedSNRatio = 0.;
  int seedCoordinate = -1;
  int nChannels = 0;

  // The seed and neighbours have a signal-noise ratio requirement
  // The SNratio*100000 should be saved as the raw pixel value in the eventLoader
  // Which is the case for the ALiBaVa eventLoader
  for(auto pixel : pixels){
    neighbourCandidates.push_back(pixel);
    double pixelSNRatio = pixel->raw()/100000.0;
    if(pixelSNRatio > m_seedThreshold && pixelSNRatio > greatestSeedSNRatio){
      greatestSeedSNRatio = pixelSNRatio;
      seedCoordinate = nChannels;
      seed = pixel;
    }
    nChannels++;
  }

  if(seedCoordinate != -1){ // i.e. if a seed is found
    cluster->addPixel(&*seed);
    bool right_neighbourFound = false;
    bool left_neighbourFound = false;

    bool clusterContact = true;
    for(int rightChan = seedCoordinate+1; rightChan < nChannels; rightChan++){
      if(!clusterContact) break;
      clusterContact = false;
      std::shared_ptr<Pixel> rNeighbourCandidate = neighbourCandidates[rightChan];
      double rNeighCandSNRatio = rNeighbourCandidate->raw()/100000.0;
      if(rNeighCandSNRatio > m_neighbourThreshold){
        right_neighbourFound = true;
        clusterContact = true;
        cluster->addPixel(&*rNeighbourCandidate);
      }
    }

    clusterContact = true;
    for(int leftChan = seedCoordinate-1; leftChan >= 0; leftChan--){
      if(!clusterContact) break;
      clusterContact = false;
      std::shared_ptr<Pixel> lNeighbourCandidate = neighbourCandidates[leftChan];
      double lNeighCandSNRatio = lNeighbourCandidate->raw()/100000.0;
      if(lNeighCandSNRatio > m_neighbourThreshold){
        left_neighbourFound = true;
        clusterContact = true;
        cluster->addPixel(&*lNeighbourCandidate);
      }
    }

    // This is for the calculation of the eta distribution
    // and to add the neighbours to the cluster
    double seedCharge=seed->charge();
    double leftNeighbourCharge=-9999999.;
    double rightNeighbourCharge=-9999999.;
    double eta=0.0;

    if(left_neighbourFound){
      leftNeighbourCharge = neighbourCandidates[seedCoordinate-1]->charge();
    }
    if(right_neighbourFound){
       rightNeighbourCharge = neighbourCandidates[seedCoordinate+1]->charge();
     }

     if(leftNeighbourCharge > rightNeighbourCharge){
       eta = leftNeighbourCharge/(leftNeighbourCharge+seedCharge);
     }
     else if(leftNeighbourCharge < rightNeighbourCharge){
       eta = seedCharge/(seedCharge+rightNeighbourCharge);
     }
     else{
       eta = -1.0;
     }

    // Fill the histograms
    clusterSize->Fill(static_cast<double>(cluster->size()));
    clusterSeedCharge->Fill(seedCharge);
    clusterCharge->Fill(static_cast<double>(cluster->charge()));
    clusterPosition->Fill(seedCoordinate);
    if(right_neighbourFound || left_neighbourFound) etaDistribution->Fill(eta);
    // If there are values of -1 in the eta distribution, this means that
    // there are neighbours on both sides with the exact same charge
    // Put the cluster on the clipboard
    deviceClusters.push_back(cluster);
    clipboard->putData(deviceClusters, m_detector->getName());
  }

  return StatusCode::Success;
}

void ClusteringSeed::finalize(const std::shared_ptr<ReadonlyClipboard>&) {

    LOG(DEBUG) << "Analysed " << m_eventNumber << " events";
}
