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

  // Prepare variables for the seed and neighbours
  std::shared_ptr<Pixel> left_neighbour ;
  std::shared_ptr<Pixel> seed;
  std::shared_ptr<Pixel> right_neighbour;

  // Prepare cluster
  ClusterVector deviceClusters;
  auto cluster = std::make_shared<Cluster>();

  // Prepare variables for finding the pixels in the cluster
  std::shared_ptr<Pixel> previousPixel;
  double greatestSeedSNRatio = 0;
  bool pixelWasSeed = false;
  bool seedFound = false;
  bool left_neighbourFound = false;
  bool right_neighbourFound = false;
  bool firstChannel = true;

  // Find the seed and direct neighbours (to the left or right)
  // We only define up to one cluster per event
  // The seed and neighbours have a signal-noise ratio requirement
  // The SNratio*100000 should be saved as the raw pixel value in the eventLoader
  // Which is the case for the ALiBaVa eventLoader
  for(auto pixel : pixels) {
    double pixelSNRatio = pixel->raw()/100000.0;
    if(pixelWasSeed){
      if(pixelSNRatio > m_neighbourThreshold){
        right_neighbour = pixel;
        right_neighbourFound = true;
      }
      pixelWasSeed = false;
    }
    if(pixelSNRatio > m_seedThreshold && pixelSNRatio > greatestSeedSNRatio){
      seed = pixel;
      seedFound = true;
      left_neighbourFound = false;
      right_neighbourFound = false;
      if(!firstChannel && previousPixel->charge() > m_neighbourThreshold){
         left_neighbour = previousPixel;
         left_neighbourFound = true;
       }
      greatestSeedSNRatio = pixelSNRatio;
      pixelWasSeed = true;
    }
    previousPixel = pixel;
    firstChannel = false;
  }

  // If a seed was found, add it to the cluster
  // Cluster timestamp = seed timestamp
  // Left or right neighbours aren't strictly necessary to form a cluster
  if(seedFound){
    cluster->addPixel(&*seed);
    cluster->setTimestamp(seed->timestamp());

    // This is for the calculation of the eta distribution
    // and to add the neighbours to the cluster
    double seedCharge=seed->charge();;
    double leftNeighbourCharge=0.0;
    double rightNeighbourCharge=0.0;
    double eta=0.0;
    if(left_neighbourFound){
      cluster->addPixel(&*left_neighbour);
      leftNeighbourCharge = left_neighbour->charge();
    };
    if(right_neighbourFound){
       cluster->addPixel(&*right_neighbour);
       rightNeighbourCharge = right_neighbour->charge();
     };
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
    clusterSeedCharge->Fill(cluster->getSeedPixel()->charge());
    clusterCharge->Fill(seedCharge+leftNeighbourCharge+rightNeighbourCharge);
    clusterPosition->Fill(seed->column());
    if(right_neighbourFound || left_neighbourFound) etaDistribution->Fill(eta);
    // Put the cluster on the clipboard
    deviceClusters.push_back(cluster);
    clipboard->putData(deviceClusters, m_detector->getName());
  }
  return StatusCode::Success;
}

void ClusteringSeed::finalize(const std::shared_ptr<ReadonlyClipboard>&) {

    LOG(DEBUG) << "Analysed " << m_eventNumber << " events";
}
