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

StatusCode ClusteringSeed::runOLD(const std::shared_ptr<Clipboard>& clipboard) {

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


StatusCode ClusteringSeed::run(const std::shared_ptr<Clipboard>& clipboard){
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

  // This algorithm assumes that there is no limit to the size of clusters and
  // No limit to the amount of clusters

  // Take all pixels and put them in an ordered (by channel) array
  for(auto pixel : pixels){
    orderedStrip[pixel->column()] = pixel;
  }
  // Create a collection of pixelchains
  // A bool to remember if the previous pixel in the loop passed the neighbourcut
  // A bool to see if a seed was found for the current chain of pixels
  // A chain of pixels that pass at least the neighbour cut
  std::vector<PixelVector> pixelChainCollection;
  bool lastPixelNeighbour = false;
  bool seedFound = false;
  PixelVector pixelChain;
  // This loops over all possible channels for one chip
  for(int ichan = 0; ichan < 128; ichan++){
    // Check if the previous pixel in the loop was NOT a neighbour
    if(!lastPixelNeighbour){
      // If a seed was found and the last pixel was NOT a neighbour
      // this means we have a cluster and have identified it from beginning to end.
      // Add it to the collection of pixelChains (i.e. essentially clusters)
      if(seedFound){
        pixelChainCollection.push_back(pixelChain);
      }
      // Reset the bool that identifies if a seed was found
      // since we're now looking for a new cluster.
      // Clear the current pixelChain variable (i.e. cluster)
      seedFound = false;
      pixelChain.clear();
    }
    // Reset the bool for identifying if the last pixel passed the neighbourcut
    // since it has served its purpose here and will now save the status of the current pixel
    // for the next part of the loop.
    lastPixelNeighbour = false;
    // Get the pixel signal-noise ratio
    double pixelSNRatio = orderedStrip[ichan]->raw()/100000.0;
    // If the ratio passes the threshold, put the pixel in the pixelchain and update
    // the bool for identifying if the pixel passes the neighbourcut
    if(pixelSNRatio > m_neighbourThreshold){
      pixelChain.push_back(orderedStrip[ichan]);
      // Additionally, check if the current pixel is a seed.
      // If so, update the corresponding bool, which is relevant not just for one loop iteration
      // but for as long as there is an unbroken chain of pixels that pass at least the
      // neighbourcut, i.e. as long as we are in a pixelChain or (if there's a seed) cluster.
      if(pixelSNRatio > m_seedThreshold){
        seedFound = true;

      }
      lastPixelNeighbour = true;
    }
  }
  // After this loop, we should end up with a vector (pixelChainCollection) containing vectors (pixelchains).
  // The latter objects contain clusters, i.e. unbroken chains of pixels that all pass the neighbour cut AND
  // have at least one seed.


  for(auto cluster : pixelChainCollection){
    int clusterSizeVal = 0;
    double clusterSeedChargeVal = 0;
    double clusterChargeVal = 0;
    double seedPositionVal = -1;
    for(auto pixel : cluster){
      double pixelSNRatio = pixel->raw()/100000.0;
      double pixelCharge = pixel->charge();
      clusterSizeVal++;
      clusterChargeVal += pixelCharge;
      if(pixelSNRatio > m_seedThreshold && pixelSNRatio > clusterSeedChargeVal){
          clusterSeedChargeVal = pixelCharge;
          seedPositionVal = pixel->column();
      }
    }
    clusterSize->Fill(clusterSizeVal);
    clusterSeedCharge->Fill(clusterSeedChargeVal);
    clusterCharge->Fill(clusterChargeVal);
    clusterPosition->Fill(seedPositionVal);
    etaDistribution->Fill(0);
  }

/*
  // This part deals with clusters that have multiple seeds.
  // Generally, I wouldn't expect there to be more than 2 seeds per cluster.
  // If there are, then you probably haven't set reasonable neighbour and seed cuts.
  // It could probably be integrated into the above algorithm, but it would be messy.
  // Maybe do that later for optimisation/speed-up
  for(auto cluster : pixelChainCollection){
    int seedCount = 0;
    double firstSeedSNR = 0;
    double secondSeedSNR = 0;
    int squeezedNeighbourCount = 0;
    bool firstSeedAhead = true;
    bool secondSeedAhead = true;
    for(auto pixel : cluster){
      double pixelSNRatio = pixel->raw()/100000.0;
      if(pixelSNRatio > m_seedThreshold){
        if(firstSeedAhead){
          firstSeedSNR = pixelSNRatio;
          firstSeedAhead = false;
        }
        if(secondSeedAhead && !firstSeedAhead){
          secondSeedSNR = pixelSNRatio;
          secondSeedAhead = false;
        }
        seedCount++
      }
      else if(pixelSNRatio > m_neighbourThreshold){
        if(!firstSeedAhead && secondSeedAhead) squeezedNeighbourCount++;
      }
    }

    if(seedCount == 2){
      if(squeezedNeighbourCount == 0){
        for(auto pixel : cluster){
          // Make highest SNratio pixel the seed
          // Keep the cluster, just turn the other seed into a neighbour
          // This should happen automatically in Corryvreckan though...
        }
      }
      else if(squeezedNeighbourCount == 1){
        PixelVector newCluster_one;
        PixelVector newCluster_two;
        for(auto pixel : cluster){
          // Give the neighbour to the highest SNratio seed
          bool firstSeedAhead = true;
          bool secondSeedAhead = true;
          if(firstSeedSNR > secondSeedSNR){
            double pixelSNRatio = pixel->raw()/100000.0;

            if(firstSeedAhead && secondSeedAhead){
              newCluster_one.push_back(pixel)
            }
            if(pixelSNRatio > m_seedThreshold){

            }
            newCluster_one.push_back(pixel)
          }
          else{}
        }
        cluster = pixelChainCollection.erase(cluster);
      }
      else if(squeezedNeighbourCount == 2){
        for(auto pixel : cluster){
          // Split into two clusters
        }
        cluster = pixelChainCollection.erase(cluster);
      }
      else if(squeezedNeighbourCount > 2){
        LOG(WARNING) << "There are more than 2 neighbours in between 2 seeds in the same cluster, change your neighbour- and seedcuts";
      }
    }
    else if(seedCount > 2){
      LOG(WARNING) << "There are clusters with more than 2 seeds, change your neighbour- and seedcuts";
      return StatusCode::Failure;
    }
  }
*/
}



void ClusteringSeed::finalize(const std::shared_ptr<ReadonlyClipboard>&) {
    if(m_calculate_crosstalk) LOG(WARNING) << "Crosstalk correction calculation results: b1 = " << QLeftOverQSeed-QRightOverQSeed << " and b2 = " << QLeftBOOverQSeed-QRightBOOverQSeed;
    LOG(DEBUG) << "Analysed " << m_eventNumber << " events";
}
