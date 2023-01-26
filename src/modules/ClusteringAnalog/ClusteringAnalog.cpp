/**
 * @file
 * @brief Implementation of module ClusteringAnalog
 *
 * @copyright Copyright (c) 2017-2022 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "ClusteringAnalog.h"
#include "objects/Pixel.hpp"

#include <TFile.h>

using namespace corryvreckan;
using namespace std;

ClusteringAnalog::ClusteringAnalog(Configuration& config, std::shared_ptr<Detector> detector)
    : Module(config, detector), m_detector(detector) {

    rejectByROI = config_.get<bool>("reject_by_roi", false);
    thresholdSeed = config_.get<float>("threshold_seed");
    thresholdNeighbor = config_.get<float>("threshold_neighbor", thresholdSeed);
    thresholdIteration = config_.get<float>("threshold_iteration", thresholdNeighbor);
    thresholdSeedSNR = config_.get<float>("thresholdSNR_seed", thresholdSeed);
    thresholdNeighborSNR = config_.get<float>("thresholdSNR_neighbor", thresholdNeighbor);
    thresholdIterationSNR = config_.get<float>("thresholdSNR_iteration", thresholdNeighborSNR);
    thresholdClusterCharge = config_.get<float>("threshold_cluster", thresholdSeed);
    estimationMethod = config_.get<EstimationMethod>("method", EstimationMethod::CLUSTER);
    seedingMethod = config_.get<SeedingMethod>("seeding_method", SeedingMethod::MULTI);
    thresholdType = config_.get<ThresholdType>("threshold_type", ThresholdType::FIX);
    windowSize = config_.get<int>("window_size", 1);
    if(windowSize < 1) {
        throw InvalidValueError(config_, "window_size", "Invalid window size - value should be >= 1.");
    }
    includeCorners = config_.get<bool>("include_corners", false);
    useTriggerTimestamp = config_.get<bool>("use_trigger_timestamp", false);
    flagAnalysisShape = config_.get<bool>("analysis_shape", false);
    flagAnalysisSNR = thresholdType == ThresholdType::SNR || thresholdType == ThresholdType::MIX;

    auto detConf = m_detector->getConfiguration();

    if(flagAnalysisShape) {
        auto coordinates = detConf.get<std::string>("coordinates", "cartesian");
        std::transform(coordinates.begin(), coordinates.end(), coordinates.begin(), ::tolower);
        if(coordinates != "cartesian") {
            throw InvalidCombinationError(
                detConf, {"coordinates", "analysis_shape"}, "Shape analysis implemented only for cartesian coordinates");
        }
    }

    // Read calibration file
    isCalibrated = false;
    if(thresholdType == ThresholdType::SNR || thresholdType == ThresholdType::MIX) {
        if(detConf.has("calibration_file")) {
            auto calibFilePath = detConf.getPath("calibration_file", true); // Return absolute path
            if(readCalibrationFileROOT(calibFilePath)) {
                LOG(INFO) << "Calibration file " << calibFilePath << " loaded successfully";
                isCalibrated = true;
            } else {
                throw InvalidValueError(detConf, "calibration_file", "Invalid calibration file");
            }
        } else {
            throw InvalidCombinationError(
                detConf, {"calibration_file", "threshold_type"}, "Missing calibration file, required by S/N ratio analysis");
        }
    }

    neighborsSizeCentral =
        m_detector
            ->getNeighbors(
                m_detector->nPixels().X() / 2, m_detector->nPixels().Y() / 2, static_cast<size_t>(windowSize), true)
            .size();
}

bool ClusteringAnalog::readCalibrationFileROOT(const std::filesystem::path fileName) {
    auto f = new TFile(fileName.c_str(), "READ");
    if(!f->IsOpen()) {
        return false;
    }
    // Read histogram name from conf.
    string hTemp = config_.get<string>("calibration_pedestal");
    TH2F* hSensorPedestal = dynamic_cast<TH2F*>(f->Get(hTemp.c_str())->Clone("sensorPedestal"));
    hTemp = config_.get<string>("calibration_noise");
    TH2F* hSensorNoise = dynamic_cast<TH2F*>(f->Get(hTemp.c_str())->Clone("sensorNoise"));
    hSensorPedestal->SetDirectory(nullptr);
    hSensorNoise->SetDirectory(nullptr);

    if(m_detector->nPixels().X() != hSensorNoise->GetNbinsX() || m_detector->nPixels().Y() != hSensorNoise->GetNbinsY()) {
        LOG(ERROR) << "Number of columns and/or rows in calibration file does not match those of the detector!";
        return false;
    }
    int nCols = m_detector->nPixels().X();
    int nRows = m_detector->nPixels().Y();

    noisemap.resize(static_cast<size_t>(nCols), vector<double>(static_cast<size_t>(nRows)));
    for(int col = 0; col < nCols; col++) {
        for(int row = 0; row < nRows; row++) {
            noisemap[static_cast<size_t>(col)][static_cast<size_t>(row)] = hSensorNoise->GetBinContent(col + 1, row + 1);
        }
    }
    f->Close();
    return true;
}

void ClusteringAnalog::initialize() {
    std::string title;

    title = m_detector->getName() + " Cluster selection statistics";
    hCutHisto = new TH1F("hCutHisto", title.c_str(), RejectionType::kAllTypes, 0, RejectionType::kAllTypes);
    hCutHisto->GetXaxis()->SetBinLabel(RejectionType::kOutsideROI + 1, "Outside of ROI");
    hCutHisto->GetXaxis()->SetBinLabel(RejectionType::kIncompleteEdgeNxN + 1, "Edge (size < NxN)");
    hCutHisto->GetXaxis()->SetBinLabel(RejectionType::kLowCharge + 1, "Cluster charge < thr.");
    // TO DO: make ranges configurable
    // Detector alignment/shift
    double dx = m_detector->displacement().x();
    double dy = m_detector->displacement().y();
    // Cluster plots
    title = m_detector->getName() + " Cluster size;cluster size;events";
    clusterSize = new TH1F("clusterSize", title.c_str(), 100, -0.5, 99.5);
    title = m_detector->getName() + " Cluster charge;cluster charge (ADCu);events";
    clusterCharge = new TH1F("clusterCharge", title.c_str(), 2500, -4999.5, 20000.5);
    title = m_detector->getName() + " Cluster NxN charge;cluster NxN charge (ADCu);events";
    clusterNxNCharge = new TH1F("clusterNxNCharge", title.c_str(), 2100, -999.5, 20000.5);
    title = m_detector->getName() + " Seed charge;seed charge (ADCu);events";
    clusterSeedCharge = new TH1F("clusterSeedCharge", title.c_str(), 2100, -999.5, 20000.5);
    title = m_detector->getName() + " Neighbors charge;neighbors/pixel charge (ADCu);#clusters";
    clusterNeighborsCharge = new TH1F("clusterNeighborsCharge", title.c_str(), 2000, -9999.5, 10000.5);
    title = m_detector->getName() + " Sum of neighbors charge;Charge outside seed (ADCu);#pixels";
    clusterNeighborsChargeSum = new TH1F("clusterNeighborsChargeSum", title.c_str(), 2000, -9999.5, 10000.5);

    // Cluster (central)
    title = m_detector->getName() + " Cluster size (central);cluster size;events";
    clusterSizeCentral = new TH1F("clusterSizeCentral", title.c_str(), 100, -0.5, 99.5);
    title = m_detector->getName() + " Cluster charge (central);cluster charge (ADCu);events";
    clusterChargeCentral = new TH1F("clusterChargeCentral", title.c_str(), 2500, -4999.5, 20000.5);
    title = m_detector->getName() + " Cluster NxN charge (central);cluster NxN charge (ADCu);events";
    clusterNxNChargeCentral = new TH1F("clusterNxNChargeCentral", title.c_str(), 2100, -999.5, 20000.5);
    title = m_detector->getName() + " Seed charge (central);seed charge (ADCu);events";
    clusterSeedChargeCentral = new TH1F("clusterSeedChargeCentral", title.c_str(), 2100, -999.5, 20000.5);

    // Timestamps
    title = m_detector->getName() + " Cluster timestamps;cluster timestamp [ns]; # events";
    clusterTimes = new TH1F("clusterTimes", title.c_str(), 400, 0, 100e3);

    // Correlations
    title = m_detector->getName() + " Seed charge vs cluster;seed charge (ADCu);cluster charge (ADCu);events";
    clusterCharge_SeedvsCluster =
        new TH2F("clusterCharge_SeedvsCluster", title.c_str(), 110, -999.5, 10000.5, 210, -999.5, 20000.5);
    title = m_detector->getName() + " Seed charge vs sum of neighbors;seed charge (ADCu);charge outside seed (ADCu);events";
    clusterCharge_SeedvsNeighborsSum =
        new TH2F("clusterCharge_SeedvsNeighborsSum", title.c_str(), 110, -999.5, 10000.5, 200, -9999.5, 10000.5);
    title = m_detector->getName() + " Seed charge vs neighbors;seed charge (ADCu);neighbors charge (ADCu);events";
    clusterCharge_SeedvsNeighbors =
        new TH2F("clusterCharge_SeedvsNeighbors", title.c_str(), 110, -999.5, 10000.5, 200, -9999.5, 10000.5);

    // SNR
    if(flagAnalysisSNR) {
        title = m_detector->getName() + " Cluster seed S/N;S/N ratio;events";
        clusterSeedSNR = new TH1F("clusterSeedSNR", title.c_str(), 1000, -0.5, 99.5);
        title = m_detector->getName() + " Cluster neighbor S/N;S/N ratio;events";
        clusterNeighborsSNR = new TH1F("clusterNeighborsSNR", title.c_str(), 300, -10.5, 20.5);
        title = m_detector->getName() + " Seed SNR vs cluster charge;seed S/N;cluster charge (ADCu);events";
        clusterSeedSNRvsClusterCharge =
            new TH2F("clusterSeedSNRvsClusterCharge", title.c_str(), 1000, -0.5, 99.5, 210, -999.5, 20000.5);
        title = m_detector->getName() + " Seed SNR vs neighbors;seed S/N;neighbors S/N;events";
        clusterSNR_SeedvsNeighbors =
            new TH2F("clusterSNR_SeedvsNeighbors", title.c_str(), 1000, -0.5, 99.5, 300, -10.5, 20.5);
    }

    // Positions
    title = m_detector->getName() + " Cluster position (global);x [mm];y [mm];events";
    clusterPositionGlobal = new TH2F("clusterPositionGlobal",
                                     title.c_str(),
                                     min(500, 5 * m_detector->nPixels().X()),
                                     -m_detector->getSize().X() / 1.5 + dx,
                                     m_detector->getSize().X() / 1.5 + dx,
                                     min(500, 5 * m_detector->nPixels().Y()),
                                     -m_detector->getSize().Y() / 1.5 + dy,
                                     m_detector->getSize().Y() / 1.5 + dy);
    title = m_detector->getName() + " Cluster position (local);x [px];y [px];events";
    clusterPositionLocal = new TH2F("clusterPositionLocal",
                                    title.c_str(),
                                    m_detector->nPixels().X(),
                                    -0.5,
                                    m_detector->nPixels().X() - 0.5,
                                    m_detector->nPixels().Y(),
                                    -0.5,
                                    m_detector->nPixels().Y() - 0.5);
    title = m_detector->getName() + " Cluster seed position (global);x [px];y [px];events";
    clusterSeedPositionGlobal = new TH2F("clusterSeedPositionGlobal",
                                         title.c_str(),
                                         min(500, 5 * m_detector->nPixels().X()),
                                         -m_detector->getSize().X() / 1.5 + dx,
                                         m_detector->getSize().X() / 1.5 + dx,
                                         min(500, 5 * m_detector->nPixels().Y()),
                                         -m_detector->getSize().Y() / 1.5 + dy,
                                         m_detector->getSize().Y() / 1.5 + dy);
    title = m_detector->getName() + " Cluster seed position (local);x [px];y [px];events";
    clusterSeedPositionLocal = new TH2F("clusterSeedPositionLocal",
                                        title.c_str(),
                                        m_detector->nPixels().X(),
                                        -0.5,
                                        m_detector->nPixels().X() - 0.5,
                                        m_detector->nPixels().Y(),
                                        -0.5,
                                        m_detector->nPixels().Y() - 0.5);
    // Cluster shape
    if(flagAnalysisShape)
        initShapeAnalysis();
}

/** Cluster shape - charge sharing in cluster or window
 *
 * Type:
 * - pixel in local index w.r.t. seed
 * - pixel sorted by decreasing order
 * - accumulated value of highest N pixels
 *
 * Variables:
 * - charge, charge ratio, SNR
 **/
void ClusteringAnalog::initShapeAnalysis() {
    std::string title;

    // Cluster size by seeding
    title = m_detector->getName() + " Cluster size - Number of pixels matched seeding criteria;cluster size;# clusters";
    clusterShape_SeedCut = new TH1F("clusterShape_SeedCut", title.c_str(), 100, -0.5, 99.5);

    // Charge
    title =
        m_detector->getName() + " Cluster charge distribution in pixel (local index);pixel no.;pixel charge (ADCu);#cluster";
    clusterShape_Charge_LocalIndex =
        new TH2F("clusterShape_Charge_LocalIndex", title.c_str(), 30, -14.5, 15.5, 2000, -4995, 15005);
    title = m_detector->getName() +
            " Cluster charge in pixel (sorted order by decreasing charge);pixel no.;pixel charge (ADCu);#cluster";
    clusterShape_Charge_SortedIndex =
        new TH2F("clusterShape_Charge_SortedIndex", title.c_str(), 30, 0.5, 30.5, 2000, -4995, 15005);
    title = m_detector->getName() + " Accumulated charge of highest N pixels;pixel no.;pixel charge (ADCu);#cluster";
    clusterShape_Charge_Accumulated =
        new TH2F("clusterShape_Charge_Accumulated", title.c_str(), 30, 0.5, 30.5, 2000, -4995, 15005);

    // Charge Ratio
    title = m_detector->getName() + " Charge ratio over cluster in pixel (local index);pixel no.;Charge ratio;#cluster";
    clusterShape_ChargeRatio_LocalIndex =
        new TH2F("clusterShape_ChargeRatio_LocalIndex", title.c_str(), 30, -14.5, 15.5, 2000, -0.4995, 1.5005);
    title = m_detector->getName() +
            " Charge ratio over cluster in pixel (sorted order by decreasing charge);pixel no.;Charge ratio;#cluster";
    clusterShape_ChargeRatio_SortedIndex =
        new TH2F("clusterShape_ChargeRatio_SortedIndex", title.c_str(), 30, 0.5, 30.5, 2000, -0.4995, 1.5005);
    title =
        m_detector->getName() + " Accumulated charge ratio of highest N pixels over cluster;pixel no.;Charge ratio;#cluster";
    clusterShape_ChargeRatio_Accumulated =
        new TH2F("clusterShape_ChargeRatio_Accumulated", title.c_str(), 30, 0.5, 30.5, 2000, 0.0005, 2.0005);

    // SNR
    if(flagAnalysisSNR) {
        title = m_detector->getName() + " SNR of pixel (local index);pixel no.;Signal/Noise ratio;#cluster";
        clusterShape_SNR_LocalIndex =
            new TH2F("clusterShape_SNR_LocalIndex", title.c_str(), 30, -14.5, 15.5, 2500, -49.95, 200.05);
        title = m_detector->getName() +
                " SNR of pixel (sorted order by decreasing charge);pixel no.;Signal/Noise ratio;#cluster";
        clusterShape_SNR_SortedIndex =
            new TH2F("clusterShape_SNR_SortedIndex", title.c_str(), 30, 0.5, 30.5, 2500, -49.95, 200.05);
    }
}
// Cluster shape - analysis histograms for charge, charge ratio and SNR
void ClusteringAnalog::fillHistogramsShapeAnalysis(const std::shared_ptr<Cluster>& cluster) {
    auto seed = cluster->getSeedPixel();

    // Sort pixels by decreasing charge
    auto pixelVector = cluster->pixels();
    std::sort(pixelVector.begin(), pixelVector.end(), [](const Pixel* pa, const Pixel* pb) {
        return (pa->charge() > pb->charge());
    });
    // Loop all pixels in cluster
    double chargeSum = 0.;
    double ratioSum = 0.;
    int counter = 0;
    int counterSeed = 0;
    for(auto px : pixelVector) {
        // number of added pixel ordered by decreasing charge
        counter++;
        double chargePixel = px->charge();
        chargeSum += chargePixel;
        double ratioPixel = chargePixel / cluster->charge();
        ratioSum += ratioPixel;
        // "Seed" count - mimic binary output
        if(isAboveSeedThreshold(px))
            counterSeed++;
        // Define index in seeding window
        int index = (windowSize * 2 + 1) * (px->row() - seed->row()) + px->column() - seed->column();
        // Charge
        clusterShape_Charge_LocalIndex->Fill(index, chargePixel);
        clusterShape_Charge_SortedIndex->Fill(counter, chargePixel);
        clusterShape_Charge_Accumulated->Fill(counter, chargeSum);
        // Ratio
        clusterShape_ChargeRatio_LocalIndex->Fill(index, ratioPixel);
        clusterShape_ChargeRatio_SortedIndex->Fill(counter, ratioPixel);
        clusterShape_ChargeRatio_Accumulated->Fill(counter, ratioSum);
        // SNR
        if(flagAnalysisSNR) {
            double snrPixel = SNR(px);
            clusterShape_SNR_LocalIndex->Fill(index, snrPixel);
            clusterShape_SNR_SortedIndex->Fill(counter, snrPixel);
        }
    }
    clusterShape_SeedCut->Fill(counterSeed);
}

// Signal/Noise ratio
// return charge, if calibration file is not available.
float ClusteringAnalog::SNR(const Pixel* px) {
    if(!isCalibrated) {
        LOG_ONCE(WARNING) << "Calibration file NOT initialized - return raw charge of (" << px->column() << "," << px->row()
                          << ")";
        return float(px->charge());
    }
    double pNoise = noisemap[static_cast<size_t>(px->column())][static_cast<size_t>(px->row())];
    if(pNoise > 0.)
        return float(px->charge() / pNoise);
    else {
        LOG_ONCE(WARNING) << "Invalid noise value <" << pNoise << "> - return raw charge of (" << px->column() << ","
                          << px->row() << ")";
        return float(px->charge());
    }
}

bool ClusteringAnalog::isAboveSeedThreshold(const Pixel* px) {
    switch(thresholdType) {
    case ThresholdType::SNR:
        return (SNR(px) > thresholdSeedSNR);
        break;
    case ThresholdType::MIX:
        return (SNR(px) > thresholdSeedSNR && px->charge() > thresholdSeed);
        break;
    case ThresholdType::FIX:
    default:
        return (px->charge() > thresholdSeed);
        break;
    }
}

bool ClusteringAnalog::isAboveNeighborThreshold(const Pixel* px) {
    switch(thresholdType) {
    case ThresholdType::SNR:
        return (SNR(px) > thresholdNeighborSNR);
        break;
    case ThresholdType::MIX:
        return (SNR(px) > thresholdNeighborSNR && px->charge() > thresholdNeighbor);
        break;
    case ThresholdType::FIX:
    default:
        return (px->charge() > thresholdNeighbor);
        break;
    }
}

bool ClusteringAnalog::isAboveIterationThreshold(const Pixel* px) {
    switch(thresholdType) {
    case ThresholdType::SNR:
        return (SNR(px) > thresholdIterationSNR);
        break;
    case ThresholdType::MIX:
        return (SNR(px) > thresholdIterationSNR && px->charge() > thresholdIteration);
        break;
    case ThresholdType::FIX:
    default:
        return (px->charge() > thresholdIteration);
        break;
    }
}

// Fill analog clusters - charge, SNR, shape and correlations
void ClusteringAnalog::fillHistograms(const std::shared_ptr<Cluster>& cluster, double chargeSum = 0.) {
    auto seed = cluster->getSeedPixel();
    // Basic information
    clusterSize->Fill(static_cast<double>(cluster->size()));
    clusterCharge->Fill(cluster->charge());
    clusterSeedCharge->Fill(seed->charge());
    clusterNxNCharge->Fill(chargeSum);

    double neighborsChargeSum = cluster->charge() - seed->charge();
    clusterNeighborsChargeSum->Fill(neighborsChargeSum);
    clusterCharge_SeedvsNeighborsSum->Fill(seed->charge(), neighborsChargeSum);
    clusterCharge_SeedvsCluster->Fill(seed->charge(), cluster->charge());

    auto seedSNR = 0.;
    if(flagAnalysisSNR) {
        seedSNR = SNR(seed);
        clusterSeedSNR->Fill(seedSNR);
        clusterSeedSNRvsClusterCharge->Fill(seedSNR, cluster->charge());
    }

    // Central seeds
    if(m_detector->getNeighbors(seed->column(), seed->row(), static_cast<size_t>(windowSize), true).size() ==
       neighborsSizeCentral) {
        clusterSizeCentral->Fill(static_cast<double>(cluster->size()));
        clusterChargeCentral->Fill(cluster->charge());
        clusterSeedChargeCentral->Fill(seed->charge());
        clusterNxNChargeCentral->Fill(chargeSum);
    }

    // Fill position
    clusterPositionGlobal->Fill(cluster->global().x(), cluster->global().y());
    clusterPositionLocal->Fill(cluster->column(), cluster->row());

    auto seedLocal = m_detector->getLocalPosition(seed->column(), seed->row());
    auto seedGlobal = m_detector->localToGlobal(seedLocal);
    clusterSeedPositionGlobal->Fill(seedGlobal.x(), seedGlobal.y());
    clusterSeedPositionLocal->Fill(seed->column(), seed->row());

    for(auto px : cluster->pixels()) {
        if(px == seed)
            continue;
        clusterNeighborsCharge->Fill(px->charge());
        clusterCharge_SeedvsNeighbors->Fill(seed->charge(), px->charge());
        if(flagAnalysisSNR) {
            clusterNeighborsSNR->Fill(SNR(px));
            clusterSNR_SeedvsNeighbors->Fill(seedSNR, SNR(px));
        }
    }

    clusterTimes->Fill(static_cast<double>(Units::convert(cluster->timestamp(), "ns")));

    // Cluster shape
    if(flagAnalysisShape)
        fillHistogramsShapeAnalysis(cluster);
}

// Cluster criteria
// - ROI rejection
// - Edge detection for sumNxN
// - Charge threshold
bool ClusteringAnalog::acceptCluster(const std::shared_ptr<Cluster>& cluster) {
    // check if the cluster is within ROI
    if(rejectByROI && !m_detector->isWithinROI(cluster.get())) {
        LOG(DEBUG) << "Rejecting cluster at (" << cluster->column() << "," << cluster->row() << ") outside of "
                   << m_detector->getName() << " ROI";
        hCutHisto->Fill(RejectionType::kOutsideROI);
        return false;
    }
    // Edge detection of ROI and masked pixels
    int cluster_size = static_cast<int>(cluster->size());
    if(estimationMethod == EstimationMethod::WINDOW && cluster_size < (windowSize * 2 + 1) * (windowSize * 2 + 1)) {
        LOG(DEBUG) << "Rejecting incomplete cluster at (" << cluster->column() << "," << cluster->row()
                   << ") with size = " << cluster->size() << " and window requirement "
                   << (windowSize * 2 + 1) * (windowSize * 2 + 1);
        hCutHisto->Fill(RejectionType::kIncompleteEdgeNxN);
        return false;
    }
    // Reject noise by cluster charge
    if(cluster->charge() < thresholdClusterCharge) {
        LOG(DEBUG) << "Rejecting small cluster at (" << cluster->column() << "," << cluster->row() << ") with charge "
                   << cluster->charge() << " lower than threshold " << thresholdClusterCharge;
        hCutHisto->Fill(RejectionType::kLowCharge);
        return false;
    }
    return true;
}

StatusCode ClusteringAnalog::run(const std::shared_ptr<Clipboard>& clipboard) {

    // Get the pixels
    auto pixels = clipboard->getData<Pixel>(m_detector->getName());
    if(pixels.empty()) {
        LOG(DEBUG) << "Detector " << m_detector->getName() << " does not have any pixels on the clipboard";
        return StatusCode::Success;
    }

    // Make the cluster container and the maps for clustering
    ClusterVector deviceClusters;
    map<std::shared_ptr<Pixel>, bool> used;
    map<int, map<int, std::shared_ptr<Pixel>>> hitmap;

    // Seeding
    PixelVector seedCandidates;
    for(auto pixel : pixels) {
        // Pre-fill the hitmap with pixels
        hitmap[pixel->column()][pixel->row()] = pixel;
        // Select seeds by threshold
        if(isAboveSeedThreshold(pixel.get())) {
            switch(seedingMethod) {
            case SeedingMethod::MAX:
                if(seedCandidates.empty() || pixel->charge() > seedCandidates[0]->charge()) {
                    seedCandidates.clear();
                    seedCandidates.push_back(pixel);
                }
                break;
            case SeedingMethod::MULTI:
            default:
                seedCandidates.push_back(pixel);
                break;
            }
        }
    } // Loop pixels pre-clustering

    // Reconstruct clusters from seeds
    // sorted seedCandidates ordered by decreasing charge
    std::sort(
        seedCandidates.begin(), seedCandidates.end(), [](const std::shared_ptr<Pixel> pa, const std::shared_ptr<Pixel> pb) {
            return (pa->charge() > pb->charge());
        });
    for(auto seed : seedCandidates) {
        if(used[seed])
            continue;

        LOG(DEBUG) << "== Making cluster";

        // New seed/pixel => new cluster
        auto cluster = std::make_shared<Cluster>();
        cluster->addPixel(&*seed);
        used[seed] = true;

        if(useTriggerTimestamp) {
            if(!clipboard->getEvent()->triggerList().empty()) {
                double trigger_ts = clipboard->getEvent()->triggerList().begin()->second;
                LOG(DEBUG) << "Using trigger timestamp " << Units::display(trigger_ts, "us") << " as cluster timestamp.";
                cluster->setTimestamp(trigger_ts);
            } else {
                LOG(WARNING) << "No trigger available. Use pixel timestamp " << Units::display(seed->timestamp(), "us")
                             << " as cluster timestamp.";
                cluster->setTimestamp(seed->timestamp());
            }
        } else {
            // assign pixel timestamp
            LOG(DEBUG) << "Pixel has timestamp " << Units::display(seed->timestamp(), "us")
                       << ", set as cluster timestamp. ";
            cluster->setTimestamp(seed->timestamp());
        }

        // Search neighbors around seed
        PixelVector neighbors;
        // for cluster
        double chargeSum = seed->charge();
        double colChargeWeighted = double(seed->column()) * seed->charge();
        double rowChargeWeighted = double(seed->row()) * seed->charge();

        // Nested lambda function for recursion
        std::function<void(std::shared_ptr<Pixel>&)> searchNeighbors;
        searchNeighbors = [&](std::shared_ptr<Pixel>& pxCenter) -> void {
            for(auto neighborCandidate : m_detector->getNeighbors(pxCenter, 1, includeCorners)) {
                // Traverse adjacent pixels around seed (central pixel)
                auto pixel = hitmap[neighborCandidate.first][neighborCandidate.second];
                // No pixel or already be accepted in cluster
                if(!pixel || used[pixel])
                    continue;
                // Criteria for neighboring pixels
                if(isAboveNeighborThreshold(pixel.get())) {
                    cluster->addPixel(&*pixel);
                    used[pixel] = true;
                    neighbors.push_back(pixel);
                    LOG(DEBUG) << "- pixel col, row, charge : " << pixel->column() << "," << pixel->row() << ","
                               << pixel->charge();
                    // Position by centre-of-gravity
                    chargeSum += pixel->charge();
                    colChargeWeighted += double(pixel->column()) * pixel->charge();
                    rowChargeWeighted += double(pixel->row()) * pixel->charge();
                    // Recursive searching
                    if(isAboveIterationThreshold(pixel.get())) {
                        searchNeighbors(pixel);
                    } // Neighbors found over threshold
                }
            } // Loop neighbors
        };
        searchNeighbors(seed);

        // Window estimation
        double chargeTotalWindow = 0.;
        double colWeightTotalWindow = 0.;
        double rowWeightTotalWindow = 0.;
        // cluster as the whole window
        auto clusterWindow = std::make_shared<Cluster>();
        // Traverse all pixels in the NxN window around seed
        for(auto neighbor : m_detector->getNeighbors(seed, static_cast<size_t>(windowSize), true)) {
            auto pixel = hitmap[neighbor.first][neighbor.second];
            if(!pixel)
                continue;
            clusterWindow->addPixel(&*pixel);
            chargeTotalWindow += pixel->charge();
            colWeightTotalWindow += static_cast<double>(pixel->column()) * pixel->charge();
            rowWeightTotalWindow += static_cast<double>(pixel->row()) * pixel->charge();
        }
        // End - window estimation

        // Finalize the cluster and save it
        LOG(DEBUG) << "- cluster has " << cluster->pixels().size() << " pixels";

        switch(estimationMethod) {
        case EstimationMethod::SEED:
            cluster->setRow(seed->row());
            cluster->setColumn(seed->column());
            cluster->setCharge(seed->charge());
            break;
        case EstimationMethod::WINDOW:
            cluster = clusterWindow;
            cluster->setCharge(chargeTotalWindow);
            if(chargeTotalWindow > 0) {
                cluster->setRow(rowWeightTotalWindow / chargeTotalWindow);
                cluster->setColumn(colWeightTotalWindow / chargeTotalWindow);
            } else {
                // From low signal and noise fluctuation in unfired pixels
                // stats. recorded as RejectionType::kLowCharge
                LOG(DEBUG) << "Zero charge found in clustering (sumNxN method) - "
                           << " cluster charge = " << chargeTotalWindow << ", seed charge = " << seed->charge() << endl;
                cluster->setRow(seed->row());
                cluster->setColumn(seed->column());
            }
            break;
        case EstimationMethod::BINARY:
            rowChargeWeighted = 0.;
            colChargeWeighted = 0.;
            for(auto px : cluster->pixels()) {
                rowChargeWeighted += px->row();
                colChargeWeighted += px->column();
            }
            cluster->setRow(rowChargeWeighted / static_cast<double>(cluster->size()));
            cluster->setColumn(colChargeWeighted / static_cast<double>(cluster->size()));
            cluster->setCharge(chargeSum);
            break;
        case EstimationMethod::CLUSTER:
        default:
            assert(chargeSum > 0.);
            cluster->setRow(rowChargeWeighted / chargeSum);
            cluster->setColumn(colChargeWeighted / chargeSum);
            cluster->setCharge(chargeSum);
            break;
        }

        // Cluster criteria
        if(!acceptCluster(cluster))
            continue;

        LOG(DEBUG) << "- cluster col, row: " << cluster->column() << "," << cluster->row() << " with charge "
                   << cluster->charge();

        // Set uncertainty on position from intrinsic detector spatial resolution:
        cluster->setError(m_detector->getSpatialResolution());

        // Create object with local cluster position
        auto positionLocal = m_detector->getLocalPosition(cluster->column(), cluster->row());
        // Calculate global cluster position
        auto positionGlobal = m_detector->localToGlobal(positionLocal);

        cluster->setDetectorID(pixels.front()->detectorID());
        cluster->setClusterCentre(positionGlobal);
        cluster->setClusterCentreLocal(positionLocal);

        deviceClusters.push_back(cluster);

        LOG(DEBUG) << m_detector->getName() << " - cluster local: (" << cluster->column() << "," << cluster->row()
                   << ") - cluster global: " << cluster->global();

        // Output
        fillHistograms(cluster, chargeTotalWindow);
    } // Loop - seedCandidates

    clipboard->putData(deviceClusters, m_detector->getName());
    LOG(DEBUG) << "Put " << deviceClusters.size() << " clusters on the clipboard for detector " << m_detector->getName()
               << ". From " << pixels.size() << " pixels";

    // Return value telling analysis to keep running
    return StatusCode::Success;
}
