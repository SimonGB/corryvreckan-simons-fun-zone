/**
 * @file
 * @brief Definition of module ClusteringAnalog
 *
 * @copyright Copyright (c) 2017-2022 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#ifndef ClusteringAnalog_H
#define ClusteringAnalog_H 1

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>
#include "core/module/Module.hpp"
#include "objects/Cluster.hpp"

namespace corryvreckan {
    /** @ingroup Modules
     */
    class ClusteringAnalog : public Module {

    public:
        // Constructors and destructors
        ClusteringAnalog(Configuration& config, std::shared_ptr<Detector> detector);
        ~ClusteringAnalog() {}

        // Functions
        void initialize() override;
        void initShapeAnalysis();
        StatusCode run(const std::shared_ptr<Clipboard>& clipboard) override;

    private:
        // Methods
        float SNR(const Pixel* px); // Signal/Noise ratio
        bool isAboveSeedThreshold(const Pixel* px);
        bool isAboveNeighborThreshold(const Pixel* px);
        bool isAboveIterationThreshold(const Pixel* px);
        bool acceptCluster(const std::shared_ptr<Cluster>& cluster);
        void fillHistogramsShapeAnalysis(const std::shared_ptr<Cluster>& cluster);
        void fillHistograms(const std::shared_ptr<Cluster>& cluster, double chargeTotal);
        bool readCalibrationFileROOT(const std::filesystem::path fileName);

        std::shared_ptr<Detector> m_detector;

        // Cluster histograms
        TH1F* clusterSize;
        TH1F* clusterSeedCharge;
        TH1F* clusterCharge;
        TH1F* clusterNxNCharge;
        TH1F* clusterNeighborsCharge;
        TH1F* clusterNeighborsChargeSum;

        TH1F* clusterSeedSNR;
        TH1F* clusterNeighborsSNR;

        TH1F* clusterTimes;

        // Seeding - 2D correlation
        TH2F* clusterCharge_SeedvsNeighbors;
        TH2F* clusterSNR_SeedvsNeighbors;
        TH2F* clusterCharge_SeedvsNeighborsSum;
        TH2F* clusterCharge_SeedvsCluster;
        TH2F* clusterSeedSNRvsClusterCharge;

        // Cluster shape (charge sharing)
        TH1F* clusterShape_SeedCut; // Number of pixels matched seeding criteria
        TH2F* clusterShape_Charge_LocalIndex;
        TH2F* clusterShape_Charge_SortedIndex;
        TH2F* clusterShape_Charge_Accumulated;
        TH2F* clusterShape_ChargeRatio_LocalIndex;
        TH2F* clusterShape_ChargeRatio_SortedIndex;
        TH2F* clusterShape_ChargeRatio_Accumulated;
        TH2F* clusterShape_SNR_LocalIndex;
        TH2F* clusterShape_SNR_SortedIndex;
        // TH2F* clusterShape_SNR_Accumulated; // SNR is meaningless to cluster or more than 1 pixel

        TH1F* clusterSizeCentral;
        TH1F* clusterSeedChargeCentral;
        TH1F* clusterChargeCentral;
        TH1F* clusterNxNChargeCentral;

        TH2F* clusterPositionGlobal;
        TH2F* clusterPositionLocal;

        TH2F* clusterSeedPositionGlobal;
        TH2F* clusterSeedPositionLocal;

        // module parameters
        enum RejectionType {
            kOutsideROI,
            kIncompleteEdgeNxN,
            kLowCharge,
            kAllTypes,
        };

        enum class EstimationMethod {
            SEED,
            CLUSTER,
            WINDOW,
            BINARY,
        } estimationMethod;

        enum class SeedingMethod {
            MAX,
            MULTI,
        } seedingMethod;

        enum class ThresholdType {
            FIX,
            SNR,
            MIX,
        } thresholdType;

        int windowSize; // Cluster matrix to search neighbors
        size_t neighborsSizeCentral;
        bool includeCorners;
        bool rejectByROI;
        // Threshold - raw value (ADC unit, TOT, ...)
        float thresholdSeed;
        float thresholdNeighbor;
        float thresholdIteration;
        float thresholdClusterCharge;
        // Threshold - Signal/Noise ratio (require calibration)
        float thresholdSeedSNR;
        float thresholdNeighborSNR;
        float thresholdIterationSNR;
        // Calibration file
        std::vector<std::vector<double>> noisemap;
        // Configure associated cluster time
        bool useTriggerTimestamp;
        // Analysis functionality
        TH1F* hCutHisto;
        bool flagAnalysisSNR;   // Enable SNR estimation and histograms for analysis
        bool flagAnalysisShape; // Enable analysis for charge sharing
        bool isCalibrated;
    };
} // namespace corryvreckan
#endif // ClusteringAnalog_H
