/**
 * @file
 * @brief Definition of module ClusteringSeed
 *
 * @copyright Copyright (c) 2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */
//test
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>
#include "core/module/Module.hpp"
#include "objects/Cluster.hpp"
#include "objects/Pixel.hpp"
#include "objects/Track.hpp"

namespace corryvreckan {
    /** @ingroup Modules
     * @brief Module to do function
     *
     * More detailed explanation of module
     */
    class ClusteringSeed : public Module {

    public:
        /**
         * @brief Constructor for this unique module
         * @param config Configuration object for this module as retrieved from the steering file
         * @param detector Pointer to the detector for this module instance
         */
        ClusteringSeed(Configuration& config, std::shared_ptr<Detector> detector);
        ~ClusteringSeed() {}

        /**
         * @brief [Initialise this module]
         */
        void initialize() override;

        /**
         * @brief [Run the function of this module]
         */
        StatusCode runOLD(const std::shared_ptr<Clipboard>& clipboard);
        StatusCode run(const std::shared_ptr<Clipboard>& clipboard) override;

        /**
         * @brief [Finalise module]
         */
        void finalize(const std::shared_ptr<ReadonlyClipboard>& clipboard) override;

    private:
        int m_eventNumber;
        double m_neighbourThreshold;
        double m_seedThreshold;
        int m_lower_channel;
        int m_upper_channel;
        bool m_calculate_crosstalk;

        double QLeftOverQSeed = 0;
        double QRightOverQSeed = 0;
        double QLeftBOOverQSeed = 0;
        double QRightBOOverQSeed = 0;
        int debugVar = 0;

        std::shared_ptr<Detector> m_detector;

        TH1F* clusterSize;
        TH1F* clusterSeedCharge;
        TH1F* clusterCharge;
        TH1F* clusterPosition;
        TH1F* etaDistribution;
        TH1F* leftSignal;
        TH1F* rightSignal;
        TH1F* leftBOSignal;
        TH1F* rightBOSignal;

        TH1F* DEBUG_event_SNR;
        TH1F* DEBUG_event_charge;
    };

} // namespace corryvreckan
