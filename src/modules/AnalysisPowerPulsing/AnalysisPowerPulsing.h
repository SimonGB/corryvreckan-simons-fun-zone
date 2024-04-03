/**
 * @file
 * @brief Definition of module AnalysisPowerPulsing
 *
 * @copyright Copyright (c) 2017-2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#ifndef AnalysisPowerPulsing_H
#define AnalysisPowerPulsing_H 1

#include <iostream>

#include "core/module/Module.hpp"

#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"

namespace corryvreckan {
    /** @ingroup Modules
     */
    class AnalysisPowerPulsing : public Module {

    public:
        /**
         * @brief Constructor for this unique module
         * @param config Configuration object for this module as retrieved from the steering file
         * @param detector Pointer to the detector for this module instance
         */
        AnalysisPowerPulsing(Configuration& config, std::shared_ptr<Detector> detector);

        // Functions
        void initialize() override;
        StatusCode run(const std::shared_ptr<Clipboard>& clipboard) override;

    private:
        std::shared_ptr<Detector> m_detector;

        TH1D* openSignalVersusTime;
        TH1D* closeSignalVersusTime;
        TH1D* clustersVersusTime;
        TH1D* clustersVersusPowerOnTime;
        TH2D* timeSincePowerOnVersusClusterTime;
        TH1D* minTimeSincePowerOnVersusClusterTime;
        TH1D* diffMinTimeSincePowerOnVersusClusterTime;

        long long int m_powerOnTime;
        long long int m_powerOffTime;
        double m_shutterOpenTime;
        double m_shutterCloseTime;

        std::vector<double> v_minTime;
    };

} // namespace corryvreckan
#endif // AnalysisPowerPulsing_H
