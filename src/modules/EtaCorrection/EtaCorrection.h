/**
 * @file
 * @brief Definition of module EtaCorrection
 *
 * @copyright Copyright (c) 2017-2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#ifndef EtaCorrection_H
#define EtaCorrection_H 1

#include <TCanvas.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <iostream>

#include "core/module/Module.hpp"
#include "objects/Cluster.hpp"
#include "objects/Pixel.hpp"
#include "objects/Track.hpp"

namespace corryvreckan {
    /** @ingroup Modules
     */
    class EtaCorrection : public Module {

    public:
        // Constructors and destructors
        EtaCorrection(Configuration& config, std::shared_ptr<Detector> detector);
        ~EtaCorrection() {}

        // Functions
        void initialize() override;
        StatusCode run(const std::shared_ptr<Clipboard>& clipboard) override;

    private:
        void applyEta(Cluster* cluster);

        std::shared_ptr<Detector> detector_;
        std::string etaFormulaX_;
        TF1* etaCorrectorX_;
        bool correctX_;
        std::string etaFormulaY_;
        TF1* etaCorrectorY_;
        bool correctY_;

        // Histograms
        TProfile* etaDistributionXprofile_;
        TProfile* etaDistributionYprofile_;
    };
} // namespace corryvreckan
#endif // EtaCorrection_H
