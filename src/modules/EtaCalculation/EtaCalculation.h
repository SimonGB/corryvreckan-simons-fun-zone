/**
 * @file
 * @brief Definition of module EtaCalculation
 *
 * @copyright Copyright (c) 2017-2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#ifndef EtaCalculation_H
#define EtaCalculation_H 1

#include <iostream>

#include <TH2F.h>
#include <TProfile.h>

#include "core/module/Module.hpp"
#include "objects/Cluster.hpp"
#include "objects/Track.hpp"

namespace corryvreckan {
    /** @ingroup Modules
     */
    class EtaCalculation : public Module {

    public:
        // Constructors and destructors
        EtaCalculation(Configuration& config, std::shared_ptr<Detector> detector);

        // Functions
        void initialize() override;
        StatusCode run(const std::shared_ptr<Clipboard>& clipboard) override;
        void finalize(const std::shared_ptr<ReadonlyClipboard>& clipboard) override;

    private:
        void calculate_eta(const Track* track, const Cluster* cluster);
        std::string fit(const std::string& fname, double pitch, TProfile* profile) const;

        std::shared_ptr<Detector> detector_;
        double chi2ndof_cut_;

        // Histograms
        TH2F* etaDistributionX_;
        TH2F* etaDistributionY_;
        TProfile* etaDistributionXprofile_;
        TProfile* etaDistributionYprofile_;
    };
} // namespace corryvreckan
#endif // EtaCalculation_H
