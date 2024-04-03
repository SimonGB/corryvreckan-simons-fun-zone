/**
 * @file
 * @brief Definition of module AnalysisTracks
 *
 * @copyright Copyright (c) 2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */
#ifndef ANALYSISTRACKS_H
#define ANALYSISTRACKS_H 1
#include <TCanvas.h>
#include <TDirectory.h>
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
    class AnalysisTracks : public Module {

    public:
        /**
         * @brief Constructor for this unique module
         * @param config Configuration object for this module as retrieved from the steering file
         * @param detectors Vector of pointers to the detectors
         */
        AnalysisTracks(Configuration& config, std::vector<std::shared_ptr<Detector>> detectors);
        ~AnalysisTracks() {}
        /**
         * @brief [Initialise this module]
         */
        void initialize() override;

        /**
         * @brief [Run the function of this module]
         */
        StatusCode run(const std::shared_ptr<Clipboard>& clipboard) override;

    private:
        std::map<std::string, TH1F*> _distance_between_tracks_{};
        std::map<std::string, TH1F*> _tracks_per_hit_{};
        std::map<std::string, TH2F*> clusters_vs_tracks_{};
    };

} // namespace corryvreckan
#endif // ANALYSISTRACKS_H
