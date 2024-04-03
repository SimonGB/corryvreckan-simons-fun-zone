/**
 * @file
 * @brief Definition of module Prealignment
 *
 * @copyright Copyright (c) 2018-2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#ifndef PREALIGNMENT_H
#define PREALIGNMENT_H 1

#include <TCanvas.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>
#include "core/module/Module.hpp"
#include "objects/Cluster.hpp"
#include "objects/Pixel.hpp"

namespace corryvreckan {
    enum class PrealignMethod {
        MEAN = 0,
        MAXIMUM,
        MAXIMUM2D,
        GAUSS_FIT,
    };

    /** @ingroup Modules
     */
    class Prealignment : public Module {

    public:
        // Constructors and destructors
        Prealignment(Configuration& config, std::shared_ptr<Detector> detector);
        ~Prealignment() {}

        // Functions
        void initialize() override;
        StatusCode run(const std::shared_ptr<Clipboard>& clipboard) override;
        void finalize(const std::shared_ptr<ReadonlyClipboard>& clipboard) override;

    private:
        std::shared_ptr<Detector> m_detector;

        // Correlation plots
        TH1F* correlationX;
        TH1F* correlationY;
        TH2F* correlationXY;
        TH2F* correlationX2Dlocal;
        TH2F* correlationY2Dlocal;
        TH2F* correlationX2D;
        TH2F* correlationY2D;

        // Parameters which can be set by user
        double max_correlation_rms;
        double damping_factor;
        double timeCut;
        double range_abs;
        PrealignMethod method;
        int fit_range_rel;
        std::vector<std::string> fixed_planes_;
    };
} // namespace corryvreckan
#endif // PREALIGNMENT_H
