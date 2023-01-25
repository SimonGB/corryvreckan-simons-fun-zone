/**
 * @file
 * @brief Definition of AnalysisTrackIntercept module
 *
 * @copyright Copyright (c) 2022 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#ifndef CORRYVRECKAN_TRACKINTERCEPT_ANALYSIS_H
#define CORRYVRECKAN_TRACKINTERCEPT_ANALYSIS_H

#include <TDirectory.h>
#include <TH2F.h>
#include <vector>
#include "core/module/Module.hpp"

namespace corryvreckan {
    /** @ingroup Modules
     */
    class AnalysisTrackIntercept : public Module {

    public:
        // Constructors and destructors
        AnalysisTrackIntercept(Configuration& config, std::vector<std::shared_ptr<Detector>> detector);
        ~AnalysisTrackIntercept() {}

        // Functions
        void initialize() override;
        StatusCode run(const std::shared_ptr<Clipboard>& clipboard) override;
        void finalize(const std::shared_ptr<ReadonlyClipboard>& clipboard) override;

    private:
        int n_bins_x, n_bins_y;
        double xmin, xmax, ymin, ymax;

        std::vector<double> m_planes_z;
        std::vector<TH2F*> m_intercepts;
    };
} // namespace corryvreckan

#endif // CORRYVRECKAN_TRACKINTERCEPT_ANALYSIS_H
