/**
 * @file
 * @brief Definition of module AnalysisTiming
 *
 * @copyright Copyright (c) 2023 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include <memory>
#include <string>

#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile2D.h>

#include "core/module/Module.hpp"

namespace corryvreckan {
    /** @ingroup Modules
     * @brief Module to analysise the timing between two detectors
     */
    class AnalysisTiming : public Module {

    public:
        /**
         * @brief Constructor of this module
         * @param config Configuration object for this module as retrieved from the steering file
         * @param detector Pointer to the detector for this module instance
         */
        AnalysisTiming(Configuration& config, std::shared_ptr<Detector> detector);

        /**
         * @brief Initialise function of this module
         */
        void initialize() override;

        /**
         * @brief Run function of this module
         */
        StatusCode run(const std::shared_ptr<Clipboard>& clipboard) override;

        /**
         * @brief Finalise function of this module
         */
        void finalize(const std::shared_ptr<ReadonlyClipboard>& clipboard) override;

        // For track selection cut
        enum ETrackSelection : int {
            kAllTrack = 0,
            kPassedChi2Ndf,
            kClusterOnRef,
            kClusterOnDUT,
            kNSelection,
        };

    private:
        std::shared_ptr<Detector> detector_;

        // detector settings
        std::string reference_name_;
        bool reference_associated_clusters_;
        bool dut_associated_clusters_;

        // cuts
        double chi2_ndof_cut_;

        // histograms
        TH1D* hTimeResidual_;
        TH2F* hTimeResidualOverTime_;
        TProfile2D* hResidualMeanSensor_;
        TProfile2D* hResidualStdDevSensor_;
        TProfile2D* hResidualMeanInpix_;
        TProfile2D* hResidualStdDevInpix_;
        TH1F* hCutHisto_;

        // histogram settings
        double time_range_;
        double time_binning_;
        double time_offset_;
        ROOT::Math::XYPoint inpixel_bin_size_;
    };

} // namespace corryvreckan
