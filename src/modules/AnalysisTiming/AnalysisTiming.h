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
#include <optional>
#include <string>

#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile2D.h>

#include "core/module/Module.hpp"

namespace corryvreckan {
    class Cluster;
    class Track;

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

        class InvalidDetectorRoleError : public ConfigurationError {
        public:
            InvalidDetectorRoleError(const Detector* detector) {
                error_message_ = "Role \"" + std::string(magic_enum::enum_name(detector->getRoles())) + "\" of detector \"" +
                                 detector->getName() + "\" is not supported in AnalysisTiming module";
            }
        };

    private:
        // Enum for cut histogram
        enum ETrackSelection : int {
            kAllTrack = 0,
            kPassedChi2Ndf,
            kTimestampOnRef,
            kClusterOnDUT,
            kNSelection,
        };

        // Enum for defining how to extract reference timestamp
        enum class TimestampReferenceType {
            DUT,   // Use associated cluster from different detector
            PLANE, // Use track cluster from different detector
            TRACK, // Use track timestamp
        };

        /**
         * @brief Get closest associated cluster for DUT
         */
        Cluster* get_dut_cluster_associated(const Track* track) const;

        /**
         * @brief Get cluster used in track for DUT
         */
        Cluster* get_dut_cluster_tracking(const Track* track) const;

        // Function pointer type for get_dut_cluster_associated and get_dut_cluster_tracking
        using DutClusterFunc = decltype(&AnalysisTiming::get_dut_cluster_associated);

        /**
         * @brief Get reference timestamp from closest associated cluster
         */
        std::optional<double> get_ref_timestamp_dut(const Track* track) const;

        /**
         * @brief Get reference timestamp from cluster used in track
         */
        std::optional<double> get_ref_timestamp_plane(const Track* track) const;

        /**
         * @brief Get reference timestamp from track timestamp
         */
        std::optional<double> get_ref_timestamp_track(const Track* track) const;

        // Function pointer type for get_ref_timestamp_dut, get_ref_timestamp_plane and get_ref_timestamp_track
        using RefTimestampFunc = decltype(&AnalysisTiming::get_ref_timestamp_dut);

    private:
        // Detector variables
        std::shared_ptr<Detector> detector_;
        std::string detector_name_;
        DutClusterFunc dut_cluster_func_;
        std::string reference_name_;
        RefTimestampFunc ref_timestamp_func_;

        // Cuts
        double chi2_ndof_cut_;

        // Histograms
        TH1D* hTimeResidual_;
        TH2F* hTimeResidualOverTime_;
        TProfile2D* hResidualMeanSensor_;
        TProfile2D* hResidualStdDevSensor_;
        TProfile2D* hResidualMeanInpix_;
        TProfile2D* hResidualStdDevInpix_;
        TH1F* hCutHisto_;

        // Histogram settings
        double time_range_;
        double time_binning_;
        double time_offset_;
        ROOT::Math::XYPoint inpixel_bin_size_;
    };

} // namespace corryvreckan
