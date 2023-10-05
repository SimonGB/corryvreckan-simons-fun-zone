/**
 * @file
 * @brief Definition of module TrackingMultiplet
 *
 * @copyright Copyright (c) 2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#ifndef CORRYVRECKAN_TRACKINGMULTIPLET_H
#define CORRYVRECKAN_TRACKINGMULTIPLET_H 1

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>

#include "core/module/Module.hpp"
#include "objects/Cluster.hpp"
#include "objects/Multiplet.hpp"
#include "objects/Pixel.hpp"
#include "objects/Track.hpp"
#include "tools/kdtree.h"

namespace corryvreckan {
    /** @ingroup Modules
     */

    // enum to differentiate between up- and downstream arm in functions
    enum streams { upstream, downstream };
    // enum to differentiate between all and chosen tracks in functions
    enum selection { all, chosen };

    class TrackingMultiplet : public Module {

    public:
        // Constructors and destructors
        TrackingMultiplet(Configuration& config, std::vector<std::shared_ptr<Detector>> detectors);

        // Init, run and finalise functions
        void initialize();
        StatusCode run(const std::shared_ptr<Clipboard>& clipboard);

        /**
         * @brief Find tracklets for upstream or downstream tracklets
         */
        TrackVector find_multiplet_tracklets(const streams& stream,
                                             std::map<std::shared_ptr<Detector>, KDTree<Cluster>>& cluster_trees,
                                             std::shared_ptr<Detector> reference_first,
                                             std::shared_ptr<Detector> reference_last);

        /**
         * @brief Fill histograms for upstream or downstream tracklets
         */
        void fill_tracklet_histograms(const streams& stream, const selection& selected, TrackVector tracklets);

    private:
        // Configuration members
        std::map<std::shared_ptr<Detector>, double> time_cuts_;
        std::map<std::shared_ptr<Detector>, XYVector> spatial_cuts_;

        std::vector<std::shared_ptr<Detector>> m_upstream_detectors;
        std::vector<std::shared_ptr<Detector>> m_downstream_detectors;
        std::vector<std::shared_ptr<Detector>> m_require_detectors;

        double scatterer_position_;
        double scatterer_matching_cut_;
        double isolation_cut_;
        double momentum_;
        double beta_;
        int charge_;
        size_t min_hits_upstream_;
        size_t min_hits_downstream_;
        bool refit_gbl_{};
        bool unique_cluster_usage_{};

        // track model for up/downstream fit
        std::string track_model_;
        std::string timestamp_from_;
        std::vector<std::string> require_detectors_;
        std::vector<std::string> exclude_from_seed_;

        // Member histograms
        std::map<std::string, TH1F*> trackletMultiplicity;
        std::map<std::string, TH1F*> clustersPerTracklet;

        std::map<std::string, TH1F*> trackletAngleX;
        std::map<std::string, TH1F*> trackletAngleY;
        std::map<std::string, TH1F*> trackletPositionAtScattererX;
        std::map<std::string, TH1F*> trackletPositionAtScattererY;

        std::map<std::string, TH1F*> residualsX_local;
        std::map<std::string, TH1F*> residualsY_local;
        std::map<std::string, TH1F*> residualsX_global;
        std::map<std::string, TH1F*> residualsY_global;

        TH1F* multipletMultiplicity;
        TH1F* trackChi2;
        TH1F* trackChi2ndof;
        TH1F* trackChi2_refit;
        TH1F* trackChi2ndof_refit;

        TH1F* matchingDistanceAtScattererX;
        TH1F* matchingDistanceAtScattererY;

        TH1F* multipletOffsetAtScattererX;
        TH1F* multipletOffsetAtScattererY;

        TH1F* multipletKinkAtScattererX;
        TH1F* multipletKinkAtScattererY;

        // Function to calculate the weighted average timestamp from the clusters of a track
        double calculate_average_timestamp(const Track* track);
        // Function to refit the multiplet tracks at the end, using GBL
        TrackVector refit(MultipletVector multiplets);

        bool duplicated_hit(const Track* a, const Track* b);
        template <class T> T remove_duplicate(T tracks);
    };

} // namespace corryvreckan

#endif // CORRYVRECKAN_TRACKINGMULTIPLET_H
