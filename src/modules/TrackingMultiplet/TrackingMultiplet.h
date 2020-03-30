/**
 * @file
 * @brief Definition of module TrackingMultiplet
 *
 * @copyright Copyright (c) 2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>
#include "core/module/Module.hpp"
#include "objects/Cluster.hpp"
#include "objects/Multiplet.hpp"
#include "objects/Pixel.hpp"
#include "objects/Track.hpp"

namespace corryvreckan {
    /** @ingroup Modules
     * @brief Module to do function
     *
     * More detailed explanation of module
     */
    class TrackingMultiplet : public Module {

    public:
        /**
         * @brief Constructor for this unique module
         * @param config Configuration object for this module as retrieved from the steering file
         * @param detectors Vector of pointers to the detectors
         */
        TrackingMultiplet(Configuration config, std::vector<std::shared_ptr<Detector>> detectors);

        /**
         * @brief [Initialise this module]
         */
        void initialise();

        /**
         * @brief [Run the function of this module]
         */
        StatusCode run(std::shared_ptr<Clipboard> clipboard);

        /**
         * @brief [Finalise module]
         */
        void finalise();

    private:
        double time_cut_reference_;
        std::map<std::shared_ptr<Detector>, double> time_cuts_;
        std::map<std::shared_ptr<Detector>, XYVector> spatial_cuts_;

        double scatterer_position_;
        double scatterer_matching_cut_;

        std::vector<std::string> upstream_detectors_;
        std::vector<std::string> downstream_detectors_;

        TrackVector m_upstreamTracks;
        TrackVector m_downstreamTracks;

        MultipletVector m_multiplets;

        TH1F* upstreamMultiplicity;
        TH1F* downstreamMultiplicity;
        TH1F* multipletMultiplicity;

        TH1F* upstreamAngleX;
        TH1F* upstreamAngleY;
        TH1F* downstreamAngleX;
        TH1F* downstreamAngleY;

        TH1F* upstreamPositionXAtScatterer;
        TH1F* upstreamPositionYAtScatterer;
        TH1F* downstreamPositionXAtScatterer;
        TH1F* downstreamPositionYAtScatterer;

        TH1F* matchingDistanceXAtScatterer;
        TH1F* matchingDistanceYAtScatterer;

        TH1F* multipletOffsetXAtScatterer;
        TH1F* multipletOffsetYAtScatterer;

        TH1F* multipletKinkXAtScatterer;
        TH1F* multipletKinkYAtScatterer;
    };

} // namespace corryvreckan
