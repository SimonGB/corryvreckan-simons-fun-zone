/**
 * @file
 * @brief Definition of module AnalysisFASTPIX
 *
 * @copyright Copyright (c) 2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH2Poly.h>
#include <TProfile2D.h>
#include <TProfile2Poly.h>
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
    class AnalysisFASTPIX : public Module {

    public:
        /**
         * @brief Constructor for this unique module
         * @param config Configuration object for this module as retrieved from the steering file
         * @param detector Pointer to the detector for this module instance
         */
        AnalysisFASTPIX(Configuration& config, std::shared_ptr<Detector> detector);

        /**
         * @brief [Initialise this module]
         */
        void initialize() override;

        /**
         * @brief [Run the function of this module]
         */
        StatusCode run(const std::shared_ptr<Clipboard>& clipboard) override;

        /**
         * @brief [Finalise module]
         */
        void finalize(const std::shared_ptr<ReadonlyClipboard>& clipboard) override;

    private:
        bool inRoi(PositionVector3D<Cartesian3D<double>>);

        template <typename T> Int_t fillTriangle(T* hist, double x, double y, double val = 1);

        TH2F *hitmap, *hitmapIntercept, *hitmapNoIntercept, *hitmapTrigger, *hitmapTimecuts, *hitmapAssoc;
        TH2F* hitmapTriggerAssoc;

        TProfile2D *clusterSizeMap, *clusterChargeMap, *seedChargeMap;
        TProfile2D *clusterSizeMap_inpix, *clusterChargeMap_inpix, *seedChargeMap_inpix;
        TProfile2D* clusterSizeMap_intercept;
        TProfile2Poly *clusterSizeMap_inpix3, *clusterChargeMap_inpix3, *seedChargeMap_inpix3;
        TH2Poly* hitmapTriggerAssoc_inpix3;

        TH2Poly *hitmapTrigger_inpix3, *hitmapTimecuts_inpix3, *hitmapAssoc_inpix3;

        TH1F *clusterSize, *clusterSizeROI;

        TH2F* binningIneff_inpix;

        double chi2_ndof_cut_, time_cut_frameedge_, time_cut_deadtime_, time_cut_trigger_, time_cut_trigger_assoc_;
        bool use_closest_cluster_;

        ROOT::Math::XYVector roi_min, roi_max;
        double roi_margin_;
        int triangle_bins_;

        double last_timestamp = 0;
        double pitch, height;

        std::shared_ptr<Detector> m_detector;
    };

} // namespace corryvreckan
