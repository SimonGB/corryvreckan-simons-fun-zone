/**
 * @file
 * @brief Definition of module AlignmentDUTResidual
 *
 * @copyright Copyright (c) 2017-2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include <TCanvas.h>
#include <TFormula.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>

#include "core/module/Module.hpp"
#include "core/utils/ThreadPool.hpp"
#include "objects/Cluster.hpp"
#include "objects/Pixel.hpp"
#include "objects/Track.hpp"

namespace corryvreckan {
    /**
     * @brief Version of std::make_shared that does not delete the pointer
     *
     * This version is needed because some pointers are deleted by ROOT's MINUIT2 internally, but they are stored as
     * std::shared_ptr in this framework.
     */
    template <typename T, typename... Args> static std::shared_ptr<T> make_shared_no_delete(Args... args) {
        return std::shared_ptr<T>(new T(args...), [](T*) {});
    }

    /** @ingroup Modules
     * @brief Module to do function
     *
     * More detailed explanation of module
     */
    class AlignmentDUTResidual : public Module {

    public:
        /**
         * @brief Constructor for this unique module
         * @param config Configuration object for this module as retrieved from the steering file
         * @param detector Pointer to the detector for this module instance
         */
        AlignmentDUTResidual(Configuration& config, std::shared_ptr<Detector> detector);

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
        static void MinimiseResiduals(Int_t& npar, Double_t* grad, Double_t& result, Double_t* par, Int_t flag);
        void SetResidualsFunctions();

        std::shared_ptr<Detector> m_detector;
        int m_discardedtracks{};

        // Global container declarations
        static TrackVector globalTracks;
        static std::shared_ptr<Detector> globalDetector;
        static ThreadPool* thread_pool;

        unsigned int m_workers;
        size_t nIterations;
        bool m_pruneTracks;
        bool m_alignPosition;
        bool m_alignOrientation;
        std::string m_alignPosition_axes;
        std::string m_alignOrientation_axes;
        size_t m_maxAssocClusters;
        double m_maxTrackChi2;
        double m_spatial_cut_sensoredge;
        TH1F* residualsXPlot;
        TH1F* residualsYPlot;
        TH1F* residualsRPlot{};
        TH1F* residualsPhiPlot{};

        static std::shared_ptr<TFormula> formula_residual_x;
        static std::shared_ptr<TFormula> formula_residual_y;

        TProfile* profile_dY_X;
        TProfile* profile_dY_Y;
        TProfile* profile_dX_X;
        TProfile* profile_dX_Y;
        TGraph* align_correction_shiftX;
        TGraph* align_correction_shiftY;
        TGraph* align_correction_rot0;
        TGraph* align_correction_rot1;
        TGraph* align_correction_rot2;
    };

} // namespace corryvreckan
