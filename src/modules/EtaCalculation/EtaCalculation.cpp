/**
 * @file
 * @brief Implementation of module EtaCalculation
 *
 * @copyright Copyright (c) 2017-2022 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "EtaCalculation.h"
#include "objects/Pixel.hpp"

#include <TF1.h>

using namespace corryvreckan;
using namespace std;

EtaCalculation::EtaCalculation(Configuration& config, std::shared_ptr<Detector> detector)
    : Module(config, detector), detector_(detector) {

    config_.setDefault<double>("chi2ndof_cut", 100.);
    config_.setDefault<bool>("calculate_x", true);
    config_.setDefault<bool>("calculate_y", true);

    chi2ndof_cut_ = config_.get<double>("chi2ndof_cut");
    calculate_x_ = config_.get<bool>("calculate_x");
    calculate_y_ = config_.get<bool>("calculate_y");

    if(calculate_x_) {
        config_.setDefault<std::string>("eta_formula_x", "[0] + [1]*x + [2]*x^2 + [3]*x^3 + [4]*x^4 + [5]*x^5");
    }
    if(calculate_y_) {
        config_.setDefault<std::string>("eta_formula_y", "[0] + [1]*x + [2]*x^2 + [3]*x^3 + [4]*x^4 + [5]*x^5");
    }
}

void EtaCalculation::initialize() {

    // Initialise histograms
    auto pitch_x = detector_->getPitch().X();
    auto pitch_y = detector_->getPitch().Y();
    std::string mod_axes_x = "in-2pixel x_{cluster} [mm];in-2pixel x_{track} [mm];";
    std::string mod_axes_y = "in-2pixel y_{cluster} [mm];in-2pixel y_{track} [mm];";

    if(calculate_x_) {
        auto bins_x = std::min(static_cast<int>(Units::convert(pitch_x, "um") * 2), 1000);

        std::string title = "2D #eta distribution X;" + mod_axes_x + "No. entries";
        etaDistributionX_ = new TH2F(
            "etaDistributionX", title.c_str(), bins_x, -pitch_x / 2, pitch_x / 2, bins_x, -pitch_x / 2, pitch_x / 2);

        title = "#eta distribution X;" + mod_axes_x;
        etaDistributionXprofile_ = new TProfile(
            "etaDistributionXprofile", title.c_str(), bins_x, -pitch_x / 2, pitch_x / 2, -pitch_x / 2, pitch_x / 2);
    }

    if(calculate_y_) {
        auto bins_y = std::min(static_cast<int>(Units::convert(pitch_y, "um") * 2), 1000);

        std::string title = "2D #eta distribution Y;" + mod_axes_y + "No. entries";
        etaDistributionY_ = new TH2F(
            "etaDistributionY", title.c_str(), bins_y, -pitch_y / 2, pitch_y / 2, bins_y, -pitch_y / 2, pitch_y / 2);

        title = "#eta distribution Y;" + mod_axes_y;
        etaDistributionYprofile_ = new TProfile(
            "etaDistributionYprofile", title.c_str(), bins_y, -pitch_y / 2, pitch_y / 2, -pitch_y / 2, pitch_y / 2);
    }
}

void EtaCalculation::calculate_eta(const Track* track, const Cluster* cluster) {
    // Ignore single pixel clusters
    if(cluster->size() == 1) {
        return;
    }
    auto localIntercept = detector_->getLocalIntercept(track);

    if(cluster->columnWidth() == 2 && calculate_x_) {
        LOG(DEBUG) << "Calculating correction in local X";
        auto reference_col = 0;
        for(auto& pixel : cluster->pixels()) {
            if(pixel->column() > reference_col) {
                reference_col = pixel->column();
            }
        }
        // Map residual onto range -pitch/2 to pitch/2
        auto reference_X = detector_->getPitch().X() * (reference_col - 0.5 * detector_->nPixels().X());
        auto xmod_cluster = cluster->local().X() - reference_X;
        auto xmod_track = localIntercept.X() - reference_X;
        etaDistributionX_->Fill(xmod_cluster, xmod_track);
        etaDistributionXprofile_->Fill(xmod_cluster, xmod_track);
    }
    if(cluster->rowWidth() == 2 && calculate_y_) {
        LOG(DEBUG) << "Calculating correction in local Y";
        auto reference_row = 0;
        for(auto& pixel : cluster->pixels()) {
            if(pixel->row() > reference_row) {
                reference_row = pixel->row();
            }
        }

        // Map residual onto range -pitch/2 to pitch/2
        auto reference_Y = detector_->getPitch().Y() * (reference_row - 0.5 * detector_->nPixels().Y());
        auto ymod_cluster = cluster->local().Y() - reference_Y;
        auto ymod_track = localIntercept.Y() - reference_Y;

        etaDistributionY_->Fill(ymod_cluster, ymod_track);
        etaDistributionYprofile_->Fill(ymod_cluster, ymod_track);
    }
}

StatusCode EtaCalculation::run(const std::shared_ptr<Clipboard>& clipboard) {

    // Loop over all tracks and look at the associated clusters to plot the eta distribution
    auto tracks = clipboard->getData<Track>();
    for(const auto& track : tracks) {

        // Cut on the chi2/ndof
        if(track->getChi2ndof() > chi2ndof_cut_) {
            continue;
            LOG(DEBUG) << "Skipping track with chi2 = " << track->getChi2ndof() << " which is above cut of "
                       << chi2ndof_cut_;
        }

        // Look at the associated clusters and plot the eta function
        for(auto& dutCluster : track->getAssociatedClusters(detector_->getName())) {
            calculate_eta(track.get(), dutCluster);
        }

        // Do the same for all clusters of the track:
        for(auto& cluster : track->getClusters()) {
            if(cluster->detectorID() != detector_->getName()) {
                continue;
            }
            calculate_eta(track.get(), cluster);
        }
    }

    // Return value telling analysis to keep running
    return StatusCode::Success;
}

std::string EtaCalculation::fit(const std::string& fname, double pitch, TProfile* profile) const {

    auto formula = config_.get<std::string>(fname);
    auto function = new TF1(fname.c_str(), formula.c_str(), -pitch, pitch);
    std::stringstream parameters;

    // Get the eta distribution profiles and fit them to extract the correction parameters
    auto fit_result = profile->Fit(function, "q");
    if(!fit_result) {
        LOG(ERROR) << "Fit for " << fname << " failed!";
        return {};
    }

    // Retrieve fit parameters:
    TF1* fit = profile->GetFunction(fname.c_str());
    if(!fit) {
        LOG(ERROR) << "Could not obtain fit function for " << fname << "!";
        return {};
    }

    for(int i = 0; i < fit->GetNumberFreeParameters(); i++) {
        parameters << " " << fit->GetParameter(i);
    }
    return parameters.str();
}

void EtaCalculation::finalize(const std::shared_ptr<ReadonlyClipboard>&) {

    std::stringstream config;
    if(calculate_x_) {
        LOG(INFO) << "Calculating correction in local X";
        config << "eta_formula_x = \"" << config_.get<std::string>("eta_formula_x") << "\"" << std::endl
               << "eta_constants_x"
               << " = " << fit("eta_formula_x", detector_->getPitch().X(), etaDistributionXprofile_) << std::endl;
    }

    if(calculate_y_) {
        LOG(INFO) << "Calculating correction in local Y";
        config << "eta_formula_y = \"" << config_.get<std::string>("eta_formula_y") << "\"" << std::endl
               << "eta_constants_y"
               << " = " << fit("eta_formula_y", detector_->getPitch().Y(), etaDistributionYprofile_);
    }
    LOG(INFO) << "To apply this correction, place the following in the configuration:" << std::endl
              << "[EtaCorrection]" << std::endl
              << "name = " << detector_->getName() << std::endl
              << config.str();
}
