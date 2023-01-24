/**
 * @file
 * @brief Implementation of module EtaCorrection
 *
 * @copyright Copyright (c) 2017-2022 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "EtaCorrection.h"

using namespace corryvreckan;
using namespace std;

EtaCorrection::EtaCorrection(Configuration& config, std::shared_ptr<Detector> detector)
    : Module(config, detector), detector_(detector) {

    config_.setDefault<std::string>("eta_formula_x", "[0] + [1]*x + [2]*x^2 + [3]*x^3 + [4]*x^4 + [5]*x^5");
    config_.setDefault<std::string>("eta_formula_y", "[0] + [1]*x + [2]*x^2 + [3]*x^3 + [4]*x^4 + [5]*x^5");

    etaFormulaX_ = config_.get<std::string>("eta_formula_x");
    etaFormulaY_ = config_.get<std::string>("eta_formula_y");
}

void EtaCorrection::initialize() {

    // Initialise histograms
    auto pitch_x = detector_->getPitch().X();
    auto pitch_y = detector_->getPitch().Y();
    std::string mod_axes_x = "in-2pixel x_{cluster} [mm];in-2pixel x_{corrected} [mm];";
    std::string mod_axes_y = "in-2pixel y_{cluster} [mm];in-2pixel y_{corrected} [mm];";

    std::string title = "#eta distribution X;" + mod_axes_x;
    etaDistributionXprofile_ = new TProfile("etaDistributionXprofile",
                                            title.c_str(),
                                            static_cast<int>(Units::convert(pitch_x, "um") * 2),
                                            -pitch_x / 2,
                                            pitch_x / 2,
                                            -pitch_x / 2,
                                            pitch_x / 2);

    title = "#eta distribution Y;" + mod_axes_y;
    etaDistributionYprofile_ = new TProfile("etaDistributionYprofile",
                                            title.c_str(),
                                            static_cast<int>(Units::convert(pitch_y, "um") * 2),
                                            -pitch_y / 2,
                                            pitch_y / 2,
                                            -pitch_y / 2,
                                            pitch_y / 2);

    // Get info from configuration:
    std::vector<double> etaConstantsX_ = config_.getArray<double>("eta_constants_x", {});
    std::vector<double> etaConstantsY_ = config_.getArray<double>("eta_constants_y", {});
    if(!etaConstantsX_.empty() || !etaConstantsY_.empty()) {
        LOG(INFO) << "Found Eta correction factors for detector \"" << detector_->getName()
                  << "\": " << (etaConstantsX_.empty() ? "" : "X ") << (etaConstantsY_.empty() ? "" : "Y ");
    }

    if(!etaConstantsX_.empty()) {
        correctX_ = true;
        etaCorrectorX_ =
            new TF1("etaCorrectorX", etaFormulaX_.c_str(), -1 * detector_->getPitch().X(), detector_->getPitch().X());
        for(size_t x = 0; x < etaConstantsX_.size(); x++) {
            etaCorrectorX_->SetParameter(static_cast<int>(x), etaConstantsX_[x]);
        }
        // add the loaded function to the plot
        etaDistributionXprofile_->GetListOfFunctions()->Add(etaCorrectorX_);
    } else {
        correctX_ = false;
    }

    if(!etaConstantsY_.empty()) {
        correctY_ = true;
        etaCorrectorY_ =
            new TF1("etaCorrectorY", etaFormulaY_.c_str(), -1 * detector_->getPitch().Y(), detector_->getPitch().Y());
        for(size_t y = 0; y < etaConstantsY_.size(); y++) {
            etaCorrectorY_->SetParameter(static_cast<int>(y), etaConstantsY_[y]);
        }
        // add the loaded function to the plot
        etaDistributionYprofile_->GetListOfFunctions()->Add(etaCorrectorY_);
    } else {
        correctY_ = false;
    }
}

void EtaCorrection::applyEta(Cluster* cluster) {
    // Ignore single pixel clusters
    if(cluster->size() == 1) {
        return;
    }
    double newX = cluster->local().x();
    double newY = cluster->local().y();

    if(cluster->columnWidth() == 2) {
        if(correctX_) {
            auto reference_col = 0;
            for(auto& pixel : cluster->pixels()) {
                if(pixel->column() > reference_col) {
                    reference_col = pixel->column();
                }
            }
            auto reference_X = detector_->getPitch().X() * (reference_col - 0.5 * detector_->nPixels().X());
            auto xmod_cluster = cluster->local().X() - reference_X;
            newX = etaCorrectorX_->Eval(xmod_cluster);
            etaDistributionXprofile_->Fill(xmod_cluster, newX);
            newX += reference_X;
        }
    }

    if(cluster->rowWidth() == 2) {
        if(correctY_) {
            auto reference_row = 0;
            for(auto& pixel : cluster->pixels()) {
                if(pixel->row() > reference_row) {
                    reference_row = pixel->row();
                }
            }
            auto reference_Y = detector_->getPitch().Y() * (reference_row - 0.5 * detector_->nPixels().Y());
            auto ymod_cluster = cluster->local().Y() - reference_Y;
            newY = etaCorrectorY_->Eval(ymod_cluster);
            etaDistributionYprofile_->Fill(ymod_cluster, newY);
            newY += reference_Y;
        }
    }

    PositionVector3D<Cartesian3D<double>> positionLocal(newX, newY, 0);
    PositionVector3D<Cartesian3D<double>> positionGlobal = detector_->localToGlobal(positionLocal);
    cluster->setClusterCentre(positionGlobal);
    cluster->setClusterCentreLocal(positionLocal);
}

StatusCode EtaCorrection::run(const std::shared_ptr<Clipboard>& clipboard) {

    // Get the clusters
    auto clusters = clipboard->getData<Cluster>(detector_->getName());
    for(auto& cluster : clusters) {
        applyEta(cluster.get());
    }

    // Return value telling analysis to keep running
    return StatusCode::Success;
}
