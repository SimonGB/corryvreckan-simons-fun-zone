/**
 * @file
 * @brief Implementation of module AlignmentDUTResidual
 *
 * @copyright Copyright (c) 2017-2022 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "AlignmentDUTResidual.h"

#include <TFormula.h>
#include <TMath.h>
#include <TProfile.h>
#include <TVirtualFitter.h>

using namespace corryvreckan;

// Global container declarations
TrackVector AlignmentDUTResidual::globalTracks;
std::shared_ptr<Detector> AlignmentDUTResidual::globalDetector;
ThreadPool* AlignmentDUTResidual::thread_pool;
std::shared_ptr<TFormula> AlignmentDUTResidual::formula_residual_x;
std::shared_ptr<TFormula> AlignmentDUTResidual::formula_residual_y;

AlignmentDUTResidual::AlignmentDUTResidual(Configuration& config, std::shared_ptr<Detector> detector)
    : Module(config, detector), m_detector(detector) {

    config_.setDefault<size_t>("iterations", 3);
    config_.setDefault<bool>("prune_tracks", false);
    config_.setDefault<bool>("align_position", true);
    config_.setDefault<bool>("align_orientation", true);
    config_.setDefault<std::string>("align_position_axes", "xy");
    config_.setDefault<std::string>("align_orientation_axes", "012");
    config_.setDefault<size_t>("max_associated_clusters", 1);
    config_.setDefault<double>("max_track_chi2ndof", 10.);
    config_.setDefault<double>("spatial_cut_sensoredge", 0.);
    config_.setDefault<unsigned int>("workers", std::max(std::thread::hardware_concurrency() - 1, 1u));
    config_.setDefaultArray<std::string>("residuals", {"x - y", "x - y"});

    m_workers = config.get<unsigned int>("workers");
    nIterations = config_.get<size_t>("iterations");
    m_pruneTracks = config_.get<bool>("prune_tracks");
    m_spatial_cut_sensoredge = config_.get<bool>("spatial_cut_sensoredge");

    m_alignPosition = config_.get<bool>("align_position");
    m_alignOrientation = config_.get<bool>("align_orientation");
    m_alignPosition_axes = config_.get<std::string>("align_position_axes");
    m_alignOrientation_axes = config_.get<std::string>("align_orientation_axes");

    if(!std::all_of(m_alignOrientation_axes.begin(), m_alignOrientation_axes.end(), ::isdigit)) {
        LOG(WARNING) << "Found non-digit characters in orientation axis designation, defaulting to all axes (\"012\")";
        m_alignOrientation_axes = "012";
    }

    std::transform(m_alignPosition_axes.begin(), m_alignPosition_axes.end(), m_alignPosition_axes.begin(), ::tolower);
    std::transform(
        m_alignOrientation_axes.begin(), m_alignOrientation_axes.end(), m_alignOrientation_axes.begin(), ::tolower);

    if(m_alignPosition) {
        LOG(INFO) << "Aligning positions";
    }

    if(m_alignOrientation) {
        LOG(INFO) << "Aligning orientations";
    }

    m_maxAssocClusters = config_.get<size_t>("max_associated_clusters");
    m_maxTrackChi2 = config_.get<double>("max_track_chi2ndof");

    // Check that we're not in a variable-alignment situation:
    if(m_detector->hasVariableAlignment()) {
        throw ModuleError("Cannot perform alignment procedure with variable alignment of detector \"" +
                          m_detector->getName() + "\"");
    }

    LOG(INFO) << "Aligning detector \"" << m_detector->getName() << "\"";
}

void AlignmentDUTResidual::initialize() {

    auto detname = m_detector->getName();
    std::string title = detname + " Residuals X;x_{track}-x [#mum];events";
    residualsXPlot = new TH1F("residualsX", title.c_str(), 1000, -500, 500);
    title = detname + " Residuals Y;y_{track}-y [#mum];events";
    residualsYPlot = new TH1F("residualsY", title.c_str(), 1000, -500, 500);
    title = detname + " Residual profile dY/X;column;y_{track}-y [#mum]";
    profile_dY_X =
        new TProfile("profile_dY_X", title.c_str(), m_detector->nPixels().x(), -0.5, m_detector->nPixels().x() - 0.5);
    title = detname + " Residual profile dY/Y;row;y_{track}-y [#mum]";
    profile_dY_Y =
        new TProfile("profile_dY_Y", title.c_str(), m_detector->nPixels().y(), -0.5, m_detector->nPixels().y() - 0.5);
    title = detname + " Residual profile dX/X;column;x_{track}-x [#mum]";
    profile_dX_X =
        new TProfile("profile_dX_X", title.c_str(), m_detector->nPixels().x(), -0.5, m_detector->nPixels().x() - 0.5);
    title = detname + " Residual profile dX/y;row;x_{track}-x [#mum]";
    profile_dX_Y =
        new TProfile("profile_dX_Y", title.c_str(), m_detector->nPixels().y(), -0.5, m_detector->nPixels().y() - 0.5);

    // Add residuals in R and Phi if detector is polar
    if(m_detector->is<PolarDetector>()) {
        title = detname + " Residuals r;r_{track}-r [um];events";
        residualsRPlot = new TH1F("residualsR", title.c_str(), 1000, -50000, 50000);
        title = detname + " Residuals #phi;#phi_{track}-#phi [#murad];events";
        residualsPhiPlot = new TH1F("residualsPhi", title.c_str(), 800, -2000, 2000);
    }

    SetResidualsFunctions();
}

StatusCode AlignmentDUTResidual::run(const std::shared_ptr<Clipboard>& clipboard) {

    // Get the tracks
    auto tracks = clipboard->getData<Track>();

    TrackVector alignmenttracks;
    std::map<std::string, std::vector<Cluster*>> alignmentclusters;

    // Make a local copy and store it
    for(auto& track : tracks) {
        auto associated_clusters = track->getAssociatedClusters(m_detector->getName());
        // Do not put tracks without clusters on the DUT to the persistent storage
        if(associated_clusters.empty()) {
            LOG(TRACE) << "Discarding track for DUT alignment since no cluster associated";
            continue;
        }
        // remove tracks at the sensor edge
        if(!m_detector->hasIntercept(track.get(), m_spatial_cut_sensoredge)) {
            continue;
        }

        // Apply selection to tracks for alignment
        if(m_pruneTracks) {
            // Only allow one associated cluster:
            if(associated_clusters.size() > m_maxAssocClusters) {
                LOG(DEBUG) << "Discarded track with " << associated_clusters.size() << " associated clusters";
                m_discardedtracks++;
                continue;
            }

            // Only allow tracks with certain Chi2/NDoF:
            if(track->getChi2ndof() > m_maxTrackChi2) {
                LOG(DEBUG) << "Discarded track with Chi2/NDoF - " << track->getChi2ndof();
                m_discardedtracks++;
                continue;
            }
        }
        LOG(TRACE) << "Cloning track with track model \"" << track->getType() << "\" for alignment";

        // Keep this track on persistent storage for alignment:
        alignmenttracks.push_back(track);
        // Append associated clusters to the list we want to keep:
        for(const auto& cluster : associated_clusters) {
            alignmentclusters[m_detector->getName()].push_back(cluster);
        }

        // Find the cluster that needs to have its position recalculated
        for(auto& associated_cluster : associated_clusters) {
            // Local position of the cluster
            auto position = associated_cluster->local();
            auto column = associated_cluster->column();
            auto row = associated_cluster->row();

            // Get the track intercept with the detector
            auto trackIntercept = m_detector->getIntercept(track.get());
            auto intercept = m_detector->globalToLocal(trackIntercept);

            // Calculate the local residuals
            double residualX = formula_residual_x->Eval(intercept.X(), position.X());
            double residualY = formula_residual_y->Eval(intercept.Y(), position.Y());

            // Recalculate residuals for polar detectors
            if(m_detector->is<PolarDetector>()) {
                auto polar_det = std::dynamic_pointer_cast<PolarDetector>(m_detector);
                // Convert cluster and intercept positions to polar coordinates
                auto cluster_polar = polar_det->getPolarPosition(column, row);
                auto intercept_polar = polar_det->getPolarPosition(intercept);

                // Calculate polar residuals
                auto residualPhi = intercept_polar.phi() - cluster_polar.phi();
                auto residualR = intercept_polar.r() - cluster_polar.r();

                // Fill polar residual histograms
                residualsRPlot->Fill(static_cast<double>(Units::convert(residualR, "um")));
                residualsPhiPlot->Fill(static_cast<double>(Units::convert(residualPhi, "urad")));
            }

            // Fill the alignment residual profile plots
            residualsXPlot->Fill(static_cast<double>(Units::convert(residualX, "um")));
            residualsYPlot->Fill(static_cast<double>(Units::convert(residualY, "um")));
            profile_dY_X->Fill(column, static_cast<double>(Units::convert(residualY, "um")), 1);
            profile_dY_Y->Fill(row, static_cast<double>(Units::convert(residualY, "um")), 1);
            profile_dX_X->Fill(column, static_cast<double>(Units::convert(residualX, "um")), 1);
            profile_dX_Y->Fill(row, static_cast<double>(Units::convert(residualX, "um")), 1);
        }

        // Since we need to refit the full track, also store the track clusters:
        for(const auto& cluster : track->getClusters()) {
            alignmentclusters[cluster->detectorID()].push_back(cluster);
        }
    }

    // Store all tracks we want for alignment on the permanent storage:
    clipboard->putPersistentData(alignmenttracks, m_detector->getName());

    // Copy the objects of all track clusters on the clipboard to persistent storage:
    for(auto& clusters : alignmentclusters) {
        clipboard->copyToPersistentData(clusters.second, clusters.first);
    }

    // Otherwise keep going
    return StatusCode::Success;
}

// METHOD 1
// This method will move the detector in question and try to minimise the
// (unbiased) residuals. It uses
// the associated cluster container on the track (no refitting of the track)
void AlignmentDUTResidual::MinimiseResiduals(Int_t&, Double_t*, Double_t& result, Double_t* par, Int_t) {

    static size_t fitIterations = 0;

    // Apply new alignment conditions
    AlignmentDUTResidual::globalDetector->update(XYZPoint(par[0], par[1], par[2]), XYZVector(par[3], par[4], par[5]));
    LOG(DEBUG) << "Updated parameters for " << AlignmentDUTResidual::globalDetector->getName();

    // The chi2 value to be returned
    result = 0.;

    LOG(DEBUG) << "Looping over " << AlignmentDUTResidual::globalTracks.size() << " tracks";

    std::vector<std::shared_future<double>> result_futures;
    auto track_refit = [&](auto& track) {
        LOG(TRACE) << "track has chi2 " << track->getChi2();
        // Update geometry of plane with new detector geometry and refit to obtain new track state, need to check if the fit
        // has failed in previous iteration
        if(track->isFitted()) {
            track->updatePlane(AlignmentDUTResidual::globalDetector->getName(),
                               AlignmentDUTResidual::globalDetector->origin().z(),
                               AlignmentDUTResidual::globalDetector->materialBudget(),
                               AlignmentDUTResidual::globalDetector->toLocal());
        } else {
            track->registerPlane(AlignmentDUTResidual::globalDetector->getName(),
                                 AlignmentDUTResidual::globalDetector->origin().z(),
                                 AlignmentDUTResidual::globalDetector->materialBudget(),
                                 AlignmentDUTResidual::globalDetector->toLocal());
            // and fit again
            track->fit();
        }
        if(!track->isFitted()) {
            LOG(WARNING) << "Refit failed - track will be discarded for this alignment step ";
            return 0.0;
        }

        double track_result = 0.;

        // Find the cluster that needs to have its position recalculated
        for(auto& associatedCluster : track->getAssociatedClusters(AlignmentDUTResidual::globalDetector->getName())) {

            // Get the track intercept with the detector
            auto position = associatedCluster->local();
            auto intercept = AlignmentDUTResidual::globalDetector->getLocalIntercept(track.get());

            // Calculate the residuals in local coordinates
            double residualX = formula_residual_x->Eval(intercept.X(), position.X());
            double residualY = formula_residual_y->Eval(intercept.Y(), position.Y());

            double errorX = associatedCluster->errorX();
            double errorY = associatedCluster->errorY();

            // Recalculate for polar detectors
            if(AlignmentDUTResidual::globalDetector->is<PolarDetector>()) {
                auto polar_det = std::dynamic_pointer_cast<PolarDetector>(AlignmentDUTResidual::globalDetector);
                // Convert cluster and intercept positions to polar coordinates
                auto cluster_polar = polar_det->getPolarPosition(associatedCluster->column(), associatedCluster->row());
                auto intercept_polar = polar_det->getPolarPosition(intercept);

                // Interpreting (Phi,R) as (X,Y)
                residualX = intercept_polar.phi() - cluster_polar.phi();
                residualY = intercept_polar.r() - cluster_polar.r();
            }

            LOG(TRACE) << "- track has intercept (" << intercept.X() << "," << intercept.Y() << ")";
            LOG(DEBUG) << "- cluster has position (" << position.X() << "," << position.Y() << ")";

            double deltachi2 = (residualX * residualX) / (errorX * errorX) + (residualY * residualY) / (errorY * errorY);
            LOG(TRACE) << "- delta chi2 = " << deltachi2;
            // Add the new residual2
            track_result += deltachi2;
            LOG(TRACE) << "- result is now " << result;
        }
        return track_result;
    };

    // Loop over all tracks
    for(auto& track : AlignmentDUTResidual::globalTracks) {
        result_futures.push_back(AlignmentDUTResidual::thread_pool->submit(track_refit, track));
    }

    for(auto& result_future : result_futures) {
        result += result_future.get();
    }

    LOG_PROGRESS(INFO, "t") << "Refit of " << result_futures.size() << " track, MINUIT iteration " << fitIterations;
    fitIterations++;
    AlignmentDUTResidual::thread_pool->wait();
}

void AlignmentDUTResidual::SetResidualsFunctions() {
    // Get definition of residuals, default x-y
    auto m_residuals = config_.getArray<std::string>("residuals");
    // Check size of array
    if(m_residuals.size() != 2) {
        throw InvalidValueError(
            config_, "residuals", "Both and only the residual_x and residual_y functions must be defined");
    }
    LOG(DEBUG) << "Definition of residual_x: " << m_residuals.at(0).c_str()
               << ", definition of residual_y: " << m_residuals.at(1).c_str()
               << " [x = track intercept, y = cluster position]";
    // Get parameters for the new definition of residuals
    auto m_parameters_residuals = config_.getArray<double>("parameters_residuals", {});
    // Define residual
    AlignmentDUTResidual::formula_residual_x =
        make_shared_no_delete<TFormula>("formula_residual_x", m_residuals.at(0).c_str(), false);
    AlignmentDUTResidual::formula_residual_y =
        make_shared_no_delete<TFormula>("formula_residual_y", m_residuals.at(1).c_str(), false);
    // Check formulas
    if(!formula_residual_x->IsValid() || !formula_residual_y->IsValid()) {
        throw InvalidValueError(config_, "residuals", "Expression is not a valid function");
    }
    // Check number of parameters
    const auto N_params_x = formula_residual_x->GetNpar();
    const auto N_params_y = formula_residual_y->GetNpar();
    if(static_cast<size_t>(N_params_x + N_params_y) != m_parameters_residuals.size()) {
        throw InvalidValueError(
            config_,
            "parameters_residuals",
            "The number of function parameters does not line up with the amount of parameters in the functions.");
    }

    // Apply parameters to the functions
    for(auto n = 0; n < N_params_x; ++n) {
        formula_residual_x->SetParameter(n, m_parameters_residuals.at(static_cast<size_t>(n)));
        LOG(DEBUG) << "residual_x: Parameter [" << n << "] = " << m_parameters_residuals.at(static_cast<size_t>(n));
    }
    for(auto n = 0; n < N_params_y; ++n) {
        formula_residual_y->SetParameter(n, m_parameters_residuals.at(static_cast<size_t>(n + N_params_x)));
        LOG(DEBUG) << "residual_y: Parameter [" << n
                   << "] = " << m_parameters_residuals.at(static_cast<size_t>(n + N_params_x));
    }
}

void AlignmentDUTResidual::finalize(const std::shared_ptr<ReadonlyClipboard>& clipboard) {

    if(m_discardedtracks > 0) {
        LOG(STATUS) << "Discarded " << m_discardedtracks << " input tracks.";
    }

    // Make the fitting object
    TVirtualFitter* residualFitter = TVirtualFitter::Fitter(nullptr, 50);
    residualFitter->SetFCN(MinimiseResiduals);

    // Set the global parameters
    AlignmentDUTResidual::globalTracks = clipboard->getPersistentData<Track>(m_detector->getName());

    // Create thread pool:
    ThreadPool::registerThreadCount(m_workers);
    AlignmentDUTResidual::thread_pool =
        new ThreadPool(m_workers,
                       m_workers * 1024,
                       [log_level = corryvreckan::Log::getReportingLevel(), log_format = corryvreckan::Log::getFormat()]() {
                           // clang-format on
                           // Initialize the threads to the same log level and format as the master setting
                           corryvreckan::Log::setReportingLevel(log_level);
                           corryvreckan::Log::setFormat(log_format);
                       });

    // Set the printout arguments of the fitter
    Double_t arglist[10];
    arglist[0] = -1;
    residualFitter->ExecuteCommand("SET PRINT", arglist, 1);

    // Set some fitter parameters
    arglist[0] = 1000;  // number of function calls
    arglist[1] = 0.001; // tolerance

    // Store the alignment shifts per detector:
    std::vector<double> shiftsX;
    std::vector<double> shiftsY;
    std::vector<double> rot0;
    std::vector<double> rot1;
    std::vector<double> rot2;

    AlignmentDUTResidual::globalDetector = m_detector;
    auto name = m_detector->getName();

    size_t n_getAssociatedClusters = 0;
    // count associated clusters:
    for(auto& track : AlignmentDUTResidual::globalTracks) {
        auto associatedClusters = track->getAssociatedClusters(name);
        for(auto& associatedCluster : associatedClusters) {
            std::string detectorID = associatedCluster->detectorID();
            if(detectorID != name) {
                continue;
            }
            n_getAssociatedClusters++;
            break;
        }
    }
    if(n_getAssociatedClusters < AlignmentDUTResidual::globalTracks.size() / 2) {
        LOG(WARNING) << "Only "
                     << 100 * static_cast<double>(n_getAssociatedClusters) /
                            static_cast<double>(AlignmentDUTResidual::globalTracks.size())
                     << "% of all tracks have associated clusters on detector " << name;
    } else {
        LOG(INFO) << 100 * static_cast<double>(n_getAssociatedClusters) /
                         static_cast<double>(AlignmentDUTResidual::globalTracks.size())
                  << "% of all tracks have associated clusters on detector " << name;
    }

    LOG(STATUS) << name << " initial alignment: " << std::endl
                << "T" << Units::display(m_detector->displacement(), {"mm", "um"}) << " R"
                << Units::display(m_detector->rotation(), {"deg"});

    // Add the parameters to the fitter (z displacement not allowed to move!)
    if(m_alignPosition && m_alignPosition_axes.find('x') != std::string::npos) {
        residualFitter->SetParameter(0, (name + "_displacementX").c_str(), m_detector->displacement().X(), 0.01, -50, 50);
    } else {
        residualFitter->SetParameter(0, (name + "_displacementX").c_str(), m_detector->displacement().X(), 0, -50, 50);
    }
    if(m_alignPosition && m_alignPosition_axes.find('y') != std::string::npos) {
        residualFitter->SetParameter(1, (name + "_displacementY").c_str(), m_detector->displacement().Y(), 0.01, -50, 50);
    } else {
        residualFitter->SetParameter(1, (name + "_displacementY").c_str(), m_detector->displacement().Y(), 0, -50, 50);
    }

    // Z is never changed:
    residualFitter->SetParameter(2, (name + "_displacementZ").c_str(), m_detector->displacement().Z(), 0, -10, 500);

    if(m_alignOrientation && m_alignOrientation_axes.find('0') != std::string::npos) {
        residualFitter->SetParameter(3, (name + "_rotation0").c_str(), m_detector->rotation().X(), 0.001, -6.30, 6.30);
    } else {
        residualFitter->SetParameter(3, (name + "_rotation0").c_str(), m_detector->rotation().X(), 0, -6.30, 6.30);
    }
    if(m_alignOrientation && m_alignOrientation_axes.find('1') != std::string::npos) {
        residualFitter->SetParameter(4, (name + "_rotation1").c_str(), m_detector->rotation().Y(), 0.001, -6.30, 6.30);
    } else {
        residualFitter->SetParameter(4, (name + "_rotation1").c_str(), m_detector->rotation().Y(), 0, -6.30, 6.30);
    }
    if(m_alignOrientation && m_alignOrientation_axes.find('2') != std::string::npos) {
        residualFitter->SetParameter(5, (name + "_rotation2").c_str(), m_detector->rotation().Z(), 0.001, -6.30, 6.30);
    } else {
        residualFitter->SetParameter(5, (name + "_rotation2").c_str(), m_detector->rotation().Z(), 0, -6.30, 6.30);
    }

    for(size_t iteration = 0; iteration < nIterations; iteration++) {

        auto old_position = m_detector->displacement();
        auto old_orientation = m_detector->rotation();

        // Fit this plane (minimising global track chi2)
        residualFitter->ExecuteCommand("MIGRAD", arglist, 2);

        // Set the alignment parameters of this plane to be the optimised values from the alignment
        m_detector->update(
            XYZPoint(residualFitter->GetParameter(0), residualFitter->GetParameter(1), residualFitter->GetParameter(2)),
            XYZVector(residualFitter->GetParameter(3), residualFitter->GetParameter(4), residualFitter->GetParameter(5)));

        // Store corrections:
        shiftsX.push_back(static_cast<double>(Units::convert(m_detector->displacement().X() - old_position.X(), "um")));
        shiftsY.push_back(static_cast<double>(Units::convert(m_detector->displacement().Y() - old_position.Y(), "um")));
        rot0.push_back(static_cast<double>(Units::convert(m_detector->rotation().X() - old_orientation.X(), "deg")));
        rot1.push_back(static_cast<double>(Units::convert(m_detector->rotation().Y() - old_orientation.Y(), "deg")));
        rot2.push_back(static_cast<double>(Units::convert(m_detector->rotation().Z() - old_orientation.Z(), "deg")));

        LOG(INFO) << m_detector->getName() << "/" << iteration << " dT"
                  << Units::display(m_detector->displacement() - old_position, {"mm", "um"}) << " dR"
                  << Units::display(m_detector->rotation() - old_orientation, {"deg"});
    }

    LOG(STATUS) << m_detector->getName() << " new alignment: " << std::endl
                << "T" << Units::display(m_detector->displacement(), {"mm", "um"}) << " R"
                << Units::display(m_detector->rotation(), {"deg"});

    std::vector<double> iterations(nIterations);
    std::iota(std::begin(iterations), std::end(iterations), 0);

    std::string graph_name = "alignment_correction_displacementX_" + m_detector->getName();
    align_correction_shiftX = new TGraph(static_cast<int>(shiftsX.size()), &iterations[0], &shiftsX[0]);
    align_correction_shiftX->GetXaxis()->SetTitle("# iteration");
    align_correction_shiftX->GetYaxis()->SetTitle("correction [#mum]");
    align_correction_shiftX->Write(graph_name.c_str());

    graph_name = "alignment_correction_displacementY_" + m_detector->getName();
    align_correction_shiftY = new TGraph(static_cast<int>(shiftsY.size()), &iterations[0], &shiftsY[0]);
    align_correction_shiftY->GetXaxis()->SetTitle("# iteration");
    align_correction_shiftY->GetYaxis()->SetTitle("correction [#mum]");
    align_correction_shiftY->Write(graph_name.c_str());

    graph_name = "alignment_correction_rotation0_" + m_detector->getName() + "_in_mode_" + m_detector->orientation_mode();
    align_correction_rot0 = new TGraph(static_cast<int>(rot0.size()), &iterations[0], &rot0[0]);
    align_correction_rot0->GetXaxis()->SetTitle("# iteration");
    align_correction_rot0->GetYaxis()->SetTitle("correction [#mum]");
    align_correction_rot0->Write(graph_name.c_str());

    graph_name = "alignment_correction_rotation1_" + m_detector->getName() + "_in_mode_" + m_detector->orientation_mode();
    align_correction_rot1 = new TGraph(static_cast<int>(rot1.size()), &iterations[0], &rot1[0]);
    align_correction_rot1->GetXaxis()->SetTitle("# iteration");
    align_correction_rot1->GetYaxis()->SetTitle("correction [#mum]");
    align_correction_rot1->Write(graph_name.c_str());

    graph_name = "alignment_correction_rotation2_" + m_detector->getName() + "_in_mode_" + m_detector->orientation_mode();
    align_correction_rot2 = new TGraph(static_cast<int>(rot2.size()), &iterations[0], &rot2[0]);
    align_correction_rot2->GetXaxis()->SetTitle("# iteration");
    align_correction_rot2->GetYaxis()->SetTitle("correction [#mum]");
    align_correction_rot2->Write(graph_name.c_str());

    // Clean up local track storage
    AlignmentDUTResidual::globalTracks.clear();
    AlignmentDUTResidual::globalDetector.reset();
}
