/**
 * @file
 * @brief Implementation of [AnalysisEfficiency] module
 *
 * @copyright Copyright (c) 2017-2022 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "AnalysisEfficiency.h"

#include "objects/Cluster.hpp"
#include "objects/Pixel.hpp"
#include "objects/Track.hpp"

using namespace corryvreckan;

AnalysisEfficiency::AnalysisEfficiency(Configuration& config, std::shared_ptr<Detector> detector)
    : Module(config, detector) {
    m_detector = detector;

    config_.setDefault<double>("time_cut_frameedge", Units::get<double>(20, "ns"));
    config_.setDefault<double>("chi2ndof_cut", 3.);
    config_.setDefault<ROOT::Math::XYPoint>("inpixel_bin_size", {Units::get(1.0, "um"), Units::get(1.0, "um")});
    config_.setDefault<ROOT::Math::XYPoint>("inpixel_cut_edge", {Units::get(5.0, "um"), Units::get(5.0, "um")});
    config_.setDefault<double>("masked_pixel_distance_cut", 1.);
    config_.setDefault<double>("spatial_cut_sensoredge", 1.);
    config_.setDefault<FakeRateMethod>("fake_rate_method", FakeRateMethod::RADIUS);
    config_.setDefault<double>("fake_rate_distance", 2.);
    config_.setDefault<int>("n_charge_bins", 1000);
    config_.setDefault<double>("charge_histo_range", 1000.0);

    m_timeCutFrameEdge = config_.get<double>("time_cut_frameedge");
    m_chi2ndofCut = config_.get<double>("chi2ndof_cut");
    require_associated_cluster_on_ = config_.getArray<std::string>("require_associated_cluster_on", {});
    m_maskedPixelDistanceCut = config_.get<int>("masked_pixel_distance_cut");
    spatial_cut_sensoredge = config_.get<double>("spatial_cut_sensoredge");
    m_fake_rate_method = config_.get<FakeRateMethod>("fake_rate_method");
    m_fake_rate_distance = config_.get<double>("fake_rate_distance");
    m_n_charge_bins = config_.get<int>("n_charge_bins");
    m_charge_histo_range = config_.get<double>("charge_histo_range");

    if(config_.getArray<double>("inpixel_bin_size").size() == 2) {
        m_inpixelBinSize = config_.get<ROOT::Math::XYPoint>("inpixel_bin_size");
    } else {
        auto binsize = config_.get<double>("inpixel_bin_size");
        m_inpixelBinSize = ROOT::Math::XYPoint(binsize, binsize);
    }

    if(config_.getArray<double>("inpixel_cut_edge").size() == 2) {
        m_inpixelEdgeCut = config_.get<ROOT::Math::XYPoint>("inpixel_cut_edge");
    } else {
        auto edgecut = config_.get<double>("inpixel_cut_edge");
        m_inpixelEdgeCut = ROOT::Math::XYPoint(edgecut, edgecut);
    }
}
void AnalysisEfficiency::initialize() {

    hPixelEfficiency = new TH1D(
        "hPixelEfficiency", "hPixelEfficiency; single pixel efficiency; # entries", 201, 0, 1.005); // get 0.5%-wide bins

    hPixelEfficiencyMatrix = new TH1D("hPixelEfficiencyMatrix",
                                      "hPixelEfficiencyMatrix; single pixel efficiency; # entries",
                                      201,
                                      0,
                                      1.005); // get 0.5%-wide bins

    auto pitch_x = static_cast<double>(Units::convert(m_detector->getPitch().X(), "um"));
    auto pitch_y = static_cast<double>(Units::convert(m_detector->getPitch().Y(), "um"));

    auto nbins_x = static_cast<int>(std::ceil(m_detector->getPitch().X() / m_inpixelBinSize.x()));
    auto nbins_y = static_cast<int>(std::ceil(m_detector->getPitch().Y() / m_inpixelBinSize.y()));
    if(nbins_x > 1e4 || nbins_y > 1e4) {
        throw InvalidValueError(config_, "inpixel_bin_size", "Too many bins for in-pixel histograms.");
    }
    std::string title =
        m_detector->getName() + " Pixel efficiency map;in-pixel x_{track} [#mum];in-pixel y_{track} #mum;#epsilon";
    hPixelEfficiencyMap_trackPos_TProfile = new TProfile2D("pixelEfficiencyMap_trackPos_TProfile",
                                                           title.c_str(),
                                                           nbins_x,
                                                           -pitch_x / 2.,
                                                           pitch_x / 2.,
                                                           nbins_y,
                                                           -pitch_y / 2.,
                                                           pitch_y / 2.,
                                                           0,
                                                           1);
    hPixelEfficiencyMap_trackPos = new TEfficiency("pixelEfficiencyMap_trackPos",
                                                   title.c_str(),
                                                   nbins_x,
                                                   -pitch_x / 2.,
                                                   pitch_x / 2.,
                                                   nbins_y,
                                                   -pitch_y / 2.,
                                                   pitch_y / 2.);
    hPixelEfficiencyMap_trackPos->SetDirectory(this->getROOTDirectory());

    title = m_detector->getName() + " Chip efficiency map;x [px];y [px];#epsilon";
    hChipEfficiencyMap_trackPos_TProfile = new TProfile2D("chipEfficiencyMap_trackPos_TProfile",
                                                          title.c_str(),
                                                          m_detector->nPixels().X(),
                                                          -0.5,
                                                          m_detector->nPixels().X() - 0.5,
                                                          m_detector->nPixels().Y(),
                                                          -0.5,
                                                          m_detector->nPixels().Y() - 0.5,
                                                          0,
                                                          1);
    hChipEfficiencyMap_trackPos = new TEfficiency("chipEfficiencyMap_trackPos",
                                                  title.c_str(),
                                                  m_detector->nPixels().X(),
                                                  -0.5,
                                                  m_detector->nPixels().X() - 0.5,
                                                  m_detector->nPixels().Y(),
                                                  -0.5,
                                                  m_detector->nPixels().Y() - 0.5);
    hChipEfficiencyMap_trackPos->SetDirectory(this->getROOTDirectory());
    title = m_detector->getName() + " Pixel efficiency matrix;x [px];y [px];#epsilon";
    hPixelEfficiencyMatrix_TProfile = new TProfile2D("hPixelEfficiencyMatrixTProfile",
                                                     title.c_str(),
                                                     m_detector->nPixels().X(),
                                                     -0.5,
                                                     m_detector->nPixels().X() - 0.5,
                                                     m_detector->nPixels().Y(),
                                                     -0.5,
                                                     m_detector->nPixels().Y() - 0.5,
                                                     0,
                                                     1);

    title = m_detector->getName() + " Global efficiency map;x [mm];y [mm];#epsilon";
    hGlobalEfficiencyMap_trackPos_TProfile = new TProfile2D("globalEfficiencyMap_trackPos_TProfile",
                                                            title.c_str(),
                                                            300,
                                                            m_detector->displacement().X() - 1.5 * m_detector->getSize().X(),
                                                            m_detector->displacement().X() + 1.5 * m_detector->getSize().X(),
                                                            300,
                                                            m_detector->displacement().Y() - 1.5 * m_detector->getSize().Y(),
                                                            m_detector->displacement().Y() + 1.5 * m_detector->getSize().Y(),
                                                            0,
                                                            1);
    hGlobalEfficiencyMap_trackPos = new TEfficiency("globalEfficiencyMap_trackPos",
                                                    title.c_str(),
                                                    300,
                                                    m_detector->displacement().X() - 1.5 * m_detector->getSize().X(),
                                                    m_detector->displacement().X() + 1.5 * m_detector->getSize().X(),
                                                    300,
                                                    m_detector->displacement().Y() - 1.5 * m_detector->getSize().Y(),
                                                    m_detector->displacement().Y() + 1.5 * m_detector->getSize().Y());
    hGlobalEfficiencyMap_trackPos->SetDirectory(this->getROOTDirectory());

    hDistanceCluster = new TH1D("distanceTrackHit",
                                "distance between track and hit; | #vec{track} - #vec{dut} | [mm]",
                                static_cast<int>(std::sqrt(m_detector->getPitch().x() * m_detector->getPitch().y())),
                                0,
                                std::sqrt(m_detector->getPitch().x() * m_detector->getPitch().y()));
    hDistanceCluster_track = new TH2D("distanceTrackHit2D",
                                      "distance between track and hit; track_x - dut_x [mm]; track_y - dut_y [mm] ",
                                      150,
                                      -1.5 * m_detector->getPitch().x(),
                                      1.5 * m_detector->getPitch().x(),
                                      150,
                                      -1.5 * m_detector->getPitch().y(),
                                      1.5 * m_detector->getPitch().y());
    eTotalEfficiency = new TEfficiency("eTotalEfficiency", "totalEfficiency;;#epsilon", 1, 0, 1);
    eTotalEfficiency->SetDirectory(this->getROOTDirectory());

    efficiencyColumns = new TEfficiency("efficiencyColumns",
                                        "Efficiency vs. column number; column; #epsilon",
                                        m_detector->nPixels().X(),
                                        -0.5,
                                        m_detector->nPixels().X() - 0.5);
    efficiencyColumns->SetDirectory(this->getROOTDirectory());
    efficiencyRows = new TEfficiency("efficiencyRows",
                                     "Efficiency vs. row number; row; #epsilon",
                                     m_detector->nPixels().Y(),
                                     -0.5,
                                     m_detector->nPixels().Y() - 0.5);
    efficiencyRows->SetDirectory(this->getROOTDirectory());
    efficiencyVsTime = new TEfficiency("efficiencyVsTime", "Efficiency vs. time; time [s]; #epsilon", 3000, 0, 3000);
    efficiencyVsTime->SetDirectory(this->getROOTDirectory());
    efficiencyVsTimeLong =
        new TEfficiency("efficiencyVsTimeLong", "Efficiency vs. time; time [s]; #epsilon", 3000, 0, 30000);
    efficiencyVsTimeLong->SetDirectory(this->getROOTDirectory());

    // initialize matrix with hit timestamps to all 0:
    auto nRows = static_cast<size_t>(m_detector->nPixels().Y());
    auto nCols = static_cast<size_t>(m_detector->nPixels().X());
    std::vector<double> v_row(nRows, 0.); // create vector will zeros of length <nRows>
    prev_hit_ts.assign(nCols, v_row);     // use vector v_row to construct matrix

    createInPixelRoiPlots();
    createFakeRatePlots();
    createTrackTimePlots();
}

void AnalysisEfficiency::createInPixelRoiPlots() {
    TDirectory* directory = getROOTDirectory();
    TDirectory* in_pixel_roi = directory->mkdir("inpixelROI");
    if(in_pixel_roi == nullptr) {
        throw RuntimeError("Cannot create or access fake rate ROOT directory for module " + this->getUniqueName());
    }
    in_pixel_roi->cd();
    auto pitch_x = static_cast<double>(Units::convert(m_detector->getPitch().X(), "um"));
    auto pitch_y = static_cast<double>(Units::convert(m_detector->getPitch().Y(), "um"));

    auto nbins_x = static_cast<int>(std::ceil(m_detector->getPitch().X() / m_inpixelBinSize.X()));
    auto nbins_y = static_cast<int>(std::ceil(m_detector->getPitch().Y() / m_inpixelBinSize.Y()));

    std::string title = m_detector->getName() +
                        " Pixel efficiency map (in-pixel ROI);in-pixel x_{track} [#mum];in-pixel y_{track} #mum;#epsilon";
    hPixelEfficiencyMap_inPixelROI_trackPos_TProfile = new TProfile2D("pixelEfficiencyMap_inPixelROI_trackPos_TProfile",
                                                                      title.c_str(),
                                                                      nbins_x,
                                                                      -pitch_x / 2.,
                                                                      pitch_x / 2.,
                                                                      nbins_y,
                                                                      -pitch_y / 2.,
                                                                      pitch_y / 2.,
                                                                      0,
                                                                      1);
    eTotalEfficiency_inPixelROI = new TEfficiency(
        "eTotalEfficiency_inPixelROI", "eTotalEfficiency_inPixelROI;;#epsilon (within in-pixel ROI)", 1, 0, 1);
    eTotalEfficiency_inPixelROI->SetDirectory(this->getROOTDirectory()->GetDirectory("inpixelROI"));
}

void AnalysisEfficiency::createTrackTimePlots() {
    TDirectory* directory = getROOTDirectory();
    TDirectory* correlationsToPrevTrack = directory->mkdir("correlationsToPrevTrack");
    if(correlationsToPrevTrack == nullptr) {
        throw RuntimeError("Cannot create or access fake rate ROOT directory for module " + this->getUniqueName());
    }
    correlationsToPrevTrack->cd();

    hTrackTimeToPrevHit_matched =
        new TH1D("trackTimeToPrevHit_matched", "trackTimeToPrevHit_matched;time to prev hit [us];# events", 1e6, 0, 1e6);
    hTrackTimeToPrevHit_notmatched = new TH1D(
        "trackTimeToPrevHit_notmatched", "trackTimeToPrevHit_notmatched;time to prev hit [us];# events", 1e6, 0, 1e6);

    std::string title = m_detector->getName() + "time difference to previous track (if this has assoc cluster)";
    hTimeDiffPrevTrack_assocCluster = new TH1D("timeDiffPrevTrack_assocCluster", title.c_str(), 11000, -1000, 10000);
    hTimeDiffPrevTrack_assocCluster->GetXaxis()->SetTitle("time diff [#mus]");
    hTimeDiffPrevTrack_assocCluster->GetYaxis()->SetTitle("events");
    title = m_detector->getName() + "time difference to previous track (if this has no assoc cluster)";
    hTimeDiffPrevTrack_noAssocCluster = new TH1D("timeDiffPrevTrack_noAssocCluster", title.c_str(), 11000, -1000, 10000);
    hTimeDiffPrevTrack_noAssocCluster->GetXaxis()->SetTitle("time diff [#mus]");
    hTimeDiffPrevTrack_noAssocCluster->GetYaxis()->SetTitle("events");

    hRowDiffPrevTrack_assocCluster =
        new TH1D("rowDiffPrevTrack_assocCluster",
                 "rowDiffPrevTrack_assocCluster; row difference (matched track to prev track) [px];# events",
                 2 * m_detector->nPixels().Y(),
                 -m_detector->nPixels().Y() - 0.5,
                 m_detector->nPixels().Y() - 0.5);

    hColDiffPrevTrack_assocCluster =
        new TH1D("colDiffPrevTrack_assocCluster",
                 "colDiffPrevTrack_assocCluster;column difference (matched track to prev track) [px];# events",
                 2 * m_detector->nPixels().X(),
                 -m_detector->nPixels().X() - 0.5,
                 m_detector->nPixels().X() - 0.5);

    hRowDiffPrevTrack_noAssocCluster =
        new TH1D("rowDiffPrevTrack_noAssocCluster",
                 "rowDiffPrevTrack_noAssocCluster;row difference (non-matched track - prev track) [px];# events",
                 2 * m_detector->nPixels().Y(),
                 -m_detector->nPixels().Y() - 0.5,
                 m_detector->nPixels().Y() - 0.5);

    hColDiffPrevTrack_noAssocCluster =
        new TH1D("colDiffPrevTrack_noAassocCluster",
                 "colDiffPrevTrack_noAssocCluster;column difference (non-matched track - prev track) [px];# events",
                 2 * m_detector->nPixels().X(),
                 -m_detector->nPixels().X() - 0.5,
                 m_detector->nPixels().X() - 0.5);

    hPosDiffPrevTrack_assocCluster = new TH2D("posDiffPrevTrack_assocCluster",
                                              "posDiffPrevTrack_assocCluster;column difference (matched track - prev track) "
                                              "[px];row difference (matched track - prev track) [px];# events",
                                              2 * m_detector->nPixels().X(),
                                              -m_detector->nPixels().X() - 0.5,
                                              m_detector->nPixels().X() - 0.5,
                                              2 * m_detector->nPixels().Y(),
                                              -m_detector->nPixels().Y() - 0.5,
                                              m_detector->nPixels().Y() - 0.5);

    hPosDiffPrevTrack_noAssocCluster =
        new TH2D("posDiffPrevTrack_noAssocCluster",
                 "posDiffPrevTrack_noAssocCluster;column difference (non-matched track - prev track) [px];row difference "
                 "(non-matched track - prev track) [px];# events",
                 2 * m_detector->nPixels().X(),
                 -m_detector->nPixels().X() - 0.5,
                 m_detector->nPixels().X() - 0.5,
                 2 * m_detector->nPixels().Y(),
                 -m_detector->nPixels().Y() - 0.5,
                 m_detector->nPixels().Y() - 0.5);
    htimeRes_cluster_size =
        new TH2D("timeDiff_hit-track_vs_clustersize",
                 "Time difference between track and clusters vs cluster size; time difference [ns];cluster size [px]",
                 1000,
                 -499.5,
                 500.5,
                 6,
                 -.5,
                 5.5);
}

void AnalysisEfficiency::createFakeRatePlots() {
    TDirectory* directory = getROOTDirectory();
    TDirectory* fake_rate_directory = directory->mkdir("fake_rate");
    if(fake_rate_directory == nullptr) {
        throw RuntimeError("Cannot create or access fake rate ROOT directory for module " + this->getUniqueName());
    }
    fake_rate_directory->cd();

    std::string title = m_detector->getName() + " number of fake hits per event; hits; events";
    hFakePixelPerEvent = new TH1D("hFakePixelPerEvent", title.c_str(), 25, 0 - 0.5, 25 - 0.5);

    title = m_detector->getName() + " pixel fake hits per event;x [px];y [px]; hits";
    fakePixelPerEventMap = new TH2D("fakePixelPerEventMap",
                                    title.c_str(),
                                    m_detector->nPixels().X(),
                                    -0.5,
                                    m_detector->nPixels().X() - 0.5,
                                    m_detector->nPixels().Y(),
                                    -0.5,
                                    m_detector->nPixels().Y() - 0.5);

    title = m_detector->getName() + " pixel fake hits per event vs. time; time [s]; hits";
    fakePixelPerEventVsTime = new TProfile("fakePixelPerEventVsTime", title.c_str(), 3000, 0, 3000);

    title = m_detector->getName() + " pixel fake hits per event vs. time; time [s]; hits";
    fakePixelPerEventVsTimeLong = new TProfile("efficiencyVsTimeLong", title.c_str(), 3000, 0, 30000);

    title = m_detector->getName() + " charge distribution for fake pixels; charge [a.u.]; entries";
    hFakePixelCharge = new TH1D("hFakePixelCharge", title.c_str(), m_n_charge_bins, 0.0, m_charge_histo_range);

    title = m_detector->getName() + " charge distribution for fake clusters; charge [a.u.]; entries";
    hFakeClusterCharge = new TH1D("hFakeClusterCharge", title.c_str(), m_n_charge_bins, 0.0, m_charge_histo_range);

    title = m_detector->getName() + " number of fake clusters per event; clusters; events";
    hFakeClusterPerEvent = new TH1D("hFakeClusterPerEvent", title.c_str(), 25, 0 - 0.5, 25 - 0.5);

    title = m_detector->getName() + " cluster size of fake clusters; cluster size; events";
    hFakeClusterSize = new TH1D("hFakeClusterSize", title.c_str(), 25, 0 - 0.5, 25 - 0.5);
}

StatusCode AnalysisEfficiency::run(const std::shared_ptr<Clipboard>& clipboard) {

    // Get the telescope tracks from the clipboard
    auto tracks = clipboard->getData<Track>();

    auto pitch_x = m_detector->getPitch().X();
    auto pitch_y = m_detector->getPitch().Y();
    // Get the event:
    auto event = clipboard->getEvent();

    // Loop over all tracks
    for(auto& track : tracks) {
        n_track++;
        bool has_associated_cluster = false;
        bool is_within_roi = true;
        LOG(DEBUG) << "Looking at next track";

        // Cut on the chi2/ndof
        if(track->getChi2ndof() > m_chi2ndofCut) {
            LOG(DEBUG) << " - track discarded due to Chi2/ndof";
            n_chi2++;
            continue;
        }

        // Check if it intercepts the DUT
        auto globalIntercept = m_detector->getIntercept(track.get());
        auto localIntercept = m_detector->globalToLocal(globalIntercept);

        LOG(TRACE) << " Checking if track is outside DUT area";
        if(!m_detector->hasIntercept(track.get(), spatial_cut_sensoredge)) {
            LOG(DEBUG) << " - track outside DUT area: " << localIntercept;
            n_dut++;
            continue;
        }

        // Check that track is within region of interest using winding number algorithm
        LOG(TRACE) << " Checking if track is outside ROI";
        if(!m_detector->isWithinROI(track.get())) {
            LOG(DEBUG) << " - track outside ROI";
            n_roi++;
            is_within_roi = false;
            // here we don't continue because only some particular histograms shall be effected
        }

        // Check that it doesn't go through/near a masked pixel
        LOG(TRACE) << " Checking if track is close to masked pixel";
        if(m_detector->hitMasked(track.get(), m_maskedPixelDistanceCut)) {
            n_masked++;
            LOG(DEBUG) << " - track close to masked pixel";
            continue;
        }

        // Discard tracks which are very close to the frame edges
        if(fabs(track->timestamp() - event->end()) < m_timeCutFrameEdge) {
            // Late edge - eventEnd points to the end of the frame`
            LOG(DEBUG) << " - track close to end of readout frame: "
                       << Units::display(fabs(track->timestamp() - event->end()), {"us", "ns"}) << " at "
                       << Units::display(track->timestamp(), {"us"});
            n_frameedge++;
            continue;
        } else if(fabs(track->timestamp() - event->start()) < m_timeCutFrameEdge) {
            // Early edge - eventStart points to the beginning of the frame
            LOG(DEBUG) << " - track close to start of readout frame: "
                       << Units::display(fabs(track->timestamp() - event->start()), {"us", "ns"}) << " at "
                       << Units::display(track->timestamp(), {"us"});
            n_frameedge++;
            continue;
        }

        // check if track has an associated cluster on required detector(s):
        auto foundRequiredAssocCluster = [this](Track* t) {
            for(auto& requireAssocCluster : require_associated_cluster_on_) {
                if(!requireAssocCluster.empty() && t->getAssociatedClusters(requireAssocCluster).size() == 0) {
                    LOG(DEBUG) << "No associated cluster from required detector " << requireAssocCluster << " on the track.";
                    return false;
                }
            }
            return true;
        };
        if(!foundRequiredAssocCluster(track.get())) {
            n_requirecluster++;
            continue;
        }

        // Count this as reference track:
        total_tracks++;

        // Calculate in-pixel position of track in microns
        auto inpixel = m_detector->inPixel(localIntercept);
        auto xmod = inpixel.X();
        auto ymod = inpixel.Y();
        auto xmod_um = xmod * 1000.; // mm->um (for plotting)
        auto ymod_um = ymod * 1000.; // mm->um (for plotting)

        bool isWithinInPixelROI =
            (pitch_x - fabs(xmod * 2.) > m_inpixelEdgeCut.x()) && (pitch_y - fabs(ymod * 2.) > m_inpixelEdgeCut.y());

        // Get the DUT clusters from the clipboard, that are assigned to the track
        auto associated_clusters = track->getAssociatedClusters(m_detector->getName());
        if(associated_clusters.size() > 0) {
            auto cluster = track->getClosestCluster(m_detector->getName());
            has_associated_cluster = true;
            matched_tracks++;
            auto pixels = cluster->pixels();
            for(auto& pixel : pixels) {
                if((pixel->column() == static_cast<int>(m_detector->getColumn(localIntercept)) &&
                    pixel->row() == static_cast<int>(m_detector->getRow(localIntercept))) &&
                   isWithinInPixelROI) {
                    hPixelEfficiencyMatrix_TProfile->Fill(pixel->column(), pixel->row(), 1);
                    break; // There cannot be a second pixel within the cluster through which the track goes.
                }
            }

            auto clusterLocal = m_detector->globalToLocal(cluster->global());

            auto distance =
                ROOT::Math::XYZVector(localIntercept.x() - clusterLocal.x(), localIntercept.y() - clusterLocal.y(), 0);
            hDistanceCluster_track->Fill(distance.X(), distance.Y());
            hDistanceCluster->Fill(std::sqrt(distance.Mag2()));
        }

        if(!has_associated_cluster && isWithinInPixelROI) {
            hPixelEfficiencyMatrix_TProfile->Fill(
                m_detector->getColumn(localIntercept), m_detector->getRow(localIntercept), 0);
        }

        hGlobalEfficiencyMap_trackPos_TProfile->Fill(globalIntercept.X(), globalIntercept.Y(), has_associated_cluster);
        hGlobalEfficiencyMap_trackPos->Fill(has_associated_cluster, globalIntercept.X(), globalIntercept.Y());

        hChipEfficiencyMap_trackPos_TProfile->Fill(
            m_detector->getColumn(localIntercept), m_detector->getRow(localIntercept), has_associated_cluster);
        hChipEfficiencyMap_trackPos->Fill(
            has_associated_cluster, m_detector->getColumn(localIntercept), m_detector->getRow(localIntercept));

        // For pixels, only look at the ROI:
        if(is_within_roi) {
            hPixelEfficiencyMap_trackPos_TProfile->Fill(xmod_um, ymod_um, has_associated_cluster);
            hPixelEfficiencyMap_trackPos->Fill(has_associated_cluster, xmod_um, ymod_um);
            eTotalEfficiency->Fill(has_associated_cluster, 0); // use 0th bin for total efficiency
            efficiencyColumns->Fill(has_associated_cluster, m_detector->getColumn(localIntercept));
            efficiencyRows->Fill(has_associated_cluster, m_detector->getRow(localIntercept));
            efficiencyVsTime->Fill(has_associated_cluster, track->timestamp() / 1e9);     // convert nanoseconds to seconds
            efficiencyVsTimeLong->Fill(has_associated_cluster, track->timestamp() / 1e9); // convert nanoseconds to seconds
            if(isWithinInPixelROI) {
                hPixelEfficiencyMap_inPixelROI_trackPos_TProfile->Fill(xmod_um, ymod_um, has_associated_cluster);
                eTotalEfficiency_inPixelROI->Fill(has_associated_cluster, 0); // use 0th bin for total efficiency
            }
        }

        auto intercept_col = static_cast<size_t>(m_detector->getColumn(localIntercept));
        auto intercept_row = static_cast<size_t>(m_detector->getRow(localIntercept));

        if(has_associated_cluster) {
            for(auto c : associated_clusters) {
                htimeRes_cluster_size->Fill(track->timestamp() - c->timestamp(),
                                            ((c->pixels().size() > 4) ? 5.0 : static_cast<double>(c->pixels().size())));
            }
            hTimeDiffPrevTrack_assocCluster->Fill(
                static_cast<double>(Units::convert(track->timestamp() - last_track_timestamp, "us")));
            hRowDiffPrevTrack_assocCluster->Fill(m_detector->getRow(localIntercept) - last_track_row);
            hColDiffPrevTrack_assocCluster->Fill(m_detector->getColumn(localIntercept) - last_track_col);
            hPosDiffPrevTrack_assocCluster->Fill(m_detector->getColumn(localIntercept) - last_track_col,
                                                 m_detector->getRow(localIntercept) - last_track_row);
            if((prev_hit_ts.at(intercept_col)).at(intercept_row) != 0) {
                hTrackTimeToPrevHit_matched->Fill(static_cast<double>(
                    Units::convert(track->timestamp() - prev_hit_ts.at(intercept_col).at(intercept_row), "us")));
            }
        } else {
            hTimeDiffPrevTrack_noAssocCluster->Fill(
                static_cast<double>(Units::convert(track->timestamp() - last_track_timestamp, "us")));
            hRowDiffPrevTrack_noAssocCluster->Fill(m_detector->getRow(localIntercept) - last_track_row);
            hColDiffPrevTrack_noAssocCluster->Fill(m_detector->getColumn(localIntercept) - last_track_col);
            hPosDiffPrevTrack_noAssocCluster->Fill(m_detector->getColumn(localIntercept) - last_track_col,
                                                   m_detector->getRow(localIntercept) - last_track_row);
            if((prev_hit_ts.at(intercept_col)).at(intercept_row) != 0) {
                LOG(DEBUG) << "Found a time difference of "
                           << Units::display(track->timestamp() - prev_hit_ts.at(intercept_col).at(intercept_row), "us");
                hTrackTimeToPrevHit_notmatched->Fill(static_cast<double>(
                    Units::convert(track->timestamp() - prev_hit_ts.at(intercept_col).at(intercept_row), "us")));
            }
        }
        last_track_timestamp = track->timestamp();
        last_track_col = m_detector->getColumn(localIntercept);
        last_track_row = m_detector->getRow(localIntercept);
    } // end loop over tracks

    // Before going to the next event, loop over all pixels (all hits incl. noise)
    // and fill matrix with timestamps of previous pixels.
    auto pixels = clipboard->getData<Pixel>(m_detector->getName());
    if(pixels.empty()) {
        LOG(DEBUG) << "Detector " << m_detector->getName() << " does not have any pixels on the clipboard";
    }

    for(auto& pixel : pixels) {
        if(pixel->column() > m_detector->nPixels().X() || pixel->row() > m_detector->nPixels().Y()) {
            continue;
        }
        try {
            prev_hit_ts.at(static_cast<size_t>(pixel->column())).at(static_cast<size_t>(pixel->row())) = pixel->timestamp();
        } catch(std::out_of_range&) {
        }
    }

    // fake rate analysis
    if(m_fake_rate_method == FakeRateMethod::RADIUS) {
        LOG_ONCE(STATUS) << "Estimating fake rate based on radial cut around pixels (RADIUS method).";

        int fake_hits = 0;
        int fake_clusters = 0;

        // get and iterate dut clusters from clipboard
        auto clusters = clipboard->getData<Cluster>(m_detector->getName());
        for(auto& cluster : clusters) {

            bool track_too_close = false;
            // iterate the tracks from the clipboard
            for(auto& track : tracks) {

                // discard tracks without intercept, using a tolerance defined
                // by the radial cut which we will use later.
                if(m_detector->hasIntercept(track.get(), -m_fake_rate_distance)) {
                    continue;
                }

                // now check if the track is too close to the considered cluster,
                // using the distance in units of the pitch.
                auto interceptLocal = m_detector->getLocalIntercept(track.get());
                double norm_xdistance = (interceptLocal.X() - cluster->local().x()) / m_detector->getPitch().X();
                double norm_ydistance = (interceptLocal.Y() - cluster->local().y()) / m_detector->getPitch().Y();
                double norm_rdistance = sqrt(norm_xdistance * norm_xdistance + norm_ydistance * norm_ydistance);
                if(norm_rdistance < m_fake_rate_distance) {
                    track_too_close = true;
                    break;
                }
            }

            // now study cluster and pixel properties, which seems to be a fake
            // cluster, as there is no track nearby.
            if(!track_too_close) {

                fake_clusters++;
                hFakeClusterCharge->Fill(cluster->charge());
                hFakeClusterSize->Fill(static_cast<double>(cluster->size()));

                for(auto& pixel : cluster->pixels()) {
                    fake_hits++;
                    hFakePixelCharge->Fill(pixel->charge());
                    fakePixelPerEventMap->Fill(pixel->column(), pixel->row(), 1);
                }
            }
        }

        hFakePixelPerEvent->Fill(fake_hits);
        fakePixelPerEventVsTime->Fill(static_cast<double>(Units::convert(event->start(), "s")), fake_hits);
        fakePixelPerEventVsTimeLong->Fill(static_cast<double>(Units::convert(event->start(), "s")), fake_hits);
        hFakeClusterPerEvent->Fill(fake_clusters);
    }
    if(m_fake_rate_method == FakeRateMethod::EDGE) {
        LOG_ONCE(STATUS) << "Estimating fake rate based on events without DUT intercepting tracks (EDGE method).";

        bool track_in_active = false;
        int fake_hits = 0;
        int fake_clusters = 0;

        // iterate the tracks from the clipboard
        for(auto& track : tracks) {
            // check if one of the tracks intercepts the dut
            // A positive value of m_fake_rate_radius means that we want to
            // check an area that is larger than the sensor matrix.
            if(m_detector->hasIntercept(track.get(), -m_fake_rate_distance)) {
                track_in_active = true;
                break;
            }
        }

        // Study the dut response if there is no track pointing to the active
        // area of the dut. There might still be particles, that we failed to
        // track!
        if(!track_in_active) {

            // iterate the dut pixels from clipboard
            for(auto& pixel : pixels) {
                fake_hits++;
                hFakePixelCharge->Fill(pixel->charge());
                fakePixelPerEventMap->Fill(pixel->column(), pixel->row(), 1);
            }
            hFakePixelPerEvent->Fill(fake_hits);
            fakePixelPerEventVsTime->Fill(static_cast<double>(Units::convert(event->start(), "s")), fake_hits);
            fakePixelPerEventVsTimeLong->Fill(static_cast<double>(Units::convert(event->start(), "s")), fake_hits);

            // get and iterate dut clusters from clipboard
            auto clusters = clipboard->getData<Cluster>(m_detector->getName());
            for(auto& cluster : clusters) {
                fake_clusters++;
                hFakeClusterCharge->Fill(cluster->charge());
                hFakeClusterSize->Fill(static_cast<double>(cluster->size()));
            }
            hFakeClusterPerEvent->Fill(fake_clusters);
        }
    }

    return StatusCode::Success;
}

void AnalysisEfficiency::finalize(const std::shared_ptr<ReadonlyClipboard>&) {
    // Track selection flow:
    LOG(STATUS) << "Track selection flow:       " << n_track << std::endl
                << "* rejected by chi2          -" << n_chi2 << std::endl
                << "* track outside ROI         -" << n_roi << std::endl
                << "* track outside DUT         -" << n_dut << std::endl
                << "* track close to masked px  -" << n_masked << std::endl
                << "* track close to frame edge -" << n_frameedge << std::endl
                << "* track without an associated cluster on required detector - " << n_requirecluster << std::endl
                << "Accepted tracks:            " << total_tracks;

    double totalEff = 100 * static_cast<double>(matched_tracks) / (total_tracks > 0 ? total_tracks : 1);
    double lowerEffError = totalEff - 100 * (TEfficiency::ClopperPearson(total_tracks, matched_tracks, 0.683, false));
    double upperEffError = 100 * (TEfficiency::ClopperPearson(total_tracks, matched_tracks, 0.683, true)) - totalEff;
    LOG(STATUS) << "Total efficiency of detector " << m_detector->getName() << ": " << totalEff << "(+" << upperEffError
                << " -" << lowerEffError << ")%, measured with " << matched_tracks << "/" << total_tracks
                << " matched/total tracks";

    for(int icol = 1; icol < m_detector->nPixels().X() + 1; icol++) {
        for(int irow = 1; irow < m_detector->nPixels().Y() + 1; irow++) {
            // calculate total efficiency: (just to double check the other calculation)
            const int bin = hChipEfficiencyMap_trackPos->GetGlobalBin(icol, irow);
            double eff = hChipEfficiencyMap_trackPos->GetEfficiency(bin);
            if(eff > 0) {
                LOG(TRACE) << "col/row = " << icol << "/" << irow << ", binContent = " << eff;
                hPixelEfficiency->Fill(eff);
            }
            eff = hPixelEfficiencyMatrix_TProfile->GetBinContent(bin);
            if(eff > 0) {
                LOG(TRACE) << "col/row = " << icol << "/" << irow << ", binContent = " << eff;
                hPixelEfficiencyMatrix->Fill(eff);
            }
        }
    }
    // normalize fake rate map, if it exists
    if(hFakeClusterPerEvent->GetEntries() > 0) {
        fakePixelPerEventMap->Scale(1. / hFakePixelPerEvent->GetEntries());
    }
}
