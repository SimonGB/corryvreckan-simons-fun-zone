/**
 * @file
 * @brief Implementation of module DUTAssociation
 *
 * @copyright Copyright (c) 2017-2022 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "DUTAssociation.h"
#include "tools/cuts.h"

using namespace corryvreckan;
using namespace std;

DUTAssociation::DUTAssociation(Configuration& config, std::shared_ptr<Detector> detector)
    : Module(config, detector), m_detector(detector) {

    // Backwards compatibility: also allow timing_cut to be used for time_cut_abs and spatial_cut for spatial_cut_abs
    config_.setAlias("time_cut_abs", "timing_cut", true);
    config_.setAlias("spatial_cut_abs", "spatial_cut", true);

    config_.setDefault<bool>("use_cluster_centre", false);

    if(config_.count({"time_cut_rel", "time_cut_abs"}) == 0) {
        config_.setDefault("time_cut_rel", 3.0);
    }
    if(config_.count({"spatial_cut_rel", "spatial_cut_abs"}) == 0) {
        config_.setDefault("spatial_cut_rel", 3.0);
    }

    // timing cut, relative (x * time_resolution) or absolute:
    time_cut_ = corryvreckan::calculate_cut<double>("time_cut", config_, m_detector);

    // spatial cut, relative (x * spatial_resolution) or absolute:
    spatial_cut_ = corryvreckan::calculate_cut<XYVector>("spatial_cut", config_, m_detector);
    use_cluster_centre_ = config_.get<bool>("use_cluster_centre");
    elliptic_cut_ = config_.get<bool>("elliptic_cut", true);

    LOG(DEBUG) << "time_cut = " << Units::display(time_cut_, {"ms", "us", "ns"});
    LOG(DEBUG) << "spatial_cut = " << Units::display(spatial_cut_, {"um", "mm"});
    LOG(DEBUG) << "use_cluster_centre = " << use_cluster_centre_;
}

void DUTAssociation::initialize() {
    // Cut flow histogram
    std::string title = m_detector->getName() + ": number of tracks discarded by different cuts;cut type;clusters";
    hCutHisto = new TH1F("hCutHisto", title.c_str(), 2, 1, 3);
    hCutHisto->GetXaxis()->SetBinLabel(1, "Spatial");
    hCutHisto->GetXaxis()->SetBinLabel(2, "Timing");

    hDistX = new TH1D("hDistXClusterClosestPx",
                      "Distance cluster center to pixel closest to track; x_{cluster} - x_{closest pixel} [um]; # events",
                      2000,
                      -1000,
                      1000);
    hDistY = new TH1D("hDistYClusterClosestPx",
                      "Distance cluster center to pixel closest to track; y_{cluster} - y_{closest pixel} [um]; # events",
                      2000,
                      -1000,
                      1000);
    hDistX_1px =
        new TH1D("hDistXClusterClosestPx_1px",
                 "Distance 1px-cluster center to pixel closest to track; x_{cluster} - x_{closest pixel} [um]; # events",
                 2000,
                 -1000,
                 1000);
    hDistY_1px =
        new TH1D("hDistYClusterClosestPx_1px",
                 "Distance 1px-cluster center to pixel closest to track; y_{cluster} - y_{closest pixel} [um]; # events",
                 2000,
                 -1000,
                 1000);
    hDistX_2px =
        new TH1D("hDistXClusterClosestPx_2px",
                 "Distance 2px-cluster center to pixel closest to track; x_{cluster} - x_{closest pixel} [um]; # events",
                 2000,
                 -1000,
                 1000);
    hDistY_2px =
        new TH1D("hDistYClusterClosestPx_2px",
                 "Distance 2px-cluster center to pixel closest to track; y_{cluster} - y_{closest pixel} [um]; # events",
                 2000,
                 -1000,
                 1000);
    hDistX_3px =
        new TH1D("hDistXClusterClosestPx_3px",
                 "Distance 3px-cluster center to pixel closest to track; x_{cluster} - x_{closest pixel} [um]; # events",
                 2000,
                 -1000,
                 1000);
    hDistY_3px =
        new TH1D("hDistYClusterClosestPx_3px",
                 "Distance 3px-cluster center to pixel closest to track; y_{cluster} - y_{closest pixel} [um]; # events",
                 2000,
                 -1000,
                 1000);

    hDist2D = new TH2D("hDist_trackCluster_2D",
                       "2D distance between cluster center and track; dx [mm]; dy [mm] ",
                       500,
                       -m_detector->getSize().X(),
                       m_detector->getSize().X(),
                       500,
                       -m_detector->getSize().X(),
                       m_detector->getSize().X());
    // Nr of associated clusters per track
    title = m_detector->getName() + ": number of associated clusters per track;associated clusters;events";
    hNoAssocCls = new TH1F("no_assoc_cls", title.c_str(), 10, -0.5, 9.5);
    LOG(DEBUG) << "DUT association time cut = " << Units::display(time_cut_, {"ms", "ns"});

    // Add additional histograms for polar detectors
    if(m_detector->is<PolarDetector>()) {
        hDistR =
            new TH1D("hDistRClusterClosest",
                     "Distance cluster center to pixel closest to track; r_{cluster} - r_{closest pixel} [mm]; # events",
                     500,
                     -50,
                     50);

        hDistPhi = new TH1D("hDistPhiClusterClosest",
                            "Distance cluster center to pixel closest to track; #varphi_{cluster} - #varphi_{closest pixel} "
                            "[mrad]; # events",
                            500,
                            -150,
                            150);

        hAssocDistR = new TH1D("hAssocDistR",
                               "Distance to associated cluster; r_{cluster} - r_{pixel} [mm]; # events",
                               500,
                               -spatial_cut_.y(),
                               spatial_cut_.y());

        hAssocDistPhi = new TH1D("hAssocDistPhi",
                                 "Distance to associated cluster; #varphi_{cluster} - #varphi_{pixel} [#murad]; # events",
                                 500,
                                 -static_cast<double>(Units::convert(spatial_cut_.x(), "urad")),
                                 static_cast<double>(Units::convert(spatial_cut_.x(), "urad")));
    }
}

StatusCode DUTAssociation::run(const std::shared_ptr<Clipboard>& clipboard) {

    // Get the tracks from the clipboard
    auto tracks = clipboard->getData<Track>();
    // Get the DUT clusters from the clipboard
    auto clusters = clipboard->getData<Cluster>(m_detector->getName());

    // Loop over all tracks
    for(auto& track : tracks) {
        total_tracks_++;
        LOG(TRACE) << "Processing track with model " << track->getType() << ", chi2 of " << track->getChi2();
        int assoc_cls_per_track = 0;
        auto min_distance = std::numeric_limits<double>::max();

        if(clusters.empty()) {
            hNoAssocCls->Fill(0);
            LOG(DEBUG) << "No DUT clusters on the clipboard";
            continue;
        } else {
            LOG(TRACE) << "DUT clusters found on clipboard";
        }

        // Check distance between track and cluster
        auto interceptLocal = m_detector->getLocalIntercept(track.get());
        LOG(TRACE) << "Starting loop over all clusters";
        // Loop over all DUT clusters
        for(auto& cluster : clusters) {

            // distance of track to cluster centre
            double xdistance_centre = std::abs(interceptLocal.X() - cluster->local().x());
            double ydistance_centre = std::abs(interceptLocal.Y() - cluster->local().y());

            // distance of track to nearest pixel: initialise to maximal possible value
            auto xdistance_nearest = std::numeric_limits<double>::max();
            auto ydistance_nearest = std::numeric_limits<double>::max();
            LOG(TRACE) << "Cluster distance to track intercept " << xdistance_centre << "," << ydistance_centre;

            for(auto& pixel : cluster->pixels()) {
                // convert pixel address to local coordinates:
                auto pixelPositionLocal =
                    m_detector->getLocalPosition(static_cast<double>(pixel->column()), static_cast<double>(pixel->row()));

                xdistance_nearest = std::min(xdistance_nearest, std::abs(interceptLocal.X() - pixelPositionLocal.x()));
                ydistance_nearest = std::min(ydistance_nearest, std::abs(interceptLocal.Y() - pixelPositionLocal.y()));

                // Recalculate distances for polar detectors
                if(m_detector->is<PolarDetector>()) {
                    auto polar_det = std::dynamic_pointer_cast<PolarDetector>(m_detector);

                    // Get polar coordinates of cluster, intercept and strip
                    auto cluster_polar = polar_det->getPolarPosition(cluster->local());
                    auto intercept_polar = polar_det->getPolarPosition(interceptLocal);
                    auto strip_polar = polar_det->getPolarPosition(pixelPositionLocal);

                    // Recalculate distance to cluster centre
                    xdistance_centre = intercept_polar.phi() - cluster_polar.phi();
                    ydistance_centre = intercept_polar.r() - cluster_polar.r();

                    // Recalculate distance to nearest strip
                    xdistance_nearest = std::numeric_limits<double>::max();
                    ydistance_nearest = std::numeric_limits<double>::max();
                    xdistance_nearest = std::min(xdistance_nearest, std::abs(intercept_polar.phi() - strip_polar.phi()));
                    ydistance_nearest = std::min(ydistance_nearest, std::abs(intercept_polar.r() - strip_polar.r()));
                }
            }

            hDistX->Fill(static_cast<double>(Units::convert(xdistance_centre - xdistance_nearest, "um")));
            hDistY->Fill(static_cast<double>(Units::convert(ydistance_centre - ydistance_nearest, "um")));
            hDist2D->Fill(interceptLocal.X() - cluster->local().x(), interceptLocal.y() - cluster->local().y());
            if(cluster->columnWidth() == 1) {
                hDistX_1px->Fill(static_cast<double>(Units::convert(xdistance_centre - xdistance_nearest, "um")));
            }
            if(cluster->rowWidth() == 1) {
                hDistY_1px->Fill(static_cast<double>(Units::convert(ydistance_centre - ydistance_nearest, "um")));
            }
            if(cluster->columnWidth() == 2) {
                hDistX_2px->Fill(static_cast<double>(Units::convert(xdistance_centre - xdistance_nearest, "um")));
            }
            if(cluster->rowWidth() == 2) {
                hDistY_2px->Fill(static_cast<double>(Units::convert(ydistance_centre - ydistance_nearest, "um")));
            }
            if(cluster->columnWidth() == 3) {
                hDistX_3px->Fill(static_cast<double>(Units::convert(xdistance_centre - xdistance_nearest, "um")));
            }
            if(cluster->rowWidth() == 3) {
                hDistY_3px->Fill(static_cast<double>(Units::convert(ydistance_centre - ydistance_nearest, "um")));
            }

            // Check if the cluster is close in space (either use cluster centre of closest pixel to track)
            auto xdistance = (use_cluster_centre_ ? xdistance_centre : xdistance_nearest);
            auto ydistance = (use_cluster_centre_ ? ydistance_centre : ydistance_nearest);
            auto distance = sqrt(xdistance * xdistance + ydistance * ydistance);

            // Evaluate distances based on requested type of cut (elliptic or rectangular)
            double norm;
            if(elliptic_cut_) {
                // Check if track-cluster distance lies within ellipse defined by spatial cuts, following this example:
                // https://www.geeksforgeeks.org/check-if-a-point-is-inside-outside-or-on-the-ellipse/
                //
                // ellipse defined by: x^2/a^2 + y^2/b^2 = 1: on ellipse,
                //                                       > 1: outside,
                //                                       < 1: inside
                // Discard track if outside of ellipse
                norm = (xdistance * xdistance) / (spatial_cut_.x() * spatial_cut_.x()) +
                       (ydistance * ydistance) / (spatial_cut_.y() * spatial_cut_.y());
            } else {
                // Use rectangular cut
                norm = std::max((xdistance * xdistance) / (spatial_cut_.x() * spatial_cut_.x()),
                                (ydistance * ydistance) / (spatial_cut_.y() * spatial_cut_.y()));
            }

            // Fill distance histograms for polar detectors
            if(m_detector->is<PolarDetector>()) {
                hDistPhi->Fill(static_cast<double>(Units::convert(xdistance, "mrad")));
                hDistR->Fill(static_cast<double>(Units::convert(ydistance, "mm")));
            }

            if(norm > 1) {
                LOG(DEBUG) << "Discarding DUT cluster with distance (" << Units::display(std::abs(xdistance), {"um", "mm"})
                           << "," << Units::display(std::abs(ydistance), {"um", "mm"}) << ")"
                           << " with local track intersection at " << Units::display(interceptLocal, {"um", "mm"});
                hCutHisto->Fill(1);
                num_cluster++;
                continue;
            }

            // Check if the cluster is close in time
            if(std::abs(cluster->timestamp() - track->timestamp()) > time_cut_) {
                LOG(DEBUG) << "Discarding DUT cluster with time difference "
                           << Units::display(std::abs(cluster->timestamp() - track->timestamp()), {"ms", "s"});
                hCutHisto->Fill(2);
                num_cluster++;
                continue;
            }

            LOG(DEBUG) << "Found associated cluster with distance (" << Units::display(abs(xdistance), {"um", "mm"}) << ","
                       << Units::display(abs(ydistance), {"um", "mm"}) << ")";
            track->addAssociatedCluster(cluster.get());
            assoc_cls_per_track++;
            assoc_cluster_counter++;
            num_cluster++;

            // Fill distance to assoc. cluster for polar detectors
            if(m_detector->is<PolarDetector>()) {
                hAssocDistPhi->Fill(static_cast<double>(Units::convert(xdistance, "urad")));
                hAssocDistR->Fill(static_cast<double>(Units::convert(ydistance, "mm")));
            }

            // check if cluster is closest to track
            if(distance < min_distance) {
                min_distance = distance;
                track->setClosestCluster(cluster.get());
            }
        }
        hNoAssocCls->Fill(assoc_cls_per_track);
        if(assoc_cls_per_track > 0) {
            track_w_assoc_cls++;
        }
    }

    // Return value telling analysis to keep running
    return StatusCode::Success;
}

void DUTAssociation::finalize(const std::shared_ptr<ReadonlyClipboard>&) {
    hCutHisto->Scale(1 / double(num_cluster));
    LOG(STATUS) << "In total, " << assoc_cluster_counter << " clusters are associated to " << track_w_assoc_cls
                << " tracks.";
    LOG(INFO) << "Number of tracks with at least one associated cluster: " << track_w_assoc_cls
              << " vs total number of tracks: " << total_tracks_;
    return;
}
