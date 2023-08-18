/**
 * @file
 * @brief Implementation of module AnalysisTiming
 *
 * @copyright Copyright (c) 2023 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "AnalysisTiming.h"

#include "core/config/exceptions.h"
#include "objects/Cluster.hpp"
#include "objects/Event.hpp"
#include "objects/Pixel.hpp"
#include "objects/Track.hpp"

using namespace corryvreckan;

AnalysisTiming::AnalysisTiming(Configuration& config, std::shared_ptr<Detector> detector)
    : Module(config, detector), detector_(std::move(detector)), detector_name_(detector_->getName()) {

    // Determine whether to use associated clusters or not for DUT, set function pointer
    const auto role = detector_->getRoles();
    switch(role) {
    case DetectorRole::DUT:
        dut_cluster_func_ = &AnalysisTiming::get_dut_cluster_associated;
        LOG(INFO) << "Using closest associated clusters from " << detector_name_;
        break;
    case DetectorRole::NONE:
    case DetectorRole::REFERENCE:
        dut_cluster_func_ = &AnalysisTiming::get_dut_cluster_tracking;
        LOG(INFO) << "Using tracking clusters from " << detector_name_;
        break;
    default:
        throw InvalidDetectorRoleError(detector_.get());
        break;
    }

    // Get reference type, defaults to DUT type
    config_.setDefault<TimestampReferenceType>("reference_type", TimestampReferenceType::DUT);
    auto reference_type = config_.get<TimestampReferenceType>("reference_type");

    // Get reference name from config if not timestamps from track
    if(reference_type == TimestampReferenceType::TRACK) {
        if(config_.has("reference_name")) {
            LOG(WARNING) << "Parameter \"reference_name\" is specified but unused";
        }
    } else {
        reference_name_ = config_.get<std::string>("reference_name");
    }

    // Set function pointer for reference timestamps
    switch(reference_type) {
    case TimestampReferenceType::DUT:
        ref_timestamp_func_ = &AnalysisTiming::get_ref_timestamp_dut;
        LOG(INFO) << "Using closest associated clusters from " << reference_name_;
        break;
    case TimestampReferenceType::PLANE:
        ref_timestamp_func_ = &AnalysisTiming::get_ref_timestamp_plane;
        LOG(INFO) << "Using tracking clusters from " << reference_name_;
        break;
    case TimestampReferenceType::TRACK:
        ref_timestamp_func_ = &AnalysisTiming::get_ref_timestamp_track;
        LOG(INFO) << "Using track timestamps as reference";
        break;
    default:
        // Logic Error
        break;
    }

    // Configure cuts
    config_.setDefault<double>("chi2ndof_cut", 3.);
    chi2_ndof_cut_ = config_.get<double>("chi2ndof_cut");

    // Timing histograms
    time_range_ = config_.get<double>("time_range");
    time_binning_ = config_.get<double>("time_binning");
    config_.setDefault("time_offset", 0.);
    time_offset_ = config_.get<double>("time_offset");
    LOG(INFO) << "Histograms have a range of " << Units::display(time_range_, {"ps", "ns", "us"}) << " with offset "
              << Units::display(time_offset_, {"ps", "ns", "us"}) << " and binning "
              << Units::display(time_binning_, {"ps", "ns", "us"});

    // 2D histograms
    config_.setDefault<ROOT::Math::XYPoint>("inpixel_bin_size", {Units::get(1.0, "um"), Units::get(1.0, "um")});
    if(config_.getArray<double>("inpixel_bin_size").size() == 2) {
        inpixel_bin_size_ = config_.get<ROOT::Math::XYPoint>("inpixel_bin_size");
    } else {
        const auto binsize = config_.get<double>("inpixel_bin_size");
        inpixel_bin_size_ = ROOT::Math::XYPoint(binsize, binsize);
    }
}

void AnalysisTiming::initialize() {
    // Basic histograms
    const auto time_bins = static_cast<int>(time_range_ / time_binning_);
    const auto time_xlow = time_offset_ - 0.5 * time_range_;
    const auto time_xup = time_offset_ + 0.5 * time_range_;
    hTimeResidual_ =
        new TH1D("residualTime", "Time residual;time_{ref}-time_{dut} [ns];# entries", time_bins, time_xlow, time_xup);
    hTimeResidualOverTime_ = new TH2F("residualTimeOverTime",
                                      "Time residual over time;time_{ref} [s];time_{ref}-time_{dut} [ns];# entries",
                                      600,
                                      -3.,
                                      3600. - 3.,
                                      time_bins,
                                      time_xlow,
                                      time_xup);

    // 2D sensor histograms
    const auto npix_x = detector_->nPixels().X();
    const auto npix_y = detector_->nPixels().Y();
    hResidualMeanSensor_ = new TProfile2D("residualTime_mean_sensor",
                                          "Mean time residual sensor map;x [px];y [px];time [ns]",
                                          npix_x,
                                          -0.5,
                                          npix_x - 0.5,
                                          npix_y,
                                          -0.5,
                                          npix_y - 0.5,
                                          "s");
    hResidualStdDevSensor_ = new TProfile2D("residualTime_stddev_sensor",
                                            "Standard deviation time residual sensor map;x [px];y [px];time [ns]",
                                            npix_x,
                                            -0.5,
                                            npix_x - 0.5,
                                            npix_y,
                                            -0.5,
                                            npix_y - 0.5);

    // 2D inpixel histograms
    const auto pitch_x = detector_->getPitch().X();
    const auto pitch_y = detector_->getPitch().Y();
    const auto nbins_inpix_x = static_cast<int>(std::ceil(pitch_x / inpixel_bin_size_.X()));
    const auto nbins_inpix_y = static_cast<int>(std::ceil(pitch_y / inpixel_bin_size_.Y()));
    if(nbins_inpix_x > 1e4 || nbins_inpix_y > 1e4) {
        throw InvalidValueError(config_, "inpixel_bin_size", "Too many bins for in-pixel histograms.");
    }
    hResidualMeanInpix_ = new TProfile2D("residualTime_mean_inpix",
                                         "Mean time residual in-pixel map;x [mm];y [mm];time [ns]",
                                         nbins_inpix_x,
                                         -0.5 * pitch_x,
                                         0.5 * pitch_x,
                                         nbins_inpix_y,
                                         -0.5 * pitch_y,
                                         0.5 * pitch_y,
                                         "s");
    hResidualStdDevInpix_ = new TProfile2D("residualTime_stddev_inpix",
                                           "Standard deviation time residual in-pixel map;x [mm];y [mm];time [ns]",
                                           nbins_inpix_x,
                                           -0.5 * pitch_x,
                                           0.5 * pitch_x,
                                           nbins_inpix_y,
                                           -0.5 * pitch_y,
                                           0.5 * pitch_y);

    // Cut histogram
    hCutHisto_ = new TH1F("cutHisto",
                          "number of tracks discarded by different cuts;cut type;tracks",
                          ETrackSelection::kNSelection,
                          0,
                          ETrackSelection::kNSelection);
    hCutHisto_->GetXaxis()->SetBinLabel(1 + ETrackSelection::kAllTrack, "All tracks");
    hCutHisto_->GetXaxis()->SetBinLabel(1 + ETrackSelection::kPassedChi2Ndf, "Passed Chi2/ndof");
    hCutHisto_->GetXaxis()->SetBinLabel(1 + ETrackSelection::kTimestampOnRef, "Timestamp on ref");
    hCutHisto_->GetXaxis()->SetBinLabel(1 + ETrackSelection::kClusterOnDUT, "Cluster on DUT");
}

StatusCode AnalysisTiming::run(const std::shared_ptr<Clipboard>& clipboard) {
    // Get the telescope tracks from the clipboard
    const auto& tracks = clipboard->getData<Track>();

    // Loop over all tracks
    for(const auto& track : tracks) {
        hCutHisto_->Fill(ETrackSelection::kAllTrack);

        // Cut on Chi2/ndof
        if(track->getChi2ndof() > chi2_ndof_cut_) {
            LOG(TRACE) << "Track discarded due to Chi2/ndof";
            continue;
        }
        hCutHisto_->Fill(ETrackSelection::kPassedChi2Ndf);

        // Find cluster for the reference
        const auto timestamp_ref_opt = (this->*ref_timestamp_func_)(track.get());
        if(!timestamp_ref_opt.has_value()) {
            LOG(TRACE) << "Reference has no timestamp, skipping track";
            continue;
        }
        hCutHisto_->Fill(ETrackSelection::kTimestampOnRef);

        // Find associated cluster on the DUT
        const auto* cluster_dut = (this->*dut_cluster_func_)(track.get());
        if(cluster_dut == nullptr) {
            LOG(TRACE) << "DUT has no cluster, skipping track";
            continue;
        }
        hCutHisto_->Fill(ETrackSelection::kClusterOnDUT);

        // Get time residual
        const auto timestamp_ref = timestamp_ref_opt.value();
        const auto time_residual = timestamp_ref - cluster_dut->timestamp();
        LOG(DEBUG) << "Found time residual " << Units::display(time_residual, {"ps", "ns", "us", "ms"})
                   << " at reference time " << Units::display(timestamp_ref, {"ps", "ns", "us", "ms", "s"});

        // Fill basic histograms
        hTimeResidual_->Fill(time_residual);
        hTimeResidualOverTime_->Fill(static_cast<double>(Units::convert(timestamp_ref, "s")), time_residual);

        // Fill 2D histograms
        auto intercept_local = detector_->getLocalIntercept(track.get());
        auto intercept_pixel = detector_->getInterceptPixel(intercept_local);
        auto intercept_inpix = detector_->inPixel(intercept_local);
        hResidualMeanSensor_->Fill(intercept_pixel.first, intercept_pixel.second, time_residual);
        hResidualMeanInpix_->Fill(intercept_inpix.X(), intercept_inpix.Y(), time_residual);
    }

    // Return value telling analysis to keep running
    return StatusCode::Success;
}

void AnalysisTiming::finalize(const std::shared_ptr<ReadonlyClipboard>&) {
    LOG(INFO) << "Found time residuals in " << hCutHisto_->GetBinContent(1 + ETrackSelection::kClusterOnDUT) << " out of "
              << hCutHisto_->GetBinContent(1 + ETrackSelection::kAllTrack) << " tracks";

    // Fill 2D sensor standard deviation histogram
    for(Int_t nx = 1; nx <= hResidualMeanSensor_->GetNbinsX(); ++nx) {
        for(Int_t ny = 1; ny <= hResidualMeanSensor_->GetNbinsY(); ++ny) {
            const auto std_dev = hResidualMeanSensor_->GetBinError(nx, ny);
            if(std_dev > 0.) {
                hResidualStdDevSensor_->Fill(static_cast<double>(nx - 1), static_cast<double>(ny - 1), std_dev);
            }
        }
    }
    // Fill 2D inpixel standard deviation histogram
    const auto pitch_x = detector_->getPitch().X();
    const auto pitch_y = detector_->getPitch().Y();
    const auto bin_width_x = pitch_x / hResidualMeanInpix_->GetNbinsX();
    const auto bin_width_y = pitch_y / hResidualMeanInpix_->GetNbinsY();
    for(Int_t nx = 1; nx <= hResidualMeanInpix_->GetNbinsX(); ++nx) {
        for(Int_t ny = 1; ny <= hResidualMeanInpix_->GetNbinsY(); ++ny) {
            const auto std_dev = hResidualMeanInpix_->GetBinError(nx, ny);
            if(std_dev > 0.) {
                const auto x = -0.5 * pitch_x + (nx - 0.5) * bin_width_x;
                const auto y = -0.5 * pitch_y + (ny - 0.5) * bin_width_y;
                hResidualStdDevInpix_->Fill(x, y, std_dev);
            }
        }
    }
}

Cluster* AnalysisTiming::get_dut_cluster_associated(const Track* track) const {
    if(track->getAssociatedClusters(detector_name_).size() > 0) {
        return track->getClosestCluster(detector_name_);
    }
    return nullptr;
}

Cluster* AnalysisTiming::get_dut_cluster_tracking(const Track* track) const {
    return track->getClusterFromDetector(detector_name_);
}

std::optional<double> AnalysisTiming::get_ref_timestamp_dut(const Track* track) const {
    if(track->getAssociatedClusters(reference_name_).size() > 0) {
        // getClosestCluster never returns nullptr
        return track->getClosestCluster(reference_name_)->timestamp();
    }
    return std::nullopt;
}

std::optional<double> AnalysisTiming::get_ref_timestamp_plane(const Track* track) const {
    const auto* cluster_ref = track->getClusterFromDetector(reference_name_);
    if(cluster_ref != nullptr) {
        return cluster_ref->timestamp();
    }
    return std::nullopt;
}

std::optional<double> AnalysisTiming::get_ref_timestamp_track(const Track* track) const {
    const auto track_timestamp = track->timestamp();
    return (track_timestamp == 0.) ? std::nullopt : std::optional(track_timestamp);
}
