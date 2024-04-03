/**
 * @file
 * @brief Implementation of Multiplet base object
 *
 * @copyright Copyright (c) 2017-2022 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "Multiplet.hpp"
#include "TMath.h"
#include "core/utils/log.h"
#include "core/utils/unit.h"
#include "exceptions.h"

using namespace corryvreckan;
Multiplet::Multiplet(std::shared_ptr<Track> upstream, std::shared_ptr<Track> downstream) {
    m_upstream = std::move(upstream);
    m_downstream = std::move(downstream);

    // All clusters from up- and downstream should be referenced from this track:
    for(auto& cluster : m_upstream->getClusters()) {
        this->addCluster(cluster);
    }
    for(auto& cluster : m_downstream->getClusters()) {
        this->addCluster(cluster);
    }
    // if it's listed as plane for upstream or downstream, it should be a plane for the multiplet too
    for(auto& upp : m_upstream->getPlanes()) {
        planes_.push_back(std::move(upp));
    }
    for(auto& dop : m_downstream->getPlanes()) {
        auto pl = std::find_if(
            planes_.begin(), planes_.end(), [&dop](const Plane& plane) { return plane.getName() == dop.getName(); });
        if(pl == planes_.end()) {
            planes_.push_back(std::move(dop));
        } else {
            *pl = std::move(dop);
        }
    }
}

ROOT::Math::XYPoint Multiplet::getKinkAt(const std::string&) const { return ROOT::Math::XYPoint(0, 0); }

void Multiplet::calculateChi2() {

    chi2_ = m_upstream->getChi2() + m_downstream->getChi2();
    ndof_ = m_upstream->getNdof() + m_downstream->getNdof();
    chi2ndof_ = (ndof_ <= 0) ? -1 : (chi2_ / static_cast<double>(ndof_));
}

void Multiplet::calculateResiduals() {
    for(const auto& c : track_clusters_) {
        auto* cluster = c.get();
        residual_global_[cluster->detectorID()] = cluster->global() - getIntercept(cluster->global().z());
        if(get_plane(cluster->detectorID()) != nullptr) {
            residual_local_[cluster->detectorID()] =
                cluster->local() - get_plane(cluster->detectorID())->getToLocal() * getIntercept(cluster->global().z());
        }
    }
}

void Multiplet::fit() {

    // tracks
    m_positionAtScatterer = ((m_downstream->getIntercept(m_scattererPosition) -
                              (ROOT::Math::XYZPoint(0, 0, 0) - m_upstream->getIntercept(m_scattererPosition))) /
                             2.);
    m_offsetAtScatterer = m_downstream->getIntercept(m_scattererPosition) - m_upstream->getIntercept(m_scattererPosition);

    // Calculate the angle
    ROOT::Math::XYZVector slopeUp = m_upstream->getDirection(m_scattererPosition);
    ROOT::Math::XYZVector slopeDown = m_downstream->getDirection(m_scattererPosition);
    //
    ROOT::Math::XYZVector kinks = (slopeDown /= slopeDown.z()) - (slopeUp /= slopeUp.z());
    m_kinkAtScatterer = ROOT::Math::XYVector(kinks.x(), kinks.y());

    this->calculateChi2();
    this->calculateResiduals();
    isFitted_ = true;
}

ROOT::Math::XYZPoint Multiplet::getIntercept(double z) const {
    return z == m_scattererPosition
               ? m_positionAtScatterer
               : (z < m_scattererPosition ? m_upstream->getIntercept(z) : m_downstream->getIntercept(z));
}

ROOT::Math::XYZPoint Multiplet::getState(const std::string& detectorID) const {
    if(!isFitted_) {
        throw TrackError(typeid(*this), " not fitted");
    }

    LOG(TRACE) << "Planes known to this track: ";
    for(auto pKnown : planes_) {
        LOG(TRACE) << " - " << pKnown.getName();
    }

    auto plane =
        std::find_if(planes_.begin(), planes_.end(), [&detectorID](Plane const& p) { return p.getName() == detectorID; });
    if(plane == planes_.end()) {
        throw TrackError(typeid(*this), " does not have any entry for plane " + detectorID);
    }
    auto zpos = plane->getPosition();
    LOG(TRACE) << "Plane z position of " << detectorID << " is " << zpos << ", scatterer position is "
               << m_scattererPosition;

    LOG(TRACE) << "upstream track type " << m_upstream->getType() << ", downstream track type " << m_downstream->getType();
    if(zpos <= m_scattererPosition) {
        return m_upstream->getState(detectorID);
    } else {
        return m_downstream->getState(detectorID);
    }
}

ROOT::Math::XYZVector Multiplet::getDirection(const std::string& detectorID) const {
    return getClusterFromDetector(detectorID)->global().z() <= m_scattererPosition ? m_upstream->getDirection(detectorID)
                                                                                   : m_downstream->getDirection(detectorID);
}

XYZVector Multiplet::getDirection(const double& z) const {
    return (z <= m_scattererPosition ? m_upstream->getDirection(z) : m_downstream->getDirection(z));
}

void Multiplet::print(std::ostream& out) const {
    out << "Multiplet " << this->m_scattererPosition << ", " << this->m_positionAtScatterer << ", "
        << this->m_offsetAtScatterer << ", " << this->m_kinkAtScatterer << ", " << this->chi2_ << ", " << this->ndof_ << ", "
        << this->chi2ndof_ << ", " << this->timestamp();
}
