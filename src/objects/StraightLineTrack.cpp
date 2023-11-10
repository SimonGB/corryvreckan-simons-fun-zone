/**
 * @file
 * @brief Implementation of StraightLine track object
 *
 * @copyright Copyright (c) 2017-2022 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "StraightLineTrack.hpp"
#include "Eigen/Dense"
#include "Track.hpp"
#include "core/utils/log.h"
#include "exceptions.h"

using namespace corryvreckan;

ROOT::Math::XYPoint StraightLineTrack::distance(const Cluster* cluster) const {

    if(get_plane(cluster->detectorID()) == nullptr) {
        throw MissingReferenceException(typeid(*this), typeid(Plane));
    }
    auto trackIntercept = get_plane(cluster->detectorID())->getToLocal() * getState(cluster->detectorID());

    auto dist = cluster->local() - trackIntercept;

    // Return the distance^2
    return ROOT::Math::XYPoint(dist.x(), dist.y());
}

ROOT::Math::XYPoint StraightLineTrack::getKinkAt(const std::string&) const {
    return ROOT::Math::XYPoint(0, 0);
}

ROOT::Math::XYZPoint StraightLineTrack::getState(const std::string& detectorID) const {
    LOG(TRACE) << "Requesting state at: " << detectorID;
    auto plane =
        std::find_if(planes_.begin(), planes_.end(), [&detectorID](Plane const& p) { return p.getName() == detectorID; });
    if(plane == planes_.end()) {
        throw MissingReferenceException(typeid(*this), typeid(Plane));
    }
    auto toGlobal = plane->getToGlobal();
    ROOT::Math::XYZVector planeU, planeV, planeN;
    toGlobal.Rotation().GetComponents(planeU, planeV, planeN);
    ROOT::Math::XYZPoint origin = toGlobal.Translation() * ROOT::Math::XYZPoint(0, 0, 0);

    ROOT::Math::XYZVector distance = m_state - origin;
    double pathLength = -distance.Dot(planeN) / m_direction.Dot(planeN);
    ROOT::Math::XYZPoint position = m_state + pathLength * m_direction;

    return position;
}

ROOT::Math::XYZVector StraightLineTrack::getDirection(const std::string&) const {
    return m_direction;
}

ROOT::Math::XYZVector StraightLineTrack::getDirection(const double&) const {
    return m_direction;
}

void StraightLineTrack::calculateChi2() {

    // Get the number of clusters
    // We do have a 2-dimensional offset(x_0,y_0) and slope (dx,dy). Each hit provides two measurements.
    // ndof_ = 2*num_planes - 4 = 2 * (num_planes -2)
    ndof_ = (track_clusters_.size() - 2) * 2;
    chi2_ = 0.;
    chi2ndof_ = 0.;

    // Loop over all clusters
    for(auto& cl : track_clusters_) {
        auto* cluster = cl.get();
        if(cluster == nullptr) {
            throw MissingReferenceException(typeid(*this), typeid(Cluster));
        }
        if(get_plane(cluster->detectorID()) == nullptr) {
            throw MissingReferenceException(typeid(*this), typeid(Plane));
        }

        // Get the distance and the error
        auto intercept = get_plane(cluster->detectorID())->getToLocal() * getState(cluster->detectorID());
        auto dist = cluster->local() - intercept;
        double ex2 = cluster->errorX() * cluster->errorX();
        double ey2 = cluster->errorY() * cluster->errorY();
        chi2_ += ((dist.x() * dist.x() / ex2) + (dist.y() * dist.y() / ey2));
    }

    // Store also the chi2/degrees of freedom
    chi2ndof_ = (ndof_ <= 0) ? -1 : (chi2_ / static_cast<double>(ndof_));
}

void StraightLineTrack::calculateResiduals() {
    for(const auto& c : track_clusters_) {
        auto* cluster = c.get();
        if(get_plane(cluster->detectorID()) == nullptr) {
            throw MissingReferenceException(typeid(*this), typeid(Plane));
        }
        residual_global_[cluster->detectorID()] = cluster->global() - getState(cluster->detectorID());
        residual_local_[cluster->detectorID()] =
            cluster->local() - get_plane(cluster->detectorID())->getToLocal() * getState(cluster->detectorID());
    }
}

double StraightLineTrack::operator()(const double* parameters) {

    // Update the StraightLineTrack gradient and intercept
    this->m_direction.SetX(parameters[0]);
    this->m_state.SetX(parameters[1]);
    this->m_direction.SetY(parameters[2]);
    this->m_state.SetY(parameters[3]);

    // Calculate the chi2
    this->calculateChi2();

    // Return this to minuit
    return chi2_;
}

void StraightLineTrack::fit() {

    isFitted_ = false;
    Eigen::Matrix4d mat(Eigen::Matrix4d::Zero());
    Eigen::Vector4d vec(Eigen::Vector4d::Zero());

    // Loop over all clusters and fill the matrices
    for(auto& cl : track_clusters_) {
        auto* cluster = cl.get();
        if(cluster == nullptr) {
            throw MissingReferenceException(typeid(*this), typeid(Cluster));
        }

        // Get the measurement and its errors
        double x = cluster->global().x();
        double y = cluster->global().y();
        double z = cluster->global().z();
        Eigen::Vector2d pos(x, y);
        auto errorMatrix = cluster->errorMatrixGlobal();
        Eigen::Matrix2d V;
        V << errorMatrix(0, 0), errorMatrix(0, 1), errorMatrix(1, 0), errorMatrix(1, 1);
        Eigen::Matrix<double, 2, 4> C;
        C << 1., z, 0., 0., 0., 0., 1., z;

        // Fill the matrices
        if(fabs(V.determinant()) < std::numeric_limits<double>::epsilon()) {
            throw TrackFitError(typeid(this), "Error matrix inversion in straight line fit failed");
        }
        vec += C.transpose() * V.inverse() * pos;
        mat += C.transpose() * V.inverse() * C;
    }

    // Check for singularities.
    if(fabs(mat.determinant()) < std::numeric_limits<double>::epsilon()) {
        throw TrackFitError(typeid(this), "Martix inversion in straight line fit failed");
    }
    // Get the StraightLineTrack parameters
    Eigen::Vector4d res = mat.inverse() * vec;

    // Set the StraightLineTrack parameters
    m_state.SetX(res(0));
    m_state.SetY(res(2));
    m_state.SetZ(0.);

    m_direction.SetX(res(1));
    m_direction.SetY(res(3));
    m_direction.SetZ(1.);

    // Calculate the chi2
    this->calculateChi2();
    this->calculateResiduals();
    isFitted_ = true;
}

ROOT::Math::XYZPoint StraightLineTrack::getIntercept(double z) const {
    return m_state + m_direction * z;
}

void StraightLineTrack::print(std::ostream& out) const {
    out << "StraightLineTrack " << this->m_state << ", " << this->m_direction << ", " << this->chi2_ << ", " << this->ndof_
        << ", " << this->chi2ndof_ << ", " << this->timestamp();
}
