/**
 * @file
 * @brief Implementation of the detector model
 *
 * @copyright Copyright (c) 2017-2022 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include <fstream>
#include <map>
#include <string>

#include "Math/RotationX.h"
#include "Math/RotationY.h"
#include "Math/RotationZ.h"
#include "Math/RotationZYX.h"

#include "Detector.hpp"
#include "HexagonalPixelDetector.hpp"
#include "core/utils/log.h"
#include "exceptions.h"

using namespace ROOT::Math;
using namespace corryvreckan;

Detector::Detector(const Configuration& config) : m_role(DetectorRole::NONE) {

    // Role of this detector:
    auto roles = config.getArray<DetectorRole>("role", {DetectorRole::NONE});
    for(auto& role : roles) {
        LOG(DEBUG) << "Adding role " << corryvreckan::to_string(role);
        m_role |= role;
    }

    // Auxiliary devices cannot hold other roles:
    if(hasRole(DetectorRole::AUXILIARY) && m_role != DetectorRole::AUXILIARY) {
        throw InvalidValueError(config, "role", "Auxiliary devices cannot hold any other detector role");
    }

    if(hasRole(DetectorRole::PASSIVE) && m_role != DetectorRole::PASSIVE) {
        throw InvalidValueError(config, "role", "Passive detector cannot hold any other role");
    }

    m_detectorName = config.getName();

    // Material budget of detector, including support material
    if(!config.has("material_budget")) {
        m_materialBudget = 0.0;
        LOG(WARNING) << "No material budget given for " << m_detectorName << ", assuming " << m_materialBudget;
    } else if(config.get<double>("material_budget") < 0) {
        throw InvalidValueError(config, "material_budget", "Material budget has to be positive");
    } else {
        m_materialBudget = config.get<double>("material_budget");
    }

    m_detectorType = config.get<std::string>("type");
    std::transform(m_detectorType.begin(), m_detectorType.end(), m_detectorType.begin(), ::tolower);
    m_detectorCoordinates = config.get<std::string>("coordinates", "cartesian");
    std::transform(m_detectorCoordinates.begin(), m_detectorCoordinates.end(), m_detectorCoordinates.begin(), ::tolower);
    m_timeOffset = config.get<double>("time_offset", 0.0);
    if(m_timeOffset > 0.) {
        LOG(TRACE) << "Time offset: " << m_timeOffset;
    }

    // Time resolution - default to negative number, i.e. unknown. This will trigger an exception
    // when calling getTimeResolution
    m_timeResolution = config.get<double>("time_resolution", -1.0);
    if(m_timeResolution > 0) {
        LOG(TRACE) << "  Time resolution: " << Units::display(m_timeResolution, {"ms", "us"});
    }

    if(config.has("calibration_file")) {
        m_calibrationfile = config.getPath("calibration_file", true);
        LOG(DEBUG) << "Found calibration file for detector " << getName() << " at " << m_calibrationfile.value_or("");
    }

    // Get alignment done:
    alignment_ = std::make_shared<Alignment>(config);
}

std::shared_ptr<Detector> corryvreckan::Detector::factory(const Configuration& config) {
    // default coordinate is cartesian coordinate
    auto coordinates = config.get<std::string>("coordinates", "cartesian");
    std::transform(coordinates.begin(), coordinates.end(), coordinates.begin(), ::tolower);
    if(coordinates == "cartesian") {
        return std::make_shared<PixelDetector>(config);
    } else if(coordinates == "hexagonal") {
        return std::make_shared<HexagonalPixelDetector>(config);
    } else {
        throw InvalidValueError(config, "coordinates", "Coordinates can only set to be cartesian now");
    }
}

Detector::Alignment::Alignment(const Configuration& config) {

    // Set update granularity for alignment transformations
    granularity_ = config.get<double>("alignment_update_granularity", 1e9);

    // Get the orientation right - we keep this constant:
    orientation_ = config.get<ROOT::Math::XYZVector>("orientation", ROOT::Math::XYZVector());
    auto mode = config.get<std::string>("orientation_mode", "xyz");

    if(mode == "xyz") {
        LOG(DEBUG) << "Interpreting Euler angles as XYZ rotation";
        // First angle given in the configuration file is around x, second around y, last around z:
        rotation_fct_ = [](const ROOT::Math::XYZVector& rot) {
            return RotationZ(rot.Z()) * RotationY(rot.Y()) * RotationX(rot.X());
        };
    } else if(mode == "zyx") {
        LOG(DEBUG) << "Interpreting Euler angles as ZYX rotation";
        // First angle given in the configuration file is around z, second around y, last around x:
        rotation_fct_ = [](const ROOT::Math::XYZVector& rot) {
            return static_cast<ROOT::Math::Rotation3D>(RotationZYX(rot.x(), rot.y(), rot.z()));
        };
    } else if(mode == "zxz") {
        LOG(DEBUG) << "Interpreting Euler angles as ZXZ rotation";
        // First angle given in the configuration file is around z, second around x, last around z:
        rotation_fct_ = [](const ROOT::Math::XYZVector& rot) {
            return static_cast<ROOT::Math::Rotation3D>(EulerAngles(rot.x(), rot.y(), rot.z()));
        };
    } else {
        throw InvalidValueError(config, "orientation_mode", "orientation_mode should be either 'zyx', xyz' or 'zxz'");
    }

    // First try to only read a regular position value from the config:
    try {
        displacement_ = config.get<ROOT::Math::XYZPoint>("position", ROOT::Math::XYZPoint());
        needs_update_ = false;

        // Force calculation with the given position and orientation
        update(displacement_, orientation_);
    } catch(InvalidKeyError&) {
    }

    // Let's get the formulae for the positions:
    auto position_functions = config.getArray<std::string>("position");
    if(position_functions.size() != 3) {
        throw InvalidValueError(config, "position", "Position needs to have three components");
    }

    px = std::make_shared<TFormula>("px", position_functions.at(0).c_str(), false);
    py = std::make_shared<TFormula>("py", position_functions.at(1).c_str(), false);
    pz = std::make_shared<TFormula>("pz", position_functions.at(2).c_str(), false);

    // Check that the formulae could correctly be compiled
    if(!(px->IsValid() && py->IsValid() && pz->IsValid())) {
        throw InvalidValueError(config, "position", "Invalid formulae");
    }

    // Check that we have the correct number of dimensions
    if(px->GetNdim() > 1 || py->GetNdim() > 1 || pz->GetNdim() > 1) {
        throw InvalidValueError(config, "position", "Invalid number of dimensions, only 1d is supported");
    }

    // We have constant values, no update needed
    if(px->GetNdim() == 0 && py->GetNdim() == 0 && pz->GetNdim() == 0) {
        LOG(DEBUG) << "Constant functions, no updates needed";
        needs_update_ = false;
    } else {
        needs_update_ = true;
    }

    // Check if we expect parameters
    auto allpars = static_cast<size_t>(px->GetNpar() + py->GetNpar() + pz->GetNpar());
    if(allpars != 0) {
        LOG(DEBUG) << "Formulae require " << allpars << " parameters.";

        // Parse parameters:
        auto position_parameters = config.getArray<double>("position_parameters");
        if(allpars != position_parameters.size()) {
            throw InvalidValueError(config,
                                    "position_parameters",
                                    "The number of position parameters does not line up with the sum of "
                                    "parameters in all functions.");
        }

        // Apply parameters to the functions
        for(auto n = 0; n < px->GetNpar(); ++n) {
            px->SetParameter(n, position_parameters.at(static_cast<size_t>(n)));
        }
        for(auto n = 0; n < py->GetNpar(); ++n) {
            py->SetParameter(n, position_parameters.at(static_cast<size_t>(n + px->GetNpar())));
        }
        for(auto n = 0; n < pz->GetNpar(); ++n) {
            pz->SetParameter(n, position_parameters.at(static_cast<size_t>(n + px->GetNpar() + py->GetNpar())));
        }
    }

    // Force first calculation at t = 0
    update(0., true);
}

void Detector::Alignment::update(double time, bool force) {
    // Check if we need to update already
    if(!force && (time < last_time_ + granularity_ || !needs_update_)) {
        return;
    }

    LOG(DEBUG) << "Calculating updated transformations at t = " << Units::display(time, {"ns", "us", "ms", "s"});

    // Calculate current translation from formulae
    displacement_ = ROOT::Math::XYZVector(px->Eval(time), py->Eval(time), pz->Eval(time));
    LOG(TRACE) << "Displacement " << displacement_;

    recalculate();

    // Update time
    last_time_ = time;
}

void Detector::Alignment::update(const ROOT::Math::XYZPoint& displacement, const ROOT::Math::XYZVector& orientation) {

    LOG(DEBUG) << "Calculating updated transformations with external displacement and orientation";
    LOG(TRACE) << "Displacement " << displacement;
    displacement_ = displacement;
    LOG(TRACE) << "Orientation " << orientation;
    orientation_ = orientation;

    recalculate();
}

void Detector::Alignment::recalculate() {

    auto translations = Translation3D(displacement_.X(), displacement_.Y(), displacement_.Z());
    auto rotations = rotation_fct_(orientation_);

    // Calculate current local-to-global transformation and its inverse:
    local2global_ = Transform3D(rotations, translations);
    global2local_ = local2global_.Inverse();

    // Find the normal to the detector surface. Build two points, the origin and a unit step in z,
    // transform these points to the global coordinate frame and then make a vector pointing between them
    origin_ = ROOT::Math::XYZVector(0., 0., 0.);
    origin_ = local2global_ * origin_;
    LOG(TRACE) << "Origin " << origin_;

    auto local_z = local2global_ * ROOT::Math::XYZPoint(0., 0., 1.);
    normal_ = ROOT::Math::XYZVector(local_z.X() - origin_.X(), local_z.Y() - origin_.Y(), local_z.Z() - origin_.Z());
    LOG(TRACE) << "Normal " << normal_;
}

double Detector::getTimeResolution() const {
    if(m_timeResolution > 0) {
        return m_timeResolution;
    } else {
        throw InvalidSettingError(this, "time_resolution", "Time resolution not set but requested");
    }
}

std::string Detector::getName() const {
    return m_detectorName;
}

std::string Detector::getType() const {
    return m_detectorType;
}

bool Detector::isReference() const {
    return static_cast<bool>(m_role & DetectorRole::REFERENCE);
}

bool Detector::isDUT() const {
    return static_cast<bool>(m_role & DetectorRole::DUT);
}

bool Detector::isAuxiliary() const {
    return static_cast<bool>(m_role & DetectorRole::AUXILIARY);
}

bool Detector::isPassive() const {
    return static_cast<bool>(m_role & DetectorRole::PASSIVE);
}

DetectorRole Detector::getRoles() const {
    return m_role;
}

bool Detector::hasRole(DetectorRole role) const {
    return static_cast<bool>(m_role & role);
}

// Function to set the channel maskfile
void Detector::maskFile(std::filesystem::path file) {
    m_maskfile = std::move(file);
}

// Function to update transforms (such as during alignment)
void Detector::update(double time) {
    alignment_->update(time);
}

void Detector::update(const ROOT::Math::XYZPoint& displacement, const ROOT::Math::XYZVector& orientation) {
    alignment_->update(displacement, orientation);
}

Configuration Detector::getConfiguration() const {

    Configuration config(getName());
    config.set("type", m_detectorType);

    if(m_detectorCoordinates != "cartesian") {
        config.set("coordinates", m_detectorCoordinates);
    }

    // Store the role of the detector
    std::vector<std::string> roles;
    if(this->isDUT()) {
        roles.emplace_back("dut");
    }
    if(this->isReference()) {
        roles.emplace_back("reference");
    }
    if(this->isAuxiliary()) {
        roles.emplace_back("auxiliary");
    }
    if(this->isPassive()) {
        roles.push_back("passive");
    }

    if(!roles.empty()) {
        config.setArray("role", roles);
    }

    if(m_timeOffset != 0.) {
        config.set("time_offset", m_timeOffset, {"ns", "us", "ms", "s"});
    }

    config.set("time_resolution", m_timeResolution, {"ns", "us", "ms", "s"});

    // different for PixelDetector and StripDetector
    this->configure_pos_and_orientation(config);

    // material budget
    if(m_materialBudget > std::numeric_limits<double>::epsilon()) {
        config.set("material_budget", m_materialBudget);
    }

    // only if detector is not auxiliary:
    if(!this->isAuxiliary()) {
        this->configure_detector(config);
    }

    // add detector calibration file path, if set in the main configuration
    if(this->calibrationFile() != "") {
        config.set("calibration_file", this->calibrationFile());
    }

    return config;
}
