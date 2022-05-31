/** @file
 *  @brief Implementation of the detector model
 *  @copyright Copyright (c) 2017-2022 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include <fstream>
#include <map>
#include <string>

#include "Math/RotationX.h"
#include "Math/RotationY.h"
#include "Math/RotationZ.h"
#include "Math/RotationZYX.h"

#include "PolarDetector.hpp"
#include "core/utils/log.h"
#include "exceptions.h"

using namespace ROOT::Math;
using namespace corryvreckan;

PolarDetector::PolarDetector(const Configuration& config) : Detector(config) {

    // Set detector position and direction from configuration file
    SetPositionAndOrientation(config);

    // Auxiliary devices don't have: number_of_pixels, pixel_pitch, spatial_resolution, mask_file, region-of-interest
    if(!isAuxiliary()) {
        build_axes(config);
    }

    // initialize transform
    this->initialise();

    // Control printouts for testing purposes
    LOG(TRACE) << "== CONTROL PRINTOUTS FOR RADIAL DETECTORS:";
    LOG(TRACE) << "====> Numbers of strips:";
    for(const auto val : number_of_strips) {
        LOG(TRACE) << "====> " << val;
    }
    LOG(TRACE) << "====> Row radii:";
    for(const auto val : row_radius) {
        LOG(TRACE) << "====> " << val;
    }

    LOG(TRACE) << "====> Angular pitch:";
    for(const auto val : angular_pitch) {
        LOG(TRACE) << "====> " << val;
    }
    LOG(TRACE) << "===> Spatial res.: " << m_spatial_resolution;
    LOG(TRACE) << "===> Stereo angle: " << stereo_angle;
    LOG(TRACE) << "===> Focus shift: " << focus_translation;
}

void PolarDetector::build_axes(const Configuration& config) {

    // Get numbers of strips
    number_of_strips = config.getArray<unsigned int>("number_of_strips");

    // Get angular pitches
    angular_pitch = config.getArray<double>("angular_pitch");

    // Get row radii
    row_radius = config.getArray<double>("row_radius");

    // Get stereo angle
    stereo_angle = config.get<double>("stereo_angle", 0.0);

    // Calculate translation from sensor origin to its focal point
    focus_translation = {getCenterRadius() * sin(stereo_angle), getCenterRadius() * (1 - cos(stereo_angle)), 0};

    // Set reasonable pixel pitch placeholders - length of the strip edge and strip length
    m_pitch = {row_radius.at(3) * angular_pitch.at(3), row_radius.at(3) - row_radius.at(2)};

    LOG(TRACE) << "Initialized \"" << m_detectorType;

    // Intrinsic spatial resolution, defaults to pitch/sqrt(12):
    // Resolution in X direction: maximum strip pitch / sqrt(12)
    // Resolution in Y direction: length of the first strip row / sqrt(12)
    m_spatial_resolution =
        config.get<ROOT::Math::XYVector>("spatial_resolution",
                                         {*std::max_element(angular_pitch.begin(), angular_pitch.end()) / sqrt(12),
                                          (row_radius.at(1) - row_radius.at(0)) / sqrt(12)});
    if(!config.has("spatial_resolution")) {
        LOG(WARNING) << "Spatial resolution for detector '" << m_detectorName << "' not set." << std::endl
                     << "Using pitch/sqrt(12) as default";
    }

    // region of interest:
    m_roi = config.getMatrix<int>("roi", std::vector<std::vector<int>>());

    if(config.has("mask_file")) {
        auto mask_file = config.getPath("mask_file", true);
        LOG(DEBUG) << "Adding mask to detector \"" << config.getName() << "\", reading from " << mask_file;
        maskFile(mask_file);
        process_mask_file();
    }
}

void PolarDetector::SetPositionAndOrientation(const Configuration& config) {
    // Detector position and orientation
    m_displacement = config.get<ROOT::Math::XYZPoint>("position", ROOT::Math::XYZPoint());
    m_orientation = config.get<ROOT::Math::XYZVector>("orientation", ROOT::Math::XYZVector());
    m_orientation_mode = config.get<std::string>("orientation_mode", "xyz");

    if(m_orientation_mode != "xyz" && m_orientation_mode != "zyx" && m_orientation_mode != "zxz") {
        throw InvalidValueError(config, "orientation_mode", "orientation_mode should be either 'zyx', xyz' or 'zxz'");
    }

    LOG(TRACE) << "  Position:    " << Units::display(m_displacement, {"mm", "um"});
    LOG(TRACE) << "  Orientation: " << Units::display(m_orientation, {"deg"}) << " (" << m_orientation_mode << ")";
}

void PolarDetector::process_mask_file() {
    // Open the file with masked pixels
    std::ifstream inputMaskFile(m_maskfile, std::ios::in);
    if(!inputMaskFile.is_open()) {
        LOG(WARNING) << "Could not open mask file " << m_maskfile;
    } else {
        int row = 0, col = 0;
        std::string id;
        // loop over all lines and apply masks
        while(inputMaskFile >> id) {
            if(id == "c") {
                inputMaskFile >> col;
                if(col > nPixels().X() - 1) {
                    LOG(WARNING) << "Column " << col << " outside of pixel matrix, chip has only " << nPixels().X()
                                 << " columns!";
                }
                LOG(TRACE) << "Masking column " << col;
                for(int r = 0; r < nPixels().Y(); r++) {
                    maskChannel(col, r);
                }
            } else if(id == "r") {
                inputMaskFile >> row;
                if(row > nPixels().Y() - 1) {
                    LOG(WARNING) << "Row " << col << " outside of pixel matrix, chip has only " << nPixels().Y() << " rows!";
                }
                LOG(TRACE) << "Masking row " << row;
                for(int c = 0; c < nPixels().X(); c++) {
                    maskChannel(c, row);
                }
            } else if(id == "p") {
                inputMaskFile >> col >> row;
                if(col > nPixels().X() - 1 || row > nPixels().Y() - 1) {
                    LOG(WARNING) << "Pixel " << col << " " << row << " outside of pixel matrix, chip has only "
                                 << nPixels().X() << " x " << nPixels().Y() << " pixels!";
                }
                LOG(TRACE) << "Masking pixel " << col << " " << row;
                maskChannel(col, row); // Flag to mask a pixel
            } else {
                LOG(WARNING) << "Could not parse mask entry (id \"" << id << "\")";
            }
        }
        LOG(INFO) << m_masked.size() << " masked pixels";
    }
}

void PolarDetector::maskChannel(int chX, int chY) {
    int channelID = chX + chY;
    m_masked[channelID] = true;
}

bool PolarDetector::masked(int chX, int chY) const {
    int channelID = chX + chY;
    if(m_masked.count(channelID) > 0)
        return true;
    return false;
}

// Function to initialise transforms
void PolarDetector::initialise() {

    // Make the local to global transform, built from a displacement and rotation
    Translation3D translations = Translation3D(m_displacement.X(), m_displacement.Y(), m_displacement.Z());

    Rotation3D rotations;
    if(m_orientation_mode == "xyz") {
        LOG(DEBUG) << "Interpreting Euler angles as XYZ rotation";
        // First angle given in the configuration file is around x, second around y, last around z:
        rotations = RotationZ(m_orientation.Z()) * RotationY(m_orientation.Y()) * RotationX(m_orientation.X());
    } else if(m_orientation_mode == "zyx") {
        LOG(DEBUG) << "Interpreting Euler angles as ZYX rotation";
        // First angle given in the configuration file is around z, second around y, last around x:
        rotations = RotationZYX(m_orientation.x(), m_orientation.y(), m_orientation.z());
    } else if(m_orientation_mode == "zxz") {
        LOG(DEBUG) << "Interpreting Euler angles as ZXZ rotation";
        // First angle given in the configuration file is around z, second around x, last around z:
        rotations = EulerAngles(m_orientation.x(), m_orientation.y(), m_orientation.z());
    } else {
        throw InvalidSettingError(this, "orientation_mode", "orientation_mode should be either 'zyx', xyz' or 'zxz'");
    }

    // Additional translation to have local coordinates calculated from the sensor origin
    auto origin_trf = Translation3D(0, -getCenterRadius(), 0);
    m_localToGlobal = Transform3D(rotations, translations * origin_trf);
    m_globalToLocal = m_localToGlobal.Inverse();

    // Find the normal to the detector surface. Build two points, the origin and a unit step in z,
    // transform these points to the global coordinate frame and then make a vector pointing between them
    m_origin = PositionVector3D<Cartesian3D<double>>(0., 0., 0.);
    m_origin = m_localToGlobal * m_origin;
    PositionVector3D<Cartesian3D<double>> localZ(0., 0., 1.);
    localZ = m_localToGlobal * localZ;
    m_normal = PositionVector3D<Cartesian3D<double>>(
        localZ.X() - m_origin.X(), localZ.Y() - m_origin.Y(), localZ.Z() - m_origin.Z());
}

// Only if detector is not auxiliary
void PolarDetector::configure_detector(Configuration& config) const {

    // Number of strips
    config.setArray("number_of_strips", number_of_strips);

    // Strip angular pitches
    config.setArray("angular_pitch", angular_pitch);

    // Row radii
    config.setArray("row_radius", row_radius);

    // Stereo angle
    config.set("stereo_angle", stereo_angle);

    // Intrinsic resolution:
    config.set("spatial_resolution", m_spatial_resolution, {{"um"}});

    // Pixel mask file:
    if(!m_maskfile.empty()) {
        config.set("mask_file", m_maskfile.string());
    }

    // Region-of-interest:
    config.setMatrix("roi", m_roi);
}

void PolarDetector::configure_pos_and_orientation(Configuration& config) const {
    config.set("position", m_displacement, {"um", "mm"});
    config.set("orientation_mode", m_orientation_mode);
    config.set("orientation", m_orientation, {{"deg"}});
}

// Function to get global intercept with a track
PositionVector3D<Cartesian3D<double>> PolarDetector::getIntercept(const Track* track) const {

    // FIXME: this is else statement can only be temporary
    if(track->getType() == "GblTrack") {
        return track->getState(getName());
    } else {
        // Get the distance from the plane to the track initial state
        double distance = (m_origin.X() - track->getState(m_detectorName).X()) * m_normal.X();
        distance += (m_origin.Y() - track->getState(m_detectorName).Y()) * m_normal.Y();
        distance += (m_origin.Z() - track->getState(m_detectorName).Z()) * m_normal.Z();
        distance /= (track->getDirection(m_detectorName).X() * m_normal.X() +
                     track->getDirection(m_detectorName).Y() * m_normal.Y() +
                     track->getDirection(m_detectorName).Z() * m_normal.Z());

        // Propagate the track
        PositionVector3D<Cartesian3D<double>> globalIntercept(
            track->getState(m_detectorName).X() + distance * track->getDirection(m_detectorName).X(),
            track->getState(m_detectorName).Y() + distance * track->getDirection(m_detectorName).Y(),
            track->getState(m_detectorName).Z() + distance * track->getDirection(m_detectorName).Z());
        return globalIntercept;
    }
}

PositionVector3D<Cartesian3D<double>> PolarDetector::getLocalIntercept(const Track* track) const {
    return globalToLocal(getIntercept(track));
}

// Function to check if a track intercepts with a plane
bool PolarDetector::hasIntercept(const Track* track, double pixelTolerance) const {

    // First, get the track intercept in global coordinates with the plane
    PositionVector3D<Cartesian3D<double>> globalIntercept = this->getIntercept(track);
    // Convert to local coordinates
    PositionVector3D<Cartesian3D<double>> localIntercept = this->m_globalToLocal * globalIntercept;

    // Get the row and column numbers
    double row = this->getRow(localIntercept);
    double column = this->getColumn(localIntercept);

    // Check if the row and column are outside of the chip
    // Chip reaches from -0.5 to nPixels-0.5
    bool intercept = true;
    if(row < pixelTolerance - 0.5 || row > (static_cast<double>(number_of_strips.size()) - pixelTolerance - 0.5) ||
       column < pixelTolerance - 0.5 ||
       column > (number_of_strips.at(static_cast<unsigned int>(floor(row + 0.5))) - pixelTolerance - 0.5)) {
        intercept = false;
    }

    return intercept;
}

// Function to check if a track goes through/near a masked pixel
bool PolarDetector::hitMasked(const Track* track, int tolerance) const {

    // First, get the track intercept in global coordinates with the plane
    PositionVector3D<Cartesian3D<double>> globalIntercept = this->getIntercept(track);

    // Convert to local coordinates
    PositionVector3D<Cartesian3D<double>> localIntercept = this->m_globalToLocal * globalIntercept;

    // Get the row and column numbers
    int row = static_cast<int>(floor(this->getRow(localIntercept) + 0.5));
    int column = static_cast<int>(floor(this->getColumn(localIntercept) + 0.5));

    // Check if the pixels around this pixel are masked
    bool hitmasked = false;
    for(int r = (row - tolerance); r <= (row + tolerance); r++) {
        for(int c = (column - tolerance); c <= (column + tolerance); c++) {
            if(this->masked(c, r)) {
                hitmasked = true;
            }
        }
    }

    return hitmasked;
}

// Functions to get row and column from local position
double PolarDetector::getRow(const PositionVector3D<Cartesian3D<double>> localPosition) const {
    // Convert local position to polar coordinates
    auto polar_pos = getPositionPolar(localPosition);

    // Assign to a strip row
    unsigned int strip_y{};
    for(unsigned int row = 0; row < number_of_strips.size(); row++) {
        // Find the correct strip row by comparing to inner and outer row radii
        if(polar_pos.r() > row_radius.at(row) && polar_pos.r() <= row_radius.at(row + 1)) {
            strip_y = row;
            break;
        }
    }

    // Calculate the fraction of the strip row
    auto fract = (polar_pos.r() - row_radius.at(strip_y) - 0.5) / (row_radius.at(strip_y + 1) + row_radius.at(strip_y));

    return static_cast<double>(strip_y) + fract;
}

double PolarDetector::getColumn(const PositionVector3D<Cartesian3D<double>> localPosition) const {
    // Convert local position to polar coordinates
    auto polar_pos = getPositionPolar(localPosition);

    auto strip_y = static_cast<unsigned int>(getRow(localPosition));

    // Get the strip pitch in the correct strip row
    auto pitch = angular_pitch.at(strip_y);
    // Calculate the strip x-index
    auto strip_x = (polar_pos.phi() + stereo_angle + pitch * number_of_strips.at(strip_y) / 2) / pitch;

    return strip_x;
}

PositionVector3D<Polar3D<double>>
PolarDetector::getPositionPolar(const PositionVector3D<Cartesian3D<double>> localPosition) const {
    // Calculate the radial component
    auto r = sqrt(localPosition.x() * localPosition.x() + localPosition.y() * localPosition.y());
    // Shift the coordinate origin to the strip focal point
    auto focus_pos = localPosition - focus_translation;
    // Calculate the angular component obtained from the corrected position
    auto phi = atan2(focus_pos.x(), focus_pos.y());

    return {r, 0, phi};
}

PositionVector3D<Cartesian3D<double>>
PolarDetector::getPositionCartesian(const PositionVector3D<Polar3D<double>> localPosition) const {
    // Length of the translation vector from the local center to the focal point
    auto len_foc = std::sqrt(focus_translation.mag2());
    // Calculate two relevant angles needed for the transformation of the angular
    // component to be measured from the local
    // coordinate center instead of the strip focal point
    auto alpha = std::acos(len_foc / (2 * getCenterRadius()));
    auto gamma = asin(len_foc * sin(alpha + localPosition.phi() + stereo_angle) / localPosition.r());
    // Transform the angle
    auto phi = 2 * alpha + gamma + localPosition.phi() + stereo_angle - ROOT::Math::Pi();

    return {localPosition.r() * sin(phi), localPosition.r() * cos(phi), 0.0};
}

// Function to get local position from row and column
PositionVector3D<Cartesian3D<double>> PolarDetector::getLocalPosition(double column, double row) const {
    // Whole and decimal part of the column and row
    unsigned int column_base = static_cast<unsigned int>(floor(column + 0.5));
    auto column_dec = column - column_base;
    unsigned int row_base = static_cast<unsigned int>(floor(row + 0.5));
    auto row_dec = row - row_base;

    // Calculate the radial coordinate of the strip center
    auto local_r = (row_radius.at(row_base) + row_radius.at(row_base + 1)) / 2;
    // Adjust for in-strip position
    local_r += row_dec * (row_radius.at(row_base + 1) - row_radius.at(row_base));

    // Calculate the angular coordinate of the strip center
    auto local_phi = -angular_pitch.at(row_base) * number_of_strips.at(row_base) / 2 +
                     (column_base + 0.5) * angular_pitch.at(row_base) - stereo_angle;
    local_phi += column_dec * angular_pitch.at(row_base);

    // Convert polar coordinates to cartesian
    auto pos = getPositionCartesian({local_r, 0, local_phi});

    return PositionVector3D<Cartesian3D<double>>(pos.x(), pos.y(), 0.0);
}

// Function to get in-pixel position
ROOT::Math::XYVector PolarDetector::inPixel(const double column, const double row) const {
    // Transform received position to polar coordinates and get the coordinates of the strip center
    auto local_polar = getPositionPolar(getLocalPosition(column, row));
    auto strip_polar = getPositionPolar(getLocalPosition(floor(column + 0.5), floor(row + 0.5)));

    auto delta_phi = local_polar.phi() - strip_polar.phi();

    return {local_polar.r() * sin(delta_phi), local_polar.r() * cos(delta_phi) - strip_polar.r()};
}

ROOT::Math::XYVector PolarDetector::inPixel(const PositionVector3D<Cartesian3D<double>> localPosition) const {
    double column = getColumn(localPosition);
    double row = getRow(localPosition);
    return inPixel(column, row);
}

// Check if track position is within ROI:
bool PolarDetector::isWithinROI(const Track* track) const {

    // Empty region of interest:
    if(m_roi.empty()) {
        return true;
    }

    // Check that track is within region of interest using winding number algorithm
    auto localIntercept = this->getLocalIntercept(track);
    auto coordinates = std::make_pair(this->getColumn(localIntercept), this->getRow(localIntercept));
    if(winding_number(coordinates, m_roi) != 0) {
        return true;
    }

    // Outside ROI:
    return false;
}

// Check if cluster is within ROI and/or touches ROI border:
bool PolarDetector::isWithinROI(Cluster* cluster) const {

    // Empty region of interest:
    if(m_roi.empty()) {
        return true;
    }

    // Loop over all pixels of the cluster
    for(auto& pixel : cluster->pixels()) {
        if(winding_number(pixel->coordinates(), m_roi) == 0) {
            return false;
        }
    }
    return true;
}

XYVector PolarDetector::getSize() const {
    auto max_row = number_of_strips.size() - 1;
    return {angular_pitch.at(max_row) * number_of_strips.at(max_row) * row_radius.at(max_row),
            row_radius.at(max_row) - row_radius.at(0)};
}

/* Winding number test for a point in a polygon
 * via: http://geomalgorithms.com/a03-_inclusion.html
 *      Input:   x, y = a point,
 *               polygon = vector of vertex points of a polygon V[n+1] with V[n]=V[0]
 *      Return:  wn = the winding number (=0 only when P is outside)
 */
int PolarDetector::winding_number(std::pair<int, int> probe, std::vector<std::vector<int>> polygon) {
    // Two points don't make an area
    if(polygon.size() < 3) {
        LOG(DEBUG) << "No ROI given.";
        return 0;
    }

    int wn = 0; // the  winding number counter

    // loop through all edges of the polygon

    // edge from V[i] to  V[i+1]
    for(size_t i = 0; i < polygon.size(); i++) {
        auto point_this = std::make_pair(polygon.at(i).at(0), polygon.at(i).at(1));
        auto point_next = (i + 1 < polygon.size() ? std::make_pair(polygon.at(i + 1).at(0), polygon.at(i + 1).at(1))
                                                  : std::make_pair(polygon.at(0).at(0), polygon.at(0).at(1)));

        // start y <= P.y
        if(point_this.second <= probe.second) {
            // an upward crossing
            if(point_next.second > probe.second) {
                // P left of  edge
                if(isLeft(point_this, point_next, probe) > 0) {
                    // have  a valid up intersect
                    ++wn;
                }
            }
        } else {
            // start y > P.y (no test needed)

            // a downward crossing
            if(point_next.second <= probe.second) {
                // P right of  edge
                if(isLeft(point_this, point_next, probe) < 0) {
                    // have  a valid down intersect
                    --wn;
                }
            }
        }
    }
    return wn;
}
/* isLeft(): tests if a point is Left|On|Right of an infinite line.
 * via: http://geomalgorithms.com/a03-_inclusion.html
 *    Input:  three points P0, P1, and P2
 *    Return: >0 for P2 left of the line through P0 and P1
 *            =0 for P2  on the line
 *            <0 for P2  right of the line
 *    See: Algorithm 1 "Area of Triangles and Polygons"
 */
int PolarDetector::isLeft(std::pair<int, int> pt0, std::pair<int, int> pt1, std::pair<int, int> pt2) {
    return ((pt1.first - pt0.first) * (pt2.second - pt0.second) - (pt2.first - pt0.first) * (pt1.second - pt0.second));
}

// Check if a pixel touches any of the pixels in a cluster
bool PolarDetector::isNeighbor(const std::shared_ptr<Pixel>& neighbor,
                               const std::shared_ptr<Cluster>& cluster,
                               const int neighbor_radius_row,
                               const int neighbor_radius_col) {
    for(const auto* pixel : cluster->pixels()) {
        int row_distance = abs(pixel->row() - neighbor->row());
        int col_distance = abs(pixel->column() - neighbor->column());

        if(row_distance <= neighbor_radius_row && col_distance <= neighbor_radius_col) {
            if(row_distance > 1 || col_distance > 1) {
                cluster->setSplit(true);
            }
            return true;
        }
    }
    return false;
}
