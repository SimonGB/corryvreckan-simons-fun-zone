/**
 * @file
 * @brief Implementation of the polar detector model
 *
 * @copyright Copyright (c) 2017-2023 CERN and the Corryvreckan authors.
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

#include "PolarDetector.hpp"
#include "core/utils/log.h"
#include "exceptions.h"

using namespace ROOT::Math;
using namespace corryvreckan;

PolarDetector::PolarDetector(const Configuration& config) : Detector(config) {
    // Auxiliary devices don't have: number_of_pixels, pixel_pitch, spatial_resolution, mask_file, region-of-interest
    if(!isAuxiliary()) {
        build_axes(config);
    }

    // Compute the spatial resolution in the global coordinates by rotating the error ellipsis
    TMatrixD errorMatrix(3, 3);
    TMatrixD locToGlob(3, 3), globToLoc(3, 3);
    errorMatrix(0, 0) = getSpatialResolution().x() * getSpatialResolution().x();
    errorMatrix(1, 1) = getSpatialResolution().y() * getSpatialResolution().y();
    alignment_->local2global().Rotation().GetRotationMatrix(locToGlob);
    alignment_->global2local().Rotation().GetRotationMatrix(globToLoc);
    m_spatial_resolution_matrix_global = locToGlob * errorMatrix * globToLoc;

    // Print detector parameters
    LOG(DEBUG) << "Polar detector \"" << m_detectorName << "\" parameters:";
    LOG(DEBUG) << "  Number of strips:";
    for(const auto val : number_of_strips)
        LOG(DEBUG) << "    " << val;

    LOG(DEBUG) << "  Row radii:";
    for(const auto val : row_radius)
        LOG(DEBUG) << "    " << Units::display(val, "mm");

    LOG(DEBUG) << "  Strip lengths:";
    for(const auto val : strip_length)
        LOG(DEBUG) << "    " << Units::display(val, "mm");

    LOG(DEBUG) << "  Angular pitch:";
    for(const auto val : angular_pitch)
        LOG(DEBUG) << "    " << Units::display(val, "urad");

    LOG(DEBUG) << "  Spatial resolution: " << Units::display(getSpatialResolution(), {"urad", "mm"});
    LOG(DEBUG) << "  Stereo angle: " << Units::display(stereo_angle, "mrad");
    LOG(DEBUG) << "  Center radius: " << Units::display(center_radius, "mm");
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

    // If central radius was not provided, calculate as average radius
    if(config.has("center_radius")) {
        center_radius = config.get<double>("center_radius");
    } else {
        LOG(WARNING) << "Center radius not provided, calculating as radius average.";
        center_radius = (row_radius.at(0) + row_radius.at(number_of_strips.size())) / 2;
    }

    // If strip lengths were not provided, calculate from row radii
    if(config.has("strip_length")) {
        strip_length = config.getArray<double>("strip_length");
    } else {
        LOG(WARNING) << "Strip lengths not provided, calculating from row radii. "
                     << "This yields approximate results and is inaccurate for larger stereo angles." << std::endl
                     << "Consider defining the lengths using the 'strip_length' parameter.";
        for(unsigned int i = 1; i < row_radius.size(); i++) {
            strip_length.push_back(row_radius.at(i) - row_radius.at(i - 1));
        }
    }

    LOG(TRACE) << "Initialized \"" << m_detectorType << "\"";

    // region of interest:
    m_roi = config.getMatrix<int>("roi", std::vector<std::vector<int>>());

    if(config.has("mask_file")) {
        auto mask_file = config.getPath("mask_file", true);
        LOG(DEBUG) << "Adding mask to detector \"" << config.getName() << "\", reading from " << mask_file;
        maskFile(mask_file);
        process_mask_file();
    }
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

// Only if detector is not auxiliary
void PolarDetector::configure_detector(Configuration& config) const {

    // Number of strips
    config.setArray("number_of_strips", number_of_strips);

    // Strip angular pitches
    config.setArray("angular_pitch", angular_pitch);

    // Row radii
    config.setArray("row_radius", row_radius);

    // Stereo angle
    config.set("stereo_angle", stereo_angle, {{"mrad"}});

    // Intrinsic resolution:
    config.set("spatial_resolution", getSpatialResolution(), {"urad", "mm"});

    // Pixel mask file:
    if(!m_maskfile.empty()) {
        config.set("mask_file", m_maskfile.string());
    }

    // Region-of-interest:
    config.setMatrix("roi", m_roi);
}

void PolarDetector::configure_pos_and_orientation(Configuration& config) const {
    config.set("position", alignment_->displacement(), {"um", "mm"});
    config.set("orientation", alignment_->orientation(), {{"deg"}});
    config.set("coordinates", "polar");
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

bool PolarDetector::hitMasked(const Track* track, int tolerance) const {
    // First, get the track intercept in global coordinates with the plane
    PositionVector3D<Cartesian3D<double>> globalIntercept = this->getIntercept(track);

    // Convert to local coordinates
    PositionVector3D<Cartesian3D<double>> localIntercept = alignment_->global2local() * globalIntercept;

    // Get the row and column numbers
    int row = static_cast<int>(floor(this->getRow(localIntercept) + 0.5));
    int column = static_cast<int>(floor(this->getColumn(localIntercept) + 0.5));

    // Check if the strips around this strip are masked
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

PositionVector3D<Cartesian3D<double>> PolarDetector::getIntercept(const Track* track) const {
    return track->getState(getName());
}

PositionVector3D<Cartesian3D<double>> PolarDetector::getLocalIntercept(const Track* track) const {
    return globalToLocal(getIntercept(track));
}

bool PolarDetector::hasIntercept(const Track* track, double stripTolerance) const {

    // First, get the track intercept in global coordinates with the plane
    PositionVector3D<Cartesian3D<double>> globalIntercept = this->getIntercept(track);
    // Convert to local coordinates
    PositionVector3D<Cartesian3D<double>> localIntercept = alignment_->global2local() * globalIntercept;

    // Get the row and column numbers
    double row = this->getRow(localIntercept);
    double column = this->getColumn(localIntercept);

    // Check if the row and column are outside of the chip
    // Chip reaches from -0.5 to nPixels-0.5
    bool intercept = true;
    if(row < stripTolerance - 0.5 || row > (static_cast<double>(number_of_strips.size()) - stripTolerance - 0.5) ||
       column < stripTolerance - 0.5 ||
       column > (number_of_strips.at(static_cast<unsigned int>(floor(row + 0.5))) - stripTolerance - 0.5)) {
        intercept = false;
    }

    return intercept;
}

double PolarDetector::getRow(const PositionVector3D<Cartesian3D<double>> localPosition) const {
    // Convert local position to polar coordinates
    auto polar_pos = getPolarPosition(localPosition);

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
    auto polar_pos = getPolarPosition(localPosition);

    // Assign to a strip row
    auto strip_y = static_cast<unsigned int>(getRow(localPosition));

    // Get the strip pitch in the correct strip row
    auto pitch = angular_pitch.at(strip_y);

    // Calculate the strip column
    return (number_of_strips.at(strip_y) - 1) / 2 - polar_pos.phi() / pitch;
}

PositionVector3D<Cartesian3D<double>> PolarDetector::getLocalPosition(double column, double row) const {
    // Get position in polar coordinates
    auto local_polar = getPolarPosition(column, row);

    // Convert from local polar to local cartesian coordinates
    return getLocalPosition(local_polar);
}

PositionVector3D<Cartesian3D<double>>
PolarDetector::getLocalPosition(const PositionVector3D<Polar3D<double>> polarPosition) const {
    auto local_x = polarPosition.R() * sin(polarPosition.Phi() + stereo_angle) - center_radius * sin(stereo_angle);
    auto local_y = polarPosition.R() * cos(polarPosition.Phi() + stereo_angle) - center_radius * cos(stereo_angle);

    return {local_x, local_y, 0};
}

PositionVector3D<Polar3D<double>> PolarDetector::getPolarPosition(double column, double row) const {
    // Assign to a strip row
    unsigned int strip_y = (row < 0) ? 0 : static_cast<unsigned int>(row);

    // Get the strip pitch and number of strips in the strip row
    auto pitch = angular_pitch.at(strip_y);
    auto n_strips = number_of_strips.at(strip_y);

    // Calculate polar angle
    auto phi = pitch * ((n_strips - 1) / 2 - column);

    // Split row index into integer and fractional part
    unsigned int row_int = (row < 0) ? 0 : static_cast<unsigned int>(floor(row));
    auto row_fract = row - row_int + 0.5;

    // Get inner and outer row radii
    auto r1 = row_radius.at(row_int);
    auto r2 = row_radius.at(row_int + 1);
    // Calculate radius
    auto r = r1 + (r2 - r1) * row_fract;

    // Convert radius from beam frame to strip frame
    // Reference transformation from ATLAS12ECTechnicalSpecs_v2.3
    auto b = -4.0 * center_radius * sin(stereo_angle / 2.0) * sin(stereo_angle / 2.0 + phi);
    auto c = std::pow((2.0 * center_radius * sin(stereo_angle / 2.0)), 2) - r * r;
    auto r_conv = 0.5 * (-b + sqrt(b * b - 4 * c));

    return {r_conv, 0, phi};
}

PositionVector3D<Polar3D<double>>
PolarDetector::getPolarPosition(const PositionVector3D<Cartesian3D<double>> localPosition) const {
    // Calculate polar angle
    auto phi = (atan2(center_radius * sin(stereo_angle) + localPosition.X(),
                      center_radius * cos(stereo_angle) + localPosition.Y()) -
                stereo_angle);

    // Calculate radius
    auto r = (center_radius * cos(stereo_angle) + localPosition.Y()) / cos(stereo_angle + phi);

    return {r, 0, phi};
}

ROOT::Math::XYVector PolarDetector::inPixel(const double column, const double row) const {
    // Transform received position to polar coordinates and get the coordinates of the strip center
    auto local_polar = getPolarPosition(getLocalPosition(column, row));
    auto strip_polar = getPolarPosition(getLocalPosition(floor(column + 0.5), floor(row + 0.5)));

    auto delta_phi = local_polar.phi() - strip_polar.phi();

    return {local_polar.r() * sin(delta_phi), local_polar.r() * cos(delta_phi) - strip_polar.r()};
}

ROOT::Math::XYVector PolarDetector::inPixel(const PositionVector3D<Cartesian3D<double>> localPosition) const {
    double column = getColumn(localPosition);
    double row = getRow(localPosition);
    return inPixel(column, row);
}

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
    /* The size of a polar detector is approximated as:
     * - X direction: Length of the arc defining the final (longest)
     *                strip row, as [number of strips]*[angular pitch]*[radius]
     * - Y direction: Sum of strip lengths + 20%, to account for the radial shape
     */

    // Final strip row
    auto max_row = number_of_strips.size() - 1;

    // Calculate approximate sizes
    auto size_x = number_of_strips.at(max_row) * angular_pitch.at(max_row) * row_radius.at(max_row + 1);
    auto size_y = (row_radius.at(max_row + 1) - row_radius.at(0)) * 1.2;

    return {size_x, size_y};
}

XYVector PolarDetector::getPitch() const {
    /* The strip pitch is approximated as:
     * - X direction: Length of the arc defining the longer edge of the biggest
     *                possible strip, as [max angular pitch] * [max radius]
     * - Y direction: Length of the longest strip
     */

    auto max_pitch = *std::max_element(angular_pitch.begin(), angular_pitch.end());
    auto pitch_x = max_pitch * row_radius.at(number_of_strips.size());

    auto pitch_y = *std::max_element(strip_length.begin(), strip_length.end());

    return {pitch_x, pitch_y};
}

XYVector PolarDetector::getSpatialResolution(double, double row) const {
    // Get integer row
    auto row_int = static_cast<unsigned int>(floor(row + 0.5));

    // Get strip length and pitch for the given row
    auto strip_r = strip_length.at(row_int);
    auto strip_phi = angular_pitch.at(row_int);

    // Resolution is angular pitch and strip length / sqrt(12)
    return {strip_phi / sqrt(12), strip_r / sqrt(12)};
}

bool PolarDetector::isNeighbor(const std::shared_ptr<Pixel>& neighbor,
                               const std::shared_ptr<Cluster>& cluster,
                               const int neighbor_radius_row,
                               const int neighbor_radius_col) const {
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

std::set<std::pair<int, int>>
PolarDetector::getNeighbors(const int col, const int row, const size_t distance, const bool) const {
    // Vector to hold the neighbor indices
    std::set<std::pair<int, int>> neighbors;

    // Position of the global seed in polar coordinates
    auto seed_pol = getPolarPosition(col, row);

    // Iterate over eligible strip rows
    for(int y = static_cast<int>(-distance); y <= static_cast<int>(distance); y++) {
        // Skip row if outside of strip matrix
        if(!isWithinMatrix(0, row + y)) {
            continue;
        }

        // Set starting position of a row seed to the global seed position
        auto row_seed_r = seed_pol.r();

        // Move row seed position to the center of a requested row
        for(unsigned int shift_y = 1; shift_y <= std::labs(y); shift_y++) {
            // Add or subtract position based on whether given row is below or above global seed
            row_seed_r += (y < 0) ? -strip_length.at(static_cast<unsigned int>(row) - shift_y + 1) / 2 -
                                        strip_length.at(static_cast<unsigned int>(row) - shift_y) / 2
                                  : strip_length.at(static_cast<unsigned int>(row) + shift_y - 1) / 2 +
                                        strip_length.at(static_cast<unsigned int>(row) + shift_y) / 2;
        }

        // Get cartesian position and pixel indices of the row seed
        auto row_seed = getLocalPosition({row_seed_r, 0, seed_pol.phi()});
        auto row_seed_x = static_cast<int>(getColumn(row_seed));
        auto row_seed_y = static_cast<int>(getRow(row_seed));

        // Iterate over potential neighbors of the row seed
        for(int j = static_cast<int>(-distance); j <= static_cast<int>(distance); j++) {
            // Add to final neighbors if strip is within the pixel matrix
            if(isWithinMatrix(row_seed_x + j, row_seed_y)) {
                neighbors.insert({row_seed_x + j, row_seed_y});
            }
        }
    }

    return {neighbors.begin(), neighbors.end()};
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
