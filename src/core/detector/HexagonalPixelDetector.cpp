/**
 * @file
 * @brief Implementation of hexagonal pixel detector model
 *
 * @copyright Copyright (c) 2021-2022 CERN and the Allpix Squared authors.
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

#include "HexagonalPixelDetector.hpp"
#include "core/utils/log.h"
#include "exceptions.h"

using namespace ROOT::Math;
using namespace corryvreckan;

HexagonalPixelDetector::HexagonalPixelDetector(const Configuration& config) : PixelDetector(config) {
    // pitch_y is not along the y axis but tilted by 30°
    // only supports identical pitch values for now
    if(m_pitch.X() != m_pitch.Y()) {
        throw InvalidValueError(config, "pixel_pitch", "pitch_x != pitch_y is not supported");
    }

    m_height = 2.0 / std::sqrt(3) * m_pitch.X();
}

// Function to check if a track intercepts with a plane
bool HexagonalPixelDetector::hasIntercept(const Track* track, double /*pixelTolerance*/) const {

    // First, get the track intercept in global coordinates with the plane
    PositionVector3D<Cartesian3D<double>> globalIntercept = this->getIntercept(track);

    // Convert to local coordinates
    PositionVector3D<Cartesian3D<double>> localIntercept = alignment_->global2local() * globalIntercept;

    // Get the row and column numbers
    auto hex = getInterceptPixel(localIntercept);

    bool intercept = true;

    if(hex.second < 0 || hex.second >= m_nPixels.Y() || hex.first + hex.second / 2 < 0 ||
       hex.first + hex.second / 2 >= m_nPixels.X()) {
        intercept = false;
    }

    return intercept;
}

// Function to check if a track goes through/near a masked pixel
bool HexagonalPixelDetector::hitMasked(const Track* track, int tolerance) const {

    // First, get the track intercept in global coordinates with the plane
    PositionVector3D<Cartesian3D<double>> globalIntercept = this->getIntercept(track);

    // Convert to local coordinates
    PositionVector3D<Cartesian3D<double>> localIntercept = alignment_->global2local() * globalIntercept;

    // Get the row and column numbers
    auto pos = getInterceptPixel(localIntercept);

    int column = pos.first;
    int row = pos.second;

    // Check if the pixels around this pixel are masked
    bool hitmasked = false;
    for(int r = (row - tolerance); r <= (row + tolerance); r++) {
        for(int c = (column - tolerance); c <= (column + tolerance); c++) {
            if(this->masked(c, r))
                hitmasked = true;
        }
    }

    return hitmasked;
}

// Functions to get row and column from local position
double HexagonalPixelDetector::getRow(const PositionVector3D<Cartesian3D<double>> localPosition) const {
    auto size = getSize();

    double y = localPosition.Y() + 0.5 * (size.Y() - m_height);
    double row = (y * 2.0 / std::sqrt(3)) / m_pitch.Y();

    return row;
}

double HexagonalPixelDetector::getColumn(const PositionVector3D<Cartesian3D<double>> localPosition) const {
    auto size = getSize();

    double x = localPosition.X() + 0.5 * (size.X() - m_pitch.X());
    double y = localPosition.Y() + 0.5 * (size.Y() - m_height);
    double column = (x - y / std::sqrt(3)) / m_pitch.X();

    return column;
}

// Function to get local position from row and column
PositionVector3D<Cartesian3D<double>> HexagonalPixelDetector::getLocalPosition(double column, double row) const {

    auto size = getSize();

    // Offset by 1/2 width/height of matrix and 1/2 pixel to properly center matrix on 0,0
    double x = 0.5 * (size.X() - m_pitch.X());
    double y = 0.5 * (size.Y() - m_height);

    return PositionVector3D<Cartesian3D<double>>(
        (1.0 * column + 0.5 * row) * m_pitch.X() - x, (0.0 * column + std::sqrt(3) * 0.5 * row) * m_pitch.Y() - y, 0.);
}

// Function to get row and column of pixel
std::pair<int, int> HexagonalPixelDetector::getInterceptPixel(PositionVector3D<Cartesian3D<double>> localPosition) const {
    return round_to_nearest_hex(getColumn(localPosition), getRow(localPosition));
}

// Function to get in-pixel position
ROOT::Math::XYVector HexagonalPixelDetector::inPixel(const double column, const double row) const {
    auto hex = round_to_nearest_hex(column, row);
    double c = column - hex.first;
    double r = row - hex.second;
    return XYVector((1.0 * c + 0.5 * r) * m_pitch.X(), (0.0 * c + std::sqrt(3) * 0.5 * r) * m_pitch.Y());
}

ROOT::Math::XYVector HexagonalPixelDetector::inPixel(const PositionVector3D<Cartesian3D<double>> localPosition) const {
    double column = getColumn(localPosition);
    double row = getRow(localPosition);
    return inPixel(column, row);
}

// Check if track position is within ROI:
bool HexagonalPixelDetector::isWithinROI(const Track* track) const {

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
bool HexagonalPixelDetector::isWithinROI(Cluster* cluster) const {

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

XYVector HexagonalPixelDetector::getSize() const {
    double matrix_x = (m_nPixels.X() + 0.5) * m_pitch.X();
    double matrix_y = (3.0 / 4.0 * m_nPixels.Y() + 1.0 / 4.0) * m_height;

    return XYVector(matrix_x, matrix_y);
}

// Check if a pixel touches any of the pixels in a cluster
bool HexagonalPixelDetector::isNeighbor(const std::shared_ptr<Pixel>& neighbor,
                                        const std::shared_ptr<Cluster>& cluster,
                                        const int /*neighbor_radius_row*/,
                                        const int neighbor_radius_col) const {
    for(auto pixel : cluster->pixels()) {
        // fixme: take column and row radius into account
        if(hex_distance(pixel->row(), pixel->column(), neighbor->row(), neighbor->column()) <=
           static_cast<size_t>(neighbor_radius_col)) {
            return true;
        }
    }
    return false;
}

// Rounding is more easy in cubic coordinates, so we need to reconstruct the third coordinate from the other two as z = - x -
// y:
std::pair<int, int> HexagonalPixelDetector::round_to_nearest_hex(double x, double y) const {
    auto q = static_cast<int>(std::round(x));
    auto r = static_cast<int>(std::round(y));
    auto s = static_cast<int>(std::round(-x - y));
    double q_diff = std::abs(q - x);
    double r_diff = std::abs(r - y);
    double s_diff = std::abs(s - (-x - y));
    if(q_diff > r_diff and q_diff > s_diff) {
        q = -r - s;
    } else if(r_diff > s_diff) {
        r = -q - s;
    }
    return {q, r};
}

/*
 * In an axial-coordinates hexagon grid, simply checking for x and y to be between 0 and number_of_pixels will create
 * a rhombus which does lack the upper-left pixels and which has surplus pixels at the upper-right corner. We
 * therefore need to check the allowed range along x as a function of the y coordinate. The integer division by two
 * ensures we allow for one more x coordinate every other row in y.
 */
bool HexagonalPixelDetector::isWithinMatrix(const int x, const int y) const {
    // Check the valid pixel indices - this depends on the orientation of the axial index coordinate system with respect to
    // the cartesian local coordinate system, so we need to allow different indices depending on the hexagon orientation:
    return !(y < 0 || y >= m_nPixels.y() || x < 0 - y / 2 || x >= m_nPixels.x() - y / 2);
}

// The distance between two hexagons in cubic coordinates is half the Manhattan distance. To use axial coordinates, we have
// to reconstruct the third coordinate z = - x - y:
size_t HexagonalPixelDetector::hex_distance(double x1, double y1, double x2, double y2) const {
    return static_cast<size_t>(std::abs(x1 - x2) + std::abs(y1 - y2) + std::abs(-x1 - y1 + x2 + y2)) / 2;
}

std::set<std::pair<int, int>>
HexagonalPixelDetector::getNeighbors(const int col, const int row, const size_t distance, const bool) const {
    std::set<std::pair<int, int>> neighbors;

    for(int x = col - static_cast<int>(distance); x <= col + static_cast<int>(distance); x++) {
        for(int y = row - static_cast<int>(distance); y <= row + static_cast<int>(distance); y++) {
            if(hex_distance(col, row, x, y) <= distance) {

                if(!isWithinMatrix(x, y)) {
                    continue;
                }
                neighbors.insert({x, y});
            }
        }
    }
    return neighbors;
}
