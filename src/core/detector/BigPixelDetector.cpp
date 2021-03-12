/** @file
 *  @brief Detector model class
 *  @copyright Copyright (c) 2017-2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "BigPixelDetector.hpp"
#include "core/utils/log.h"
#include "exceptions.h"

using namespace ROOT::Math;
using namespace corryvreckan;

BigPixelDetector::BigPixelDetector(const Configuration& config) : PixelDetector(config) {
    LOG(DEBUG) << "Constructing a BIG Pixel Detector";
    // Get the information on big pixels here?
}

// Function to get local position from row and column
PositionVector3D<Cartesian3D<double>> BigPixelDetector::getLocalPosition(double column, double row) const {

    // FIXME: Replace with new coordinate transformation
    return PositionVector3D<Cartesian3D<double>>(m_pitch.X() * (column - static_cast<double>(m_nPixels.X() - 1) / 2.),
                                                 m_pitch.Y() * (row - static_cast<double>(m_nPixels.Y() - 1) / 2.),
                                                 0.);
}
