/** @file
 *  @brief Detector model class
 *  @copyright Copyright (c) 2017-2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "PixelModuleDetector.hpp"
#include "core/utils/log.h"
#include "exceptions.h"

using namespace ROOT::Math;
using namespace corryvreckan;

PixelModuleDetector::PixelModuleDetector(const Configuration& config) : PixelDetector(config) {

    // Auxiliary devices don't have: number_of_pixels, pixel_pitch, spatial_resolution, mask_file, region-of-interest
    if(!isAuxiliary()) {
        config_bigpixel(config);
        build_axes(config);
    }

    if(isDUT()) {
        LOG(ERROR) << "A PixelModuleDetector was configured to be a DUT (" << this->getName()
                   << "). Please be aware that regions with large pixels might be treated unexpectedly.";
    }
}

void PixelModuleDetector::config_bigpixel(const Configuration& config) {

    m_big_pixel = config.getMatrix<int>("big_pixels");
    big_pixel_x.assign(m_big_pixel.at(0).begin(), m_big_pixel.at(0).end());
    big_pixel_y.assign(m_big_pixel.at(1).begin(), m_big_pixel.at(1).end());

    m_big_pixel_spatial_resolution =
        config.get<ROOT::Math::XYVector>("big_pixel_spatial_resolution", 2. * m_spatial_resolution);

    LOG(INFO) << "Numbers of Big Rows (X) : " << big_pixel_x.size();
    LOG(INFO) << "Numbers of Big Columns (Y) : " << big_pixel_y.size();

    // sort big_pixel
    sort(big_pixel_x.begin(), big_pixel_x.end());
    sort(big_pixel_y.begin(), big_pixel_y.end());

    // transformed big pixel : treating big pixel as 2 regular pixels
    for(unsigned int i = 0; i < big_pixel_x.size(); i++) {
        transformed_big_pixel_x.push_back(big_pixel_x[i] + i);
        transformed_big_pixel_x.push_back(big_pixel_x[i] + i + 1);
    }
    for(unsigned int i = 0; i < big_pixel_y.size(); i++) {
        transformed_big_pixel_y.push_back(big_pixel_y[i] + i);
        transformed_big_pixel_y.push_back(big_pixel_y[i] + i + 1);
    }
    LOG(DEBUG) << "Numbers of transformed Big Rows (X) : " << transformed_big_pixel_x.size();
    LOG(DEBUG) << "Numbers of transformed Big Columns (Y) : " << transformed_big_pixel_y.size();

    for(auto i = transformed_big_pixel_x.begin(); i != transformed_big_pixel_x.end(); ++i) {
        LOG(DEBUG) << "Transform big pixel vector in X : " << *i;
    }

    for(auto i = transformed_big_pixel_y.begin(); i != transformed_big_pixel_y.end(); ++i) {
        LOG(DEBUG) << "Transform big pixel vector in Y : " << *i;
    }
}

// Functions to get row and column from local position// FIXME: Replace with new coordinate transformation
double PixelModuleDetector::getRow(const PositionVector3D<Cartesian3D<double>> localPosition) const {

    int n_big_y_left = 0;
    bool is_big_y_pixel = 0;
    double row = 0;

    double tempPosition = ((localPosition.Y() + getSize().Y() / 2.) / m_pitch.Y()) - 0.5;

    for(unsigned int i = 0; i < transformed_big_pixel_y.size(); i++) {
        if(transformed_big_pixel_y[i] <= tempPosition) {
            n_big_y_left += 1;
            if(fabs(tempPosition - static_cast<double>(transformed_big_pixel_y[i])) < 0.5) {
                is_big_y_pixel = true;
            }
        } else {
            if(fabs(tempPosition - static_cast<double>(transformed_big_pixel_y[i])) < 0.5) {
                is_big_y_pixel = true;
            }
            break;
        }
    }

    if(is_big_y_pixel == true) {
        for(unsigned int i = 0; i < transformed_big_pixel_y.size(); i = i + 2) {
            if(fabs(tempPosition - transformed_big_pixel_y[i]) <= 2) {
                row = (tempPosition - transformed_big_pixel_y[i]) / 2. + big_pixel_y[i / 2] - 0.25;
            }
        }
    } else {
        row = tempPosition - n_big_y_left / 2.;
    }

    return row;
}

double PixelModuleDetector::getColumn(const PositionVector3D<Cartesian3D<double>> localPosition) const {

    int n_big_x_left = 0;
    bool is_big_x_pixel = 0;
    double column = 0;

    double tempPosition = ((localPosition.X() + getSize().X() / 2.) / m_pitch.X()) - 0.5;

    for(unsigned int i = 0; i < transformed_big_pixel_x.size(); i++) {
        if(transformed_big_pixel_x[i] <= tempPosition) {
            n_big_x_left += 1;
            if(fabs(tempPosition - static_cast<double>(transformed_big_pixel_x[i])) < 0.5) {
                is_big_x_pixel = true;
            }
        } else {
            if(fabs(tempPosition - static_cast<double>(transformed_big_pixel_x[i])) < 0.5) {
                is_big_x_pixel = true;
            }
            break;
        }
    }

    if(is_big_x_pixel == true) {
        for(unsigned int i = 0; i < transformed_big_pixel_x.size(); i = i + 2) {
            if(abs(tempPosition - transformed_big_pixel_x[i] - 0.5) <= 2) {
                column = (tempPosition - transformed_big_pixel_x[i]) / 2. + big_pixel_x[i / 2] - 0.25;
            }
        }
    } else {
        column = tempPosition - n_big_x_left / 2.;
    }

    return column;
}

// Function to get local position from row and column
PositionVector3D<Cartesian3D<double>> PixelModuleDetector::getLocalPosition(double column, double row) const {

    int n_big_x_left = 0;
    bool is_big_x_pixel = 0;
    int n_big_y_left = 0;
    bool is_big_y_pixel = 0;
    double col_integer, row_integer;

    for(unsigned int i = 0; i < big_pixel_x.size(); i++) {
        if(big_pixel_x[i] <= column - 0.5) {
            n_big_x_left += 1;
            if(fabs(column - big_pixel_x[i]) < 0.5) {
                is_big_x_pixel = true;
            }
        } else {
            if(fabs(column - big_pixel_x[i]) < 0.5) {
                is_big_x_pixel = true;
            }
            break;
        }
    }

    for(unsigned int i = 0; i < big_pixel_y.size(); i++) {
        if(big_pixel_y[i] <= row - 0.5) {
            n_big_y_left += 1;
            if(fabs(row - big_pixel_y[i]) < 0.5) {
                is_big_y_pixel = true;
            }
        } else {
            if(fabs(row - big_pixel_y[i]) < 0.5) {
                is_big_y_pixel = true;
            }
            break;
        }
    }

    return PositionVector3D<Cartesian3D<double>>(
        m_pitch.X() * (column + 0.5 + n_big_x_left + (is_big_x_pixel ? std::modf((column + 0.5), &col_integer) : 0)) -
            getSize().X() / 2.,
        m_pitch.Y() * (row + 0.5 + n_big_y_left + (is_big_y_pixel ? std::modf((row + 0.5), &row_integer) : 0)) -
            getSize().Y() / 2.,
        0.);
}

XYVector PixelModuleDetector::getSize() const {
    return XYVector(m_pitch.X() * (m_nPixels.X() + static_cast<double>(big_pixel_x.size())),
                    m_pitch.Y() * (m_nPixels.Y() + static_cast<double>(big_pixel_y.size())));
}

XYVector PixelModuleDetector::getSpatialResolution(double column = 0, double row = 0) const {
    bool is_big_x_pixel = 0;
    bool is_big_y_pixel = 0;

    for(unsigned int i = 0; i < big_pixel_x.size(); i++) {
        if(fabs(column - big_pixel_x[i]) < 0.5) {
            is_big_x_pixel = true;
            break;
        }
    }

    for(unsigned int i = 0; i < big_pixel_y.size(); i++) {
        if(fabs(row - big_pixel_y[i]) < 0.5) {
            is_big_y_pixel = true;
            break;
        }
    }

    double resolution_x = is_big_x_pixel ? m_big_pixel_spatial_resolution.x() : m_spatial_resolution.x();
    double resolution_y = is_big_y_pixel ? m_big_pixel_spatial_resolution.y() : m_spatial_resolution.y();
    return XYVector(resolution_x, resolution_y);
}

TMatrixD PixelModuleDetector::getSpatialResolutionMatrixGlobal(double column = 0, double row = 0) const {
    TMatrixD errorMatrix(3, 3);
    TMatrixD locToGlob(3, 3), globToLoc(3, 3);
    auto spatial_resolution = getSpatialResolution(column, row);
    errorMatrix(0, 0) = spatial_resolution.x() * spatial_resolution.x();
    errorMatrix(1, 1) = spatial_resolution.y() * spatial_resolution.y();
    alignment_->local2global().Rotation().GetRotationMatrix(locToGlob);
    alignment_->global2local().Rotation().GetRotationMatrix(globToLoc);
    return (locToGlob * errorMatrix * globToLoc);
}

Configuration PixelModuleDetector::getConfiguration() const {
    auto config = PixelDetector::getConfiguration();

    config.setMatrix("big_pixels", m_big_pixel);
    config.set<XYVector>("big_pixel_spatial_resolution", m_big_pixel_spatial_resolution);

    return config;
}
