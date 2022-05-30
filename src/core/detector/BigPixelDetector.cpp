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

    // Set detector position and direction from configuration file
    SetPostionAndOrientation(config);

    // initialize transform
    this->initialise();

    // Auxiliary devices don't have: number_of_pixels, pixel_pitch, spatial_resolution, mask_file, region-of-interest
    if(!isAuxiliary()) {
        config_bigpixel(config);
        build_axes(config);
    }

    // Printing out some test point for debugging BigPixel conversion

    m_coordinates = "cartesian_big";
}

void BigPixelDetector::config_bigpixel(const Configuration& config) {

    m_big_pixel = config.getMatrix<int>("big_pixel");
    big_pixel_x.assign(m_big_pixel.at(0).begin(), m_big_pixel.at(0).end());
    big_pixel_y.assign(m_big_pixel.at(1).begin(), m_big_pixel.at(1).end());

    LOG(INFO) << "Numbers of Big Pixels in X : " << big_pixel_x.size();
    LOG(INFO) << "Numbers of Big Pixels in Y : " << big_pixel_y.size();

    // transformed big pixel : treating big pixel as 2 regular pixels
    for(unsigned int i = 0; i < big_pixel_x.size(); i++) {
        transformed_big_pixel_x.push_back(big_pixel_x[i] + i);
        transformed_big_pixel_x.push_back(big_pixel_x[i] + i + 1);
    }
    for(unsigned int i = 0; i < big_pixel_y.size(); i++) {
        transformed_big_pixel_y.push_back(big_pixel_y[i] + i);
        transformed_big_pixel_y.push_back(big_pixel_y[i] + i + 1);
    }
    LOG(INFO) << "Numbers of transformed Big Pixels in X : " << transformed_big_pixel_x.size();
    LOG(INFO) << "Numbers of transformed Big Pixels in Y : " << transformed_big_pixel_y.size();

    for(auto i = transformed_big_pixel_x.begin(); i != transformed_big_pixel_x.end(); ++i) {
        LOG(DEBUG) << "Transform big pixel vector in X : " << *i << ' ';
    }

    for(auto i = transformed_big_pixel_y.begin(); i != transformed_big_pixel_y.end(); ++i) {
        LOG(DEBUG) << "Transform big pixel vector in Y : " << *i << ' ';
    }
}

// Functions to get row and column from local position// FIXME: Replace with new coordinate transformation
double BigPixelDetector::getRow(const PositionVector3D<Cartesian3D<double>> localPosition) const {

    int n_big_y_left = 0;
    bool is_big_y_pixel = 0;
    double row = 0;

    double tempPosition = ((localPosition.Y() + getSize().Y() / 2.) / m_pitch.Y()) - 0.5;

    for(unsigned int i = 0; i < transformed_big_pixel_y.size(); i++) {
        if(transformed_big_pixel_y[i] <= tempPosition) {
            n_big_y_left += 1;
        }
        if(fabs(tempPosition - static_cast<double>(transformed_big_pixel_y[i])) < 0.5) {
            is_big_y_pixel += 1;
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
    // LOG(INFO) << "Pos: " << localPosition << "; row : " << row ;

    return row;
}

double BigPixelDetector::getColumn(const PositionVector3D<Cartesian3D<double>> localPosition) const {

    int n_big_x_left = 0;
    bool is_big_x_pixel = 0;
    double column = 0;

    double tempPosition = ((localPosition.X() + getSize().X() / 2.) / m_pitch.X()) - 0.5;

    for(unsigned int i = 0; i < transformed_big_pixel_x.size(); i++) {
        if(transformed_big_pixel_x[i] <= tempPosition) {
            n_big_x_left += 1;
        }
        if(fabs(tempPosition - static_cast<double>(transformed_big_pixel_x[i])) < 0.5) {
            is_big_x_pixel += 1;
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
PositionVector3D<Cartesian3D<double>> BigPixelDetector::getLocalPosition(double column, double row) const {

    int n_big_x_left = 0;
    bool is_big_x_pixel = 0;
    int n_big_y_left = 0;
    bool is_big_y_pixel = 0;
    double col_integer, row_integer;

    for(unsigned int i = 0; i < big_pixel_x.size(); i++) {
        if(big_pixel_x[i] + 1 <= column + 0.5) {
            n_big_x_left += 1;
        }
        if(fabs(column - big_pixel_x[i]) < 0.5) {
            is_big_x_pixel += 1;
        }
    }

    for(unsigned int i = 0; i < big_pixel_y.size(); i++) {
        if(big_pixel_y[i] + 1 <= row + 0.5) {
            n_big_y_left += 1;
        }
        if(fabs(row - big_pixel_y[i]) < 0.5) {
            is_big_y_pixel += 1;
        }
    }

    // LOG(INFO) << "n_big_left: "<<n_big_y_left;
    // LOG(INFO) << "is big pixel?: " << is_big_y_pixel;

    return PositionVector3D<Cartesian3D<double>>(
        m_pitch.X() * (column + 0.5 + n_big_x_left +
                       static_cast<double>(is_big_x_pixel ? std::modf((column + 0.5), &col_integer) : 0)) -
            getSize().X() / 2.,
        m_pitch.Y() *
                (row + 0.5 + n_big_y_left + static_cast<double>(is_big_y_pixel ? std::modf((row + 0.5), &row_integer) : 0)) -
            getSize().Y() / 2.,
        0.);
}

// Function to get in-pixel position
ROOT::Math::XYVector BigPixelDetector::inPixel(const double column, const double row) const {
    // FIXME: Replace with new coordinate transformation
    // a pixel ranges from (col-0.5) to (col+0.5)
    return XYVector(m_pitch.X() * (column - floor(column + 0.5)), m_pitch.Y() * (row - floor(row + 0.5)));
}

ROOT::Math::XYVector BigPixelDetector::getSize() const {
    return XYVector(m_pitch.X() * (m_nPixels.X() + static_cast<double>(big_pixel_x.size())),
                    m_pitch.Y() * (m_nPixels.Y() + static_cast<double>(big_pixel_y.size())));
}
