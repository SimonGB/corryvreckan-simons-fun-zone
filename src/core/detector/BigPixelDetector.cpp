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

    LOG(INFO) << "Constructing a BIG Pixel Detector";

    m_coordinates = "cartesian_big";
}

void BigPixelDetector::config_bigpixel(const Configuration& config) {

    m_big_pixel = config.getMatrix<int>("big_pixel");
    big_pixel_x.assign(m_big_pixel.at(0).begin(), m_big_pixel.at(0).end());
    big_pixel_y.assign(m_big_pixel.at(1).begin(), m_big_pixel.at(1).end());

    LOG(INFO) << "Numbers of Big Pixels in X : " << big_pixel_x.size();
    LOG(INFO) << "Numbers of Big Pixels in Y : " << big_pixel_y.size();

    // transformed big pixel : treating big pixel as 2 regular pixels
    for(int i = 0; i < big_pixel_x.size(); i++) {
        transformed_big_pixel_x.push_back(big_pixel_x[i] + i);
        transformed_big_pixel_x.push_back(big_pixel_x[i] + i + 1);
    }
    for(int i = 0; i < big_pixel_y.size(); i++) {
        transformed_big_pixel_y.push_back(big_pixel_y[i] + i);
        transformed_big_pixel_y.push_back(big_pixel_y[i] + i + 1);
    }
    LOG(INFO) << "Numbers of transformed Big Pixels in X : " << transformed_big_pixel_x.size();
    LOG(INFO) << "Numbers of transformed Big Pixels in Y : " << transformed_big_pixel_y.size();
}

// Functions to get row and column from local position// FIXME: Replace with new coordinate transformation
double BigPixelDetector::getRow(const PositionVector3D<Cartesian3D<double>> localPosition) const {

    int n_big_y_left = 0;
    bool is_big_y_pixel = 0;
    double row = 0;

    double tempPosition = (localPosition.Y() + getSize().Y() / 2.) / m_pitch.Y();

    for(int i = 0; i < transformed_big_pixel_y.size(); i++) {
        if(transformed_big_pixel_y[i] <= tempPosition) {
            n_big_y_left += 1;
        }
        if(std::count(transformed_big_pixel_y.begin(), transformed_big_pixel_y.end(), (floor(row + 0.5)))) {
            is_big_y_pixel = true;
        }
    }

    if(is_big_y_pixel == true) {
        for(int i = 0; i < transformed_big_pixel_y.size(); i = i + 2) {
            if(abs(tempPosition - transformed_big_pixel_y[i] - 0.5) <= 2) {
                row = (tempPosition - transformed_big_pixel_y[i] - 0.5) / 2. + big_pixel_y[i / 2] - 0.5;
            }
        }
    } else {
        row = tempPosition - 0.5 - n_big_y_left / 2.;
    }

    return row;
}

double BigPixelDetector::getColumn(const PositionVector3D<Cartesian3D<double>> localPosition) const {

    int n_big_x_left = 0;
    bool is_big_x_pixel = 0;
    double column = 0;

    double tempPosition = (localPosition.X() + getSize().X() / 2.) / m_pitch.X();

    for(int i = 0; i < transformed_big_pixel_x.size(); i++) {
        if(transformed_big_pixel_x[i] <= tempPosition) {
            n_big_x_left += 1;
        }
        if(std::count(transformed_big_pixel_x.begin(), transformed_big_pixel_x.end(), (floor(column + 0.5)))) {
            is_big_x_pixel = true;
        }
    }

    if(is_big_x_pixel == true) {
        for(int i = 0; i < transformed_big_pixel_x.size(); i = i + 2) {
            if(abs(tempPosition - transformed_big_pixel_x[i] - 0.5) <= 2) {
                column = (tempPosition - transformed_big_pixel_x[i] - 0.5) / 2. + big_pixel_x[i / 2] - 0.5;
            }
        }
    } else {
        column = tempPosition - 0.5 - n_big_x_left / 2.;
    }

    return column;
}

// Function to get local position from row and column
PositionVector3D<Cartesian3D<double>> BigPixelDetector::getLocalPosition(double column, double row) const {

    int n_big_x_left = 0;
    bool is_big_x_pixel = 0;
    int n_big_y_left = 0;
    bool is_big_y_pixel = 0;

    for(int i = 0; i < big_pixel_x.size(); i++) {
        if(big_pixel_x[i] <= column) {
            n_big_x_left += 1;
        }
        if(std::count(big_pixel_x.begin(), big_pixel_x.end(), (floor(column + 0.5)))) {
            is_big_x_pixel = true;
        }
    }

    for(int i = 0; i < big_pixel_y.size(); i++) {
        if(big_pixel_y[i] <= row) {
            n_big_y_left += 1;
        }
        if(std::count(big_pixel_y.begin(), big_pixel_y.end(), (floor(row + 0.5)))) {
            is_big_y_pixel = true;
        }
    }

    return PositionVector3D<Cartesian3D<double>>(
        m_pitch.X() * (column + 1 + n_big_x_left - static_cast<double>(is_big_x_pixel ? 0 : 1) / 2.) - getSize().X() / 2.,
        m_pitch.Y() * (row + 1 + n_big_y_left - static_cast<double>(is_big_y_pixel ? 0 : 1) / 2.) - getSize().Y() / 2.,
        0.);

    // return PositionVector3D<Cartesian3D<double>>(m_pitch.X() * (column - static_cast<double>(m_nPixels.X() - 1) / 2.),
    //                                             m_pitch.Y() * (row - static_cast<double>(m_nPixels.Y() - 1) / 2.),
    //                                             0.);
}

// Function to get in-pixel position
ROOT::Math::XYVector BigPixelDetector::inPixel(const double column, const double row) const {
    // FIXME: Replace with new coordinate transformation
    // a pixel ranges from (col-0.5) to (col+0.5)
    return XYVector(m_pitch.X() * (column - floor(column + 0.5)), m_pitch.Y() * (row - floor(row + 0.5)));
}

ROOT::Math::XYVector BigPixelDetector::getSize() const {
    return XYVector(m_pitch.X() * (m_nPixels.X() + big_pixel_x.size()), m_pitch.Y() * (m_nPixels.Y() + big_pixel_y.size()));
}
