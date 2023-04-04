/** @file
 *  @brief Detector model class
 *  @copyright Copyright (c) 2017-2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#ifndef CORRYVRECKAN_BIGPLANARDETECTOR_H
#define CORRYVRECKAN_BIGPLANARDETECTOR_H

#include "Detector.hpp"
#include "core/config/Configuration.hpp"
#include "core/utils/ROOT.h"
#include "core/utils/log.h"

namespace corryvreckan {

    class PixelModuleDetector : public PixelDetector {
    public:
        /**
         * Delete default constructor
         */
        PixelModuleDetector() = delete;

        /**
         * Default destructor
         */
        ~PixelModuleDetector() = default;

        /**
         * @brief Constructs a detector in the geometry
         * @param config Configuration object describing the detector
         */
        PixelModuleDetector(const Configuration& config);

        void config_bigpixel(const Configuration& config);

        double getRow(const PositionVector3D<Cartesian3D<double>> localPosition) const override;
        double getColumn(const PositionVector3D<Cartesian3D<double>> localPosition) const override;

        // Function to get local position from column (x) and row (y) coordinates
        PositionVector3D<Cartesian3D<double>> getLocalPosition(double column, double row) const override;

        /**
         * @brief Get intrinsic spatial resolution of the detector
         * @return Intrinsic spatial resolution in X and Y
         *
         * @note For a detector with variable pixel sizes this declaration could be changed to take column and row pixel
         * indices to calculate the resolution for a specific pixel
         */
        XYVector getSpatialResolution(double, double) const override;

        /**
         * @brief Get intrinsic spatial resolution in global coordinates of the detector
         * @return Intrinsic spatial resolution in global X and Y
         */
        TMatrixD getSpatialResolutionMatrixGlobal(double, double) const override;

        XYVector getSize() const override;

        /**
         * @brief Retrieve configuration object from detector, containing all (potentially updated) parameters
         * @return Configuration object for this detector
         */
        Configuration getConfiguration() const override;

    private:
        std::vector<std::vector<int>> m_big_pixel;
        std::vector<unsigned int> big_pixel_x{};
        std::vector<unsigned int> big_pixel_y{};
        std::vector<unsigned int> transformed_big_pixel_x{};
        std::vector<unsigned int> transformed_big_pixel_y{};
        XYVector m_big_pixel_spatial_resolution{};
    };

} // namespace corryvreckan

#endif // CORRYVRECKAN_BIGPLANARDETECTOR_H
