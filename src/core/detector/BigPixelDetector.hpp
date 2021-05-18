/** @file
 *  @brief Detector model class
 *  @copyright Copyright (c) 2017-2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#ifndef CORRYVRECKAN_BIGPLANARDETECTOR_H
#define CORRYVRECKAN_BIGPLANARDETECTOR_H

#include "Detector.hpp"
#include "core/config/Configuration.hpp"
#include "core/utils/ROOT.h"
#include "core/utils/log.h"

namespace corryvreckan {

    class BigPixelDetector : public PixelDetector {
    public:
        /**
         * Delete default constructor
         */
        BigPixelDetector() = delete;

        /**
         * Default destructor
         */
        ~BigPixelDetector() = default;

        /**
         * @brief Constructs a detector in the geometry
         * @param config Configuration object describing the detector
         */
        BigPixelDetector(const Configuration& config);

        void config_bigpixel(const Configuration& config);

        double getRow(const PositionVector3D<Cartesian3D<double>> localPosition) const override;
        double getColumn(const PositionVector3D<Cartesian3D<double>> localPosition) const override;

        // Function to get local position from column (x) and row (y) coordinates
        PositionVector3D<Cartesian3D<double>> getLocalPosition(double column, double row) const override;

        ROOT::Math::XYVector inPixel(const double column, const double row) const override;

        ROOT::Math::XYVector getSize() const override;

    private:
        std::vector<int> big_pixel_x{};
        std::vector<int> big_pixel_y{};
        std::vector<int> transformed_big_pixel_x{};
        std::vector<int> transformed_big_pixel_y{};
    };

} // namespace corryvreckan

#endif // CORRYVRECKAN_BIGPLANARDETECTOR_H
