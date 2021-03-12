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

        // Function to get local position from column (x) and row (y) coordinates
        PositionVector3D<Cartesian3D<double>> getLocalPosition(double column, double row) const override;
    };

} // namespace corryvreckan

#endif // CORRYVRECKAN_BIGPLANARDETECTOR_H
