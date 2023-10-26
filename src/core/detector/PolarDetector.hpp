/**
 * @file
 * @brief Polar detector model class
 *
 * @copyright Copyright (c) 2017-2023 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#ifndef CORRYVRECKAN_POLARDETECTOR_H
#define CORRYVRECKAN_POLARDETECTOR_H

#include <fstream>
#include <map>
#include <string>

#include <Math/DisplacementVector2D.h>
#include <Math/Vector2D.h>
#include <Math/Vector3D.h>
#include "Math/Transform3D.h"
#include "Math/Vector3D.h"

#include "Detector.hpp"
#include "core/config/Configuration.hpp"
#include "core/utils/ROOT.h"
#include "core/utils/log.h"
#include "objects/Cluster.hpp"
#include "objects/Pixel.hpp"
#include "objects/Track.hpp"

namespace corryvreckan {

    /**
     * @brief PolarDetector representation derived from Detector interface in the reconstruction chain
     *
     * Contains the PolarDetector with all its specific properties such as angular pitches, stereo angle and polar coordinate
     * system
     */
    class PolarDetector : public Detector {
    public:
        /**
         * Delete default constructor
         */
        PolarDetector() = delete;

        /**
         * Default destructor
         */
        ~PolarDetector() override = default;

        /**
         * @brief Constructs a detector in the geometry
         * @param config Configuration object describing the detector
         */
        explicit PolarDetector(const Configuration& config);

        /**
         * @brief Mark a detector channel as masked
         * @param chX X coordinate of the pixel to be masked
         * @param chY Y coordinate of the pixel to be masked
         */
        void maskChannel(int chX, int chY) override;

        /**
         * @brief Check if a detector channel is masked
         * @param chX X coordinate of the pixel to check
         * @param chY Y coordinate of the pixel to check
         * @return Mask status of the pixel in question
         */
        bool masked(int chX, int chY) const override;

        /**
         * @brief Check if a track goes through/near a masked strip
         * @param track Track to check
         * @param tolerance Number of surrounding strips around a masked one that are also considered masked
         */
        bool hitMasked(const Track* track, int tolerance = 0.) const override;

        /**
         * @brief Get global intercept of a track with the detector
         * @param track Track that intercepts the detector
         * @return Global position of the intercept
         */
        PositionVector3D<Cartesian3D<double>> getIntercept(const Track* track) const override;

        /**
         * @brief Get the local intercept of a track with the detector
         * @param track Track that intercepts the detector
         * @return Local position of the track intercept
         */
        PositionVector3D<Cartesian3D<double>> getLocalIntercept(const Track* track) const override;

        /**
         * @brief Get the row and column of a strip
         * @param localPosition Local position in cartesian coordinates
         * @return Strip indices
         */
        std::pair<int, int> getInterceptPixel(PositionVector3D<Cartesian3D<double>> localPosition) const override {
            return {floor(getColumn(localPosition)), floor(getRow(localPosition))};
        }

        /**
         * @brief Check if a track intercepts with a detector plane
         * @param track Track to check
         * @param stripTolerance The strip matrix is virtually extended in both directions by this number for intercept
         * evaluation
         * @return True if track intercepts the detector plane, false otherwise
         */
        bool hasIntercept(const Track* track, double stripTolerance = 0.0) const override;

        /**
         * @brief Get the strip row from local cartesian position
         * @param localPosition Local position in cartesian coordinates
         * @return Strip row (y-index)
         */
        double getRow(PositionVector3D<Cartesian3D<double>> localPosition) const override;

        /**
         * @brief Get the strip column from local cartesian position
         * @param localPosition Local position in cartesian coordinates
         * @return Strip column (x-index)
         */
        double getColumn(PositionVector3D<Cartesian3D<double>> localPosition) const override;

        /**
         * @brief Get the local position in cartesian coordinates
         * @param column Strip column (x-index)
         * @param row Strip row (y-index)
         * @return Local position in cartesian coordinates
         */
        PositionVector3D<Cartesian3D<double>> getLocalPosition(double column, double row) const override;

        /**
         * @brief Convert the local position in polar coordinates to cartesian coordinates
         * @param polarPosition Local position in polar coordinates
         * @return Local position in cartesian coordinates
         */
        PositionVector3D<Cartesian3D<double>> getLocalPosition(const PositionVector3D<Polar3D<double>> polarPosition) const;

        /**
         * @brief Get the local position in polar coordinates
         * @param column Strip column (x-index)
         * @param row Strip row (y-index)
         * @return Local position in polar coordinates
         */
        PositionVector3D<Polar3D<double>> getPolarPosition(double column, double row) const;

        /**
         * @brief Converts the local position in cartesian coordinates to polar coordinates
         * @param localPosition Local position in cartesian coordinates
         * @return Local position in polar coordinates
         *
         * @note The polar coordinates are defined in the strip frame where:
         *  - R is measured from the polar origin
         *  - Phi is measured from the strip focal point
         */
        PositionVector3D<Polar3D<double>> getPolarPosition(const PositionVector3D<Cartesian3D<double>> localPosition) const;

        /**
         * @brief Transform from local (sensor) coordinates to in-pixel coordinates
         * @param column Column address ranging from int_column-0.5*pitch to int_column+0.5*pitch
         * @param row Row address ranging from int_column-0.5*pitch to int_column+0.5*pitch
         * @return Position within a single pixel cell, given in units of length
         */
        XYVector inPixel(const double column, const double row) const override;

        /**
         * @brief Transformation from local (sensor) coordinates to in-pixel coordinates
         * @param localPosition Local position on the sensor
         * @return Position within a single pixel cell, given in units of length
         */
        XYVector inPixel(PositionVector3D<Cartesian3D<double>> localPosition) const override;

        /**
         * @brief Check whether given track is within the detector's region-of-interest
         * @param track The track to be checked
         * @return Boolean indicating cluster affiliation with region-of-interest
         */
        bool isWithinROI(const Track* track) const override;

        /**
         * @brief Check whether given cluster is within the detector's region-of-interest
         * @param cluster The cluster to be checked
         * @return Boolean indicating cluster affiliation with region-of-interest
         */
        bool isWithinROI(Cluster* cluster) const override;

        /**
         * @brief Get the total size of the active matrix
         * @return 2D vector with the dimensions of the strip matrix in X and Y
         *
         * @note Due to the shape of the polar detectors, the returned dimensions are only approximate
         */
        XYVector getSize() const override;

        /**
         * @brief Get pitch of a single pixel
         * @return Pitch of a pixel in X and Y
         *
         * @note In polar detectors, strip pitches can vary from row to row. This function returns the largest pitches found
         * on the detector.
         */
        XYVector getPitch() const override;

        /**
         * @brief Get intrinsic spatial resolution of the detector
         * @param column Strip column (x-index) to calculate the resolution for
         * @param row Strip row (y-index) to calculate the resolution for
         * @return Intrinsic spatial resolution in polar dimensions (Phi and R)
         */
        XYVector getSpatialResolution(double column = 0, double row = 0) const override;

        /**
         * @brief Get intrinsic spatial resolution in global coordinates of the detector
         * @return Intrinsic spatial resolution in global X and Y
         */
        virtual TMatrixD getSpatialResolutionMatrixGlobal(double, double) const override {
            return m_spatial_resolution_matrix_global;
        }

        /**
         * @brief Get number of pixels in x and y
         * @return Number of two-dimensional pixels
         *
         * @note For polar detectors, the numbers are defined as
         * - x: the largest number of strips in any row
         * - y: the number of strip rows
         */
        ROOT::Math::DisplacementVector2D<ROOT::Math::Cartesian2D<int>> nPixels() const override {
            return {static_cast<int>(*std::max_element(number_of_strips.begin(), number_of_strips.end())),
                    static_cast<int>(number_of_strips.size())};
        }

        /**
         * @brief Test whether one pixel touches the cluster
         * @return true if it fulfills the condition
         * @note users should define their specific clustering method in the detector class, for pixel detector, the default
         * is 2D clustering
         */
        bool isNeighbor(const std::shared_ptr<Pixel>&, const std::shared_ptr<Cluster>&, const int, const int) const override;

        /**
         * @brief Checks if a given strip index lies within the strip matrix of the detector
         * @param col Strip column (x-index)
         * @param row Strip row (y-index)
         * @return True if strip index is within matrix bounds, false otherwise
         */
        bool isWithinMatrix(const int col, const int row) const override {
            return !(row < 0 || row >= nPixels().y() || col < 0 ||
                     col >= static_cast<int>(number_of_strips.at(static_cast<unsigned int>(row))));
        }

        /**
         * @brief Return a set containing all strips neighboring the given one with a configurable maximum distance
         * @param col       Column of strip in question
         * @param row       Row of strip in question
         * @param distance  Distance for strips to be considered neighbors
         * @param include_corners Boolean to select whether strips only touching via corners should be returned
         * @return Set of neighboring strip indices, including the initial strip
         *
         * @note The returned set should always also include the initial pixel indices the neighbors are calculated for
         */
        std::set<std::pair<int, int>>
        getNeighbors(const int col, const int row, const size_t distance, const bool include_corners) const override;

    private:
        // Build axis, for devices which are not auxiliary
        void build_axes(const Configuration& config) override;

        // Config detector, for devices which are not auxiliary
        void configure_detector(Configuration& config) const override;

        // Config position, orientation, mode of detector
        void configure_pos_and_orientation(Configuration& config) const override;

        // Functions to set and check channel masking
        void process_mask_file() override;

        // Seems to be used in other coordinate
        inline static int isLeft(std::pair<int, int> pt0, std::pair<int, int> pt1, std::pair<int, int> pt2);
        static int winding_number(std::pair<int, int> probe, std::vector<std::vector<int>> polygon);

        // For planar detector
        TMatrixD m_spatial_resolution_matrix_global{3, 3};
        std::vector<std::vector<int>> m_roi{};

        // For polar detectors
        std::vector<unsigned int> number_of_strips{};
        std::vector<double> row_radius{};
        std::vector<double> strip_length{};
        std::vector<double> angular_pitch{};
        double stereo_angle{};
        double center_radius{};
    };
} // namespace corryvreckan

#endif // CORRYVRECKAN_POLARDETECTOR_H
