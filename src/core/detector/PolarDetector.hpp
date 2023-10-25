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
     * Contains the PolarDetector with all its properties such as position and orientation, pitch, spatial resolution
     * etc.
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
         * @return    Mask status of the pixel in question
         */
        bool masked(int chX, int chY) const override;

        // Function to get global intercept with a track
        PositionVector3D<Cartesian3D<double>> getIntercept(const Track* track) const override;
        // Function to get local intercept with a track
        PositionVector3D<Cartesian3D<double>> getLocalIntercept(const Track* track) const override;

        // Function to check if a track intercepts with a plane
        bool hasIntercept(const Track* track, double pixelTolerance = 0.) const override;

        // Function to check if a track goes through/near a masked pixel
        bool hitMasked(const Track* track, int tolerance = 0.) const override;

        // Functions to get row and column from local position
        double getRow(PositionVector3D<Cartesian3D<double>> localPosition) const override;

        double getColumn(PositionVector3D<Cartesian3D<double>> localPosition) const override;

        // Function to get local position from column (x) and row (y) coordinates
        PositionVector3D<Cartesian3D<double>> getLocalPosition(double column, double row) const override;

        /**
         * Transformation from local (sensor) coordinates to in-pixel coordinates
         * @param  column Column address ranging from int_column-0.5*pitch to int_column+0.5*pitch
         * @param  row Row address ranging from int_column-0.5*pitch to int_column+0.5*pitch
         * @return               Position within a single pixel cell, given in units of length
         */
        XYVector inPixel(const double column, const double row) const override;

        /**
         * Transformation from local (sensor) coordinates to in-pixel coordinates
         * @param  localPosition Local position on the sensor
         * @return               Position within a single pixel cell, given in units of length
         */
        XYVector inPixel(PositionVector3D<Cartesian3D<double>> localPosition) const override;

        /**
         * @brief Check whether given track is within the detector's region-of-interest
         * @param  track The track to be checked
         * @return       Boolean indicating cluster affiliation with region-of-interest
         */
        bool isWithinROI(const Track* track) const override;

        /**
         * @brief Check whether given cluster is within the detector's region-of-interest
         * @param  cluster The cluster to be checked
         * @return         Boolean indicating cluster affiliation with region-of-interest
         */
        bool isWithinROI(Cluster* cluster) const override;

        /**
         * @brief Get the total size of the active matrix, i.e. pitch * number of pixels in both dimensions
         * @return 2D vector with the dimensions of the pixle matrix in X and Y
         */
        XYVector getSize() const override;

        /**
         * @brief Get pitch of a single pixel
         * @return Pitch of a pixel
         */
        XYVector getPitch() const override { return m_pitch; }

        /**
         * @brief Get intrinsic spatial resolution of the detector
         * @return Intrinsic spatial resolution in X and Y
         */
        XYVector getSpatialResolution(double column = 0, double row = 0) const override;

        /**
         * @brief Get number of pixels in x and y
         * @return Number of two dimensional pixels
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
         * @brief Converts the local position in cartesian coordinates to polar coordinates
         * @param localPosition Position in local cartesian coordinates of the detector model
         * @return Local position in polar coordinates
         *
         * @note The polar coordinates are defined in a system where:
         *  - R is measured from the local coordinate center
         *  - Phi is measured from the strip focal point
         */
        PositionVector3D<Polar3D<double>> getPositionPolar(const PositionVector3D<Cartesian3D<double>> localPosition) const;

        /**
         * @brief Converts the position in polar coordinates to cartesian coordinates in the local frame.
         * @param polarPosition Position in local polar coordinates of the detector model
         * @return Local position in cartesian coordinates
         */
        PositionVector3D<Cartesian3D<double>>
        getPositionCartesian(const PositionVector3D<Polar3D<double>> polarPosition) const;

        /**
         * @brief Get the radius of the strip sensor center
         * @return Radius of the strip sensor center
         */
        double getCenterRadius() const { return center_radius; }

        // Function to get row and column of pixel
        std::pair<int, int> getInterceptPixel(PositionVector3D<Cartesian3D<double>> localPosition) const override;

        /**
         * @brief Checks if a given strip index lies within the strip matrix of the detector
         * @return True if strip index is within matrix bounds, false otherwise
         */
        bool isWithinMatrix(const int col, const int row) const override {
            return !(row < 0 || row >= nPixels().y() || col < 0 ||
                     col >= static_cast<int>(number_of_strips.at(static_cast<unsigned int>(row))));
        }

        /**
         * @brief Get intrinsic spatial resolution in global coordinates of the detector
         * @return Intrinsic spatial resolution in global X and Y
         */
        virtual TMatrixD getSpatialResolutionMatrixGlobal(double, double) const override {
            return m_spatial_resolution_matrix_global;
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

        std::vector<double> getRowRadius() { return row_radius; };
        std::vector<double> getAngularPitch() { return angular_pitch; };
        unsigned int nStrips(const int row) { return number_of_strips.at(static_cast<unsigned int>(row)); }

        PositionVector3D<Polar3D<double>> getLocalPolarPosition(double column, double row) const;

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
        XYVector m_pitch{};
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
