/**
 * @file
 * @brief Detector model class
 *
 * @copyright Copyright (c) 2017-2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#ifndef CORRYVRECKAN_DETECTOR_H
#define CORRYVRECKAN_DETECTOR_H

#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <string>

#include <Math/DisplacementVector2D.h>
#include <Math/RotationX.h>
#include <Math/RotationY.h>
#include <Math/RotationZ.h>
#include <Math/RotationZYX.h>
#include <Math/Transform3D.h>
#include <Math/Vector2D.h>
#include <Math/Vector3D.h>
#include <TFormula.h>
#include <TMatrixD.h>

#include "core/config/Configuration.hpp"
#include "core/utils/ROOT.h"
#include "core/utils/log.h"
#include "objects/Track.hpp"

namespace corryvreckan {
    using namespace ROOT::Math;

    /**
     * @brief Role of the detector
     */
    enum class DetectorRole : int {
        NONE = 0x0,           ///< No specific detector role
        REFERENCE = (1 << 0), ///< Reference detector
        DUT = (1 << 1),       ///< Detector used as device under test
        AUXILIARY = (1 << 2), ///< Auxiliary device which should not participate in regular reconstruction but might provide
                              /// additional information
        PASSIVE = (1 << 3),   ///< Passive device which only acts as scatterer. This is ignored for detector modules
    };

    /**
     * @brief Detector coordinate system
     */
    enum class Coordinates {
        CARTESIAN = 1,    ///< Cartesian coordinates
        HEXAGONAL,        ///< Hexagonal pixel coordinates, axial coordinate system
        CARTESIAN_MODULE, ///< Cartesian coordinates for detector modules
        POLAR,            ///< Polar coordinates
    };

    inline constexpr DetectorRole operator&(DetectorRole x, DetectorRole y) {
        return static_cast<DetectorRole>(static_cast<int>(x) & static_cast<int>(y));
    }

    inline constexpr DetectorRole operator|(DetectorRole x, DetectorRole y) {
        return static_cast<DetectorRole>(static_cast<int>(x) | static_cast<int>(y));
    }

    inline DetectorRole& operator&=(DetectorRole& x, DetectorRole y) {
        x = x & y;
        return x;
    }

    inline DetectorRole& operator|=(DetectorRole& x, DetectorRole y) {
        x = x | y;
        return x;
    }

    /**
     * @brief Detector interface in the reconstruction chain
     *
     * Contains the detector with common properties such as type, name, coordinate
     * etc.
     */
    class Detector {
        class Alignment {
        public:
            explicit Alignment(const Configuration& config);

            // Transforms from local to global and back
            const Transform3D& local2global() const { return local2global_; };

            const Transform3D& global2local() const { return global2local_; };

            // Normal to the detector surface and point on the surface
            const ROOT::Math::XYZVector& normal() const { return normal_; };

            const ROOT::Math::XYZPoint& origin() const { return origin_; };

            const ROOT::Math::XYZPoint& displacement() const { return displacement_; };

            const ROOT::Math::XYZVector& orientation() const { return orientation_; };

            const std::string& mode() const { return mode_; };

            void update(double time, bool force = false);

            void update(const ROOT::Math::XYZPoint& displacement, const ROOT::Math::XYZVector& orientation);

            bool isVariable() const { return needs_update_; };

        private:
            void recalculate();

            std::array<std::shared_ptr<TFormula>, 3> parse_formulae(const Configuration& config, const std::string& key);

            // Cache for the orientation mode:
            std::string mode_{};

            // Cache for last time the transformations were renewed, in ns:
            double last_time_{};
            double granularity_{};
            bool needs_update_{false};

            // Cached position and orientation
            ROOT::Math::XYZPoint displacement_;
            ROOT::Math::XYZVector orientation_;

            // Cache for calculated transformations
            ROOT::Math::XYZPoint origin_;
            ROOT::Math::XYZVector normal_;
            ROOT::Math::Transform3D local2global_;
            ROOT::Math::Transform3D global2local_;

            // The formulae
            std::array<std::shared_ptr<TFormula>, 3> formulae_pos_;

            // Function to generate rotation matrix, depending on mode
            std::function<ROOT::Math::Rotation3D(const ROOT::Math::XYZVector& rot)> rotation_fct_;
        };

    public:
        /**
         * Delete default constructor
         */
        Detector() = delete;

        /**
         * Default destructor
         */
        virtual ~Detector() = default;

        /**
         * @brief Constructs a detector in the geometry
         * @param config Configuration object describing the detector
         */
        explicit Detector(const Configuration& config);

        /**
         * @brief Factory to dynamically create detectors
         * @param config Configuration object describing the detector
         * @return shared_ptr that contains the real detector
         */
        static std::shared_ptr<Detector> factory(const Configuration& config);

        /**
         * @brief Get type of the detector
         * @return Type of the detector model
         */
        std::string getType() const;

        /**
         * @brief Get name of the detector
         * @return Detector name
         */
        std::string getName() const;

        /**
         * @brief Check whether detector is registered as reference
         * @return Reference status
         */
        bool isReference() const;

        /**
         * @brief Check whether detector is registered as DUT
         * @return DUT status
         */
        bool isDUT() const;

        /**
         * @brief Check whether detector is registered as auxiliary device and should not parttake in the reconstruction
         * @return Auxiliary status
         */
        bool isAuxiliary() const;
        /**
         * @brief Check whether detector is registered as Passive
         * @return Passive status
         */
        bool isPassive() const;
        /**
         * @brief Obtain roles assigned to this detector
         * @return List of detector roles
         */
        DetectorRole getRoles() const;

        /**
         * @brief Check if this detector has a certain role assigned
         * @param  role Role to be checked for
         * @return True if detector holds this role
         */
        bool hasRole(DetectorRole role) const;

        /**
         * @brief Retrieve configuration object from detector, containing all (potentially updated) parameters
         * @return Configuration object for this detector
         */
        virtual Configuration getConfiguration() const;

        /**
         * @brief Get the total size of the active matrix, i.e. pitch * number of pixels in both dimensions
         * @return 2D vector with the dimensions of the pixle matrix in X and Y
         * @todo: this is designed for PixelDetector, find a proper interface for other Detector type
         */
        virtual XYVector getSize() const = 0;

        /**
         * @brief Get pitch of a single pixel
         * @return Pitch of a pixel
         * @todo: this is designed for PixelDetector, find a proper interface for other Detector type
         */
        virtual XYVector getPitch() const = 0;

        /**
         * @brief Checks if a given pixel index lies within the pixel matrix of the detector
         * @return True if pixel index is within matrix bounds, false otherwise
         */
        virtual bool isWithinMatrix(const int col, const int row) const = 0;

        /**
         * @brief Get intrinsic spatial resolution of the detector
         * @param column Column (x-) index of the pixel to calculate the spatial resolution from
         * @param row Row (y-) index of the pixel to calculate the spatial resolution from
         * @return Intrinsic spatial resolution in X and Y
         * @todo: this is designed for PixelDetector, find a proper interface for other Detector type
         */
        virtual XYVector getSpatialResolution(double column = 0, double row = 0) const = 0;

        /**
         * @brief Get intrinsic spatial resolution in global coordinates of the detector
         * @param column Column (x-) index of the pixel to calculate the spatial resolution from
         * @param row Row (y-) index of the pixel to calculate the spatial resolution from
         * @return Intrinsic spatial resolution in global X and Y
         */
        virtual TMatrixD getSpatialResolutionMatrixGlobal(double column = 0, double row = 0) const = 0;

        /**
         * @brief Get number of pixels in x and y
         * @return Number of two dimensional pixels
         * @todo: this is designed for PixelDetector, find a proper interface for other Detector type
         */
        virtual ROOT::Math::DisplacementVector2D<ROOT::Math::Cartesian2D<int>> nPixels() const = 0;

        /**
         * @brief Get detector time offset from global clock, can be used to correct for constant shifts or time of flight
         * @return Time offset of respective detector
         */
        double timeOffset() const { return m_timeOffset; }

        /**
         * @brief Set the time offset of the detector
         * @param time New time offset of the detector
         */
        void setTimeOffset(double time) { m_timeOffset = time; }

        /**
         * @brief Get detector time resolution, used for timing cuts during clustering, track formation, etc.
         * @return Time resolutiom of respective detector
         */
        double getTimeResolution() const;

        /**
         * @brief Get position in the world
         * @return Global position in Cartesian coordinates
         */
        XYZPoint displacement() const { return alignment_->displacement(); }

        /**
         * @brief Get orientation in the world
         * @return Vector with three rotation angles
         */
        XYZVector rotation() const { return alignment_->orientation(); };

        /**
         * @brief Get normal vector to sensor surface
         * @return Normal vector to sensor surface
         */
        ROOT::Math::XYZVector normal() const { return alignment_->normal(); }
        /**
         * @brief Get orientation mode
         * @return string with mode
         */
        std::string orientation_mode() const { return alignment_->mode(); }

        /**
         * @brief Get origin vector to sensor surface
         * @return Origin vector to sensor surface
         */
        ROOT::Math::XYZPoint origin() const { return alignment_->origin(); }

        /**
         * @brief Get path of the file with calibration information
         * @return Path of the calibration file.
         *
         * @note The data contained in the calibration file is detector-specific and is
         * not parsed. This is left to the individual modules decoding the detector data.
         */
        std::filesystem::path calibrationFile() const { return m_calibrationfile.value_or(""); }

        /**
         * @brief Set the file with pixel mask information
         * @param file New mask file name
         */
        void maskFile(std::filesystem::path file);

        /**
         * @brief Get path of the file with pixel mask information
         * @return Path of the pixel mask file
         */
        std::filesystem::path maskFile() const { return m_maskfile; }

        /**
         * @brief Mark a detector channel as masked
         * @param chX X coordinate of the pixel to be masked
         * @param chY Y coordinate of the pixel to be masked
         * @todo: This is designed for PixelDetector, the parameters can be different with other type of Detector
         */
        virtual void maskChannel(int chX, int chY) = 0;

        /**
         * @brief Check if a detector channel is masked
         * @param chX X coordinate of the pixel to check
         * @param chY Y coordinate of the pixel to check
         * @return    Mask status of the pixel in question
         * @todo: This is designed for PixelDetector, the parameters can be different with other type of Detector
         */
        virtual bool masked(int chX, int chY) const = 0;

        /**
         * @brief Update coordinate transformations based on position and orientation values at the given time
         * @param time Time for which the alignment should be calculated
         */
        void update(double time);

        /**
         * @brief Update coordinate transformations based on the provided displacement and orientation
         * @param displacement  New displacement of the detector
         * @param orientation   New orientation of the detector
         */
        void update(const ROOT::Math::XYZPoint& displacement, const ROOT::Math::XYZVector& orientation);

        /**
         * @brief Check whether this detector has a variable alignment configured
         * @return true if alignment changes over time, false otherwise
         */
        bool hasVariableAlignment() const { return alignment_->isVariable(); }

        // Function to get global intercept with a track
        virtual PositionVector3D<Cartesian3D<double>> getIntercept(const Track* track) const = 0;
        // Function to get local intercept with a track
        virtual PositionVector3D<Cartesian3D<double>> getLocalIntercept(const Track* track) const = 0;

        // Function to check if a track intercepts with a plane
        virtual bool hasIntercept(const Track* track, double pixelTolerance = 0.) const = 0;

        // Function to check if a track goes through/near a masked pixel
        virtual bool hitMasked(const Track* track, int tolerance = 0.) const = 0;

        // Functions to get row and column from local position
        virtual double getRow(PositionVector3D<Cartesian3D<double>> localPosition) const = 0;
        virtual double getColumn(PositionVector3D<Cartesian3D<double>> localPosition) const = 0;

        // Function to get local position from column (x) and row (y) coordinates
        virtual PositionVector3D<Cartesian3D<double>> getLocalPosition(double column, double row) const = 0;

        // Function to get row and column of pixel
        virtual std::pair<int, int> getInterceptPixel(PositionVector3D<Cartesian3D<double>> localPosition) const = 0;

        /**
         * Transformation from local (sensor) coordinates to in-pixel coordinates
         * @param  column Column address ranging from int_column-0.5*pitch to int_column+0.5*pitch
         * @param  row Row address ranging from int_column-0.5*pitch to int_column+0.5*pitch
         * @return               Position within a single pixel cell, given in units of length
         */
        virtual XYVector inPixel(const double column, const double row) const = 0;

        /**
         * Transformation from local (sensor) coordinates to in-pixel coordinates
         * @param  localPosition Local position on the sensor
         * @return               Position within a single pixel cell, given in units of length
         */
        virtual XYVector inPixel(PositionVector3D<Cartesian3D<double>> localPosition) const = 0;

        /**
         * @brief Transform local coordinates of this detector into global coordinates
         * @param  local Local coordinates in the reference frame of this detector
         * @return       Global coordinates
         */
        XYZPoint localToGlobal(const XYZPoint& local) const { return alignment_->local2global() * local; };

        /**
         * @brief Transform global coordinates into detector-local coordinates
         * @param  global Global coordinates
         * @return        Local coordinates in the reference frame of this detector
         */
        XYZPoint globalToLocal(const XYZPoint& global) const { return alignment_->global2local() * global; };

        /**
         * @brief Check whether given track is within the detector's region-of-interest
         * @param  track The track to be checked
         * @return       Boolean indicating cluster affiliation with region-of-interest
         */
        virtual bool isWithinROI(const Track* track) const = 0;

        /**
         * @brief Check whether given cluster is within the detector's region-of-interest
         * @param  cluster The cluster to be checked
         * @return         Boolean indicating cluster affiliation with region-of-interest
         */
        virtual bool isWithinROI(Cluster* cluster) const = 0;

        /**
         * @brief Return the thickness of the senosr assembly layer (sensor+support) in fractions of radiation length
         * @return thickness in fractions of radiation length
         */
        double materialBudget() const { return m_materialBudget; }

        /**
         * @brief toGlobal Get the transformation from local to global coordinates
         * @return Transform3D local to global
         */
        Transform3D toGlobal() const { return alignment_->local2global(); }

        /**
         * @brief toLocal Get the transformation from global to local coordinates
         * @return Transform3D global to local
         */
        Transform3D toLocal() const { return alignment_->global2local(); }

        /**
         * @brief Test whether one pixel touches the cluster
         * @return true if it fulfills the condition
         * @note users should define their specific clustering method in the detector class
         */
        virtual bool
        isNeighbor(const std::shared_ptr<Pixel>&, const std::shared_ptr<Cluster>&, const int, const int) const = 0;

        /**
         * @brief Return a set containing all pixels neighboring the given one with a configurable maximum distance
         * @param px        Pixel in question
         * @param distance  Distance for pixels to be considered neighbors
         * @param include_corners Boolean to select whether pixels only touching via corners should be returned
         * @return Set of neighboring pixel indices, including the initial pixel
         *
         * @note The returned set should always also include the initial pixel indices the neighbors are calculated for
         *
         * @note alias for getNeighbors(const int col, const int row, const size_t distance, const bool include_corners)
         */
        std::set<std::pair<int, int>>
        getNeighbors(const std::shared_ptr<Pixel>& px, const size_t distance, const bool include_corners) {
            return getNeighbors(px->column(), px->row(), distance, include_corners);
        }

        /**
         * @brief Return a set containing all pixels neighboring the given one with a configurable maximum distance
         * @param col       Column of pixel in question
         * @param row       Row of pixel in question
         * @param distance  Distance for pixels to be considered neighbors
         * @param include_corners Boolean to select whether pixels only touching via corners should be returned
         * @return Set of neighboring pixel indices, including the initial pixel
         *
         * @note The returned set should always also include the initial pixel indices the neighbors are calculated for
         *
         * @note This method is purely virtual and must be implemented by the respective concrete detector model classes
         */
        virtual std::set<std::pair<int, int>>
        getNeighbors(const int col, const int row, const size_t distance, const bool include_corners) const = 0;

        /**
         * @brief Helper method to determine if this detector is of a given type
         * The template parameter needs to be specified specifically, i.e.
         *     if(model->is<PolarDetector>()) { }
         * @return Boolean indication whether this detector is of the given type or not
         */
        template <class T> bool is() { return dynamic_cast<T*>(this) != nullptr; }

    protected:
        // Roles of the detector
        DetectorRole m_role;

        // Build axis, for devices which are not auxiliary
        // Different in Pixel/Strip Detector
        virtual void build_axes(const Configuration& config) = 0;

        // Config detector, for devices which are not auxiliary
        // Different in Pixel/Strip Detector
        virtual void configure_detector(Configuration& config) const = 0;
        // Set position, orientation, mode of detector
        // Different in Pixel/Strip Detector
        virtual void configure_pos_and_orientation(Configuration& config) const = 0;

        virtual void process_mask_file() = 0;

        // Detector information
        std::string m_detectorType;
        std::string m_detectorName;
        std::string m_detectorCoordinates;

        double m_timeOffset;
        double m_timeResolution;
        double m_materialBudget;

        // Alignment and coordinate transformation information:
        std::shared_ptr<Alignment> alignment_;

        // Path of calibration file
        std::optional<std::filesystem::path> m_calibrationfile;

        // List of masked channels
        std::map<int, bool> m_masked;
        std::filesystem::path m_maskfile;
    };
} // namespace corryvreckan

#include "PixelDetector.hpp"
// #include "HexagonalPixelDetector.hpp"
#include "PolarDetector.hpp"
#endif // CORRYVRECKAN_DETECTOR_H
