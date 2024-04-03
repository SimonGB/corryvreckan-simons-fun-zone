/**
 * @file
 * @brief Definition of cluster object
 *
 * @copyright Copyright (c) 2017-2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#ifndef CORRYVRECKAN_CLUSTER_H
#define CORRYVRECKAN_CLUSTER_H 1

#include <Math/Point3D.h>
#include <Math/Vector2D.h>
#include <TMatrixD.h>
#include <TRef.h>

#include <iostream>

#include "Object.hpp"
#include "Pixel.hpp"

namespace corryvreckan {
    /**
     * @ingroup Objects
     * @brief This class is a simple cluster class which is used as a base class to interface with the track class. Anything
     * which inherits from it can be placed on a track and used for fitting.
     */
    class Cluster : public Object {

    public:
        // Constructors and destructors
        Cluster() = default;

        /**
         * @brief Static member function to obtain base class for storage on the clipboard.
         * This method is used to store objects from derived classes under the typeid of their base classes
         *
         * @warning This function should not be implemented for derived object classes
         *
         * @return Class type of the base object
         */
        static std::type_index getBaseType() { return typeid(Cluster); }

        // Functions
        // Add a new pixel to the cluster
        void addPixel(const Pixel* pixel);

        // Retrieve cluster parameters
        double column() const { return m_column; }
        double row() const { return m_row; }
        double charge() const { return m_charge; }
        double error() const;
        double errorX() const { return m_error.X(); }
        double errorY() const { return m_error.Y(); }
        TMatrixD errorMatrixGlobal() const { return m_error_matrix_global; }

        bool isSplit() const { return m_split; }
        void setSplit(bool split);

        ROOT::Math::XYZPoint global() const { return m_global; }
        ROOT::Math::XYZPoint local() const { return m_local; }

        size_t size() const { return pixels_.size(); }
        size_t columnWidth() const { return m_columnWidth; }
        size_t rowWidth() const { return m_rowWidth; }
        std::vector<const Pixel*> pixels() const;

        /**
         * @brief Retrieve the seed pixel of the cluster.
         *
         * The seed pixel is defined as the one with the highest charge. In case all pixels have the same charge return a
         * random pixel.
         *
         * @return Seed pixel of the cluster
         */
        const Pixel* getSeedPixel() const;
        /**
         * @brief Retrieve the earliest pixel of the cluster.  In case all pixels have the same timstamp return a random
         * pixel.
         *
         * @return Earliest pixel of the cluster
         */
        const Pixel* getEarliestPixel() const;

        // Set cluster parameters
        void setColumn(double col) { m_column = col; }
        void setRow(double row) { m_row = row; }
        void setCharge(double charge) { m_charge = charge; }
        void setClusterCentre(ROOT::Math::XYZPoint global) { m_global = std::move(global); }
        void setClusterCentreLocal(ROOT::Math::XYZPoint local) { m_local = std::move(local); }
        void setErrorX(double error) { m_error.SetX(error); }
        void setErrorY(double error) { m_error.SetY(error); }
        void setError(ROOT::Math::XYVector error) { m_error = std::move(error); }
        void setErrorMatrixGlobal(TMatrixD errorMatrix) { m_error_matrix_global = std::move(errorMatrix); }

        /**
         * @brief Print an ASCII representation of Cluster to the given stream
         * @param out Stream to print to
         */
        void print(std::ostream& out) const override;

        void loadHistory() override;
        void petrifyHistory() override;

    private:
        // Member variables
        std::vector<PointerWrapper<Pixel>> pixels_;
        double m_column;
        double m_row;
        double m_charge;
        ROOT::Math::XYVector m_error;
        TMatrixD m_error_matrix_global{3, 3};
        size_t m_columnWidth{0};
        size_t m_rowWidth{0};
        bool m_split{false};

        ROOT::Math::XYZPoint m_local;
        ROOT::Math::XYZPoint m_global;

        std::map<int, bool> m_rowHits;
        std::map<int, bool> m_columnHits;

        // ROOT I/O class definition - update version number when you change this class!
        ClassDefOverride(Cluster, 15)
    };

    // Vector type declaration
    using ClusterVector = std::vector<std::shared_ptr<Cluster>>;
} // namespace corryvreckan

#endif // CORRYVRECKAN_CLUSTER_H
