/**
 * @file
 * @brief Definition of module AlignmentTime
 *
 * @copyright Copyright (c) 2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>
#include "core/module/Module.hpp"
#include "objects/Cluster.hpp"
#include "objects/Pixel.hpp"
#include "objects/Track.hpp"

namespace corryvreckan {
    /** @ingroup Modules
     * @brief Module to do function
     *
     * More detailed explanation of module
     */
    class AlignmentTime : public Module {

    public:
        /**
         * @brief Constructor for this unique module
         * @param config Configuration object for this module as retrieved from the steering file
         * @param detector Pointer to the detector for this module instance
         */
        AlignmentTime(Configuration& config, std::shared_ptr<Detector> detector);

        /**
         * @brief [Initialise this module]
         */
        void initialize() override;

        /**
         * @brief [Run the function of this module]
         */
        StatusCode run(const std::shared_ptr<Clipboard>& clipboard) override;

        /**
         * @brief [Finalise module]
         */
        void finalize(const std::shared_ptr<ReadonlyClipboard>& clipboard) override;

    private:
        // Handling time reference
        std::string reference_name_;
        bool reference_filled_;

        // Container for time stamps
        std::map<std::string, std::vector<double>> timestamps_;

        // Returns the array element closest to the target value.
        double findClosest(std::vector<double> const&, double);
        // Helper for findClosest.
        // Compares two values to target and returns the closer one.
        double whichCloser(double, double, double);

        // Histograms
        std::map<std::string, TH1D*> hTimeStamps;
        std::map<std::string, TH1D*> hTimeStamps_long;
        TH1D* hTimeStampsRef;
        TH1D* hTimeStampsRef_long;

        int m_eventNumber;
    };

} // namespace corryvreckan
