/**
 * @file
 * @brief Definition of EventLoaderMuPixTelescope module
 *
 * @copyright Copyright (c) 2019-2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>
#include "core/module/Module.hpp"

#include "blockfile.hpp"
#include "telescope_frame.hpp"

namespace corryvreckan {
    /** @ingroup Modules
     * @brief Module to do function
     *
     * More detailed explanation of module
     */
    class EventLoaderMuPixTelescope : public Module {

    public:
        /**
         * @brief Constructor for this unique module
         * @param config Configuration object for this module as retrieved from the steering file
         * @param detectors Vector of pointers to the detectors
         */
        EventLoaderMuPixTelescope(Configuration config, std::vector<std::shared_ptr<Detector>> detectors);

        /**
         * @brief [Initialise this module]
         */
        void initialise();

        /**
         * @brief [Run the function of this module]
         */
        StatusCode run(std::shared_ptr<Clipboard> clipboard);

    private:
        std::string m_inputDirectory;
        bool m_isSorted;
        bool m_ts2IsGray;
        int m_runNumber;
        BlockFile* m_blockFile;
        TelescopeFrame m_tf;

        // Histograms
        TH2F* hHitMap;
        TH1F* hPixelToT;
        TH1F* hTimeStamp;
    };

} // namespace corryvreckan
