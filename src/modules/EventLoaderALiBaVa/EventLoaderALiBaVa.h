/**
 * @file
 * @brief Definition of module EventLoaderALiBaVa
 *
 * @copyright Copyright (c) 2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

 #ifndef EventLoaderALiBaVa_H
 #define EventLoaderALiBaVa_H 1

#include "core/module/Module.hpp"
#include "ALiBaVa/DataFileRoot.h"

namespace corryvreckan {
    /** @ingroup Modules
     * @brief Module to do function
     *
     * More detailed explanation of module
     */
    class EventLoaderALiBaVa : public Module {

    public:
        /**
         * @brief Constructor for this unique module
         * @param config Configuration object for this module as retrieved from the steering file
         * @param detector Pointer to the detector for this module instance
         */
        EventLoaderALiBaVa(Configuration& config, std::shared_ptr<Detector> detector);

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
        std::shared_ptr<Detector> m_detector;
        DataFileRoot * ALiBaVaPointer;
        int iEvent = 0;
        int nEvents = 0;

        std::string m_datafilename;
        std::string m_pedestalfilename;
        std::string m_calibrationfilename;
    };

} // namespace corryvreckan
#endif
