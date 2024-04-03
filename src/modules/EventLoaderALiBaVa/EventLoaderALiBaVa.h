/**
 * @file
 * @brief Definition of module EventLoaderALiBaVa
 *
 * @copyright Copyright (c) 2021-2022 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#ifndef EventLoaderALiBaVa_H
#define EventLoaderALiBaVa_H 1

#include <TFormula.h>
#include "ALiBaVa/DataFileRoot.h"
#include "core/module/Module.hpp"

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
        std::shared_ptr<TFormula> parse_formula();

        std::shared_ptr<Detector> detector_;
        std::shared_ptr<DataFileRoot> alibava_;

        TH1F* hChargeSignal{};
        TH1F* hADCSignal{};
        TH1F* hSNR{};
        TH1F* hPedestal{};
        TH1F* hNoise{};
        TH1F* hPedestalCorrect{};
        TH1F* hNoiseCorrect{};
        TH1F* hNoiseElectrons{};
        TProfile* hTimeProfile{};
        TProfile* hTemperature{};
        TH2F* hPedestalCorrect2D{};
        TH2F* hNoiseCorrect2D{};
        TH2F* hNoiseElectrons2D{};

        double chargecut_{};
        std::shared_ptr<TFormula> charge_calibration_formula_;
        std::vector<unsigned int> roi_ch_{};
    };

} // namespace corryvreckan
#endif
