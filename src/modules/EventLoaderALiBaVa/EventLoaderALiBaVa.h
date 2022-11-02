/**
 * @file
 * @brief Definition of module EventLoaderALiBaVa
 *
 * @copyright Copyright (c) 2021-2022 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#ifndef EventLoaderALiBaVa_H
#define EventLoaderALiBaVa_H 1

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
        std::shared_ptr<Detector> detector_;
        DataFileRoot* ALiBaVaPointer;

        TH1F* hChargeSignal{};
        TH1F* hADCSignal{};
        TH1F* hSNR{};
        TH1F* hPedestal{};
        TH1F* hNoise{};
        TH1F* hPedestalCorrect{};
        TH1F* hNoiseCorrect{};
        TProfile* hTimeProfile{};

        int m_run{};
        double m_timecut_low{};
        double m_timecut_up{};
        int m_ignore_events{};
        int m_lower_channel{};
        int m_upper_channel{};
        double m_chargecut{};
        std::string m_inputDirectory{};
        bool m_correct_crosstalk{};
        double m_calibration_constant{};
        double m_b_one{};
        double m_b_two{};
        std::vector<unsigned int> m_roi{};
        std::vector<unsigned int> m_roi_ch;
        bool m_horizontal;
        int m_polarity;

        std::string m_datafilename;
        std::string m_pedestalfilename;
        std::string m_calibrationfilename;
    };

} // namespace corryvreckan
#endif
