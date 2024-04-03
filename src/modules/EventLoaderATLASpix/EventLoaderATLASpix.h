/**
 * @file
 * @brief Definition of module EventLoaderATLASpix
 *
 * @copyright Copyright (c) 2017-2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#ifndef EventLoaderATLASpix_H
#define EventLoaderATLASpix_H 1

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile2D.h>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <queue>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include "core/module/Module.hpp"
#include "objects/Cluster.hpp"
#include "objects/Pixel.hpp"
#include "objects/Track.hpp"

namespace corryvreckan {
    /** @ingroup Modules
     */
    class EventLoaderATLASpix : public Module {

    public:
        // Constructors and destructors
        EventLoaderATLASpix(Configuration& config, std::shared_ptr<Detector> detector);
        ~EventLoaderATLASpix() {}

        // Functions
        void initialize() override;
        StatusCode run(const std::shared_ptr<Clipboard>& clipboard) override;

    private:
        /*
         * @brief Converts gray encoded data to binary number
         */
        uint32_t gray_decode(uint32_t gray);

        /*
         * @brief Read data in the format written by the Caribou readout system and fill time-sorted buffer
         * @return Bool which is false when reaching the end-of-file and true otherwise
         */
        bool read_caribou_data();

        // custom comparator for time-sorted priority_queue
        struct CompareTimeGreater {
            bool operator()(const std::shared_ptr<Pixel> a, const std::shared_ptr<Pixel> b) {
                return a->timestamp() > b->timestamp();
            }
        };
        // Buffer of timesorted pixel hits: (need to use greater here!)
        std::priority_queue<std::shared_ptr<Pixel>, PixelVector, CompareTimeGreater> sorted_pixels_;

        std::shared_ptr<Detector> m_detector;
        std::string m_filename;
        std::ifstream m_file;

        // Resuming in next event:
        unsigned long long readout_ts_ = 0;
        unsigned long long fpga_ts_ = 0;
        unsigned long long fpga_ts1_ = 0;
        unsigned long long fpga_ts2_ = 0;
        unsigned long long fpga_ts3_ = 0;
        bool new_ts1_ = false;
        bool new_ts2_ = false;
        size_t t0_seen_ = 0;
        bool eof_reached = false;

        // int ts1Range;
        int ts2Range;

        TH1F* hMessages;
        TH2F* hHitMap;
        TH2F* hHitMap_highTot;
        TProfile2D* hHitMap_totWeighted;
        TH1F* hPixelToT;
        TH1F* hPixelToT_beforeCorrection;
        TH1F* hPixelCharge;
        TH1F* hPixelToA;
        TH1F* hPixelTimeEventBeginResidual;
        TH1F* hPixelTimeEventBeginResidual_wide;
        TH2F* hPixelTimeEventBeginResidualOverTime;

        std::map<size_t, TH1D*> hPixelTriggerTimeResidual;
        TH2D* hPixelTriggerTimeResidualOverTime;
        TH1D* hTriggersPerEvent;

        TH1F* hPixelMultiplicity;
        TH1F* hPixelTimes;
        TH1F* hPixelTimes_long;

        // TS1 and TS2 specific histograms:
        TH1F* hPixelTS1;
        TH1F* hPixelTS2;
        TH1F* hPixelTS1bits;
        TH1F* hPixelTS2bits;
        TH1F* hPixelTS1_lowToT;
        TH1F* hPixelTS2_lowToT;
        TH1F* hPixelTS1bits_lowToT;
        TH1F* hPixelTS2bits_lowToT;
        TH1F* hPixelTS1_highToT;
        TH1F* hPixelTS2_highToT;
        TH1F* hPixelTS1bits_highToT;
        TH1F* hPixelTS2bits_highToT;

        // Parameters:
        std::string m_inputDirectory;
        double m_clockCycle;
        int m_highToTCut;
        // int m_clkdivendM;
        int m_clkdivend2M;
        int m_buffer_depth;
        double m_time_offset;
    };
} // namespace corryvreckan
#endif // EventLoaderATLASpix_H
