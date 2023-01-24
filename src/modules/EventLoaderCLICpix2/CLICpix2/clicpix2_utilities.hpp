/**
 * @file
 * @brief CLICpix2 frame decorder
 * Based on System Verilog CLICpix2_Readout_scoreboard
 *
 * @copyright Copyright (c) 2022 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#ifndef CLICPIX2_UTILITIES_HPP
#define CLICPIX2_UTILITIES_HPP

#include <array>

#include "clicpix2_pixels.hpp"

namespace clicpix2_utils {

    /* Routine to read the pixel matrix configuration from file and store it
     */
    std::map<std::pair<uint8_t, uint8_t>, caribou::pixelConfig> readMatrix(std::string filename);
} // namespace clicpix2_utils
#endif
