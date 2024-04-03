/**
 * @file
 * @brief Interface to the core framework
 *
 * @copyright Copyright (c) 2017-2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

/**
 * @defgroup Managers Managers
 * @brief The global set of managers used in framework
 */

#ifndef CORRYVRECKAN_CORRYVRECKAN_H
#define CORRYVRECKAN_CORRYVRECKAN_H

#include <atomic>
#include <fstream>
#include <memory>
#include <string>

#include "config/ConfigManager.hpp"
#include "module/ModuleManager.hpp"

namespace corryvreckan {
    /**
     * @brief Provides the link between the core framework and the executable.
     *
     * Supply the path location the main configuration which should be provided to the executable. Hereafter this class
     * should be used to load, initialize, run and finalize all the modules.
     */

    class Corryvreckan {
    public:
        /**
         * @brief Constructs Corryvreckan and initialize all managers
         * @param config_file_name Path of the main configuration file
         * @param module_options List of extra configuration options for module instances
         * @param detector_options List of extra configuration options for detector instances
         */
        explicit Corryvreckan(std::string config_file_name,
                              const std::vector<std::string>& module_options = std::vector<std::string>(),
                              const std::vector<std::string>& detector_options = std::vector<std::string>());

        /**
         * @brief Load modules from the main configuration and construct them
         * @warning Should be called after the \ref Corryvreckan() "constructor"
         */
        void load();

        /**
         * @brief Initialize all modules (pre-run)
         * @warning Should be called after the \ref Corryvreckan::load "load function"
         */
        void init();

        /**
         * @brief Run all modules for the number of events (run)
         * @warning Should be called after the \ref Corryvreckan::init "init function"
         */
        void run();

        /**
         * @brief Finalize all modules (post-run)
         * @warning Should be called after the \ref Corryvreckan::run "run function"
         */
        void finalize();

        /**
         * @brief Request termination as early as possible without changing the standard flow
         */
        void terminate();

    private:
        /**
         * @brief Sets the default unit conventions
         */
        void add_units();

        /**
         * @brief Set the default ROOT plot style
         */
        void set_style();

        // Indicate the framework should terminate
        std::atomic<bool> terminate_;
        std::atomic<bool> has_run_;

        // Log file if specified
        std::ofstream log_file_;

        // All managers in the framework
        std::unique_ptr<ModuleManager> mod_mgr_;
        std::unique_ptr<ConfigManager> conf_mgr_;
    };
} // namespace corryvreckan

#endif /* CORRYVRECKAN_CORRYVRECKAN_H */
