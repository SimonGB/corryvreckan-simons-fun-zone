# SPDX-FileCopyrightText: 2017-2022 CERN and the Corryvreckan authors
# SPDX-License-Identifier: MIT

# Define module and return the generated name as MODULE_NAME
CORRYVRECKAN_GLOBAL_MODULE(MODULE_NAME)

# Add source files to library
CORRYVRECKAN_MODULE_SOURCES(${MODULE_NAME}
    AnalysisTrackIntercept.cpp
)

# Provide standard install target
CORRYVRECKAN_MODULE_INSTALL(${MODULE_NAME})
