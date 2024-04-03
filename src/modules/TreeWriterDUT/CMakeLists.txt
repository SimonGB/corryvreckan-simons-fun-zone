# SPDX-FileCopyrightText: 2017-2022 CERN and the Corryvreckan authors
# SPDX-License-Identifier: MIT

# Define module and return the generated name as MODULE_NAME
CORRYVRECKAN_DUT_MODULE(MODULE_NAME)

# Add source files to library
CORRYVRECKAN_MODULE_SOURCES(${MODULE_NAME}
    TreeWriterDUT.cpp
)

TARGET_LINK_LIBRARIES(${MODULE_NAME} ROOT::Tree)

# Provide standard install target
CORRYVRECKAN_MODULE_INSTALL(${MODULE_NAME})
