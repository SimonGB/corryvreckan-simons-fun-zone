# Define module and return the generated name as MODULE_NAME
CORRYVRECKAN_DETECTOR_MODULE(MODULE_NAME) # was originally DETECTOR not DUT
#CORRYVRECKAN_EXCLUDE_AUX(${MODULE_NAME}) # was originally not here

# Add source files to library
CORRYVRECKAN_MODULE_SOURCES(${MODULE_NAME}
    ClusteringSeed.cpp
    # ADD SOURCE FILES HERE...
)

# Provide standard install target
CORRYVRECKAN_MODULE_INSTALL(${MODULE_NAME})
