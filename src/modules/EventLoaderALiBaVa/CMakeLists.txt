# Define module and return the generated name as MODULE_NAME
CORRYVRECKAN_DETECTOR_MODULE(MODULE_NAME)

# Add source files to library
CORRYVRECKAN_MODULE_SOURCES(${MODULE_NAME}
    EventLoaderALiBaVa.cpp
    ALiBaVa/AsciiRoot.cpp
    ALiBaVa/auxfunctions.cpp
    ALiBaVa/ChanList.cpp
    ALiBaVa/DataFileRoot.cpp
    ALiBaVa/HDFRoot.cpp
    ALiBaVa/Tracer.cpp
    ALiBaVa/utils.cpp
    ALiBaVa/Hit.cpp
)

FIND_PACKAGE(HDF5 REQUIRED)
TARGET_INCLUDE_DIRECTORIES(${MODULE_NAME} SYSTEM PRIVATE ${HDF5_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(${MODULE_NAME} ${HDF5_LIBRARIES})
TARGET_COMPILE_DEFINITIONS(${MODULE_NAME} PRIVATE HAVE_HDF5=1)

# Provide standard install target
CORRYVRECKAN_MODULE_INSTALL(${MODULE_NAME})
