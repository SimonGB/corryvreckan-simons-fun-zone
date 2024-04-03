# SPDX-FileCopyrightText: 2017-2022 CERN and the Corryvreckan authors
# SPDX-License-Identifier: MIT

# Include the standard dependencies
INCLUDE_DIRECTORIES(SYSTEM ${CORRYVRECKAN_DEPS_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(SYSTEM ${GBL_INCLUDE_DIRS})

# Include the RootMacros.cmake file
INCLUDE("${ROOT_USE_FILE}")

# Generate the ROOT dictionary
ROOT_GENERATE_DICTIONARY(CorryvreckanObjectsDictionary
    ${CMAKE_SOURCE_DIR}/3rdparty/GeneralBrokenLines/include/GblPoint.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Object.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Pixel.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Cluster.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SpidrSignal.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Track.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/StraightLineTrack.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GblTrack.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/MCParticle.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Event.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Multiplet.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Waveform.hpp
    LINKDEF
    ${CMAKE_CURRENT_SOURCE_DIR}/Linkdef.h
    OPTIONS
    -inlineInputHeader
    -I${CMAKE_CURRENT_SOURCE_DIR}
    MODULE
    CorryvreckanObjects
)

SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/CorryvreckanObjectsDictionary.cxx
  PROPERTIES
  GENERATED TRUE
  COMPILE_FLAGS "-Wno-unused-function -Wno-overlength-strings -Wno-zero-as-null-pointer-constant -w"
)

# Define the library adding the object file created above
ADD_LIBRARY(CorryvreckanObjects SHARED
    Object.cpp
    Pixel.cpp
    Cluster.cpp
    Track.cpp
    StraightLineTrack.cpp
    GblTrack.cpp
    MCParticle.cpp
    Event.cpp
    Multiplet.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/CorryvreckanObjectsDictionary.cxx
)

# Link the standard dependencies
TARGET_INCLUDE_DIRECTORIES(CorryvreckanObjects INTERFACE ${GBL_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(CorryvreckanObjects ${CORRYVRECKAN_LIBRARIES} ${CORRYVRECKAN_DEPS_LIBRARIES})
TARGET_COMPILE_OPTIONS(CorryvreckanObjects PRIVATE ${CORRYVRECKAN_CXX_FLAGS})

# Specify install for the messages
INSTALL(TARGETS CorryvreckanObjects
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib)

# Also install the dictionary objects
INSTALL(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/libCorryvreckanObjects_rdict.pcm
    ${CMAKE_CURRENT_BINARY_DIR}/libCorryvreckanObjects.rootmap
    DESTINATION lib)
