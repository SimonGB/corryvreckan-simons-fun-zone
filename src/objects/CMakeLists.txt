# Include the standard dependencies
INCLUDE_DIRECTORIES(SYSTEM ${CORRYVRECKAN_DEPS_INCLUDE_DIRS})

# Find the RootNewMacros.cmake file and include it
GET_FILENAME_COMPONENT(ROOT_CMAKE_DIR ${ROOT_USE_FILE} DIRECTORY)
FIND_FILE(ROOT_MACROS_FILE
          NAMES RootNewMacros.cmake
          HINTS ${ROOT_CMAKE_DIR}/modules/ $ENV{ROOTSYS}/cmake/modules)

IF(NOT ROOT_MACROS_FILE)
    MESSAGE(WARNING "Cannot find ROOT macros, including generic file as fallback.")
    INCLUDE(${ROOT_USE_FILE})
ELSE()
    INCLUDE(${ROOT_MACROS_FILE})
ENDIF()

# Generate the ROOT dictionary
ROOT_GENERATE_DICTIONARY(CorryvreckanObjectsDictionary
    ${CMAKE_CURRENT_SOURCE_DIR}/Cluster.h
    ${CMAKE_CURRENT_SOURCE_DIR}/GuiDisplay.h
    ${CMAKE_CURRENT_SOURCE_DIR}/KDTree.h
    ${CMAKE_CURRENT_SOURCE_DIR}/KDTreeTimepix3.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Pixel.h
    ${CMAKE_CURRENT_SOURCE_DIR}/SpidrSignal.h
    ${CMAKE_CURRENT_SOURCE_DIR}/TestBeamObject.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Timepix3Cluster.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Timepix3Pixel.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Timepix3Track.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Track.h
    MODULE
    CorryvreckanObjects
)

# Explicitly add all the dependent include dirs to the manual command
FOREACH(dir IN ITEMS ${CORRYVRECKAN_DEPS_INCLUDE_DIRS})
    IF(NOT INCLUDE_DIR_ARGS)
        SET(INCLUDE_DIR_ARGS "-I${dir}")
    ELSE()
        SET(INCLUDE_DIR_ARGS "${INCLUDE_DIR_ARGS} -I${dir}")
    ENDIF()
ENDFOREACH()

# Compile the dictionary through a special target
# WARNING This hack is necessary to prevent standard warnings and those in clang-tidy for use in the CI
ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CorryvreckanObjectsDictionary.cxx.o COMMAND ${CMAKE_CXX_COMPILER} -fPIC -std=c++${CMAKE_CXX_STANDARD} -I${CMAKE_SOURCE_DIR}/src ${INCLUDE_DIR_ARGS} -o ${CMAKE_CURRENT_BINARY_DIR}/CorryvreckanObjectsDictionary.cxx.o -c ${CMAKE_CURRENT_BINARY_DIR}/CorryvreckanObjectsDictionary.cxx DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CorryvreckanObjectsDictionary.cxx)

# Define the library adding the object file created above
ADD_LIBRARY(CorryvreckanObjects SHARED
    TestBeamObject.C
    ${CMAKE_CURRENT_BINARY_DIR}/CorryvreckanObjectsDictionary.cxx.o
)

# Link the standard dependencies
TARGET_LINK_LIBRARIES(CorryvreckanObjects ${CORRYVRECKAN_DEPS_LIBRARIES})

# Specify install for the messages
INSTALL(TARGETS CorryvreckanObjects
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib)

# Also install the dictionary objects
INSTALL(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/libCorryvreckanObjects_rdict.pcm
    ${CMAKE_CURRENT_BINARY_DIR}/libCorryvreckanObjects.rootmap
    DESTINATION lib)