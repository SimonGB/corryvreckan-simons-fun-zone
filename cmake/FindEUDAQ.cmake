# - Try to find pxarCore used in native reader library for exporting RAW to pxarCore format needed for reconstruction
# Once done this will define
#  EUDAQ_FOUND - System has EUDAQ
#  EUDAQ_INCLUDE_DIR - The EUDAQ main include directories
#  EUDAQ_LIBRARY - The libraries needed to use EUDAQ

MESSAGE(STATUS "Looking for EUDAQ...")

FIND_PATH(EUDAQ_INCLUDE_DIR NAMES "eudaq/FileReader.hh" PATHS "$ENV{EUDAQPATH}/main/include")
MESSAGE(STATUS "FileReader.hh => ${EUDAQ_INCLUDE_DIR}")
IF(EUDAQ_INCLUDE_DIR)
   SET(EUDAQ_INC_FOUND TRUE)
   MESSAGE(STATUS "Found EUDAQ headers: ${EUDAQ_INCLUDE_DIR}")
ENDIF()

FIND_LIBRARY(EUDAQ_LIBRARY NAMES "EUDAQ" HINTS "$ENV{EUDAQPATH}/lib")
MESSAGE(STATUS "libEUDAQ => ${EUDAQ_LIBRARY}")
IF(EUDAQ_LIBRARY)
   SET(EUDAQ_LIB_FOUND TRUE)
   MESSAGE(STATUS "Found EUDAQ library: ${EUDAQ_LIBRARY}")
ENDIF()

IF(EUDAQ_LIB_FOUND AND EUDAQ_INC_FOUND)
   SET(EUDAQ_FOUND TRUE)
ENDIF()

mark_as_advanced(EUDAQ_LIBRARY EUDAQ_INCLUDE_DIR)
