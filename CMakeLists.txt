# Minimum CMake version
cmake_minimum_required(VERSION 3.24)

# Project name and language
project(TAF LANGUAGES CXX)

# Enable testing
enable_testing()

# Define the library name
set(libname lib${PROJECT_NAME})
message("Library name: ${libname}")

# Set directory paths
set(DTDIR ${PROJECT_SOURCE_DIR})
set(BUILD_DIR "${DTDIR}/bin/bin")
set(OBJ_DIR "${DTDIR}/bin/obj")
set(LIB_DIR "${DTDIR}/bin/lib")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BUILD_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB_DIR})
set(CMAKE_INSTALL_LIBDIR ${LIB_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${OBJ_DIR})

message("Build directories:")
message("  BUILD_DIR = ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
message("  OBJ_DIR   = ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
message("  LIB_DIR   = ${CMAKE_INSTALL_LIBDIR}")

# C++ standard and compiler flags
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fPIC -O2 -g")

# Clang-specific settings
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()

# Define source and include directories
set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/code/include)
set(SRC_DIR ${PROJECT_SOURCE_DIR}/code/src)

# Define source files for DT and M components
set(DTSRCS
  ${SRC_DIR}/MMain.cxx
  ${SRC_DIR}/TNTBoardReader.cxx
  ${SRC_DIR}/PXIBoardReader.cxx
  ${SRC_DIR}/GIGBoardReader.cxx
  ${SRC_DIR}/IMGBoardReader.cxx
  ${SRC_DIR}/BoardReader.cxx
  ${SRC_DIR}/VMEBoardReader.cxx
  ${SRC_DIR}/MCBoardReader.cxx
  ${SRC_DIR}/BoardReaderIHEP.cxx
  ${SRC_DIR}/AliMIMOSA22RawStreamVASingle.cxx
  ${SRC_DIR}/DecoderM18.cxx
  ${SRC_DIR}/DecoderGeant.cxx
  ${SRC_DIR}/BoardReaderMIMOSIS.cxx
  ${SRC_DIR}/DSession.cxx
  ${SRC_DIR}/DSetup.cxx
  ${SRC_DIR}/DAcq.cxx
  ${SRC_DIR}/DTracker.cxx
  ${SRC_DIR}/DPlane.cxx
  ${SRC_DIR}/DStrip.cxx
  ${SRC_DIR}/DHit.cxx
  ${SRC_DIR}/DTrack.cxx
  ${SRC_DIR}/DLine.cxx
  ${SRC_DIR}/DR3.cxx
  ${SRC_DIR}/DCut.cxx
  ${SRC_DIR}/DAlign.cxx
  ${SRC_DIR}/DEvent.cxx
  ${SRC_DIR}/DEventMC.cxx
  ${SRC_DIR}/DParticle.cxx
  ${SRC_DIR}/DGlobalTools.cxx
  ${SRC_DIR}/DPrecAlign.cxx
  ${SRC_DIR}/DPixel.cxx
  ${SRC_DIR}/DLadder.cxx
  ${SRC_DIR}/DMiniVector.cxx
  ${SRC_DIR}/DHelix.cxx
  ${SRC_DIR}/DHelixFitter.cxx
  ${SRC_DIR}/DTrackFitter.cxx
  ${SRC_DIR}/DBeaster.cxx
  ${SRC_DIR}/MKalmanFilter.cxx
  ${SRC_DIR}/MLeastChiSquare.cxx
)

set(DTHDRS
  ${INCLUDE_DIR}/TNTBoardReader.h
  ${INCLUDE_DIR}/PXIBoardReader.h
  ${INCLUDE_DIR}/GIGBoardReader.h
  ${INCLUDE_DIR}/IMGBoardReader.h
  ${INCLUDE_DIR}/BoardReader.h
  ${INCLUDE_DIR}/VMEBoardReader.h
  ${INCLUDE_DIR}/MCBoardReader.h
  ${INCLUDE_DIR}/BoardReaderIHEP.h
  ${INCLUDE_DIR}/AliMIMOSA22RawStreamVASingle.h
  ${INCLUDE_DIR}/DecoderM18.h
  ${INCLUDE_DIR}/DecoderGeant.h
  ${INCLUDE_DIR}/BoardReaderMIMOSIS.h
  ${INCLUDE_DIR}/DSession.h
  ${INCLUDE_DIR}/DSetup.h
  ${INCLUDE_DIR}/DAcq.h
  ${INCLUDE_DIR}/DTracker.h
  ${INCLUDE_DIR}/DPlane.h
  ${INCLUDE_DIR}/DStrip.h
  ${INCLUDE_DIR}/DHit.h
  ${INCLUDE_DIR}/DTrack.h
  ${INCLUDE_DIR}/DLine.h
  ${INCLUDE_DIR}/DR3.h
  ${INCLUDE_DIR}/DCut.h
  ${INCLUDE_DIR}/DAlign.h
  ${INCLUDE_DIR}/DEvent.h
  ${INCLUDE_DIR}/DEventMC.h
  ${INCLUDE_DIR}/DParticle.h
  ${INCLUDE_DIR}/DGlobalTools.h
  ${INCLUDE_DIR}/DPrecAlign.h
  ${INCLUDE_DIR}/DPixel.h
  ${INCLUDE_DIR}/DLadder.h
  ${INCLUDE_DIR}/DMiniVector.h
  ${INCLUDE_DIR}/DHelix.h
  ${INCLUDE_DIR}/DHelixFitter.h
  ${INCLUDE_DIR}/DTrackFitter.h
  ${INCLUDE_DIR}/DBeaster.h
  ${INCLUDE_DIR}/MKalmanFilter.h
  ${INCLUDE_DIR}/MLeastChiSquare.h
)

set(MSRCS
  ${SRC_DIR}/MPrep.cxx
  ${SRC_DIR}/MAnalysis.cxx
  ${SRC_DIR}/MPost.cxx
  ${SRC_DIR}/MCommands.cxx
  ${SRC_DIR}/MMCGeneration.cxx
  ${SRC_DIR}/MAlign.cxx
  ${SRC_DIR}/MHist.cxx
  ${SRC_DIR}/MRaw.cxx
  ${SRC_DIR}/MRax.cxx
  ${SRC_DIR}/MAlignment.cxx
  ${SRC_DIR}/MMillepede.cxx
  ${SRC_DIR}/MGlobalAlign.cxx
  ${SRC_DIR}/MKalmanFilter.cxx
  ${SRC_DIR}/MLeastChiSquare.cxx
)

set(MHDRS
  ${INCLUDE_DIR}/MAnalysis.h
  ${INCLUDE_DIR}/MAlign.h
  ${INCLUDE_DIR}/MHist.h
  ${INCLUDE_DIR}/MRaw.h
  ${INCLUDE_DIR}/MRax.h
  ${INCLUDE_DIR}/MAlignment.h
  ${INCLUDE_DIR}/MMillepede.h
  ${INCLUDE_DIR}/MGlobalAlign.h
  ${INCLUDE_DIR}/MKalmanFilter.h
  ${INCLUDE_DIR}/MLeastChiSquare.h
)

# ROOT setup
find_package(ROOT COMPONENTS Minuit Geom EG Gui Ged TreePlayer Matrix Physics MathCore Rio REQUIRED)
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR})
link_libraries(${ROOT_LIBRARIES})

# Generate dictionaries
message("Generating DT dictionary")
ROOT_GENERATE_DICTIONARY(DTdict ${DTHDRS} LINKDEF ${INCLUDE_DIR}/DTLinkDef.h)

message("Generating M dictionary")
ROOT_GENERATE_DICTIONARY(Mdict ${MHDRS} LINKDEF ${INCLUDE_DIR}/MLinkDef.h)

# Create library
include_directories(${INCLUDE_DIR})
add_library(${libname} SHARED ${DTSRCS} ${MSRCS} DTdict.cxx Mdict.cxx)
target_link_libraries(${libname} PUBLIC ${ROOT_LIBRARIES})

# Create executable
add_executable(TAF ${SRC_DIR}/MMain.cxx)
target_link_libraries(TAF PUBLIC ${libname} ${ROOT_LIBRARIES})

# Output properties
set_target_properties(${PROJECT_NAME}
  PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${OBJ_DIR}
  LIBRARY_OUTPUT_DIRECTORY ${LIB_DIR}
  RUNTIME_OUTPUT_DIRECTORY ${BUILD_DIR}
)

# # Custom commands for dictionary files
# add_custom_command(TARGET ${libname} POST_BUILD
#   COMMAND ${CMAKE_COMMAND} -E copy
#   ${OBJ_DIR}/libTAF.so ${LIB_DIR}/libTAF.so
# )
# add_custom_command(TARGET ${libname} POST_BUILD
#   COMMAND ${CMAKE_COMMAND} -E copy
#   ${OBJ_DIR}/libDTdict_rdict.pcm ${LIB_DIR}/libDTdict_rdict.pcm
# )
# add_custom_command(TARGET ${libname} POST_BUILD
#   COMMAND ${CMAKE_COMMAND} -E copy
#   ${OBJ_DIR}/libMdict_rdict.pcm ${LIB_DIR}/libMdict_rdict.pcm
# )
# add_custom_command(TARGET ${libname} POST_BUILD
#   COMMAND ${CMAKE_COMMAND} -E copy
#   ${OBJ_DIR}/libDTdict.rootmap ${LIB_DIR}/libDTdict.rootmap
# )
# add_custom_command(TARGET ${libname} POST_BUILD
#   COMMAND ${CMAKE_COMMAND} -E copy
#   ${OBJ_DIR}/libMdict.rootmap ${LIB_DIR}/libMdict.rootmap
# )

# Boost unit test setup
find_package(Boost CONFIG COMPONENTS unit_test_framework REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# Add and configure unit tests
enable_testing()
add_subdirectory(${PROJECT_SOURCE_DIR}/tests/unit)
