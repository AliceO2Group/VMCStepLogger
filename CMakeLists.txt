# Copyright 2019-2020 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

# @author Sandro Wenzel
# @brief  cmake setup for module Utilities/MCStepLogger

# Minimum version of CMake
CMAKE_MINIMUM_REQUIRED(VERSION 3.15.0 FATAL_ERROR)

project(MCStepLogger)
set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake)

option(MCStepLogger_BUILD_MCReplay "Build MCReplay engine" ON)
option(MCStepLogger_BUILD_TESTS "Control whether tests are built" OFF)

include(CTest)

# Install directories
SET(INSTALL_BIN_DIR ${CMAKE_INSTALL_PREFIX}/bin)
SET(INSTALL_MACRO_DIR ${CMAKE_INSTALL_PREFIX}/macro)
SET(INSTALL_INC_DIR ${CMAKE_INSTALL_PREFIX}/include)
SET(INSTALL_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib)
SET(INSTALL_CMAKE_DIR ${CMAKE_INSTALL_PREFIX}/cmake)

# setup for rpath and mac on the shoulders of O2 CMake
include(GNUInstallDirs)
if (APPLE)
  set(basePoint @loader_path)
else()
  set(basePoint $ORIGIN)
endif()
set(CMAKE_SKIP_BUILD_RPATH FALSE)
if (APPLE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
endif()

# add to the install RPATH the (automatically determined) parts of the RPATH
# that point to directories outside the build tree
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# specify libraries directory relative to binaries one.
file(RELATIVE_PATH relDir ${INSTALL_BIN_DIR} ${INSTALL_LIB_DIR})
set(CMAKE_INSTALL_RPATH ${basePoint} ${basePoint}/${relDir})

#################################
# Whether or not to build tests #
#################################
set(HAS_TEST_DATA FALSE)
set(REPO_NAME_TEST MCStepLoggerTestData)
set(TEST_DATA_DIR)
if(NOT MCStepLogger_BUILD_TESTS)
  message(WARNING "No tests are built.")
else()
  # Could potentially use CMake's FetchContent module but seems to be overhead since we only need this single line
  execute_process(COMMAND git clone https://gitlab.cern.ch/bvolkel/${REPO_NAME_TEST}.git WORKING_DIRECTORY ${CMAKE_BINARY_DIR} RESULT_VARIABLE res)
  if(NOT res EQUAL "0")
    message(WARNING "Test data could not be fetched")
  else()
    set(HAS_TEST_DATA TRUE)
    set(TEST_DATA_DIR ${CMAKE_BINARY_DIR}/${REPO_NAME_TEST})
  endif(NOT res EQUAL "0")
endif(NOT MCStepLogger_BUILD_TESTS)


#################
# Find packages #
#################

########
# ROOT #
########
# Find ROOT and get useful functions from ROOT_USE_FILE,
# e.g. ROOT_GENERATE_DICTIONARY
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS Core Hist Graf Gpad Tree RIO Geom EG)
include(${ROOT_USE_FILE})

#######
# VMC #
#######
find_package(VMC CONFIG)
if(NOT VMC_FOUND)
  if(NOT ROOT_vmc_FOUND)
    message(FATAL_ERROR "Could neither find VMC standalone nor ROOT's builtin")
  endif()
  add_library(VMCLibrary ALIAS ROOT::VMC)
endif()

#########
# Boost #
#########
find_package(Boost COMPONENTS program_options chrono unit_test_framework REQUIRED)

###############################
# Determine CXX STD from ROOT #
###############################

#
# SW: The following code was there but was
# problematic. In principle, we expect the C++ version
# to be set from the outside via (-DCMAKE_CXX_STANDARD).
# Please send email if this does not work for you.
#

#SET(CMAKE_CXX_STANDARD 11)
# Find ROOT CXX standard
#string(FIND ${ROOT_CXX_FLAGS} "-std=" POSITION)
#if (${POSITION} GREATER -1)
#    string(SUBSTRING ${ROOT_CXX_FLAGS} ${POSITION} 11 ROOT_CXX_STD)
#    if(${ROOT_CXX_STD} STREQUAL "-std=c++1z " OR ${ROOT_CXX_STD} STREQUAL "-std=c++17 ")
#        SET(CMAKE_CXX_STANDARD 17)
#    elseif(${ROOT_CXX_STD} STREQUAL "-std=c++1y " OR ${ROOT_CXX_STD} STREQUAL "-std=c++14 ")
#        SET(CMAKE_CXX_STANDARD 14)
#    endif()
#endif()
#message(STATUS "Build with CXX STD ${CMAKE_CXX_STANDARD}")

#######################
# Do the installation #
#######################
include(MCStepLoggerUtils)
set(CMAKE_INSTALL_LIBDIR ${INSTALL_LIB_DIR})
set(CMAKE_INSTALL_INCLUDEDIR ${INSTALL_INC_DIR})

add_subdirectory(MCStepLogger)
if(MCStepLogger_BUILD_MCReplay)
  add_subdirectory(MCReplay)
endif()

##################################################
# Do final global installation and configuration #
##################################################
install(EXPORT ${CMAKE_PROJECT_NAME}Exports FILE MCStepLoggerTargets.cmake DESTINATION ${INSTALL_CMAKE_DIR})

# Configure and install the config to be used by CMake of depending package
configure_file(
    "${PROJECT_SOURCE_DIR}/cmake/MCStepLoggerConfig.cmake.in"
    "${PROJECT_BINARY_DIR}/MCStepLoggerConfig.cmake" @ONLY)
install(FILES
    "${PROJECT_BINARY_DIR}/MCStepLoggerConfig.cmake"
    DESTINATION ${INSTALL_CMAKE_DIR})
