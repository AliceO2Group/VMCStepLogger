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

# @author Benedikt Volkel
# @brief  cmake setup for TMCReplay based on VMC

set(MODULE_NAME "MCReplay")

################################################################################
# Basic configurations
################################################################################
# Sources for built are here.
set(CXX_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(CXX_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Executables
set(EXECUTABLE_NAME mcreplay)

################################################################################
# Sources and headers
################################################################################

set(SRCS
   ${CXX_SOURCE_DIR}/MCReplayEngine.cxx
   ${CXX_SOURCE_DIR}/MCReplayGenericApplication.cxx
   ${CXX_SOURCE_DIR}/MCReplayGenericStack.cxx
   ${CXX_SOURCE_DIR}/MCReplayEvGen.cxx
   ${CXX_SOURCE_DIR}/MCReplayPhysics.cxx
)
set(HEADERS
   ${CXX_INCLUDE_DIR}/${MODULE_NAME}/MCReplayEngine.h
   ${CXX_INCLUDE_DIR}/${MODULE_NAME}/MCReplayGenericApplication.h
   ${CXX_INCLUDE_DIR}/${MODULE_NAME}/MCReplayGenericStack.h
   ${CXX_INCLUDE_DIR}/${MODULE_NAME}/MCReplayEvGen.h
   ${CXX_INCLUDE_DIR}/${MODULE_NAME}/MCReplayPhysics.h
)

set(EXECUTABLE_SRCS
   ${CXX_SOURCE_DIR}/replay.cxx
)

mcsl_add_library(
  TARGETNAME Core
  BASENAME ${MODULE_NAME}
  DEPENDENCIES ROOT::Core ROOT::Hist ROOT::Graf ROOT::Gpad ROOT::Tree ROOT::RIO ROOT::Geom ROOT::EG MCStepLoggerCore
  SOURCES ${SRCS}
  INCLUDE_DIRECTORIES ${CMAKE_CURRENT_LIST_DIR}/include ${Boost_INCLUDE_DIR}
  LINKDEFDIR ${CMAKE_CURRENT_LIST_DIR}/src
  ROOT_DICTIONARY_HEADERS ${HEADERS})


add_executable(${EXECUTABLE_NAME} ${EXECUTABLE_SRCS})
target_link_libraries(${EXECUTABLE_NAME} MCReplayCore Boost::program_options)
install(TARGETS ${EXECUTABLE_NAME} DESTINATION ${INSTALL_BIN_DIR} EXPORT ${CMAKE_PROJECT_NAME}Exports)

if(MCStepLogger_BUILD_TESTS AND HAS_TEST_DATA)
  add_subdirectory(test)
endif()
