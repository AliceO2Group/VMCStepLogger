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

find_package(Boost REQUIRED COMPONENTS unit_test_framework)
include_directories(${Boost_INCLUDE_DIR})

# The basic test only checks core functionality without any backend.
set(CORE_TARGET analysis)
set(TEST_CORE "test${CORE_TARGET}")
set(STEP_TEST_FILE MCStepLoggerOutput.root)
add_executable(${CORE_TARGET} "${CORE_TARGET}.cxx")
target_link_libraries(${CORE_TARGET} MCStepLoggerCore MCStepLoggerAnalysis ${Boost_LIBRARIES})
add_test(NAME ${TEST_CORE} COMMAND ${CORE_TARGET} -- ${TEST_DATA_DIR}/${CORE_TARGET}/${STEP_TEST_FILE} "./")
