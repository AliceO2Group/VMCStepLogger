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

# The basic test only checks core functionality without any backend.
if(HAS_TEST_DATA)
  find_package(Boost REQUIRED COMPONENTS unit_test_framework)
  include_directories(${Boost_INCLUDE_DIR})

  set(CORE_TARGET analysis)
  set(TEST_CORE "test${CORE_TARGET}")
  set(STEP_TEST_FILE MCStepLoggerOutput.root)
  add_executable(${CORE_TARGET} "${CORE_TARGET}.cxx")
  target_link_libraries(${CORE_TARGET} MCStepLoggerCore MCStepLoggerAnalysis ${Boost_LIBRARIES})
  add_test(NAME ${TEST_CORE} COMMAND ${CORE_TARGET} -- ${TEST_DATA_DIR}/${CORE_TARGET}/${STEP_TEST_FILE} "./")
endif()

# The interception test runs fake applications under the preloaded logger and needs no test data.
# It needs a shared interceptor (-DBUILD_SHARED_LIBS=ON) to preload.
get_target_property(INTERCEPT_TYPE MCStepLoggerInterceptSteps TYPE)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND INTERCEPT_TYPE STREQUAL "SHARED_LIBRARY")
  # named like the O2 libraries so that the by-name fallback for RTLD_LOCAL loads finds them
  add_library(O2FairRootBase SHARED fakeApplication.cxx)
  add_library(O2Field SHARED fakeField.cxx)
  foreach(fake O2FairRootBase O2Field)
    target_link_libraries(${fake} PRIVATE VMCLibrary)
    set_target_properties(${fake} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fake)
  endforeach()

  add_executable(interceptDriverLinked interceptDriver.cxx)
  target_compile_definitions(interceptDriverLinked PRIVATE FAKE_LINKED)
  target_link_libraries(interceptDriverLinked O2FairRootBase O2Field VMCLibrary ROOT::RIO ROOT::Tree)
  # DT_RPATH so the fakes win over real O2 libraries on LD_LIBRARY_PATH
  target_link_options(interceptDriverLinked PRIVATE -Wl,--disable-new-dtags)
  add_executable(interceptDriverDlopen interceptDriver.cxx)
  target_link_libraries(interceptDriverDlopen VMCLibrary ROOT::RIO ROOT::Tree ${CMAKE_DL_LIBS})

  set(PRELOAD LD_PRELOAD=$<TARGET_FILE:MCStepLoggerInterceptSteps> MCSTEPLOG_TTREE=1)
  add_test(NAME testInterceptLinked
           COMMAND ${CMAKE_COMMAND} -E env ${PRELOAD} MCSTEPLOG_OUTFILE=interceptLinked.root
                   $<TARGET_FILE:interceptDriverLinked>)
  foreach(scope global local)
    add_test(NAME testInterceptDlopen_${scope}
             COMMAND ${CMAKE_COMMAND} -E env ${PRELOAD} MCSTEPLOG_OUTFILE=interceptDlopen_${scope}.root
                     $<TARGET_FILE:interceptDriverDlopen> ${scope}
                     $<TARGET_FILE:O2FairRootBase> $<TARGET_FILE:O2Field>)
  endforeach()
else()
  message(STATUS "Interception tests need Linux and -DBUILD_SHARED_LIBS=ON; not built")
endif()
