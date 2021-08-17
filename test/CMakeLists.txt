# Build tests

set(REPO_NAME_TEST MCStepLoggerTestData)

# Could potentially use CMake's FetchContent module but seems to be overhead since we only need this single line
execute_process(COMMAND git clone https://gitlab.cern.ch/bvolkel/${REPO_NAME_TEST}.git WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} RESULT_VARIABLE res)

if(NOT res EQUAL "0")
  message(WARNING "Test data could not be fetched")
  return()
endif(NOT res EQUAL "0")

# Find boost unit tests dependencies
find_package(Boost REQUIRED COMPONENTS unit_test_framework)
include_directories(${Boost_INCLUDE_DIR})

# The basic test only checks core functionality without any backend.
set(CORE_TARGET analysis)
set(TEST_CORE "test${CORE_TARGET}")
set(STEP_TEST_FILE MCStepLoggerOutput.root)
add_executable(${CORE_TARGET} "${CORE_TARGET}.cxx")
target_link_libraries(${CORE_TARGET} MCStepLoggerCore MCStepLoggerAnalysis ${Boost_LIBRARIES})
add_test(NAME ${TEST_CORE} COMMAND ${CORE_TARGET} -- ${REPO_NAME_TEST}/${CORE_TARGET}/${STEP_TEST_FILE} "./")
