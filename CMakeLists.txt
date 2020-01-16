# @author Sandro Wenzel
# @brief  cmake setup for module Utilities/MCStepLogger

set(MODULE_NAME "MCStepLogger")

# Minimum version of CMake
CMAKE_MINIMUM_REQUIRED(VERSION 3.11.0 FATAL_ERROR)

project(${MODULE_NAME})
set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake)

# Install directories
SET(INSTALL_BIN_DIR ${CMAKE_INSTALL_PREFIX}/bin)
SET(INSTALL_MACRO_DIR ${CMAKE_INSTALL_PREFIX}/macro)
SET(INSTALL_INC_DIR ${CMAKE_INSTALL_PREFIX}/include/${MODULE_NAME})
SET(INSTALL_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib)
SET(INSTALL_CMAKE_DIR ${CMAKE_INSTALL_PREFIX}/cmake)

# Source directories
SET(IMP_SRC_DIR ${CMAKE_SOURCE_DIR}/src)
SET(INC_SRC_DIR ${CMAKE_SOURCE_DIR}/include/${MODULE_NAME})
SET(MACRO_SRC_DIR ${CMAKE_SOURCE_DIR}/macro)

# Library name
SET(LIBRARY_NAME ${MODULE_NAME})
# Executable name
SET(EXECUTABLE_NAME mcStepAnalysis)

# Required source files to build the library.
set(SRCS
    ${IMP_SRC_DIR}/MCStepInterceptor.cxx
    ${IMP_SRC_DIR}/MCStepLoggerImpl.cxx
    ${IMP_SRC_DIR}/StepInfo.cxx
    ${IMP_SRC_DIR}/MCAnalysis.cxx
    ${IMP_SRC_DIR}/BasicMCAnalysis.cxx
    ${IMP_SRC_DIR}/SimpleStepAnalysis.cxx
    ${IMP_SRC_DIR}/MCAnalysisManager.cxx
    ${IMP_SRC_DIR}/MCAnalysisFileWrapper.cxx
    ${IMP_SRC_DIR}/MCAnalysisUtilities.cxx
    ${IMP_SRC_DIR}/ROOTIOUtilities.cxx
   )

# Requried headers to build the library.
set(HEADERS
   ${INC_SRC_DIR}/StepInfo.h
   ${INC_SRC_DIR}/MetaInfo.h
   ${INC_SRC_DIR}/MCAnalysis.h
   ${INC_SRC_DIR}/BasicMCAnalysis.h
   ${INC_SRC_DIR}/SimpleStepAnalysis.h
   ${INC_SRC_DIR}/MCAnalysisManager.h
   ${INC_SRC_DIR}/MCAnalysisFileWrapper.h
   ${INC_SRC_DIR}/MCAnalysisUtilities.h
   ${INC_SRC_DIR}/ROOTIOUtilities.h
  )
include_directories(include/)

# Required source for the executable
set(EXE_SRCS
    ${IMP_SRC_DIR}/analyseMCSteps.cxx
   )

# Macros to be copied
set(MACROS
   ${MACRO_SRC_DIR}/plotAnalysisHistograms.C
  )

########
# ROOT #
########
# Find ROOT and get useful functions from ROOT_USE_FILE,
# e.g. ROOT_GENERATE_DICTIONARY
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS Core Hist Graf Gpad Tree VMC RIO Geom EG)
include(${ROOT_USE_FILE})

#########
# Boost #
#########
find_package(Boost COMPONENTS program_options chrono unit_test_framework REQUIRED)

# Generate ROOT dictionary
SET(ROOT_DICT_LINKDEF_FILE ${IMP_SRC_DIR}/MCStepLoggerLinkDef.h)
SET(ROOT_DICT_NAME "G__${MODULE_NAME}")

ROOT_GENERATE_DICTIONARY(${ROOT_DICT_NAME} ${HEADERS} LINKDEF ${ROOT_DICT_LINKDEF_FILE})
# Files produced by the dictionary generation
SET(ROOT_DICT_LIB_FILES
    "${PROJECT_BINARY_DIR}/lib${MODULE_NAME}_rdict.pcm"
    "${PROJECT_BINARY_DIR}/lib${MODULE_NAME}.rootmap"
)

###############################
# Determine CXX STD from ROOT #
###############################
SET(CMAKE_CXX_STANDARD 11)
# Find ROOT CXX standard
string(FIND ${ROOT_CXX_FLAGS} "-std=" POSITION)
if (${POSITION} GREATER -1)
    string(SUBSTRING ${ROOT_CXX_FLAGS} ${POSITION} 11 ROOT_CXX_STD)
    if(${ROOT_CXX_STD} STREQUAL "-std=c++1z " OR ${ROOT_CXX_STD} STREQUAL "-std=c++17 ")
        SET(CMAKE_CXX_STANDARD 17)
    elseif(${ROOT_CXX_STD} STREQUAL "-std=c++1y " OR ${ROOT_CXX_STD} STREQUAL "-std=c++14 ")
        SET(CMAKE_CXX_STANDARD 14)
    endif()
endif()
message(STATUS "Build with CXX STD ${CMAKE_CXX_STANDARD}")

# Build a library from the sources specified above together with generated ROOT
# dictionary
add_library(${MODULE_NAME} SHARED ${SRCS} "${ROOT_DICT_NAME}.cxx" ${HEADERS})

# Link together with ROOT libs
set(LIB_DEPS ROOT::Core ROOT::Hist ROOT::Graf ROOT::Gpad ROOT::Tree ROOT::VMC ROOT::RIO ROOT::Geom ROOT::EG)
target_link_libraries(${MODULE_NAME} ${LIB_DEPS})
set_target_properties(${MODULE_NAME} PROPERTIES INTERFACE_LINK_LIBRARIES "${LIB_DEPS}")
target_include_directories(${MODULE_NAME} INTERFACE $<INSTALL_INTERFACE:include>)

# Install headers
install(FILES ${HEADERS} DESTINATION ${INSTALL_INC_DIR})

# Install libraries (and add to export such that we can easily find its dependencies later)
install(TARGETS ${MODULE_NAME} DESTINATION ${INSTALL_LIB_DIR} EXPORT ${CMAKE_PROJECT_NAME}Exports)

# Install the ROOT dictionary files
install(FILES ${ROOT_DICT_LIB_FILES} DESTINATION ${INSTALL_LIB_DIR})

# Install executables (and add to export such that we can easily find its dependencies later)
add_executable(${EXECUTABLE_NAME} ${EXE_SRCS})
set(EXE_DEPS ${MODULE_NAME} Boost::program_options Boost::chrono Boost::unit_test_framework)
target_link_libraries(${EXECUTABLE_NAME} ${EXE_DEPS})
install(TARGETS ${EXECUTABLE_NAME} DESTINATION ${INSTALL_BIN_DIR} EXPORT ${CMAKE_PROJECT_NAME}Exports)

# Install macros
install(FILES ${MACROS} DESTINATION ${INSTALL_MACRO_DIR})

# Install exports to find targets and dependencies later
install(EXPORT ${CMAKE_PROJECT_NAME}Exports NAMESPACE MCStepLogger:: FILE MCStepLoggerTargets.cmake DESTINATION ${INSTALL_CMAKE_DIR})

# Configure and install the config to be used by CMake of depending package
configure_file(
    "${PROJECT_SOURCE_DIR}/cmake/MCStepLoggerConfig.cmake.in"
    "${PROJECT_BINARY_DIR}/MCStepLoggerConfig.cmake" @ONLY)
install(FILES
    "${PROJECT_BINARY_DIR}/MCStepLoggerConfig.cmake"
    DESTINATION ${INSTALL_CMAKE_DIR})
