# @author Sandro Wenzel
# @brief  cmake setup for module Utilities/MCStepLogger

set(MODULE_NAME "MCStepLogger")

# Minimum version of CMake
CMAKE_MINIMUM_REQUIRED(VERSION 3.11.0 FATAL_ERROR)

project(${MODULE_NAME})


# Install directories
SET(INSTALL_BIN_DIR ${CMAKE_INSTALL_PREFIX}/bin)
SET(INSTALL_INC_DIR ${CMAKE_INSTALL_PREFIX}/include/${MODULE_NAME})
SET(INSTALL_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib)

# Source directories
SET(IMP_SRC_DIR ${CMAKE_SOURCE_DIR}/src)
SET(INC_SRC_DIR ${CMAKE_SOURCE_DIR}/include/${MODULE_NAME})

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

########
# ROOT #
########
# Find ROOT and get useful functions from ROOT_USE_FILE,
# e.g. ROOT_GENERATE_DICTIONARY
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})
# Find ROOT headers
include_directories(${ROOT_INCLUDE_DIRS})

#########
# Boost #
#########
#find_package(Boost COMPONENTS container thread system timer program_options random filesystem chrono exception regex serialization log log_setup unit_test_framework date_time signals iostreams REQUIRED)
find_package(Boost COMPONENTS program_options chrono unit_test_framework REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

# Generate ROOT dictionary
SET(ROOT_DICT_LINKDEF_FILE ${IMP_SRC_DIR}/MCStepLoggerLinkDef.h)
SET(ROOT_DICT_NAME "G__${MODULE_NAME}")

ROOT_GENERATE_DICTIONARY(${ROOT_DICT_NAME} ${HEADERS} LINKDEF ${ROOT_DICT_LINKDEF_FILE})
# Files produced by the dictionary generation
SET(ROOT_DICT_LIB_FILES
    "${PROJECT_BINARY_DIR}/lib${MODULE_NAME}_rdict.pcm"
    "${PROJECT_BINARY_DIR}/lib${MODULE_NAME}.rootmap"
)

# Build a library from the sources specified above together with generated ROOT
# dictionary
add_library(${MODULE_NAME} SHARED ${SRCS} "${ROOT_DICT_NAME}.cxx" ${HEADERS})

# Link together with ROOT libs
target_link_libraries(${MODULE_NAME} -lCore -lHist -lGraf -lGpad -lTree -lVMC)

# Add the executable to do analysis with the MCStepLogger output files
add_executable(${EXECUTABLE_NAME} ${EXE_SRCS})
target_link_libraries(${EXECUTABLE_NAME} ${MODULE_NAME} ${Boost_LIBRARIES})

# Install headers
install(FILES ${HEADERS} DESTINATION ${INSTALL_INC_DIR})
# Install libraries
install(TARGETS ${MODULE_NAME} DESTINATION ${INSTALL_LIB_DIR})
# Install the ROOT dictionary files
install(FILES ${ROOT_DICT_LIB_FILES} DESTINATION ${INSTALL_LIB_DIR})
# Install executables
install(TARGETS ${EXECUTABLE_NAME} DESTINATION ${INSTALL_BIN_DIR})
