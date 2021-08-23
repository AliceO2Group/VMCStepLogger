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


# On the shoulders of https://github.com/AliceO2Group/AliceO2/blob/dev/cmake/O2AddLibrary.cmake

include_guard()

# Common CMake functionality

function(mcsl_add_library)
  cmake_parse_arguments(
    PARSE_ARGV
    0
    A
    ""
    "TARGETNAME;BASENAME"
    "DEPENDENCIES;SOURCES;INCLUDE_DIRECTORIES;ROOT_DICTIONARY_HEADERS;LINKDEFDIR"
  )
  if(A_UNPARSED_ARGUMENTS)
    message(
      FATAL_ERROR "Unexpected unparsed arguments: ${A_UNPARSED_ARGUMENTS}")
  endif()

  set(basename ${A_BASENAME})
  set(target "${basename}${A_TARGETNAME}")

  add_library(${target} ${A_SOURCES})
  #set(includeDirs $<TARGET_PROPERTY:${target},INCLUDE_DIRECTORIES>)
  if(A_INCLUDE_DIRECTORIES)
    foreach(d IN LISTS A_INCLUDE_DIRECTORIES)
      get_filename_component(adir ${d} ABSOLUTE)
      if(NOT IS_DIRECTORY ${adir})
        message(
          FATAL_ERROR "Trying to append non existing include directory ${d}")
      endif()
      target_include_directories(${target} PUBLIC $<BUILD_INTERFACE:${adir}>)
    endforeach()
  endif()


  if(A_ROOT_DICTIONARY_HEADERS)

    # Build the LD_LIBRARY_PATH required to get rootcling running fine
    #
    # Need at least root core library
    get_filename_component(LD_LIBRARY_PATH ${ROOT_Core_LIBRARY} DIRECTORY)
    # and possibly toolchain libs if we are using a toolchain
    if(DEFINED ENV{GCC_TOOLCHAIN_ROOT})
      set(LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:$ENV{GCC_TOOLCHAIN_ROOT}/lib")
      set(LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:$ENV{GCC_TOOLCHAIN_ROOT}/lib64")
    endif()

    set(dictionary G__${target})
    set(pcmFile "${CMAKE_CURRENT_BINARY_DIR}/${dictionary}_rdict.pcm")
    set(rootmapFile "${CMAKE_CURRENT_BINARY_DIR}/lib${target}.rootmap")
    set(dictionaryFile "G__${target}.cxx")
    set(linkdef "${A_LINKDEFDIR}/${target}LinkDef.h")
    message("LINKDEF for target ${target} assumed at ${linkdef}")
    # Since this is actually generated before libraries are built, we have the dictionary files available even though this is inside a CMake function
    add_custom_command(
      OUTPUT ${dictionaryFile} ${pcmFile} ${rootmapFile}
      VERBATIM
      COMMAND
      ${CMAKE_COMMAND} -E env LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ${ROOT_rootcling_CMD}
      -f ${dictionaryFile}
      -rmf ${rootmapFile}
      -rml ${target}
      -noGlobalUsingStd
      -inlineInputHeader
      -I$<JOIN:$<REMOVE_DUPLICATES:$<TARGET_PROPERTY:${target},INCLUDE_DIRECTORIES>>,$<SEMICOLON>-I>
      ${A_ROOT_DICTIONARY_HEADERS} ${linkdef}
      DEPENDS
      ${A_ROOT_DICTIONARY_HEADERS} ${linkdef}
      COMMAND_EXPAND_LISTS
    )
    # Now add to target
    target_sources(${target} PRIVATE ${dictionaryFile})
  endif()

  if(A_DEPENDENCIES)
    target_link_libraries(${target} PUBLIC ${A_DEPENDENCIES})
    set_target_properties(${target} PROPERTIES INTERFACE_LINK_LIBRARIES "${A_DEPENDENCIES}")
  endif()

  install(TARGETS ${target}
          EXPORT ${CMAKE_PROJECT_NAME}Exports
          INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
          DESTINATION ${CMAKE_INSTALL_LIBDIR})
  install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/include/${basename}
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  if(A_ROOT_DICTIONARY_HEADERS)
    install(FILES ${rootmapFile} ${pcmFile} DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif()
endfunction()
