# Copyright © 2017-2018 Martin Ueding <dev@martin-ueding.de>

cmake_minimum_required(VERSION 3.9)

project(sLapH-contractions CXX C)

set(CMAKE_CXX_STANDARD 14)

add_library(lcontract
  src/Correlators.cpp
  src/Diagram.cpp
  src/dilution-iterator.cpp
  src/EigenVector.cpp
  src/Gamma.cpp
  src/GaugeField.cpp
  src/global_data.cpp
  src/global_data_build_IO_names.cpp
  src/global_data_input_handling.cpp
  src/global_data_input_handling_utils.cpp
  src/init_lookup_tables.cpp
  src/h5-wrapper.cpp
  src/io_utils.cpp
  src/OperatorFactory.cpp
  src/Perambulator.cpp
  src/DilutedFactorY.cpp
  src/DilutedTraceFactory.cpp
  src/RandomVector.cpp
  src/ranlxs.cpp
  src/StopWatch.cpp
  )

target_include_directories(lcontract PUBLIC "include")

###############################################################################
#                          Git Revision Information                           #
###############################################################################

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(GetGitRevisionDescription)

get_git_head_revision(git_refspec git_sha1)
git_local_changes(git_changes)

message(STATUS "git sha1: ${git_sha1}")
message(STATUS "git refspec: ${git_refspec}")
message(STATUS "git state: ${git_changes}")

execute_process(COMMAND
  hostname -f
  OUTPUT_VARIABLE hostname
  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

set(git_user $ENV{USER})
set(git_host ${hostname})

message(STATUS "user: ${git_user}")
message(STATUS "host: ${git_host}")

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/include/git.hpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/include/git.hpp"
  @ONLY)
list(APPEND SOURCES "${CMAKE_CURRENT_BINARY_DIR}/include/git.hpp")
target_include_directories(lcontract PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/include")

###############################################################################
#                             External Libraries                              #
###############################################################################

# HDF5
if(NOT DEFINED HDF5_INCLUDE_DIRS OR NOT DEFINED HDF5_LIBRARIES)
  find_package(HDF5 REQUIRED COMPONENTS C CXX)
endif()
target_include_directories(lcontract PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(lcontract PUBLIC ${HDF5_LIBRARIES})

# Boost
if(NOT DEFINED Boost_INCLUDE_DIRS OR NOT DEFINED Boost_LIBRARIES)
  find_package(Boost REQUIRED COMPONENTS filesystem system program_options)
endif()
target_include_directories(lcontract PUBLIC ${Boost_INCLUDE_DIRS})
target_link_libraries(lcontract PUBLIC ${Boost_LIBRARIES})

# Eigen
if(NOT DEFINED EIGEN3_INCLUDE_DIR)
  find_package(Eigen3 REQUIRED)
endif()
target_include_directories(lcontract SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIR})
target_compile_definitions(lcontract PUBLIC EIGEN_DONT_PARALLELIZE)

# Lime
if(NOT DEFINED LIME_INCLUDE_DIRS)
  message(FATAL_ERROR "Must specify -DLIME_INCLUDE_DIRS on the command line to CMake.")
endif()
target_include_directories(lcontract SYSTEM PUBLIC ${LIME_INCLUDE_DIRS})
target_link_libraries(lcontract PUBLIC ${LIME_LIBRARIES})

# Lime
if(NOT DEFINED LIME_LIBRARIES)
  message(FATAL_ERROR "Must specify -DLIME_LIBRARIES on the command line to CMake.")
endif()

# OpenMP
find_package(OpenMP)
target_compile_options(lcontract PUBLIC ${OpenMP_CXX_FLAGS})
target_link_libraries(lcontract PUBLIC ${OpenMP_CXX_LIBRARIES} ${OpenMP_CXX_FLAGS})

###############################################################################
#                                 Build Flags                                 #
###############################################################################

# We want warnings, but not all of them.
add_definitions(-Wall -pedantic)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-fdiagnostics-color has_diagnostics_color)
if(has_diagnostics_color)
  if (CMAKE_GENERATOR STREQUAL "Ninja")
    add_definitions("-fdiagnostics-color=always")
  else()
    add_definitions("-fdiagnostics-color=auto")
  endif()
endif()

option(SLAPH_TIMING_CORRELATOR "Timing of correlator building (per thread)" OFF)
option(SLAPH_TIMING_DIAGRAMS "Timing of individual diagram multiplies and traces" OFF)
option(SLAPH_TIMING_ULTRA_FINE "Timing of DilutedFactor components" OFF)
option(SLAPH_TIMING_FINE "Timing of QuarkLineBlock constructions" OFF)

if(SLAPH_TIMING_CORRELATOR)
  add_definitions("-DSLAPH_TIMING_CORRELATOR")
endif()
if(SLAPH_TIMING_DIAGRAMS)
  add_definitions("-DSLAPH_TIMING_DIAGRAMS")
endif()
if(SLAPH_TIMING_ULTRA_FINE)
  add_definitions("-DSLAPH_TIMING_ULTRA_FINE")
endif()
if(SLAPH_TIMING_FINE)
  add_definitions("-DSLAPH_TIMING_FINE")
endif()

# Default to “Release” build type.
message(STATUS "Build Type: '${CMAKE_BUILD_TYPE}'")
if(CMAKE_BUILD_TYPE STREQUAL "")
  message(STATUS "No CMAKE_BUILD_TYPE, assuming release and enabling all optimizations.")
  set(CMAKE_BUILD_TYPE Release)
endif()

# Tune for the local architecture because this code will be used where it is
# compiled, in most cases.
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

###############################################################################
#                                 Executables                                 #
###############################################################################

### contract

add_executable(contract main/contract.cpp)
target_link_libraries(contract lcontract)

install(TARGETS contract DESTINATION bin)

### test-iterator
add_executable(test-dilution-iterator main/test-iterator.cpp)
target_link_libraries(test-dilution-iterator lcontract)

### test-derivative
add_executable(test-derivative main/test-derivative.cpp)
target_link_libraries(test-derivative lcontract)

###############################################################################
#                               Subdirectories                                #
###############################################################################

add_subdirectory(benchmarks)
add_subdirectory(tests/gtest)

###############################################################################
#                                    Tests                                    #
###############################################################################

enable_testing()

add_test(NAME unit-dilution-iterator
  COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/tests/dilution-iterator/run-test" "${CMAKE_CURRENT_SOURCE_DIR}")

add_test(NAME integration-L4
  COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration-L4/run-integration-test" "${CMAKE_CURRENT_BINARY_DIR}"
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration-L4")
