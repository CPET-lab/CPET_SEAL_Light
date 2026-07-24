# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT license.
#
# Standalone project: NOT hooked into the root CMakeLists.txt on purpose.
# Build it against the already-built Microsoft SEAL in the repo's main
# `build` directory by pointing SEAL_DIR at its exported config, e.g.:
#
#   cmake -S native/experiments -B native/experiments/build -DSEAL_DIR=${repo}/build/cmake
#   cmake --build native/experiments/build --target <exe_name>
#
# Every .cpp file directly under this directory becomes its own executable,
# named after the file (e.g. boottest.cpp -> the `boottest` executable).
# After adding a NEW .cpp file, re-run the `cmake -S ... -B ...` configure
# step once so file(GLOB ...) below picks it up.

cmake_minimum_required(VERSION 3.13)

project(SEALExperiments VERSION 4.1.2 LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type" FORCE)
endif()

find_package(SEAL 4.1.2 EXACT REQUIRED)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

if(TARGET SEAL::seal)
    set(seal_target SEAL::seal)
elseif(TARGET SEAL::seal_shared)
    set(seal_target SEAL::seal_shared)
else()
    message(FATAL_ERROR "Cannot find target SEAL::seal or SEAL::seal_shared")
endif()

file(GLOB SEAL_EXPERIMENTS_SOURCES ${CMAKE_CURRENT_LIST_DIR}/*.cpp)

foreach(source_file ${SEAL_EXPERIMENTS_SOURCES})
    get_filename_component(exe_name ${source_file} NAME_WE)
    add_executable(${exe_name} ${source_file})
    target_link_libraries(${exe_name} PRIVATE ${seal_target})
endforeach()
