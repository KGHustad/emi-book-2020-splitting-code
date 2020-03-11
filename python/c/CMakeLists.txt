cmake_minimum_required (VERSION 3.9 FATAL_ERROR)

project(emi-fdm-hpc LANGUAGES C CXX)
set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_POSITION_INDEPENDENT_CODE True)

option(USE_OPENMP "Build with OpenMP" ON)

if(USE_OPENMP)
    find_package(OpenMP)
    if(NOT OpenMP_FOUND)
        message(WARNING "OpenMP was not found. Disabling USE_OPENMP")
        set(USE_OPENMP OFF CACHE BOOL "Build with OpenMP")
    endif()
endif()


if ("${PROJECT_BINARY_DIR}" STREQUAL "${PROJECT_SOURCE_DIR}")
    #message(${PROJECT_BINARY_DIR})
    #message(${PROJECT_SOURCE_DIR})
    message(FATAL_ERROR "You cannot build in the source directory. Please run cmake from a subdirectory called 'build'")
endif()

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_C_FLAGS_DEBUG "-g -Wall -DDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -DDEBUG")

set(CMAKE_C_FLAGS_RELEASE "-O3 -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

add_subdirectory(ionic_models)
add_subdirectory(linear_solvers)
