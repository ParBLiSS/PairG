cmake_minimum_required(VERSION 3.13)

# Project settings
project(pairg LANGUAGES CXX)

## Options
# Check if the user wants to only build the library
option(BUILD_LIBRARY_ONLY "Inform whether only library should be built" OFF)
# Check if the user wants to build google test applications
option(BUILD_TESTING "Inform whether test applications should be built" OFF)

# Adding CMake module path.
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
## Finding dependencies
find_package(Kokkos QUIET)
find_package(KokkosKernels QUIET)
# Include external dependencies
include(${PROJECT_SOURCE_DIR}/ext/external.cmake)

## Defining target library `libpairg`
# Setting header files
file(GLOB HEADER_FILES "${PROJECT_SOURCE_DIR}/include/pairg/*.hpp")
add_library(libpairg INTERFACE)
# Defining target sources
target_sources(libpairg INTERFACE "$<BUILD_INTERFACE:${HEADER_FILES}>")
# Defining include directories
target_include_directories(libpairg
  INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  SYSTEM INTERFACE $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>)
# Defining link libraries
target_link_libraries(libpairg
  INTERFACE Kokkos::kokkos
  INTERFACE Kokkos::kokkoskernels)
# Use C++14
target_compile_features(libpairg INTERFACE cxx_std_14)
# Defining `pairg::libpairg` alias
add_library(pairg::libpairg ALIAS libpairg)
# Install targets
install(FILES ${HEADER_FILES} DESTINATION include/pairg)

# Save libs and executables in the same place
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib CACHE PATH "Output directory for libraries" )

# Build PairG command-line tool
if(NOT BUILD_LIBRARY_ONLY)
  add_subdirectory("${PROJECT_SOURCE_DIR}/src")
endif(NOT BUILD_LIBRARY_ONLY)

# Generate test executables if requested
if(BUILD_TESTING)
  add_subdirectory("${PROJECT_SOURCE_DIR}/tests")
endif(BUILD_TESTING)
