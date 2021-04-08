cmake_minimum_required(VERSION 3.13)

# Project settings
project(pairg LANGUAGES CXX)

## Options
# Check if the user wants to only build the library
option(BUILD_LIBRARY_ONLY "Inform whether only library should be built" OFF)
# Check if the user wants to build google test applications
option(BUILD_TESTING "Inform whether test applications should be built" OFF)

# Include external modules
include(GNUInstallDirs)
# Adding CMake module path.
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
## Finding dependencies
find_package(Kokkos QUIET)
find_package(KokkosKernels QUIET)
# Handle bundled dependencies
include(${PROJECT_SOURCE_DIR}/ext/external.cmake)

## Defining target library `libpairg`
# Defining the INTERFACE library
add_library(libpairg INTERFACE)
# Defining include directories
target_include_directories(libpairg
  INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>;$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>)
# Defining link libraries
target_link_libraries(libpairg
  INTERFACE $<BUILD_INTERFACE:Kokkos::kokkos>;$<INSTALL_INTERFACE:Kokkos::kokkos>
  INTERFACE $<BUILD_INTERFACE:Kokkos::kokkoskernels>;$<INSTALL_INTERFACE:Kokkos::kokkoskernels>)
# Use C++14
target_compile_features(libpairg INTERFACE cxx_std_14)
# Setting header files
file(GLOB HEADER_FILES RELATIVE "${PROJECT_SOURCE_DIR}/include" "${PROJECT_SOURCE_DIR}/include/pairg/*.hpp")
list(TRANSFORM HEADER_FILES PREPEND "${PROJECT_SOURCE_DIR}/include/" OUTPUT_VARIABLE BUILD_HEADER_FILES)
list(TRANSFORM HEADER_FILES PREPEND "${CMAKE_INSTALL_FULL_INCLUDEDIR}/" OUTPUT_VARIABLE INSTALL_HEADER_FILES)
# Defining target sources
target_sources(libpairg
  INTERFACE "$<BUILD_INTERFACE:${BUILD_HEADER_FILES}>;$<INSTALL_INTERFACE:${INSTALL_HEADER_FILES}>")
# Defining `pairg::libpairg` alias
add_library(pairg::libpairg ALIAS libpairg)
# Install targets
install(FILES ${BUILD_HEADER_FILES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pairg)

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
