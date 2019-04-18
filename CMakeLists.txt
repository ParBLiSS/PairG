cmake_minimum_required(VERSION 3.1)
include(CheckCXXCompilerFlag)

# project settings
project(pairmap2graph)

##### General Compilation Settings
set(CMAKE_CXX_FLAGS "-w ${CMAKE_CXX_FLAGS} --std=c++11")

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")

#Set default cmake build type to RelWithDebInfo during prototyping
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo"
    FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# Save libs and executables in the same place
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib CACHE PATH "Output directory for libraries" )

#### OpenMP
find_package(OpenMP)
# FindOpenMP defines the OpenMP_C_FLAGS and OpenMP_CXX_FLAGS.
if (OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  message(STATUS "Found OpenMP")
endif()

#### kokkos stuff
include_directories("${PROJECT_SOURCE_DIR}/ext/kokkos-kernels/src")
include_directories("${PROJECT_SOURCE_DIR}/ext/kokkos-kernels/src/impl")
include_directories("${PROJECT_SOURCE_DIR}/ext/kokkos-kernels/src/impl/tpls")
include_directories("${PROJECT_SOURCE_DIR}/ext/kokkos-kernels/src/common")
include_directories("${PROJECT_SOURCE_DIR}/ext/kokkos-kernels/src/graph")
include_directories("${PROJECT_SOURCE_DIR}/ext/kokkos-kernels/src/graph/impl")
include_directories("${PROJECT_SOURCE_DIR}/ext/kokkos-kernels/src/sparse")
include_directories("${PROJECT_SOURCE_DIR}/ext/kokkos-kernels/src/sparse/impl")

add_subdirectory("${PROJECT_SOURCE_DIR}/ext/kokkos" ${PROJECT_BINARY_DIR}/kokkos)
include_directories(${Kokkos_INCLUDE_DIRS_RET})

add_executable(pairmap2graph "${PROJECT_SOURCE_DIR}/src/spgemm.cpp")
target_link_libraries(pairmap2graph kokkos)
