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

# Check if the user wants to build google test applications
OPTION(BUILD_TESTS "Inform whether test applications should be built" ON)

# Save libs and executables in the same place
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib CACHE PATH "Output directory for libraries" )

include_directories("${PROJECT_SOURCE_DIR}/ext")
include_directories("${PROJECT_SOURCE_DIR}/src/include")

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

file(WRITE "${PROJECT_BINARY_DIR}/KokkosKernels_config.h" "//dummy file for kokkos-kernels to include")
include_directories("${PROJECT_BINARY_DIR}")

### Include and compile htslib
include_directories("${PROJECT_SOURCE_DIR}/ext/htslib")
add_custom_target(LIBHTS COMMAND make WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/ext/htslib" COMMENT "Compiling htslib") 
set(HTS_LIBRARY "${PROJECT_SOURCE_DIR}/ext/htslib/libhts.so")


### Protobuf stuff
#NOTE: user is expected to set the following variable PROTOBUF_DIR 
set(PROTOBUF_DIR CACHE STRING "absolute path to protobuf installation directory")
set(PROTOBUF_INCLUDE_DIR "${PROTOBUF_DIR}/include")
set(PROTOBUF_LIBRARY "${PROTOBUF_DIR}/lib/libprotobuf.so")
set(PROTOBUF_PROTOC_EXECUTABLE "${PROTOBUF_DIR}/bin/protoc")
if(EXISTS "${PROTOBUF_INCLUDE_DIR}" AND EXISTS "${PROTOBUF_LIBRARY}" AND EXISTS "${PROTOBUF_PROTOC_EXECUTABLE}")
  message(STATUS "Protobuf libraries, header files and executables located")
else()
  message(FATAL_ERROR "Protobuf path is invalid, cannot find required libraries, header files or executables. Set the cmake variable PROTOBUF_DIR properly.")
endif()
find_package(Protobuf REQUIRED)
include_directories(${PROTOBUF_INCLUDE_DIR})

### Include libvgio
PROTOBUF_GENERATE_CPP(PROTO_SRCS PROTO_HDRS "${PROJECT_SOURCE_DIR}/ext/libvgio/deps/vg.proto")
add_custom_target(SYMLNK COMMAND ${CMAKE_COMMAND} -E create_symlink . vg  WORKING_DIRECTORY "${PROJECT_BINARY_DIR}" COMMENT "Adding symbolic link vg") 
include_directories("${PROJECT_SOURCE_DIR}/ext/libvgio/include")
file(GLOB VGIO_SRC "${PROJECT_SOURCE_DIR}/ext/libvgio/src/*.cpp")
add_library(LIBVGIO ${PROTO_SRCS} ${VGIO_SRC})
add_dependencies(LIBVGIO SYMLNK)


###generate test executables if requested
add_subdirectory("${PROJECT_SOURCE_DIR}/tests")

###generate main executable
add_executable(pairmap2graph "${PROJECT_SOURCE_DIR}/src/main.cpp")
add_dependencies(pairmap2graph LIBHTS SYMLNK)
target_link_libraries(pairmap2graph kokkos ${PROTOBUF_LIBRARY} LIBVGIO ${HTS_LIBRARY})
