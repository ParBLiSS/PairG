if(NOT TARGET Kokkos::kokkos)
  message(STATUS "Using bundled Kokkos library")
  set(Kokkos_SOURCE_DIR ${PROJECT_SOURCE_DIR}/ext/kokkos)
  execute_process(
    COMMAND git submodule update --init --recursive -- ${Kokkos_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  add_subdirectory(${Kokkos_SOURCE_DIR})
endif()

if(NOT TARGET Kokkos::kokkoskernels)
  message(STATUS "Using bundled kokkos-kernels library")
  set(KokkosKernels_SOURCE_DIR ${PROJECT_SOURCE_DIR}/ext/kokkos-kernels)
  execute_process(
    COMMAND git submodule update --init --recursive -- ${KokkosKernels_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  add_subdirectory(${KokkosKernels_SOURCE_DIR})
endif()
