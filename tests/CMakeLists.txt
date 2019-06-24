##### add definitions to find path to test data files
add_definitions(-DPROJECT_TEST_DATA_DIR=${PROJECT_SOURCE_DIR}/tests/data)

if (BUILD_TESTS)

  add_executable(test_graph_load test_main.cpp test_graphLoad.cpp)
  add_dependencies(test_graph_load LIBHTS SYMLNK)
  target_link_libraries(test_graph_load kokkos ${PROTOBUF_LIBRARY} LIBVGIO ${HTS_LIBRARY})

  add_executable(test_pairing test_main.cpp test_pairing.cpp)
  add_dependencies(test_pairing LIBHTS SYMLNK)
  target_link_libraries(test_pairing kokkos ${PROTOBUF_LIBRARY} LIBVGIO ${HTS_LIBRARY})

  add_executable(test_bfs test_main.cpp test_bfs.cpp)
  add_dependencies(test_bfs LIBHTS SYMLNK)
  target_link_libraries(test_bfs kokkos ${PROTOBUF_LIBRARY} LIBVGIO ${HTS_LIBRARY})

endif(BUILD_TESTS)
