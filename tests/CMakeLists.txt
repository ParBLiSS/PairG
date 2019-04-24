##### add definitions to find path to test data files
add_definitions(-DPROJECT_TEST_DATA_DIR=${PROJECT_SOURCE_DIR}/tests/data)

if (BUILD_TESTS)

  add_executable(test_graph_load test_main.cpp test_graphLoad.cpp)
  target_link_libraries(test_graph_load kokkos vg_protobuf ${PROTOBUF_LIBRARY})

endif(BUILD_TESTS)
