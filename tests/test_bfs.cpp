/**
 * @file    test_bfs.hpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#include "reachability.hpp"
#include "heuristics.hpp"

//External includes
#include "catch/single_include/catch2/catch.hpp"

#define QUOTE(name) #name
#define STR(macro) QUOTE(macro)
#define FOLDER STR(PROJECT_TEST_DATA_DIR)

TEST_CASE("building valid-pair matrix for a chain graph") 
{
  Kokkos::initialize();

  //get file name
  std::string file = FOLDER;
  file = file + "/chain.txt";

  std::vector<char> RFILE(file.c_str(), file.c_str() + file.size() + 1u);

  SECTION( "distance limits are 0, 50" )
  {
    char *argv[] = {"pairmap2graph", "-m", "txt", "-r", RFILE.data(), "-l", "0", "-u", "50", "-t", "4", "-c", "0", nullptr};
    int argc = 13;

    pairg::Parameters parameters;        
    pairg::parseandSave(argc, argv, parameters);

    int V = 81189;
    int E = 81188;

    pairg::matrixOps::crsMat_t A = pairg::getAdjacencyMatrix(parameters);

    int NNZ = E;

    SECTION( "evaluating adjacency matrix size" ) {
      REQUIRE(A.numRows() == V);  
      REQUIRE(A.numCols() == V); 
      REQUIRE(A.graph.row_map.extent(0) == V + 1); 
      REQUIRE(A.graph.entries.extent(0) == NNZ); 
      REQUIRE(A.values.extent(0) == NNZ); 
    }

    SECTION( "evaluating matrix content" ) {
      REQUIRE(std::accumulate(A.values.data(), A.values.data() + A.values.extent(0), 0) == NNZ);
      REQUIRE(std::accumulate(A.values.data(), A.values.data() + A.values.extent(0), 1, std::multiplies<int>()) == 1);
    }

    SECTION( "checking whether queries are answered correctly" ) {
      //Remember that vertex ids are 0-based while querying

      REQUIRE(pairg::queryReachabilityBFS (A, parameters, 0, 0) == true);
      REQUIRE(pairg::queryReachabilityBFS (A, parameters, 0, 1) == true);
      REQUIRE(pairg::queryReachabilityBFS (A, parameters, 1, 0) == false);
      REQUIRE(pairg::queryReachabilityBFS (A, parameters, 0, 50) == true);
      REQUIRE(pairg::queryReachabilityBFS (A, parameters, 0, 51) == false);

      REQUIRE(pairg::queryReachabilityBFS (A, parameters, 81137, 81188) == false);
      REQUIRE(pairg::queryReachabilityBFS (A, parameters, 81138, 81188) == true);
    }
  }

  Kokkos::finalize();
}
