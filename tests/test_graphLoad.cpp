/**
 * @file    test_graphLoad.hpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#include "reachability.hpp"

//External includes
#include "catch/single_include/catch2/catch.hpp"

#define QUOTE(name) #name
#define STR(macro) QUOTE(macro)
#define FOLDER STR(PROJECT_TEST_DATA_DIR)

/**
 * @brief   builds a graph from BRCA1 sequence
 *          This routine checks for the correctness
 *          of graph loading (format = .txt)
 **/
TEST_CASE(".txt formatted graph can be loaded") 
{
  //get file name
  std::string file = FOLDER;
  file = file + "/chain.txt";

  //load graph

  psgl::graphLoader g;
  g.loadFromTxt(file);
  auto &graph = g.diCharGraph;

  SECTION( "evaluating graph size" ) {
    REQUIRE(graph.numVertices == 81189); 
    REQUIRE(graph.numEdges == 81188); 
  }
}

