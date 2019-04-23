/**
 * @file    main.hpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#include "spgemm_utility.hpp"
#include "utility.hpp"
#include "parseCmdArgs.hpp"
#include "reachability.hpp"

//External includes
#include "clipp/include/clipp.h"

/**
 * @brief   get a random pair of integers, each value lies in [0, MAX)
 */
std::pair<int,int> getRandomPair(int MAX)
{
  int first = rand() % MAX;
  int second = rand() % MAX;

  return std::make_pair(first, second);
}

/**
 * @brief     main function
 */
int main(int argc, char* argv[]) 
{
  Kokkos::initialize();

  //parse command line arguments   
  pairg::Parameters parameters;        
  pairg::parseandSave(argc, argv, parameters);   

  {
    pairg::timer T1;

    pairg::matrixOps::crsMat_t adj_mat = pairg::getAdjacencyMatrix(parameters);
    std::cout << "INFO, pairg::main, Time to build adjacency matrix (ms): " << T1.elapsed() << "\n";
    pairg::matrixOps::printMatrix(adj_mat, 1);

    pairg::timer T2;
    pairg::matrixOps::crsMat_t valid_pairs_mat = pairg::buildValidPairsMatrix(adj_mat, parameters); 
    std::cout << "INFO, pairg::main, Time to build result matrix (ms): " << T2.elapsed() << "\n";
    pairg::matrixOps::printMatrix(valid_pairs_mat, 1);

    pairg::timer T3;
    int search_count = 1000;
    for(int i = 0; i < search_count; i++)
    {
      auto p = getRandomPair (valid_pairs_mat.numRows());
      bool exists = pairg::matrixOps::queryValue (valid_pairs_mat, p.first, p.second); 
    }
    std::cout << "INFO, pairg::main, Time to execute " << search_count << " queries (ms): " << T3.elapsed() << "\n";
  }

  Kokkos::finalize();
}
