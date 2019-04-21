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
    std::cout << "INFO, pairg::main, Time to build adjacency matrices (ms): " << T1.elapsed() << "\n";
    pairg::matrixOps::printMatrix(adj_mat, 1);

    pairg::timer T2;
    pairg::matrixOps::crsMat_t valid_pairs_mat = pairg::buildValidPairsMatrix(adj_mat, parameters); 
    std::cout << "INFO, pairg::main, Time to build result matrix (ms): " << T2.elapsed() << "\n";
    pairg::matrixOps::printMatrix(valid_pairs_mat, 1);
  }

  Kokkos::finalize();
}
