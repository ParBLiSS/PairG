/**
 * @file    main.hpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#include "spgemm_utility.hpp"
#include "utility.hpp"

//External includes
#include "clipp/include/clipp.h"


/**
 * @brief     main function
 */
int main(int argc, char* argv[]) 
{
  Kokkos::initialize();

  //parse command-line arguments
  int dim, minNNZ, maxNNZ;

  auto cli = 
    (
     clipp::required("-d") & clipp::value("dim", dim).doc("dimension of square matrix"),
     clipp::required("-l") & clipp::value("minNNZ", minNNZ).doc("minimum non-zero vals in a row"),
     clipp::required("-h") & clipp::value("maxNNZ", maxNNZ).doc("maximum non-zero vals in a row")
    );

  if(!clipp::parse(argc, argv, cli)) 
  {
    //print help page
    clipp::operator<<(std::cout, clipp::make_man_page(cli, argv[0])) << std::endl;
    exit(1);
  }

  {

    pairg::timer T1;

    // Create and fill matrices A and B to multiply
    pairg::matrixOps::crsMat_t A = pairg::matrixOps::randomMatrix (dim, minNNZ, maxNNZ, false);
    pairg::matrixOps::crsMat_t B = pairg::matrixOps::randomMatrix (dim, minNNZ, maxNNZ, false);

    pairg::matrixOps::printMatrix(A, 1);
    pairg::matrixOps::printMatrix(B, 1);

    double constructionTime = T1.elapsed();
    std::cout << "Time to build sample matrices (ms): " << constructionTime << "\n";

    pairg::timer T2;

    pairg::matrixOps::crsMat_t C = pairg::matrixOps::multiplyMatrices (A,B);

    double multiplicationTime = T2.elapsed();
    std::cout << "Time to multiply sample matrices (ms): " << multiplicationTime << "\n";

    pairg::matrixOps::printMatrix(C, 3);

  }

    Kokkos::finalize();
}
