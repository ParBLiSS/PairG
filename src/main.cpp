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
    pairg::matrixOps::crsMat_t A = pairg::matrixOps::createRandomMatrix (dim, minNNZ, maxNNZ, false);
    pairg::matrixOps::crsMat_t B = pairg::matrixOps::createRandomMatrix (dim, minNNZ, maxNNZ, false);
    //pairg::matrixOps::crsMat_t A = pairg::matrixOps::createIdentityMatrix (dim);
    //pairg::matrixOps::crsMat_t B = pairg::matrixOps::createIdentityMatrix (dim);
    std::cout << "Time to build sample matrices (ms): " << T1.elapsed() << "\n";

    pairg::matrixOps::printMatrix(A, 3);
    pairg::matrixOps::printMatrix(B, 3);

    pairg::timer T2;
    pairg::matrixOps::crsMat_t C = pairg::matrixOps::multiplyMatrices (A,B);
    std::cout << "Time to multiply sample matrices (ms): " << T2.elapsed() << "\n";

    pairg::matrixOps::printMatrix(C, 3);

    pairg::timer T3;
    pairg::matrixOps::crsMat_t D = pairg::matrixOps::addMatrices (A,B);
    std::cout << "Time to add sample matrices (ms): " << T3.elapsed() << "\n";

    pairg::matrixOps::printMatrix(D, 3);

    pairg::timer T4;
    pairg::matrixOps::crsMat_t E = pairg::matrixOps::power (A, 10);
    std::cout << "Time to raise power of sample matrix (ms): " << T4.elapsed() << "\n";

    pairg::matrixOps::printMatrix(E, 3);
  }

    Kokkos::finalize();
}
