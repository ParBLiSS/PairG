/**
 * @file    spgemm_utility.cpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef PAIR_REACHABILITY_HPP
#define PAIR_REACHABILITY_HPP

#include "spgemm_utility.hpp"
#include "parseCmdArgs.hpp"

//External includes
#include "PaSGAL/graphLoad.hpp"

namespace pairg
{
  /**
   * @brief     build adjacency matrix from variaton graph
   */
  matrixOps::crsMat_t getAdjacencyMatrix(const Parameters &parameters) 
  {
    psgl::graphLoader g;
    {
      if (parameters.gmode.compare("vg") == 0)
        g.loadFromVG(parameters.graphfile);
      else if(parameters.gmode.compare("txt") == 0)
        g.loadFromTxt(parameters.graphfile);
      else 
      {
        std::cerr << "Invalid graph format " << parameters.gmode << std::endl;
        exit(1);
      }
    }

    //Use g.diCharGraph to build crsMat_t matrix
    typename matrixOps::lno_t nrows = g.diCharGraph.numVertices;
    typename matrixOps::size_type nnz = g.diCharGraph.numEdges;

    typename matrixOps::lno_nnz_view_t entries("entries", nnz);
    typename matrixOps::scalar_view_t values("values", nnz);
    typename matrixOps::lno_view_t rowmap("rowmap", nrows + 1);

    for(matrixOps::size_type i = 0; i < nnz; i++) {
      values(i) = 1;
    }

    for(matrixOps::size_type i = 0; i < nnz; i++) {
      entries(i) = g.diCharGraph.adjcny_out[i];
    }

    for(matrixOps::lno_t i = 0; i < nrows + 1; i++) {
      rowmap(i) = g.diCharGraph.offsets_out[i];
    }

    return matrixOps::crsMat_t("adjacency matrix", nrows, nrows, nnz, values, rowmap, entries);
  }

  /**
   * @brief           build matrix associated with valid vertices that satisfy 
   *                  distance constraints
   * @param[in] A     graph adjacency matrix
   * @return          validity matrix
   *                  cell (i,j) = 1 iff there is a valid path from v_i to v_j
   */
  matrixOps::crsMat_t buildValidPairsMatrix(const matrixOps::crsMat_t &A, const Parameters &p)
  {
    pairg::timer T1;
    matrixOps::crsMat_t B = matrixOps::addMatrices(A, matrixOps::createIdentityMatrix(A.numRows()));
    std::cout << "INFO, pairg::buildValidPairsMatrix, time to add identity matrix (ms): " << T1.elapsed() << "\n";

    pairg::timer T2;
    matrixOps::crsMat_t C = matrixOps::power(A, p.d_low);
    std::cout << "INFO, pairg::buildValidPairsMatrix, time to raise adjacency matrix (ms): " << T2.elapsed() << "\n";

    pairg::timer T3;
    matrixOps::crsMat_t D = matrixOps::power (B, p.d_up - p.d_low);
    std::cout << "INFO, pairg::buildValidPairsMatrix, time to raise adjacency+identity matrix (ms): " << T3.elapsed() << "\n";

    pairg::timer T4;
    matrixOps::crsMat_t E = matrixOps::multiplyMatrices(C,D); 
    std::cout << "INFO, pairg::buildValidPairsMatrix, time to execute final multiplication (ms): " << T4.elapsed() << "\n";

    //sort entries within each row
    pairg::timer T5;
    matrixOps::indexForQuery(E); 
    std::cout << "INFO, pairg::buildValidPairsMatrix, time to index for querying (ms): " << T5.elapsed() << "\n";

    return E;
  }
}

#endif
