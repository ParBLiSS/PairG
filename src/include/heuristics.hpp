/**
 * @file    heuristics.hpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef PAIRG_HEURISTICS_HPP
#define PAIRG_HEURISTICS_HPP

#include "spgemm_utility.hpp"

namespace pairg
{
  /**
   * @brief             implement a simple heuristic to check distance constraints
   * @param[in] A       adjacency matrix (CSR formatted)
   * @param[in] p       input parameters (pull distance constraints from here)
   * @param[in] src     source vertex
   * @param[in] target  target vertex
   * @details           run BFS from source vertex up to d_2 levels and see if we reach
   *                    the target vertex
   * @return            boolean value (true if reached, false if not)
   */
  bool queryReachabilityBFS(const matrixOps::crsMat_t &A, const Parameters &p, matrixOps::lno_t src, matrixOps::lno_t target)
  {
    if (src >= A.numRows() || target >= A.numCols()) {
      std::cout << "WARNING, pairg::matrixOps::queryValue, query index out of range" << std::endl;
      return false;
    }

    assert(p.d_up > 0);

    if (src == target)
      return true;

    //Count of bfs levels 
    int level = 0;

    //dummy value to track BFS levels
    matrixOps::lno_t dummy = std::numeric_limits<matrixOps::lno_t>::max();

    //Mark all the vertices as not visited 
    std::vector<bool> visited(A.numRows(), false);  

    //queue for putting visited vertices
    std::list<matrixOps::lno_t> Q; 
    Q.push_back(src);  visited[src] = true;
    Q.push_back(dummy);   //to travel the levels

    while(!Q.empty() && level <= p.d_up) 
    {
      //pick vertex from queue
      auto cur = Q.front();
      Q.pop_front();

      if (cur == dummy)
      {
        //moving to next level
        level += 1;
        Q.push_back(dummy);

        if (Q.front() == dummy)
          break; //all vertices visited
      }
      else
      {
        //check if we hit the target vertex
        if (cur == target)
          return true;

        //Explore its neighbors
        for(auto i = A.graph.row_map(cur); i < A.graph.row_map(cur+1); i++)
        {
          auto neighbor = A.graph.entries(i);
          if (!visited[neighbor])
          {
            visited[neighbor] = true;
            Q.push_back(neighbor);
          }
        }
      }
    }
 
    return false;
  }
}

#endif
