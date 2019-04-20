/**
 * @file    spgemm_utility.cpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

namespace pairg
{
  /**
   * @brief       supports loading of sequence graphs from VG or .txt file format
   */
  class graphLoader
  {
    public:

      //count of vertices (each labeled by a character)
      //equivalent to #rows in adjacency matrix
      int numVertices;

      //count of directed edges
      //equivalent to nnz elements in adjacency matrix
      int numEdges;

      //out-neighbors of each vertex saved in CSR form
      std::vector<int> adjcny;    //size = numEdges
      std::vector<int> offsets;   //size = numVertices + 1
        
      /**
       * @brief                 load graph from VG graph format
       * @param[in]   filename
       * @details               VG (github.com/vgteam/vg) uses .vg format to save graphs
       */
      void loadFromVG(const std::string &filename)
      {

      }

      /**
       * @brief                 load graph from .txt file
       * @param[in]   filename
       * @details               file format: 
       *                        first line specifies count of vertices
       *                        following lines specify out-neighbors and label, for each vertex per line
       *                        values in each row delimited by spaces
       */
      void loadFromTxt(const std::string &filename)
      {

      }
  };
}
