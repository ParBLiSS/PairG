/**
 * @file    spgemm_utility.cpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */


//External includes
#include "Kokkos_Core.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_spgemm.hpp"

namespace pairg
{
  template <typename crsMat_t>
    crsMat_t addMatrices(crsMat_t A, crsMat_t B)
    {
    }

  template <typename crsMat_t>
    crsMat_t multiplyMatrices(crsMat_t A, crsMat_t B)
    {
    }

  template <typename crsMat_t>
    crsMat_t power(crsMat_t A, int n)
    {
    }

  template <typename crsMat_t, typename ordinal_type>
    crsMat_t queryValue(crsMat_t A, ordinal_type row, ordinal_type col)
    {
    }

  /**
   * @brief                       create a random square matrix for testing, using kokkos-kernels
   * @tparam[in] crsMat_t         CRS matrix type
   * @tparam[in] ordinal_type     ordinal type (type for storing matrix indices)
   * @param[in] nrows             count of rows
   * @param[in] minNNZ            minimum non-zero elements per row
   * @param[in] maxNNZ            maximum non-zero elements per row
   * @param[in] sortRows          sort indices within each row 
   * @return                      the generated matrix
   * @details                     - value of each non-zero element is set to 1
   *                              - not suitable for very large matrices as the randomization procedure is expensive
   *                              - modified from kokkos-kernels repo: unit_test/sparse/Test_Sparse_spadd.hpp
   */
  template <typename crsMat_t, typename ordinal_type>
    crsMat_t randomMatrix(ordinal_type nrows, ordinal_type minNNZ, ordinal_type maxNNZ, bool sortRows)
    {
      typedef typename crsMat_t::StaticCrsGraphType graph_t;
      typedef typename graph_t::row_map_type::non_const_type size_type_view_t;
      typedef typename graph_t::entries_type::non_const_type lno_view_t;
      typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
      typedef typename size_type_view_t::non_const_value_type size_type;  //rowptr type
      typedef typename lno_view_t::non_const_value_type lno_t;            //colind type

      static_assert(std::is_same<ordinal_type, lno_t>::value, "ordinal_type should be same as lno_t from crsMat_t");

      //first, populate rowmap
      size_type_view_t rowmap("rowmap", nrows + 1);
      typename size_type_view_t::HostMirror h_rowmap = Kokkos::create_mirror_view(rowmap);
      size_type nnz = 0;
      for(lno_t i = 0; i < nrows; i++)
      {
        size_type rowEntries = rand() % (maxNNZ - minNNZ + 1) + minNNZ;
        h_rowmap(i) = nnz;
        nnz += rowEntries;
      }
      h_rowmap(nrows) = nnz;
      Kokkos::deep_copy(rowmap, h_rowmap);

      //allocate values and entries
      scalar_view_t values("values", nnz);
      //populate values
      typename scalar_view_t::HostMirror h_values = Kokkos::create_mirror_view(values);
      for(size_type i = 0; i < nnz; i++)
      {
        h_values(i) = 1;
      }
      Kokkos::deep_copy(values, h_values);

      //populate entries (make sure no entry is repeated within a row)
      lno_view_t entries("entries", nnz);
      typename lno_view_t::HostMirror h_entries = Kokkos::create_mirror_view(entries);
      std::vector<lno_t> indices(nrows);
      for(lno_t i = 0; i < nrows; i++)
      {
        for(lno_t j = 0; j < nrows; j++)
        {
          indices[j] = j;
        }
        std::random_shuffle(indices.begin(), indices.end());
        size_type rowStart = h_rowmap(i);
        size_type rowCount = h_rowmap(i + 1) - rowStart;
        if(sortRows)
        {
          std::sort(indices.begin(), indices.begin() + rowCount);
        }
        for(size_type j = 0; j < rowCount; j++)
        {
          h_entries(rowStart + j) = indices[j];
        }
      }
      Kokkos::deep_copy(entries, h_entries);

      return crsMat_t("test matrix", nrows, nrows, nnz, values, rowmap, entries);
    }
}
