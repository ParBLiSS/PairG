#include<algorithm>   //for std::random_shuffle
#include<cstdlib>     //for rand
#include<type_traits> //for std::is_same

#define KOKKOSKERNELS_IMPL_COMPILE_LIBRARY false

//External includes
#include<Kokkos_Core.hpp>
#include<KokkosSparse_CrsMatrix.hpp>
#include<KokkosSparse_spgemm.hpp>

/**
 * @brief     create a random matrix for testing
 *            borrowed from kokkos-kernels repo: unit_test/sparse/Test_Sparse_spadd.hpp
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
    h_values(i) = rand() / RAND_MAX;
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

/**
 * @brief     main function
 */
int main(int argc, char* argv[]) 
{
  Kokkos::initialize();

  typedef float scalar_t;
  typedef int lno_t;
  typedef Kokkos::OpenMP Device;

  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, Device> crsMat_t;

  // Create and fill matrices to multiply, input_mat and input_mat2
  crsMat_t input_mat = randomMatrix<crsMat_t, lno_t>(100, 1, 2, true);
  crsMat_t input_mat2 = randomMatrix<crsMat_t, lno_t>(100, 1, 2, true);

  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::entries_type::non_const_type   lno_nnz_view_t;
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;

  typedef typename lno_view_t::value_type size_type;
  typedef typename lno_nnz_view_t::value_type lno_t;
  typedef typename scalar_view_t::value_type scalar_t;

  // Create KokkosKernelHandle
  typedef KokkosKernels::Experimental::KokkosKernelsHandle
    <size_type, lno_t, scalar_t,
    typename Device::execution_space, typename Device::memory_space,typename Device::memory_space > KernelHandle;

  KernelHandle kh; 
  kh.set_team_work_size(16);
  kh.set_dynamic_scheduling(true);

  // Select an spgemm algorithm, limited by configuration at compile-time and set via the handle
  // Some options: {SPGEMM_KK_MEMORY, SPGEMM_KK_SPEED, SPGEMM_KK_MEMSPEED, /*SPGEMM_CUSPARSE, */ SPGEMM_MKL}
  KokkosSparse::SPGEMMAlgorithm spgemm_algorithm = KokkosSparse::SPGEMM_KK_MEMORY;
  kh.create_spgemm_handle(spgemm_algorithm);

  const size_t num_rows_1 = input_mat.numRows();
  const size_t num_rows_2 = input_mat2.numRows();
  const size_t num_cols_2 = input_mat2.numCols();

  const size_t num_cols_1 = input_mat.numCols();
  // Require num_rows_2 == num_cols_1 is true

  // Prepare resultant matrix
  lno_view_t row_mapC ("non_const_lnow_row", num_rows_1 + 1);
  lno_nnz_view_t  entriesC;
  scalar_view_t valuesC;


  KokkosSparse::Experimental::spgemm_symbolic (
      &kh,
      num_rows_1,
      num_rows_2,
      num_cols_2,
      input_mat.graph.row_map,
      input_mat.graph.entries,
      false,
      input_mat2.graph.row_map,
      input_mat2.graph.entries,
      false,
      row_mapC
      );

  size_t c_nnz_size = kh.get_spgemm_handle()->get_c_nnz();
  if (c_nnz_size){
    entriesC = lno_nnz_view_t (Kokkos::ViewAllocateWithoutInitializing("entriesC"), c_nnz_size);
    valuesC = scalar_view_t (Kokkos::ViewAllocateWithoutInitializing("valuesC"), c_nnz_size);
  }

  KokkosSparse::Experimental::spgemm_numeric(
      &kh,
      num_rows_1,
      num_rows_2,
      num_cols_2,
      input_mat.graph.row_map,
      input_mat.graph.entries,
      input_mat.values,
      false,
      input_mat2.graph.row_map,
      input_mat2.graph.entries,
      input_mat2.values,
      false,
      row_mapC,
      entriesC,
      valuesC
      );


  graph_t static_graph (entriesC, row_mapC);
  crsMat_t crsmat("CrsMatrix", num_cols_2, valuesC, static_graph);

  
  std::cout << "row_mapA:" << input_mat.graph.row_map.extent(0) << std::endl;
  std::cout << "entriesA:" << input_mat.graph.entries.extent(0) << std::endl;
  std::cout << "valuesA:" << input_mat.values.extent(0) << std::endl;
  KokkosKernels::Impl::print_1Dview(input_mat.graph.row_map, true);
  KokkosKernels::Impl::print_1Dview(input_mat.graph.entries, true);
  KokkosKernels::Impl::print_1Dview(input_mat.values, true);

  std::cout << "row_mapB:" << input_mat2.graph.row_map.extent(0) << std::endl;
  std::cout << "entriesB:" << input_mat2.graph.entries.extent(0) << std::endl;
  std::cout << "valuesB:" << input_mat2.values.extent(0) << std::endl;
  KokkosKernels::Impl::print_1Dview(input_mat2.graph.row_map, true);
  KokkosKernels::Impl::print_1Dview(input_mat2.graph.entries, true);
  KokkosKernels::Impl::print_1Dview(input_mat2.values, true);

  std::cout << "row_mapC:" << row_mapC.extent(0) << std::endl;
  std::cout << "entriesC:" << entriesC.extent(0) << std::endl;
  std::cout << "valuesC:" << valuesC.extent(0) << std::endl;
  KokkosKernels::Impl::print_1Dview(row_mapC, true);
  KokkosKernels::Impl::print_1Dview(entriesC, true);
  KokkosKernels::Impl::print_1Dview(valuesC, true);

  kh.destroy_spgemm_handle();

  Kokkos::finalize();
}
