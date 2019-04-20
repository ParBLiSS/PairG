#include<algorithm>   //for std::random_shuffle
#include<cstdlib>     //for rand
#include<type_traits> //for std::is_same
#include <chrono>

#define KOKKOSKERNELS_IMPL_COMPILE_LIBRARY false

//External includes
#include "clipp/include/clipp.h"
#include "Kokkos_Core.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_spgemm.hpp"



/**
 * @brief     main function
 */
int main(int argc, char* argv[]) 
{
  Kokkos::initialize();
  using timer = timer_impl<std::chrono::duration<double, std::milli> >;

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

  typedef int scalar_t;
  typedef int lno_t;
  typedef Kokkos::OpenMP Device;

  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, Device> crsMat_t;

  {

    timer T1;

    // Create and fill matrices to multiply, input_mat1 and input_mat2
    crsMat_t input_mat1 = randomMatrix<crsMat_t, lno_t> (dim, minNNZ, maxNNZ, false);
    crsMat_t input_mat2 = randomMatrix<crsMat_t, lno_t> (dim, minNNZ, maxNNZ, false);

    double constructionTime = T1.elapsed();
    std::cout << "Time to build sample matrices (ms): " << constructionTime << "\n";

    typedef typename crsMat_t::StaticCrsGraphType graph_t;
    typedef typename graph_t::row_map_type::non_const_type lno_view_t;
    typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
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

    const size_t num_rows_1 = input_mat1.numRows();
    const size_t num_cols_1 = input_mat1.numCols();
    const size_t num_rows_2 = input_mat2.numRows();
    const size_t num_cols_2 = input_mat2.numCols();

    // Require num_rows_2 == num_cols_1 is true
    timer T2;

    // Prepare resultant matrix
    lno_view_t row_mapC ("non_const_lnow_row", num_rows_1 + 1);
    lno_nnz_view_t  entriesC;
    scalar_view_t valuesC;

    KokkosSparse::Experimental::spgemm_symbolic (
        &kh,
        num_rows_1,
        num_rows_2,
        num_cols_2,
        input_mat1.graph.row_map,
        input_mat1.graph.entries,
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
        input_mat1.graph.row_map,
        input_mat1.graph.entries,
        input_mat1.values,
        false,
        input_mat2.graph.row_map,
        input_mat2.graph.entries,
        input_mat2.values,
        false,
        row_mapC,
        entriesC,
        valuesC
        );

    double multiplicationTime = T2.elapsed();
    std::cout << "Time to multiply sample matrices (ms): " << multiplicationTime << "\n";

    graph_t static_graph (entriesC, row_mapC);
    crsMat_t crsmat("CrsMatrix", num_cols_2, valuesC, static_graph);

    //Print out results
    {
      std::cout << "row_mapA:" << input_mat1.graph.row_map.extent(0) << std::endl;
      std::cout << "entriesA:" << input_mat1.graph.entries.extent(0) << std::endl;
      std::cout << "valuesA:" << input_mat1.values.extent(0) << std::endl;
      //KokkosKernels::Impl::print_1Dview(input_mat1.graph.row_map, true);
      //KokkosKernels::Impl::print_1Dview(input_mat1.graph.entries, true);
      //KokkosKernels::Impl::print_1Dview(input_mat1.values, true);

      std::cout << "row_mapB:" << input_mat2.graph.row_map.extent(0) << std::endl;
      std::cout << "entriesB:" << input_mat2.graph.entries.extent(0) << std::endl;
      std::cout << "valuesB:" << input_mat2.values.extent(0) << std::endl;
      //KokkosKernels::Impl::print_1Dview(input_mat2.graph.row_map, true);
      //KokkosKernels::Impl::print_1Dview(input_mat2.graph.entries, true);
      //KokkosKernels::Impl::print_1Dview(input_mat2.values, true);

      std::cout << "row_mapC:" << row_mapC.extent(0) << std::endl;
      std::cout << "entriesC:" << entriesC.extent(0) << std::endl;
      std::cout << "valuesC:" << valuesC.extent(0) << std::endl;
      //KokkosKernels::Impl::print_1Dview(row_mapC, true);
      //KokkosKernels::Impl::print_1Dview(entriesC, true);
      //KokkosKernels::Impl::print_1Dview(valuesC, true);
    }

    kh.destroy_spgemm_handle();
  }

  Kokkos::finalize();
}
