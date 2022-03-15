#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <vector>
#include <cmath>
#include <algorithm>

#include "StrumpackSparseSolverMPIDist.hpp"
#include "sparse/CSRMatrixMPI.hpp"
#include "misc/RandomWrapper.hpp"

// const std::string mat = "/project/projectdirs/m77/nstx_large_matrices/nstxu_180degree_antenna_phasing_write_matrix_order3/matrix.";

using real_t = double;
using scalar_t = std::complex<real_t>;


int main(int argc, char* argv[]) {
  std::string mat = std::string(argv[1]);

  int thread_level;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &thread_level);
  strumpack::MPIComm comm(MPI_COMM_WORLD);
  int rank = comm.rank(), P = comm.size();

  double t_start = MPI_Wtime();
  std::int64_t m = 0;
  using Trip_t = strumpack::Triplet<scalar_t,std::int64_t>;
  std::vector<Trip_t> triplets;
  for (int i=rank; i<256; i+=P) {
    std::stringstream ss;
    ss << mat << std::setfill('0') << std::setw(6) << i;
    std::string fname = ss.str();
    std::cout << "Rank " << rank << " opening " << fname << std::endl;
    const int max_cline = 256;
    char cline[max_cline];
    FILE *fp = fopen(fname.c_str(), "r");
    if (!fp) {
      std::cout << "Unable to open file" << std::endl;
      return 1;
    }
    std::int64_t row, col;
    double real, imag;
    while (fscanf(fp, "%lld %lld %lf %lf\n", &row, &col, &real, &imag) != EOF) {
      triplets.emplace_back(row, col, scalar_t(real, imag));
      m = std::max(m, row+1);
    }
    fclose(fp);
  }
  std::int64_t nnz = comm.all_reduce(triplets.size(), MPI_SUM);
  m = comm.all_reduce(m, MPI_MAX);
  if (!rank) {
    std::cout << "Matrix has " << m << " rows and "
              << nnz << " nonzeros" << std::endl;
    std::cout << "Reading matrix took "
	      << MPI_Wtime() - t_start << " seconds" << std::endl;
  }
  t_start = MPI_Wtime();

  // assign ~m/P rows to each proc
  std::vector<std::int64_t> dist(P+1);
  for (int p=1; p<P; p++)
    dist[p] = std::floor(p * double(m) / P);
  dist[P] = m;
  std::int64_t lrows = dist[rank+1] - dist[rank];

  if (!rank) {
    std::cout << "local rows {";
    for (std::size_t p=0; p<P; p++)
      std::cout << lrows << ", ";
    std::cout << "}" << std::endl;
  }

  // redistribute the triplets according to dist
  std::vector<std::vector<Trip_t>> sbuf(P);
  for (auto& t : triplets) {
    auto dst = std::upper_bound(dist.begin(), dist.end(), t.r)
      - dist.begin() - 1;
    assert(dst >= 0);
    assert(dst < P);
    sbuf[dst].push_back(t);
  }
  triplets.clear();
  if (!rank) std::cout << "redistributing triplets" << std::endl;
  auto ltrips = comm.all_to_all_v(sbuf);
  sbuf.clear();
  std::sort(ltrips.begin(), ltrips.end(),
            [](const Trip_t& a, const Trip_t& b) {
              return (a.r < b.r) || (a.r == b.r && a.c < b.c); });

  std::int64_t lnnz = ltrips.size();
  if (!rank) {
    std::cout << "local nnz {";
    for (std::size_t p=0; p<P; p++)
      std::cout << lnnz << ", ";
    std::cout << "}" << std::endl;
    std::cout << "Redistributing matrix took "
	      << MPI_Wtime() - t_start << " seconds" << std::endl;
  }
  t_start = MPI_Wtime();
  std::vector<std::int64_t> rptr(lrows+1), cind(lnnz);
  std::vector<scalar_t> val(lnnz);
  if (!rank) std::cout << "building local CSR" << std::endl;
  for (std::size_t i=0; i<ltrips.size(); i++) {
    auto& t = ltrips[i];
    rptr[t.r-dist[rank]+1]++;
    cind[i] = t.c;
    val[i] = t.v;
  }
  for (std::size_t i=1; i<lrows; i++)
    rptr[i+1] += rptr[i];
  assert(rptr[lrows] == lnnz);
  assert(rptr[0] == 0);
  if (!rank) std::cout << "creating CSRMatrixMPI" << std::endl;
  strumpack::CSRMatrixMPI<scalar_t,std::int64_t>
    A(lrows, rptr.data(), cind.data(), val.data(), dist.data(), comm, false);
  val.clear();
  rptr.clear();
  cind.clear();
  if (!rank) {
    std::cout << "Creating CSRMatrixMPI took "
	      << MPI_Wtime() - t_start << " seconds" << std::endl;
  }

  if (!rank) std::cout << "creating STRUMPACK solver" << std::endl;
  strumpack::StrumpackSparseSolverMPIDist<scalar_t,std::int64_t>sp(MPI_COMM_WORLD);
  sp.options().set_from_command_line(argc, argv);
  sp.options().set_matching(strumpack::MatchingJob::NONE);

  std::vector<scalar_t> b(lrows), x(lrows), x_exact(lrows);
  // construct random exact solution
  auto rgen = strumpack::random::make_default_random_generator<real_t>();
  for (auto& xi : x_exact)
    xi = scalar_t(rgen->get());
  A.spmv(x_exact.data(), b.data());

  sp.set_matrix(A);

  if (sp.reorder() != strumpack::ReturnCode::SUCCESS) {
    if (!rank)
      std::cout << "problem with reordering of the matrix." << std::endl;
    return 1;
  }
  if (sp.factor() != strumpack::ReturnCode::SUCCESS) {
    if (!rank)
      std::cout << "problem during factorization of the matrix."
                << std::endl;
    return 1;
  }
  sp.solve(b.data(), x.data());

  auto scaled_res = A.max_scaled_residual(x.data(), b.data());
  if (!rank)
    std::cout << "# COMPONENTWISE SCALED RESIDUAL = "
              << scaled_res << std::endl;
  strumpack::blas::axpy(lrows, scalar_t(-1.), x_exact.data(), 1, x.data(), 1);
  auto nrm_error = strumpack::norm2(x, comm);
  auto nrm_x_exact = strumpack::norm2(x_exact, comm);
  if (!rank)
    std::cout << "# RELATIVE ERROR = " << (nrm_error/nrm_x_exact)
              << std::endl;

  return 0;
}
