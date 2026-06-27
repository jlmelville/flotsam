#ifndef FLOTSAM_LTSA_INTERNAL_H
#define FLOTSAM_LTSA_INTERNAL_H

#include <algorithm>
#include <cmath>
#include <cpp11.hpp>
#include <cpp11/doubles.hpp>
#include <cpp11/integers.hpp>
#include <cpp11/matrix.hpp>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <vector>

#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#include "pforr.h"

struct CompactEntry {
  int row;
  double value;
};

struct SparseComponents {
  std::vector<int> i;
  std::vector<int> p;
  std::vector<double> x;
};

struct LocalWeights {
  std::vector<double> weights;
  int rank = 0;
};

struct GramLocalWeightsWorkspace {
  GramLocalWeightsWorkspace(std::size_t n_nbrs, std::size_t n_dim, int ndim,
                            bool use_row_major);

  std::size_t n_nbrs_size;
  std::size_t n_dim_size;
  int n_nbrs;
  int n_dim;
  int requested;
  std::vector<int> nni;
  std::vector<double> centered;
  std::vector<double> row_buffer;
  std::vector<double> col_means;
  std::vector<double> gram;
  std::vector<double> values;
  std::vector<double> work;
  std::vector<double> weights;
  std::vector<int> keep;
};

std::size_t checked_triplet_count(std::size_t n_obs, std::size_t n_nbrs,
                                  const char* name);

int checked_neighbor_index(int idx, std::size_t n_obs);

void checked_append_output(int row, double value, std::vector<int>& out_i,
                           std::vector<double>& out_x, std::size_t max_int);

void checked_ndim(int ndim);

int checked_lapack_dim(std::size_t value, const char* name);

std::size_t checked_row_major_copy_max_bytes(double max_bytes);

std::vector<int> flat_neighbors_zero_based(const cpp11::integers& value_nnt,
                                           std::size_t offset,
                                           std::size_t n_nbrs);

void fill_flat_neighbors_zero_based(const cpp11::integers& value_nnt,
                                    std::size_t offset, std::size_t n_nbrs,
                                    std::vector<int>& out);

void fill_centered_neighborhood_ptr(const double* x_data, std::size_t n_obs,
                                    const std::vector<int>& nni,
                                    std::vector<double>& centered,
                                    std::size_t n_dim);

void fill_centered_neighborhood_row_major(const std::vector<double>& row_major,
                                          const std::vector<int>& nni,
                                          std::vector<double>& row_buffer,
                                          std::vector<double>& col_means,
                                          std::vector<double>& centered,
                                          std::size_t n_dim);

void fill_weights_from_basis(std::size_t n_nbrs, const std::vector<int>& keep,
                             const std::vector<double>& basis,
                             std::vector<double>& weights);

int select_local_basis_columns(const std::vector<double>& values, int n_values,
                               int n_nbrs, int n_dim, int requested,
                               bool values_ascending,
                               std::vector<int>& keep);

int query_dsyev_workspace(int n, std::vector<double>& gram,
                          std::vector<double>& values);

int query_dgesdd_workspace(int n_nbrs, int n_dim, int min_dim,
                           std::vector<double>& a, std::vector<double>& d,
                           std::vector<double>& u, std::vector<double>& vt,
                           std::vector<int>& iwork);

bool row_major_copy_within_limit(std::size_t n_obs, std::size_t n_dim,
                                 std::size_t max_bytes);

void make_row_major_copy(const double* x_data, std::size_t n_obs,
                         std::size_t n_dim, std::vector<double>& row_major);

int compute_local_weights_gram_workspace(const double* x_data,
                                         std::size_t n_obs,
                                         GramLocalWeightsWorkspace& workspace,
                                         const std::vector<double>* row_major);

LocalWeights compute_local_weights_shape_routed(const cpp11::doubles_matrix<>& x,
                                                const std::vector<int>& nni,
                                                int ndim);

std::size_t checked_size_add(std::size_t lhs, std::size_t rhs,
                             const char* message);

std::size_t checked_size_mul(std::size_t lhs, std::size_t rhs,
                             const char* message);

std::size_t checked_raw_staging_bytes(std::size_t canonical_count,
                                      std::size_t full_count);

std::size_t triangular_pair_count(std::size_t n_nbrs);

std::size_t triangular_pair_offset(std::size_t local_col,
                                   std::size_t local_row);

std::string row_major_fallback_reason(bool use_gram_workspace,
                                      bool row_major_used,
                                      bool row_major_within_limit);

class LtsaTripletAssemblyBuilder {
public:
  LtsaTripletAssemblyBuilder(const cpp11::integers& value_nnt,
                             std::size_t value_n_nbrs, std::size_t n_obs,
                             std::size_t max_int);

  void append_prechecked(const std::vector<int>& nni,
                         const std::vector<double>& weights);

  SparseComponents finalize_components();

  std::size_t raw_entries_estimate() const;

  std::size_t duplicate_fallback_count() const;

private:
  std::size_t n_obs_;
  std::size_t value_n_nbrs_;
  std::size_t max_int_;
  std::size_t n_appended_ = 0;
  std::size_t raw_entries_estimate_ = 0;
  std::size_t duplicate_fallback_count_ = 0;
  bool finalized_ = false;
  std::vector<std::vector<CompactEntry>> canonical_columns_;
  std::vector<std::vector<CompactEntry>> full_columns_;
  std::vector<int> append_seen_;

  bool has_duplicate_neighbors(const std::vector<int>& nni);

  void append_full_prechecked(const std::vector<int>& nni,
                              const std::vector<double>& weights);

  void append_triangular_prechecked(const std::vector<int>& nni,
                                    const std::vector<double>& weights);

  void expand_canonical_to_full(std::vector<double>& row_sums,
                                std::vector<int>& row_seen,
                                std::vector<int>& touched_rows);
};

#endif
