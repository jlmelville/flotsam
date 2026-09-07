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

#include "ltsa_parallel_reduce.h"

struct SparseComponents {
  std::vector<int> i;
  std::vector<int> p;
  std::vector<double> x;
};

enum LocalWeightsComputationStatus {
  LOCAL_WEIGHTS_COMPUTATION_OK = 0,
  LOCAL_WEIGHTS_COMPUTATION_NONFINITE = 1
};

constexpr std::size_t ROW_MAJOR_COPY_LIMIT_BYTES =
    static_cast<std::size_t>(256) * 1024 * 1024;

template <typename T>
std::size_t checked_vector_size(std::size_t count, const char *name) {
  if (count > std::vector<T>().max_size()) {
    cpp11::stop("%s is too large", name);
  }
  return count;
}

template <typename T>
std::size_t checked_vector_size_mul(std::size_t lhs, std::size_t rhs,
                                    const char *name) {
  const std::size_t max_count = std::vector<T>().max_size();
  if (lhs != 0 && rhs > max_count / lhs) {
    cpp11::stop("%s is too large", name);
  }
  return lhs * rhs;
}

struct SvdLocalWeightsWorkspace {
  SvdLocalWeightsWorkspace(std::size_t n_nbrs, std::size_t n_features,
                           int ndim);

  std::size_t n_nbrs_size;
  std::size_t n_features_size;
  int n_nbrs;
  int n_features;
  int min_dim;
  int requested_basis_size;
  std::vector<int> neighbor_indices;
  std::vector<int> basis_columns;
  std::vector<double> centered;
  std::vector<double> a;
  std::vector<double> d;
  std::vector<double> u;
  std::vector<double> vt;
  std::vector<double> work;
  std::vector<int> iwork;
  std::vector<double> weights;
};

struct GramLocalWeightsWorkspace {
  GramLocalWeightsWorkspace(std::size_t n_nbrs, std::size_t n_features,
                            int ndim, bool use_row_major);

  std::size_t n_nbrs_size;
  std::size_t n_features_size;
  int n_nbrs;
  int n_features;
  int requested_basis_size;
  std::vector<int> neighbor_indices;
  std::vector<double> centered;
  std::vector<double> row_buffer;
  std::vector<double> gram;
  std::vector<double> values;
  std::vector<double> work;
  std::vector<double> weights;
  std::vector<int> basis_columns;
};

std::size_t checked_triplet_count(std::size_t n_obs, std::size_t n_nbrs,
                                  const char *name);

int checked_zero_based_neighbor_index(int idx, std::size_t n_obs);

void checked_append_output(int row, double value, std::vector<int> &out_i,
                           std::vector<double> &out_x, std::size_t max_int);

void checked_ndim(int ndim);

int checked_lapack_dim(std::size_t value, const char *name);

void fill_flat_neighbors_zero_based(
    const cpp11::integers &transposed_neighbor_indices, std::size_t offset,
    std::size_t n_nbrs, std::vector<int> &out);

bool fill_centered_neighborhood_column_major(
    const double *x_data, std::size_t n_obs,
    const std::vector<int> &neighbor_indices, std::vector<double> &centered,
    std::size_t n_features);

bool fill_centered_neighborhood_row_major(
    const std::vector<double> &row_major,
    const std::vector<int> &neighbor_indices, std::vector<double> &row_buffer,
    std::vector<double> &centered, std::size_t n_features);

int clean_local_basis(std::size_t n_nbrs, std::vector<int> &basis_columns,
                      std::vector<double> &basis);

void fill_weights_from_basis(std::size_t n_nbrs,
                             const std::vector<int> &basis_columns,
                             const std::vector<double> &basis,
                             std::vector<double> &weights);

[[noreturn]] void stop_local_weights_computation(int status, int neighborhood);

int select_local_basis_columns(const std::vector<double> &values, int n_values,
                               int n_nbrs, int n_features,
                               int requested_basis_size, bool values_ascending,
                               std::vector<int> &basis_columns);

int query_dsyev_workspace(int n, std::vector<double> &gram,
                          std::vector<double> &values);

int query_dgesdd_workspace(int n_nbrs, int n_features, int min_dim,
                           std::vector<double> &a, std::vector<double> &d,
                           std::vector<double> &u, std::vector<double> &vt,
                           std::vector<int> &iwork);

bool row_major_copy_within_limit(std::size_t n_obs, std::size_t n_features,
                                 std::size_t max_bytes);

void make_row_major_copy(const double *x_data, std::size_t n_obs,
                         std::size_t n_features,
                         std::vector<double> &row_major);

int compute_local_weights_svd_workspace(const double *x_data, std::size_t n_obs,
                                        SvdLocalWeightsWorkspace &workspace,
                                        int &rank, int &computation_status);

int compute_local_weights_gram_workspace(const double *x_data,
                                         std::size_t n_obs,
                                         GramLocalWeightsWorkspace &workspace,
                                         const std::vector<double> *row_major,
                                         int &rank, int &computation_status);

std::size_t checked_size_add(std::size_t lhs, std::size_t rhs,
                             const char *message);

std::size_t checked_size_mul(std::size_t lhs, std::size_t rhs,
                             const char *message);

std::size_t triangular_pair_count(std::size_t n_nbrs);

std::size_t triangular_pair_offset(std::size_t local_col,
                                   std::size_t local_row);

class TripletAssemblyBuilder {
public:
  TripletAssemblyBuilder(const cpp11::integers &transposed_neighbor_indices,
                         std::size_t n_neighbors, std::size_t n_obs,
                         std::size_t max_int);

  void append_prechecked_neighborhood(const std::vector<int> &neighbor_indices,
                                      const std::vector<double> &weights);

  SparseComponents finalize_sparse_components();

private:
  std::size_t n_obs_;
  std::size_t n_neighbors_;
  std::size_t max_int_;
  std::size_t n_appended_ = 0;
  bool finalized_ = false;
  std::vector<std::vector<CompactEntry>> canonical_columns_;
  std::vector<std::vector<CompactEntry>> full_columns_;

  void append_triangular_prechecked(const std::vector<int> &neighbor_indices,
                                    const std::vector<double> &weights);

  void expand_canonical_to_full(std::vector<double> &row_sums,
                                std::vector<int> &row_seen,
                                std::vector<int> &touched_rows);
};

#endif
