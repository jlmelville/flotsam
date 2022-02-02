#include <algorithm>
#include <cpp11.hpp>
#include <cpp11/doubles.hpp>
#include <cpp11/integers.hpp>
#include <iostream>
#include <unordered_set>
#include <vector>
using namespace cpp11;

// convert the row/column indices of neighbor idx in ns to the 1D index in the
// sparse is vector, using ps to find the row range
[[cpp11::register]] integers
sparse_idxs(const integers& is, const integers& ps, const integers& ns) {
  // contents of ns are 1-indexed
  // contents of is and ps are 0-indexed

  std::size_t n_nbrs = ns.size();
  writable::integers sp_idx(n_nbrs * n_nbrs);
  auto mbegin = is.cbegin();

  // loop over every pair in n including e.g. (j, i) as well as (i, j)
  for (std::size_t i = 0; i < n_nbrs; i++) {
    auto r = ns[i]; // r is 1-indexed

    for (std::size_t j = 0; j < n_nbrs; j++) {
      auto c = ns[j]; // c is 1-indexed

      // cbegin/cend range of data at row r

      auto cbegin = ps[r - 1];
      auto cend = ps[r];
      auto pos = std::find(integers::const_iterator(mbegin + cbegin),
                           is.cend() + cend,
                           c - 1) -
                 mbegin + 1;
      sp_idx[i * n_nbrs + j] = pos;
    }
  }

  return sp_idx;
}

  // calculates M * d where M is a sparse matrix with pointer ps and values
  // xs, and d is the vector ds
  [[cpp11::register]] doubles
  spm_times_scalar(const integers& ps, const doubles& xs, const doubles& ds)
{
  auto nrow = ds.size();
  if (nrow != ps.size() - 1) {
    stop("Inconsistent diagonal and pointer lengths");
  }
  if (xs.size() != ps[nrow]) {
    stop("Inconsistent value and pointers");
  }
  writable::doubles dxs(xs.size());

  for (R_xlen_t i = 0; i < nrow; i++) {
    R_xlen_t begin = ps[i];
    R_xlen_t end = ps[i + 1];
    for (R_xlen_t j = begin; j < end; j++) {
      dxs[j] = xs[j] * ds[i];
    }
  }

  return dxs;
}

std::vector<int>
create_idx1d(const integers& nnt, std::size_t n_nbrs, std::size_t n_obs)
{
  std::unordered_set<int> all_idx;
  std::unordered_set<int> tmp;
  for (std::size_t i = 0; i < n_obs; i++) {
    auto inbrs = i * n_nbrs;
    for (std::size_t j = 0; j < n_nbrs; j++) {
      auto nbrj = nnt[inbrs + j] - 1;
      for (std::size_t k = j; k < n_nbrs; k++) {
        auto nbrk = nnt[inbrs + k] - 1;
        auto idx = nbrj < nbrk ? nbrj * n_obs + nbrk : nbrk * n_obs + nbrj;
        tmp.insert(idx);
      }
    }
    if (tmp.size() > all_idx.size()) {
      all_idx.insert(tmp.begin(), tmp.end());
      tmp.clear();
    }
  }
  if (tmp.size() > 0) {
    all_idx.insert(tmp.begin(), tmp.end());
    tmp.clear();
  }

  // create sorted vector from set
  std::vector<int> ridx;
  ridx.reserve(all_idx.size());
  ridx.insert(ridx.end(), all_idx.begin(), all_idx.end());

  return ridx;
}

[[cpp11::register]] list
nbrhood_triplets(const integers& nnt, std::size_t n_nbrs) {
  auto n_obs = nnt.size() / n_nbrs;

  // find unique 1D indices (i < j only)
  std::vector<int> ridx = create_idx1d(nnt, n_nbrs, n_obs);

  std::sort(ridx.begin(), ridx.end());

  // convert to (i, j)
  std::vector<int> cidx(ridx.size());
  for (std::size_t i = 0; i < ridx.size(); i++) {
    auto i1d = ridx[i];
    ridx[i] = i1d % n_obs;
    cidx[i] = static_cast<int>(i1d / n_obs);
  }

  return writable::list({ "i"_nm = ridx, "j"_nm = cidx });
}
