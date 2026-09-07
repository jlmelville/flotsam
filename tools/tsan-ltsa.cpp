#include "ltsa_parallel_reduce.h"

#include <cstddef>
#include <iostream>
#include <vector>

namespace {

using flotsam_detail::ParallelReduceWorkspace;

bool entries_match(const std::vector<CompactEntry> &actual,
                   const std::vector<CompactEntry> &expected) {
  if (actual.size() != expected.size()) {
    return false;
  }
  for (std::size_t entry = 0; entry < actual.size(); entry++) {
    if (actual[entry].row != expected[entry].row ||
        actual[entry].value != expected[entry].value) {
      return false;
    }
  }
  return true;
}

bool reduction_matches(std::size_t n_threads) {
  const std::vector<std::size_t> column_starts{0, 4, 7, 12, 14};
  const std::vector<std::size_t> column_counts{4, 3, 5, 2};
  const std::vector<int> raw_rows{0, 0, 1, 1, 0, 2, 0, 3, 1, 3, 1, 2, 2, 0};
  const std::vector<double> raw_values{1, 2, 3, -3, 4, 5, -1,
                                       1, 2, 3, -2, 0, 7, -2};
  const std::vector<std::vector<CompactEntry>> expected{
      {{0, 3}}, {{0, 3}, {2, 5}}, {{3, 4}}, {{0, -2}, {2, 7}}};

  std::vector<std::vector<CompactEntry>> reduced_columns(4);
  const std::vector<pforr::IndexRange> ranges = pforr::split_input_range(
      pforr::IndexRange(0, reduced_columns.size()), n_threads, 1);
  std::vector<ParallelReduceWorkspace> workspaces;
  workspaces.reserve(ranges.size());
  for (std::size_t chunk = 0; chunk < ranges.size(); chunk++) {
    workspaces.emplace_back(reduced_columns.size());
  }

  flotsam_detail::reduce_raw_columns_parallel_core(
      column_starts, column_counts, raw_rows, raw_values, reduced_columns,
      workspaces, n_threads);

  for (std::size_t col = 0; col < expected.size(); col++) {
    if (!entries_match(reduced_columns[col], expected[col])) {
      return false;
    }
  }
  return true;
}

} // namespace

int main() {
  const std::size_t thread_counts[] = {1, 2, 3, 4, 8};
  for (std::size_t repetition = 0; repetition < 200; repetition++) {
    for (const std::size_t n_threads : thread_counts) {
      if (!reduction_matches(n_threads)) {
        std::cerr << "parallel reduction mismatch at repetition " << repetition
                  << " with " << n_threads << " requested threads\n";
        return 1;
      }
    }
  }
  return 0;
}
