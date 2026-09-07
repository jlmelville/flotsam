#ifndef FLOTSAM_LTSA_PARALLEL_REDUCE_H
#define FLOTSAM_LTSA_PARALLEL_REDUCE_H

#include <algorithm>
#include <cstddef>
#include <vector>

#include "pforr.h"

struct CompactEntry {
  int row;
  double value;
};

namespace flotsam_detail {

struct ParallelReduceWorkspace {
  explicit ParallelReduceWorkspace(std::size_t n_obs)
      : row_sums(n_obs, 0.0), row_seen(n_obs, -1) {
    touched_rows.reserve(std::min(n_obs, static_cast<std::size_t>(1024)));
  }

  std::vector<double> row_sums;
  std::vector<int> row_seen;
  std::vector<int> touched_rows;
};

struct ColumnReduceWorker {
  const std::vector<std::size_t> *column_starts;
  const std::vector<std::size_t> *column_counts;
  const std::vector<int> *raw_rows;
  const std::vector<double> *raw_values;
  std::vector<std::vector<CompactEntry>> *reduced_columns;
  std::vector<ParallelReduceWorkspace> *workspaces;

  void operator()(std::size_t begin, std::size_t end, std::size_t chunk_id) {
    ParallelReduceWorkspace &workspace = (*workspaces)[chunk_id];

    for (std::size_t col = begin; col < end; col++) {
      const int marker = static_cast<int>(col);
      workspace.touched_rows.clear();
      const std::size_t start = (*column_starts)[col];
      const std::size_t count = (*column_counts)[col];
      for (std::size_t pos = start; pos < start + count; pos++) {
        const int row = (*raw_rows)[pos];
        if (workspace.row_seen[row] != marker) {
          workspace.row_seen[row] = marker;
          workspace.row_sums[row] = 0.0;
          workspace.touched_rows.push_back(row);
        }
        workspace.row_sums[row] += (*raw_values)[pos];
      }

      std::sort(workspace.touched_rows.begin(), workspace.touched_rows.end());
      std::vector<CompactEntry> &out = (*reduced_columns)[col];
      out.reserve(workspace.touched_rows.size());
      for (const int row : workspace.touched_rows) {
        const double value = workspace.row_sums[row];
        if (value != 0.0) {
          out.push_back(CompactEntry{row, value});
        }
      }
    }
  }
};

inline void reduce_raw_columns_parallel_core(
    const std::vector<std::size_t> &column_starts,
    const std::vector<std::size_t> &column_counts,
    const std::vector<int> &raw_rows, const std::vector<double> &raw_values,
    std::vector<std::vector<CompactEntry>> &reduced_columns,
    std::vector<ParallelReduceWorkspace> &workspaces, std::size_t n_threads) {
  ColumnReduceWorker worker{&column_starts, &column_counts,   &raw_rows,
                            &raw_values,    &reduced_columns, &workspaces};
  pforr::parallel_for_indexed(0, reduced_columns.size(), worker, n_threads, 1);
}

} // namespace flotsam_detail

#endif
