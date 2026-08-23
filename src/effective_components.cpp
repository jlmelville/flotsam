#include "ltsa_internal.h"

namespace {

class EffectiveComponentSet {
public:
  explicit EffectiveComponentSet(std::size_t n_obs)
      : parent_(checked_vector_size<int>(n_obs, "component parents")),
        tree_size_(
            checked_vector_size<std::size_t>(n_obs, "component tree sizes"),
            1) {
    for (std::size_t obs = 0; obs < n_obs; obs++) {
      parent_[obs] = static_cast<int>(obs);
    }
  }

  int find_root(int node) {
    int root = node;
    while (parent_[root] != root) {
      root = parent_[root];
    }
    while (parent_[node] != node) {
      const int next = parent_[node];
      parent_[node] = root;
      node = next;
    }
    return root;
  }

  void unite(int left, int right) {
    int left_root = find_root(left);
    int right_root = find_root(right);
    if (left_root == right_root) {
      return;
    }
    if (tree_size_[left_root] < tree_size_[right_root]) {
      std::swap(left_root, right_root);
    }
    parent_[right_root] = left_root;
    tree_size_[left_root] += tree_size_[right_root];
  }

private:
  std::vector<int> parent_;
  std::vector<std::size_t> tree_size_;
};

} // namespace

[[cpp11::register]] cpp11::list compute_effective_components_cpp(
    const cpp11::integers &transposed_neighbor_indices, int n_obs,
    int n_neighbors) {
  if (n_obs < 1) {
    cpp11::stop("n_obs must be positive");
  }
  if (n_neighbors < 1) {
    cpp11::stop("n_neighbors must be positive");
  }

  const std::size_t n_obs_size = static_cast<std::size_t>(n_obs);
  const std::size_t n_neighbors_size = static_cast<std::size_t>(n_neighbors);
  const std::size_t expected_size = checked_size_mul(
      n_obs_size, n_neighbors_size, "Too many component graph indices");
  if (static_cast<std::size_t>(transposed_neighbor_indices.size()) !=
      expected_size) {
    cpp11::stop("Inconsistent component graph dimensions");
  }

  EffectiveComponentSet components(n_obs_size);
  const int *indices = INTEGER(transposed_neighbor_indices.data());
  for (std::size_t obs = 0; obs < n_obs_size; obs++) {
    const std::size_t offset = obs * n_neighbors_size;
    const int representative =
        checked_zero_based_neighbor_index(indices[offset], n_obs_size);
    for (std::size_t local = 1; local < n_neighbors_size; local++) {
      const int neighbor = checked_zero_based_neighbor_index(
          indices[offset + local], n_obs_size);
      components.unite(representative, neighbor);
    }
  }

  std::vector<int> root_labels(
      checked_vector_size<int>(n_obs_size, "component root labels"), -1);
  std::vector<int> membership(
      checked_vector_size<int>(n_obs_size, "component membership"));
  std::vector<int> sizes;
  sizes.reserve(checked_vector_size<int>(n_obs_size, "component sizes"));
  for (std::size_t obs = 0; obs < n_obs_size; obs++) {
    const int root = components.find_root(static_cast<int>(obs));
    if (root_labels[root] < 0) {
      root_labels[root] = static_cast<int>(sizes.size());
      sizes.push_back(0);
    }
    const int label = root_labels[root];
    membership[obs] = label + 1;
    sizes[label]++;
  }

  cpp11::sexp r_sizes(cpp11::as_sexp(sizes));
  cpp11::sexp r_membership(cpp11::as_sexp(membership));
  return cpp11::writable::list(
      {cpp11::named_arg("component_count") = static_cast<int>(sizes.size()),
       cpp11::named_arg("component_sizes") = r_sizes,
       cpp11::named_arg("component_membership") = r_membership});
}
