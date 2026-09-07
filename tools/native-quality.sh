#!/usr/bin/env bash

set -euo pipefail

usage() {
  printf 'Usage: %s\n' "${0##*/}"
}

case "${1-}" in
  "")
    ;;
  -h|--help)
    usage
    exit 0
    ;;
  *)
    printf '%s: unexpected argument: %s\n' "${0##*/}" "$1" >&2
    usage >&2
    exit 2
    ;;
esac

script_dir=$(CDPATH='' cd -- "$(dirname -- "$0")" && pwd)
repo_root=$(CDPATH='' cd -- "$script_dir/.." && pwd)
cd "$repo_root"

required_commands=(Rscript g++ clang++ clang-format mktemp)
for required_command in "${required_commands[@]}"; do
  if ! command -v "$required_command" >/dev/null 2>&1; then
    printf '%s: required command not found: %s\n' \
      "${0##*/}" "$required_command" >&2
    exit 1
  fi
done

required_formatter_version=21.1.8
formatter_version_output=$(clang-format --version)
if [[ $formatter_version_output =~ clang-format\ version\ ([0-9]+\.[0-9]+\.[0-9]+) ]]; then
  formatter_version=${BASH_REMATCH[1]}
else
  printf '%s: unable to determine clang-format version from: %s\n' \
    "${0##*/}" "$formatter_version_output" >&2
  exit 1
fi
if [[ $formatter_version != "$required_formatter_version" ]]; then
  printf '%s: clang-format %s is required, found %s\n' \
    "${0##*/}" "$required_formatter_version" "$formatter_version" >&2
  exit 1
fi

r_include=$(Rscript --vanilla -e 'cat(R.home("include"))')
cpp11_include=$(Rscript --vanilla -e 'cat(system.file("include", package = "cpp11"))')
if [[ ! -d "$r_include" ]]; then
  printf '%s: R include directory not found: %s\n' \
    "${0##*/}" "$r_include" >&2
  exit 1
fi
if [[ ! -d "$cpp11_include" ]]; then
  printf '%s: cpp11 is not installed or its include directory is unavailable\n' \
    "${0##*/}" >&2
  exit 1
fi

maintained_sources=(
  src/effective_components.cpp
  src/ltsa_assembly_common.cpp
  src/ltsa_checks.cpp
  src/ltsa_local_weights.cpp
  src/ltsa_parallel_assembly.cpp
  src/ltsa_serial_assembly.cpp
  src/ltsa_sparse_normalization.cpp
  src/ltsa_triplet_builder.cpp
)

format_sources=(
  "${maintained_sources[@]}"
  src/ltsa_internal.h
  src/ltsa_parallel_reduce.h
  inst/include/pforr.h
  tools/tsan-ltsa.cpp
)

compiler_flags=(
  -std=gnu++17
  -pthread
  -fPIC
  -Isrc
  -Iinst/include
  -isystem "$r_include"
  -isystem "$cpp11_include"
  -Wall
  -Wextra
  -Wpedantic
  -Wformat=2
  -Wnull-dereference
  -Werror
  -fsyntax-only
)

compilers=(g++ clang++)
for compiler in "${compilers[@]}"; do
  for source_file in "${maintained_sources[@]}"; do
    "$compiler" "${compiler_flags[@]}" "$source_file"
  done
done

scratch_dir=$(mktemp -d "${TMPDIR:-/tmp}/flotsam-native-quality.XXXXXX")
cleanup() {
  rm -rf -- "$scratch_dir"
}
trap cleanup EXIT HUP INT TERM

consumer_source="$scratch_dir/pforr-consumer.cpp"
cat >"$consumer_source" <<'EOF'
#include <pforr.h>

#include <atomic>
#include <cstddef>

struct Counter {
  std::atomic<std::size_t> *count;

  void operator()(std::size_t begin, std::size_t end) {
    count->fetch_add(end - begin, std::memory_order_relaxed);
  }
};

int main() {
  std::atomic<std::size_t> count(0);
  Counter worker{&count};
  pforr::parallel_for(0, 17, worker, 3, 2);
  return count.load(std::memory_order_relaxed) == 17 ? 0 : 1;
}
EOF

consumer_flags=(
  -std=c++11
  -pthread
  -Iinst/include
  -Wall
  -Wextra
  -Wpedantic
  -Wformat=2
  -Werror
)
for compiler in "${compilers[@]}"; do
  consumer_executable="$scratch_dir/pforr-${compiler}"
  "$compiler" "${consumer_flags[@]}" \
    "$consumer_source" -o "$consumer_executable"
  "$consumer_executable"
done

clang-format --dry-run --Werror "${format_sources[@]}"

printf 'native-quality: PASS (GCC, Clang, pforr consumer, clang-format)\n'
