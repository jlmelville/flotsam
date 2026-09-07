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

required_commands=(g++ mktemp)
for required_command in "${required_commands[@]}"; do
  if ! command -v "$required_command" >/dev/null 2>&1; then
    printf '%s: required command not found: %s\n' \
      "${0##*/}" "$required_command" >&2
    exit 1
  fi
done

if [[ $(uname -s) != Linux ]]; then
  printf '%s: ThreadSanitizer check requires Linux\n' "${0##*/}" >&2
  exit 1
fi

scratch_dir=$(mktemp -d "${TMPDIR:-/tmp}/flotsam-tsan.XXXXXX")
cleanup() {
  rm -rf -- "$scratch_dir"
}
trap cleanup EXIT HUP INT TERM

test_executable="$scratch_dir/ltsa-parallel-reduce-tsan"
g++ \
  -std=gnu++17 \
  -pthread \
  -Isrc \
  -Iinst/include \
  -Wall \
  -Wextra \
  -Wpedantic \
  -Wformat=2 \
  -Werror \
  -g \
  -O1 \
  -fsanitize=thread \
  -fno-omit-frame-pointer \
  tools/tsan-ltsa.cpp \
  -o "$test_executable"

TSAN_OPTIONS="halt_on_error=1" "$test_executable"

printf 'tsan-ltsa: PASS (production parallel column reduction)\n'
