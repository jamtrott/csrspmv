#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# A clean build avoids reusing objects compiled with different CFLAGS.
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
mode=portable
jobs=${JOBS:-2}
check=0
usage() {
    cat <<'EOF'
Usage: scripts/build.sh [--portable|--native] [--jobs N] [--check]
  --portable  Optimized build without -march=native (default).
  --native    Optimize for this machine; rebuild before using another CPU.
  --check     Also run the v2 unit suite after the build.
CC, CFLAGS and LDFLAGS may be supplied in the environment.
CFLAGS, when set, replaces the optimization/ISA defaults; required include
and OpenMP flags are appended. This script always runs make clean first.
EOF
}
while (($#)); do
    case $1 in
        --portable) mode=portable ;;
        --native) mode=native ;;
        --jobs)
            if (($# < 2)); then usage >&2; exit 2; fi
            jobs=$2; shift ;;
        --check) check=1 ;;
        -h|--help) usage; exit 0 ;;
        *) printf 'Unknown option: %s\n' "$1" >&2; usage >&2; exit 2 ;;
    esac
    shift
done
if [[ ! $jobs =~ ^[1-9][0-9]*$ ]]; then
    printf 'jobs must be a positive integer\n' >&2; exit 2
fi
flags=-O3
if [[ $mode == native ]]; then flags+=' -march=native'; fi
flags="${CFLAGS-$flags} -g -Wall -iquote src -fopenmp"
cd -- "$root"
make clean
make -j "$jobs" "CC=${CC:-cc}" "CFLAGS=$flags" "LDFLAGS=${LDFLAGS-}"
if ((check)); then
    make check "CC=${CC:-cc}" "CFLAGS=$flags" "LDFLAGS=${LDFLAGS-}"
fi
