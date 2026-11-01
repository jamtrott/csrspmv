#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Unit suite plus an exact, independent CLI oracle for the supplied matrix.
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
cd -- "$root"
if [[ ! -x ./csrspmv ]]; then
    printf 'Build first: scripts/build.sh --check\n' >&2; exit 1
fi
make check
tmp=$(mktemp -d "${TMPDIR:-/tmp}/csrspmv-check.XXXXXXXX")
trap 'rm -rf -- "$tmp"' EXIT
for kernel in original v2; do
    OMP_NUM_THREADS=2 ./csrspmv test.mtx --kernel="$kernel" --repeat=3 \
        --destination-vector="$tmp/$kernel.mtx" --destination-vector-prec=9
    awk '
        /^%/ { next }
        !header++ { if (NF != 1 || $1 != 4) exit 1; next }
        { expected[++row] = $1 }
        END {
            if (row != 4 || expected[1] != 9 || expected[2] != 3 ||
                expected[3] != 9 || expected[4] != 6) exit 1
        }' "$tmp/$kernel.mtx"
done
printf 'CLI oracle passed for original and v2 (R=3).\n'
