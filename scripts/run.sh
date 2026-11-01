#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Thin v2 launcher: all input files are resolved from the caller's directory.
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
if [[ ${1-} == -h || ${1-} == --help || $# == 0 ]]; then
    cat <<'EOF'
Usage: scripts/run.sh MATRIX.mtx [csrspmv options...]
Runs the v2 auto backend with verbose timing. Later options can override
these defaults, e.g. --kernel=original, --repeat=9, --source-vector=X.mtx.
Set OMP_NUM_THREADS to choose the requested team size. Existing OpenMP
settings are respected; unset binding/waiting settings get fixed defaults.
The script never rebuilds, warms up, or repeats the program implicitly.
EOF
    exit 0
fi
if [[ ! -x "$root/csrspmv" ]]; then
    printf 'Build first: %s/scripts/build.sh --check\n' "$root" >&2
    exit 1
fi
export OMP_DYNAMIC=${OMP_DYNAMIC:-FALSE}
export OMP_PLACES=${OMP_PLACES:-cores}
export OMP_PROC_BIND=${OMP_PROC_BIND:-spread}
export OMP_WAIT_POLICY=${OMP_WAIT_POLICY:-PASSIVE}
exec "$root/csrspmv" --kernel=v2 --verbose "$@"
