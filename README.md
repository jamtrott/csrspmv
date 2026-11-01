# csrspmv

A Matrix Market benchmark for sparse matrix-vector multiplication in
compressed sparse row (CSR) format, computing **y += A*x** on each call.

This branch adds an **opt-in v2 backend** for f32 matrix/vector data. The
original CSR backend remains the command-line default, and other value
types fall back to it. No SELL copy, cached A*x, or matrix reordering is
required. Original CSR/OpenMP functionality is retained.

## Build and check

Requirements: Linux or WSL, Bash, GNU Make, a C compiler with OpenMP
(tested with GCC 11.4), and GNU awk/standard POSIX tools for the smoke check.
No external numerical library or Python package is needed to build or run.

~~~sh
./scripts/build.sh --portable --check
./scripts/check.sh

# For benchmarks on the current CPU:
./scripts/build.sh --native --jobs 12 --check
~~~

The build script cleans old objects before building, so changing ISA or
optimization flags cannot silently reuse an older object. Portable is the
default; native binaries must be rebuilt for a different CPU. Set CC,
CFLAGS or LDFLAGS if needed; CFLAGS replaces the optimization defaults.
The script appends the required include and OpenMP flags. The existing
make interface is also available:

~~~sh
make clean
make CFLAGS='-O3 -march=native -g -Wall -iquote src -fopenmp'
make check CFLAGS='-O3 -march=native -g -Wall -iquote src -fopenmp'
~~~

## Run

~~~sh
OMP_NUM_THREADS=12 ./scripts/run.sh test.mtx --repeat=9

# The unchanged default entry remains available:
OMP_NUM_THREADS=12 ./csrspmv test.mtx --verbose

# Fair, same-binary baseline and optimized runs:
OMP_NUM_THREADS=12 ./scripts/run.sh input.mtx --repeat=9 --kernel=original
OMP_NUM_THREADS=12 ./scripts/run.sh input.mtx --repeat=9 --kernel=v2

# Save a result, with a source vector if required:
./scripts/run.sh input.mtx --source-vector=x.mtx \
    --destination-vector=y.mtx --destination-vector-prec=9
~~~

scripts/run.sh selects v2 and verbose timing, but later CLI options can
override them. It respects the caller's OpenMP settings; when unset, it
uses OMP_DYNAMIC=FALSE, OMP_PLACES=cores, OMP_PROC_BIND=spread and
OMP_WAIT_POLICY=PASSIVE. Set OMP_NUM_THREADS explicitly for reproducible
comparisons. Input/output paths are relative to the caller's directory.
The launcher does not rebuild, warm up or repeat the program implicitly.

The source vector defaults to all ones and y starts at zero. R calls
really perform R accumulations; R*A*x is not substituted. To expand the
stored triangle of a symmetric Matrix Market input, explicitly pass
--include-symmetric-part to both comparison runs.

## v2 backend

| Option | Meaning |
| --- | --- |
| --kernel=original\|v2 | Original is the executable default; v2 is opt-in. |
| --v2-kernel=auto\|csr\|short8\|merge | Direct rows, validated fixed-eight unrolling, or split-row merge. |
| --partition=auto\|rows\|work | Equal rows or contiguous nnz+rows work cuts. |
| --execution=auto\|serial\|parallel | Real serial execution for small work or a shared team. |
| --serial-threshold=N | Auto serial threshold in rows+nnz; default 32768. |
| --strict-order | Disable split-row summation; incompatible with explicit merge. |
| --prefetch=0\|32\|128 | Optional lookahead; default 0. |
| --plan-memory-mib=N | Merge scratch budget; default 4 MiB. |

Auto performs full CSR validation and gathers row statistics. It chooses
short8 only when all rows have exactly eight entries; otherwise it uses
CSR. A max/mean equal-row work ratio above 1.5 selects work cuts. A row
longer than twice the mean event count may select merge when order
reassociation is allowed. These are frozen engineering rules, not
universal hardware constants.

The plan borrows CSR storage and reads current values/x/y on every call.
Rebuild it after changes to row pointers, column indices, dimensions or
types. x and y must not overlap. Prepare/free are outside the team;
execute is collective in a team, or called serially. See [interface and
algorithm details](docs/v2.md).

## Costs, numerical behavior and evidence

Verbose output reports spmv_prepare, execution, actual team size and
spmv_cleanup. Compare **prepare + all R executions + cleanup**; process
wall time also includes Matrix Market I/O and conversion. All validation
and scratch allocation are charged. Work cuts use binary search in
existing row_ptr; merge uses 64 bytes of scratch per prepared thread.

Only the v2 translation unit disables FP contraction, and its row
function locally disables loop/SLP vectorization. The original backend
retains its existing compiler settings. Merge reassociates split-row
sums; use an error bound rather than assuming universal bitwise
equivalence. Strict-order disables splitting but does not promise
bitwise identity across compilers/FMA settings.

On the recorded Ryzen 9 7845HX/WSL2 environment, the original v2 study's
12-million-row, 96-million-entry input at T24/R9 had median total cost
0.811400 s (original) versus 0.734880 s (v2), a 9.43% decrease including
preparation/cleanup. This is a measured workload result, not an all-input
speed guarantee. Some shapes regress, and short8 alone did not have a
stable independent advantage over the matched generic CSR control.

Fixed inputs, sources, full positive/negative results, paired intervals
and reproduction scripts are archived in the
[course materials repository](https://github.com/visibychen/spmv-course-materials/tree/main/experiments/2026-10-26-v2).
Later final-route verification is recorded
[separately](https://github.com/visibychen/spmv-course-materials/tree/main/experiments/2026-10-31-v3);
timings from different batches/loaders must not be pooled.
make check covers row distributions, real changing x/y, repeated calls,
team changes, malformed inputs and simulated allocation failure.
scripts/check.sh additionally checks both CLI backends against the exact
known result of test.mtx.

## Sources and license

csrspmv is GPL-3.0-or-later; see [COPYING](COPYING). Merge-path search is
adapted from Merrill/Garland's CPU implementation, with its BSD-3-Clause
notice retained in the source and [license](docs/merge-spmv-LICENSE.txt).
Method provenance and adaptation boundaries are documented in
[docs/v2.md](docs/v2.md). Third-party real matrices are not redistributed.

The original README copyright and permission notice is retained below:

  Copyright (C) 2023 James D. Trotter

  Copying and distribution of this file, with or without modification,
  are permitted in any medium without royalty provided the copyright
  notice and this notice are preserved.
