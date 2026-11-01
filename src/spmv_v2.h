/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef SPMV_V2_H
#define SPMV_V2_H
#include "csr.h"
#include "vector.h"
#include <stddef.h>

enum v2_kernel { v2_auto, v2_csr, v2_short8, v2_merge };
enum v2_partition { v2_part_auto, v2_rows, v2_work };
enum v2_execution { v2_exec_auto, v2_serial, v2_parallel };
struct spmv_v2_options {
    enum v2_kernel kernel;
    enum v2_partition partition;
    enum v2_execution execution;
    int prefetch, strict_order;
    int64_t serial_threshold;
};
struct v2_fragment { int64_t row; float sum; };
struct v2_scratch {
    struct v2_fragment fragment[2];
    int count;
    char padding[28]; /* 64-byte stride; no over-aligned allocation required */
};
struct spmv_v2_state {
    struct spmv_v2_options options;
    enum v2_kernel kernel;
    enum v2_partition partition;
    int serial, capacity, all8;
    int64_t longest;
    double row_scan_s, column_validation_s, policy_s;
    struct v2_scratch *scratch;
};
struct spmv_plan;
void spmv_v2_defaults(struct spmv_v2_options *);
const char *spmv_v2_name(enum v2_kernel);
int spmv_v2_prepare(struct spmv_plan *, const struct csr_matrix_int32 *,
                    const struct vector *, const struct vector *,
                    const struct spmv_v2_options *, size_t);
int spmv_v2_execute(const struct spmv_plan *, const float *, float *);
/* Validated CSR only. Team-independent, 64-bit, exclusive row boundaries. */
int64_t spmv_v2_row_cut(const struct csr_matrix_int32 *, int, int, int);
void spmv_v2_merge_cut(const struct csr_matrix_int32 *, int64_t, int64_t *, int64_t *);
#endif
