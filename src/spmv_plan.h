/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef SPMV_PLAN_H
#define SPMV_PLAN_H
#include "csr.h"
#include "vector.h"
#include "spmv_v2.h"
#include <stddef.h>
enum spmv_kernel { spmv_original, spmv_v2 };
/* Borrows validated CSR: rebuild after structural mutation. Direct CSR reads
 * current values; x/y may change. Non-overlapping x/y are required.
 * prepare/free are serial; execute is collective in the current team, or
 * called outside a team for the true serial path. No concurrent shared y.
 */
struct spmv_plan {
    const struct csr_matrix_int32 *matrix;
    enum spmv_kernel kernel;
    const char *reason;
    size_t extra_bytes;
    struct spmv_v2_state v2;
};
const char *spmv_kernel_name(enum spmv_kernel);
int spmv_plan_prepare(struct spmv_plan *,const struct csr_matrix_int32 *,
    const struct vector *,const struct vector *,enum spmv_kernel,
    const struct spmv_v2_options *,size_t);
int spmv_plan_execute(const struct spmv_plan *,const struct vector *,struct vector *,int64_t *);
void spmv_plan_free(struct spmv_plan *);
#endif
