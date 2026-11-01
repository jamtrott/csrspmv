/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "spmv_plan.h"
#include <errno.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
const char *spmv_kernel_name(enum spmv_kernel k)
{ return k==spmv_original ? "original":k==spmv_v2 ? "v2":"invalid"; }
void spmv_plan_free(struct spmv_plan *p)
{ free(p->v2.scratch); p->v2.scratch=NULL; p->extra_bytes=0; }
int spmv_plan_prepare(struct spmv_plan *p,const struct csr_matrix_int32 *a,
    const struct vector *x,const struct vector *y,enum spmv_kernel k,
    const struct spmv_v2_options *options,size_t budget)
{
    if(k==spmv_v2) return spmv_v2_prepare(p,a,x,y,options,budget);
    memset(p,0,sizeof(*p));p->matrix=a;p->kernel=k;p->reason="original-source";
    if(k!=spmv_original || a->num_rows<0 || a->num_columns<0 ||
       a->num_nonzeros<0 || x->num_values!=a->num_columns ||
       y->num_values!=a->num_rows) return EINVAL;
    return 0;
}
int spmv_plan_execute(const struct spmv_plan *p,const struct vector *x,
    struct vector *y,int64_t *flops)
{
    if(p->matrix->num_rows!=y->num_values ||
       p->matrix->num_columns!=x->num_values) return EINVAL;
    if(p->kernel==spmv_original) return csr_matrix_int32_spmv(p->matrix,x,y,flops);
    if(p->kernel!=spmv_v2 || x->value_format!=vector_value_f32 ||
       y->value_format!=vector_value_f32) return EINVAL;
    if(*flops>INT64_MAX-2*p->matrix->num_nonzeros) return EOVERFLOW;
    spmv_v2_execute(p,x->values,y->values);
    *flops+=2*p->matrix->num_nonzeros; return 0;
}
