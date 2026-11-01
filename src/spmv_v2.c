/*
Copyright (c) 2011-2016, NVIDIA CORPORATION.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
   *  Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
   *  Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
   *  Neither the name of the NVIDIA CORPORATION nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
/* SPDX-License-Identifier: GPL-3.0-or-later
 * Direct CSR v2. Merge-path search adapted from NVIDIA merge-spmv,
 * cpu_spmv.cpp, commit 18895571ce9af960ee207dae541b0ffc701ea4bb.
 * Copyright (c) 2011-2016, NVIDIA CORPORATION. All rights reserved.
 * The original BSD-3-Clause notice is retained in docs/merge-spmv-LICENSE.txt.
 * Adaptations: C/f32/int64, y+=Ax, actual-team storage, delayed split rows.
 */
#include "spmv_plan.h"
#include <errno.h>
#include <limits.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double now(void)
{
    struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9*t.tv_nsec;
}
void spmv_v2_defaults(struct spmv_v2_options *o)
{
    *o=(struct spmv_v2_options){v2_auto,v2_part_auto,v2_exec_auto,0,0,32768};
}
const char *spmv_v2_name(enum v2_kernel k)
{
    switch(k) {
    case v2_auto:return "auto";
    case v2_csr:return "csr";
    case v2_short8:return "short8";
    case v2_merge:return "merge";
    }
    return "invalid";
}
static int64_t target(int64_t total, int t, int nt)
{
    return (total/nt)*t + (t < total%nt ? t : total%nt);
}
int64_t spmv_v2_row_cut(const struct csr_matrix_int32 *a, int t, int nt, int work)
{
    if (!work) return target(a->num_rows,t,nt);
    int64_t goal=target(a->num_nonzeros+a->num_rows,t,nt);
    int64_t lo=0,hi=a->num_rows;
    while(lo<hi) {
        int64_t mid=lo+(hi-lo)/2;
        if(a->row_ptr[mid]+mid<goal) lo=mid+1;
        else hi=mid;
    }
    return lo;
}
/* Stable tie: row-end events precede equal nonzero-index events. */
void spmv_v2_merge_cut(const struct csr_matrix_int32 *a, int64_t diagonal,
                        int64_t *row, int64_t *nz)
{
    int64_t lo=diagonal>a->num_nonzeros ? diagonal-a->num_nonzeros : 0;
    int64_t hi=diagonal<a->num_rows ? diagonal:a->num_rows;
    while(lo<hi) {
        int64_t mid=lo+(hi-lo)/2;
        if(a->row_ptr[mid+1] <= diagonal-mid-1) lo=mid+1;
        else hi=mid;
    }
    *row=lo; *nz=diagonal-lo;
}
int spmv_v2_prepare(struct spmv_plan *p, const struct csr_matrix_int32 *a,
                    const struct vector *x, const struct vector *y,
                    const struct spmv_v2_options *o, size_t budget)
{
    memset(p,0,sizeof(*p)); p->matrix=a; p->kernel=spmv_v2; p->reason="selected";
    struct spmv_v2_state *s=&p->v2; s->options=*o;
    if(o->kernel<v2_auto || o->kernel>v2_merge ||
       o->partition<v2_part_auto || o->partition>v2_work ||
       o->execution<v2_exec_auto || o->execution>v2_parallel ||
       (o->prefetch!=0 && o->prefetch!=32 && o->prefetch!=128) ||
       o->serial_threshold<0 || (o->kernel==v2_merge &&
       (o->strict_order || o->partition!=v2_part_auto))) return EINVAL;
    if(a->num_rows<0 || a->num_columns<0 || a->num_nonzeros<0 ||
       x->num_values!=a->num_columns || y->num_values!=a->num_rows) return EINVAL;
    if(a->value_format!=csr_value_f32 || x->value_format!=vector_value_f32 ||
       y->value_format!=vector_value_f32) {
        p->kernel=spmv_original; p->reason="unsupported-types"; return 0;
    }
    if(a->num_nonzeros>INT64_MAX/2 || a->num_nonzeros>INT64_MAX-a->num_rows)
        return EOVERFLOW;
    if(!a->row_ptr || a->row_ptr[0]!=0 ||
       a->row_ptr[a->num_rows]!=a->num_nonzeros ||
       (a->num_nonzeros && (!a->column_indices || !a->values)) ||
       (x->num_values && !x->values) || (y->num_values && !y->values)) return EINVAL;
    int nt=omp_get_max_threads();
    s->serial=o->execution==v2_serial || (o->execution==v2_exec_auto &&
        (nt==1 || a->num_nonzeros+a->num_rows<=o->serial_threshold));
    /* Full validation is always charged. Parallel only when large enough.
     * The fused row scan also establishes short8 and longest-row statistics. */
    int par=!s->serial && a->num_nonzeros+a->num_rows>=100000;
    int bad=0,all8=1; int64_t longest=0;
    double t=now();
#pragma omp parallel for if(par) reduction(|:bad) reduction(&:all8) reduction(max:longest) schedule(static)
    for(int64_t r=0;r<a->num_rows;r++) {
        int64_t begin=a->row_ptr[r],end=a->row_ptr[r+1];
        if(begin<0 || end<begin || end>a->num_nonzeros) bad=1;
        else {
            int64_t len=end-begin;
            all8 &= len==8;
            if(len>longest) longest=len;
        }
    }
    s->row_scan_s=now()-t;
    if(bad) return EINVAL;
    t=now();
#pragma omp parallel for if(par) reduction(|:bad) schedule(static)
    for(int64_t k=0;k<a->num_nonzeros;k++)
        if(a->column_indices[k]<0 || a->column_indices[k]>=a->num_columns) bad=1;
    s->column_validation_s=now()-t;
    if(bad) return EINVAL;
    t=now(); s->all8=all8; s->longest=longest;
    s->kernel=o->kernel==v2_auto ? (all8 ? v2_short8:v2_csr):o->kernel;
    if(s->kernel==v2_short8 && !all8) {
        s->kernel=v2_csr; p->reason="not-all-rows-eight";
    }
    int actual=s->serial ? 1:nt;
    int64_t total=a->num_rows+a->num_nonzeros,max_work=0;
    for(int t0=0;t0<actual;t0++) {
        int64_t b=target(a->num_rows,t0,actual),e=target(a->num_rows,t0+1,actual);
        int64_t w=a->row_ptr[e]-a->row_ptr[b]+e-b;
        if(w>max_work) max_work=w;
    }
    double ideal=(double)total/actual;
    s->partition=o->partition==v2_part_auto ?
        (max_work>1.5*ideal ? v2_work:v2_rows):o->partition;
    if(o->kernel==v2_auto && !o->strict_order && o->partition==v2_part_auto &&
       actual>1 && longest>2*ideal) s->kernel=v2_merge;
    if(s->kernel==v2_merge) {
        s->capacity=actual;
        if((size_t)actual>SIZE_MAX/sizeof(*s->scratch) ||
           (size_t)actual*sizeof(*s->scratch)>budget) {
            s->kernel=v2_csr; p->reason="v2-scratch-budget";
        } else {
            s->scratch=malloc((size_t)actual*sizeof(*s->scratch));
            if(!s->scratch) {s->kernel=v2_csr;p->reason="v2-allocation-failed";}
            else p->extra_bytes=(size_t)actual*sizeof(*s->scratch);
        }
    }
    s->policy_s=now()-t;
    return 0;
}

/* Local ordered scalar operations. No global fast-math/novec. Short8 is a
 * fixed-width adaptation of Trotter's published unrolling method, not ELL. */
#if defined(__GNUC__) && !defined(__clang__)
__attribute__((optimize("no-tree-loop-vectorize,no-tree-slp-vectorize")))
#endif
static void rows_f32(const struct csr_matrix_int32 *a, const float *restrict x,
                     float *restrict y, int64_t begin, int64_t end, int short8, int d)
{
    const float *restrict v=a->values;
    const int32_t *restrict j=a->column_indices;
    int64_t stop=a->row_ptr[end];
    for(int64_t r=begin;r<end;r++) {
        int64_t k=short8 ? 8*r:a->row_ptr[r];
        float sum=0;
        if(short8) {
            if(d && d<stop-k) __builtin_prefetch(x+j[k+d],0,1);
            sum+=v[k]*x[j[k]]; sum+=v[k+1]*x[j[k+1]];
            sum+=v[k+2]*x[j[k+2]]; sum+=v[k+3]*x[j[k+3]];
            sum+=v[k+4]*x[j[k+4]]; sum+=v[k+5]*x[j[k+5]];
            sum+=v[k+6]*x[j[k+6]]; sum+=v[k+7]*x[j[k+7]];
        } else {
            for(int64_t e=a->row_ptr[r+1];k<e;k++) {
                if(d && d<stop-k) __builtin_prefetch(x+j[k+d],0,1);
                sum+=v[k]*x[j[k]];
            }
        }
        y[r]+=sum;
    }
}
static void merge_f32(const struct spmv_plan *p,const float *restrict x,
                       float *restrict y,int tid,int nt)
{
    const struct csr_matrix_int32 *a=p->matrix;
    struct v2_scratch *slots=p->v2.scratch,*s=&slots[tid];
    const float *v=a->values; const int32_t *j=a->column_indices;
    int64_t r,k,er,ek,total=a->num_rows+a->num_nonzeros;
    spmv_v2_merge_cut(a,target(total,tid,nt),&r,&k);
    spmv_v2_merge_cut(a,target(total,tid+1,nt),&er,&ek);
    s->count=0;
    for(;r<er;r++) {
        int partial=k>a->row_ptr[r];
        float sum=0;
        for(int64_t e=a->row_ptr[r+1];k<e;k++) sum+=v[k]*x[j[k]];
        if(partial) s->fragment[s->count++]=(struct v2_fragment){r,sum};
        else y[r]+=sum;
    }
    if(k<ek) {
        float sum=0;
        for(;k<ek;k++) sum+=v[k]*x[j[k]];
        s->fragment[s->count++]=(struct v2_fragment){r,sum};
    }
#pragma omp barrier
#pragma omp single
    {
        int64_t row=-1; float sum=0;
        for(int t=0;t<nt;t++) for(int f=0;f<slots[t].count;f++) {
            struct v2_fragment v0=slots[t].fragment[f];
            if(v0.row!=row) {
                if(row>=0) y[row]+=sum;
                row=v0.row; sum=0;
            }
            sum+=v0.sum;
        }
        if(row>=0) y[row]+=sum;
    } /* single's implicit barrier completes all split-row updates */
}
int spmv_v2_execute(const struct spmv_plan *p,const float *x,float *y)
{
    int tid=omp_get_thread_num(),nt=omp_get_num_threads();
    const struct spmv_v2_state *s=&p->v2;
    if(s->kernel==v2_merge && nt>1 && nt<=s->capacity) {
        merge_f32(p,x,y,tid,nt); return 0;
    }
    int work=s->partition==v2_work;
    int64_t b=spmv_v2_row_cut(p->matrix,tid,nt,work);
    int64_t e=spmv_v2_row_cut(p->matrix,tid+1,nt,work);
    rows_f32(p->matrix,x,y,b,e,s->kernel==v2_short8,s->options.prefetch);
    /* For a larger-than-prepared team, ALL threads take this safe row path. */
    if(omp_in_parallel()) {
#pragma omp barrier
    }
    return 0;
}
