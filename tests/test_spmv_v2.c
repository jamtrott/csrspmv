/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "spmv_plan.h"
#include <assert.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int fail_after=-1,checks;
void *__real_malloc(size_t);
void *__wrap_malloc(size_t n)
{
    if(fail_after==0) {fail_after=-1;return NULL;}
    if(fail_after>0) fail_after--;
    return __real_malloc(n);
}
static unsigned rng=20261026;
static unsigned next(void) {rng=1664525u*rng+1013904223u;return rng;}
static void boundaries(const struct csr_matrix_int32 *a)
{
    for(int nt=1;nt<=257;nt=nt<24 ? nt+1:nt+233) for(int work=0;work<2;work++) {
        int64_t last=0;
        for(int t=0;t<=nt;t++) {
            int64_t cut=spmv_v2_row_cut(a,t,nt,work);
            assert(cut>=last && cut<=a->num_rows);last=cut;
        }
        assert(last==a->num_rows); checks++;
    }
    /* Independent literal merge of row-end events and nonzero events. */
    int64_t row=0,nz=0,total=a->num_rows+a->num_nonzeros;
    for(int64_t d=0;d<=total;d++) {
        int64_t r,k;spmv_v2_merge_cut(a,d,&r,&k);
        assert(r==row && k==nz);
        if(d<total) {
            if(row<a->num_rows && a->row_ptr[row+1]<=nz) row++;
            else nz++;
        }
    }
    checks++;
}
static void run_case(int m,int n,int len,int profile,int special)
{
    struct csr_matrix_int32 a={0};a.num_rows=m;a.num_columns=n;a.value_format=csr_value_f32;
    a.row_ptr=calloc((size_t)m+1,sizeof(int64_t));
    for(int r=0;r<m;r++) {
        int l=len;
        if(profile==1) l=next()%(len+1);
        if(profile==2) l=r==0 ? 10001:(r==18 ? 3007:1);
        if(profile==3 && r%100) l=0;
        if(profile==4 && r==m/2) l=7;
        if(profile==5 && r==0) l=7;
        if(profile==6 && r==m-1) l=7;
        a.row_ptr[r+1]=a.row_ptr[r]+(n ? l:0);
    }
    a.num_nonzeros=a.row_ptr[m];
    a.column_indices=malloc(((size_t)a.num_nonzeros+1)*sizeof(int32_t));
    a.values=malloc(((size_t)a.num_nonzeros+1)*sizeof(float));
    for(int64_t k=0;k<a.num_nonzeros;k++) {
        a.column_indices[k]=next()%n;
        ((float *)a.values)[k]=((int)(next()%31)-15)*0.125f;
        if(special==2) ((float *)a.values)[k]=(float[]){1e20f,1,-1e20f,0x1p-140f,-0.0f}[k%5];
    }
    boundaries(&a);
    struct vector x={n,vector_value_f32,calloc((size_t)n+1,sizeof(float))};
    struct vector y={m,vector_value_f32,calloc((size_t)m+1,sizeof(float))};
    long double *ref=malloc(((size_t)m+1)*sizeof(long double));
    long double *bound=malloc(((size_t)m+1)*sizeof(long double));
    float *ordered=malloc(((size_t)m+1)*sizeof(float));
    for(int kernel=0;kernel<4;kernel++) for(int partition=0;partition<3;partition++)
    for(int di=0;di<3;di++) {
        if(kernel==v2_merge && partition) continue;
        struct spmv_v2_options o;spmv_v2_defaults(&o);
        o.kernel=kernel;o.partition=partition;o.prefetch=(int[]){0,32,128}[di];
        o.execution=v2_parallel;
        for(int r=0;r<m;r++) ref[r]=ordered[r]=((float *)y.values)[r]=(r%7-3)*0.25f;
        memset(bound,0,((size_t)m+1)*sizeof(*bound));
        struct spmv_plan p;
        assert(!spmv_v2_prepare(&p,&a,&x,&y,&o,SIZE_MAX));
        for(int rep=0;rep<3;rep++) {
            for(int j=0;j<n;j++) ((float *)x.values)[j]=(j%11-5+rep)*0.25f;
            if(n && special==1) ((float *)x.values)[0]=INFINITY;
            if(n && special==3) ((float *)x.values)[0]=NAN;
            for(int r=0;r<m;r++) {
                long double sum=0,ab=0;float s=0;
                for(int64_t k=a.row_ptr[r];k<a.row_ptr[r+1];k++) {
                    float v=((float *)a.values)[k],xx=((float *)x.values)[a.column_indices[k]];
                    long double prod=(long double)v*xx;
                    sum+=prod;ab+=fabsl(prod);s+=v*xx;
                }
                ref[r]+=sum;ordered[r]+=s;
                bound[r]+=(a.row_ptr[r+1]-a.row_ptr[r]+16)*FLT_EPSILON*ab*8;
            }
            int err=0;int64_t flops=0;
#pragma omp parallel reduction(|:err) reduction(max:flops)
            {flops=0;err=spmv_plan_execute(&p,&x,&y,&flops);}
            assert(!err && flops==2*a.num_nonzeros);
            for(int r=0;r<m;r++) {
                float got=((float *)y.values)[r];
                if(isnan(ref[r])) assert(isnan(got));
                else if(isinf(ref[r])) assert(isinf(got) && !!signbit(got)==!!signbit(ref[r]));
                else assert(isfinite(got) && fabsl(got-ref[r])<=bound[r]+1e-35L);
                if(p.v2.kernel!=v2_merge && isfinite(got))
                    assert(!memcmp(&got,&ordered[r],sizeof(float)));
            }
            checks++;
        }
        spmv_plan_free(&p);
    }
    free(ref);free(bound);free(ordered);free(x.values);free(y.values);
    free(a.row_ptr);free(a.column_indices);free(a.values);
}
static void guards(void)
{
    int64_t rows[]={0,8,16,24};int32_t col[24]={0};float v[24]={0},xx[1]={1},yy[3]={0};
    struct csr_matrix_int32 a={3,1,rows,24,col,csr_value_f32,v};
    struct vector x={1,vector_value_f32,xx},y={3,vector_value_f32,yy};
    struct spmv_v2_options o;spmv_v2_defaults(&o);struct spmv_plan p;
    col[0]=-1;assert(spmv_v2_prepare(&p,&a,&x,&y,&o,SIZE_MAX)==EINVAL);col[0]=0;checks++;
    rows[1]=25;assert(spmv_v2_prepare(&p,&a,&x,&y,&o,SIZE_MAX)==EINVAL);rows[1]=8;checks++;
    a.num_nonzeros=INT64_MAX;assert(spmv_v2_prepare(&p,&a,&x,&y,&o,SIZE_MAX)==EOVERFLOW);a.num_nonzeros=24;checks++;
    a.value_format=csr_value_f64;assert(!spmv_v2_prepare(&p,&a,&x,&y,&o,SIZE_MAX));assert(p.kernel==spmv_original);spmv_plan_free(&p);a.value_format=csr_value_f32;checks++;
    o.kernel=v2_merge;o.execution=v2_parallel;
    assert(!spmv_v2_prepare(&p,&a,&x,&y,&o,0));assert(p.v2.kernel==v2_csr);spmv_plan_free(&p);checks++;
    fail_after=0;assert(!spmv_v2_prepare(&p,&a,&x,&y,&o,SIZE_MAX));assert(p.v2.kernel==v2_csr);spmv_plan_free(&p);checks++;
    o.strict_order=1;assert(spmv_v2_prepare(&p,&a,&x,&y,&o,SIZE_MAX)==EINVAL);o.strict_order=0;checks++;
    o.partition=v2_work;assert(spmv_v2_prepare(&p,&a,&x,&y,&o,SIZE_MAX)==EINVAL);o.partition=v2_part_auto;checks++;
    /* Prepared capacity smaller than the actual team: collective row fallback. */
    int saved=omp_get_max_threads();omp_set_num_threads(1);
    assert(!spmv_v2_prepare(&p,&a,&x,&y,&o,SIZE_MAX));
    int err=0;int64_t flops=0;
#pragma omp parallel num_threads(3) reduction(|:err) reduction(max:flops)
    {flops=0;err=spmv_plan_execute(&p,&x,&y,&flops);}
    assert(!err && flops==48);spmv_plan_free(&p);omp_set_num_threads(saved);checks++;
    rows[1]=rows[2]=rows[3]=0;a.num_nonzeros=0;yy[0]=-0.0f;
    spmv_v2_defaults(&o);assert(!spmv_v2_prepare(&p,&a,&x,&y,&o,SIZE_MAX));
    flops=0;assert(!spmv_plan_execute(&p,&x,&y,&flops));assert(!signbit(yy[0]));spmv_plan_free(&p);checks++;
}
int main(void)
{
    guards();
    int lens[]={0,1,2,4,7,8,9,15,16,17,32,64,256};
    for(int i=0;i<13;i++) {run_case(19,37,lens[i],0,0);run_case(33,7,lens[i],1,0);}
    run_case(0,7,0,0,0);run_case(7,0,0,0,0);
    run_case(47,37,1,2,0);run_case(201,13,8,3,0);
    for(int profile=4;profile<=6;profile++) run_case(19,13,8,profile,0);
    run_case(7,9,8,0,0);run_case(9,9,8,0,0);
    run_case(17,13,9,1,1);run_case(19,13,17,1,2);run_case(17,13,9,1,3);
    printf("PASS v2: %d checks, max_threads=%d\n",checks,omp_get_max_threads());return 0;
}
