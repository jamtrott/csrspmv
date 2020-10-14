/**
 * Benchmark program for CSR SpMV
 * Copyright (C) 2020 James D. Trotter
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 * Authors: James D. Trotter <james@simula.no>
 *
 * Sparse matrices in the compressed sparse row (CSR) format.
 */

#include "csr.h"
#include "matrix_market.h"
#include "vector.h"

#include <errno.h>

#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * `parse_csr_value_format()` parses a string designating a format for
 * matrices.
 *
 * On success, `parse_csr_value_format()` returns `0`. If the string
 * does not correspond to a valid format, then
 * `parse_csr_value_format()` returns `EINVAL`.
 */
int parse_csr_value_format(
    const char * s,
    enum csr_value_format * format)
{
    if (strcmp("binary", s) == 0) {
        *format = csr_value_binary;
    } else if (strcmp("int32", s) == 0) {
        *format = csr_value_int32;
    } else if (strcmp("int64", s) == 0) {
        *format = csr_value_int64;
    } else if (strcmp("f32", s) == 0) {
        *format = csr_value_f32;
    } else if (strcmp("f64", s) == 0) {
        *format = csr_value_f64;
    } else if (strcmp("complex32", s) == 0) {
        *format = csr_value_complex32;
    } else {
        return EINVAL;
    }
    return 0;
}

/**
 * `csr_value_format_size()` is the size of a single matrix nonzero
 * for the given matrix format.
 */
static int csr_value_format_size(
    enum csr_value_format value_format,
    size_t * size)
{
    switch (value_format) {
    case csr_value_binary:
        *size = 0;
        return 0;
    case csr_value_int32:
        *size = sizeof(int32_t);
        return 0;
    case csr_value_int64:
        *size = sizeof(int64_t);
        return 0;
    case csr_value_f32:
        *size = sizeof(float);
        return 0;
    case csr_value_f64:
        *size = sizeof(double);
        return 0;
    case csr_value_complex32:
        *size = 2*sizeof(float);
        return 0;
    default:
        return EINVAL;
    }
}

/**
 * `csr_matrix_int32_alloc()` allocates a sparse matrix in CSR format.
 */
int csr_matrix_int32_alloc(
    struct csr_matrix_int32 * matrix,
    int32_t num_rows,
    int32_t num_columns,
    int64_t num_nonzeros,
    enum csr_value_format value_format)
{
    int err;

    /* 1. Allocate storage for row pointers. */
    int64_t * row_ptr = (int64_t *) malloc((num_rows+1) * sizeof(int64_t));
    if (!row_ptr)
        return errno;

    /* 2. Allocate storage for column indices. */
    int32_t * column_indices = (int32_t *) malloc(
        num_nonzeros * sizeof(int32_t));
    if (!column_indices) {
        free(row_ptr);
        return errno;
    }

    /* 3. Determine storage required for matrix nonzeros. */
    size_t nonzero_size;
    err = csr_value_format_size(value_format, &nonzero_size);
    if (err) {
        free(column_indices);
        free(row_ptr);
        return err;
    }

    /* 4. Allocate storage for matrix nonzeros. */
    void * values =  malloc(num_nonzeros * nonzero_size);
    if (!values) {
        free(column_indices);
        free(row_ptr);
        return errno;
    }

    matrix->num_rows = num_rows;
    matrix->num_columns = num_columns;
    matrix->row_ptr = row_ptr;
    matrix->num_nonzeros = num_nonzeros;
    matrix->column_indices = column_indices;
    matrix->value_format = value_format;
    matrix->values = values;
    return 0;
}

/**
 * `csr_matrix_int32_free()` destroys the given matrix.
 */
void csr_matrix_int32_free(
    struct csr_matrix_int32 * matrix)
{
    free(matrix->row_ptr);
    free(matrix->column_indices);
    free(matrix->values);
}

/**
 * `csr_matrix_int32_init()` creates a sparse matrix in CSR format.
 */
int csr_matrix_int32_init(
    struct csr_matrix_int32 * matrix,
    int32_t num_rows,
    int32_t num_columns,
    int64_t * row_ptr,
    int32_t * column_indices,
    enum csr_value_format value_format,
    void * values)
{
    matrix->num_rows = num_rows;
    matrix->num_columns = num_columns;
    matrix->row_ptr = row_ptr;
    matrix->num_nonzeros = num_rows > 0 ? row_ptr[num_rows] : 0;
    matrix->column_indices = column_indices;
    matrix->value_format = value_format;
    matrix->values = values;
    return 0;
}

/**
 * `csr_matrix_int32_zero()` zeros the values of a matrix.
 */
int csr_matrix_int32_zero(
    struct csr_matrix_int32 * matrix)
{
    switch (matrix->value_format) {
    case csr_value_binary:
        break;
    case csr_value_int32:
#pragma omp parallel for
        for (int64_t k = 0; k < matrix->num_nonzeros; k++)
            ((int32_t *) matrix->values)[k] = 0;
        break;
    case csr_value_int64:
#pragma omp parallel for
        for (int64_t k = 0; k < matrix->num_nonzeros; k++)
            ((int64_t *) matrix->values)[k] = 0;
        break;
    case csr_value_f32:
#pragma omp parallel for
        for (int64_t k = 0; k < matrix->num_nonzeros; k++)
            ((float *) matrix->values)[k] = 0;
        break;
    case csr_value_f64:
#pragma omp parallel for
        for (int64_t k = 0; k < matrix->num_nonzeros; k++)
            ((double *) matrix->values)[k] = 0;
        break;
    case csr_value_complex32:
#pragma omp parallel for
        for (int64_t k = 0; k < matrix->num_nonzeros; k++)
            ((float *) matrix->values)[2*k+0] = ((float *) matrix->values)[2*k+1] = 0;
        break;
    default:
        return EINVAL;
    }
    return 0;
}

/**
 * `csr_matrix_int32_print()` prints a matrix.
 */
int csr_matrix_int32_print(
    const struct csr_matrix_int32 * matrix,
    FILE * f,
    const char * row_delim,
    const char * nonzero_delim,
    const char * value_delim,
    int colidx_width,
    int value_width,
    int value_precision)
{
    switch (matrix->value_format) {
    case csr_value_binary:
        {
            const int32_t * j = matrix->column_indices;
            int32_t i = 0;
            int64_t k = 0;
            for (i = 0; i < matrix->num_rows-1; i++) {
                for (k = matrix->row_ptr[i]; k < matrix->row_ptr[i+1]-1; k++)
                    fprintf(f, "%*d%s", colidx_width, j[k], nonzero_delim);
                fprintf(f, "%*d%s", colidx_width, j[k], row_delim);
            }
            for (k = matrix->row_ptr[i]; k < matrix->row_ptr[i+1]-1; k++)
                fprintf(f, "%*d%s", colidx_width, j[k], nonzero_delim);
            fprintf(f, "%*d", colidx_width, j[k]);
            return 0;
        }

    case csr_value_int32:
        {
            const int32_t * j = matrix->column_indices;
            const int32_t * a = (const int32_t *) matrix->values;
            int32_t i = 0;
            int64_t k = 0;
            for (i = 0; i < matrix->num_rows-1; i++) {
                for (k = matrix->row_ptr[i]; k < matrix->row_ptr[i+1]-1; k++) {
                    fprintf(f, "%*d%s%*"PRId32"%s",
                            colidx_width, j[k], value_delim,
                            value_width, a[k], nonzero_delim);
                }
                fprintf(f, "%*d%s%*"PRId32"%s",
                        colidx_width, j[k], value_delim,
                        value_width, a[k], row_delim);
            }
            for (k = matrix->row_ptr[i]; k < matrix->row_ptr[i+1]-1; k++) {
                fprintf(f, "%*d%s%*"PRId32"%s",
                        colidx_width, j[k], value_delim,
                        value_width, a[k], nonzero_delim);
            }
            fprintf(f, "%*d%s%*"PRId32,
                    colidx_width, j[k], value_delim,
                    value_width, a[k]);
            return 0;
        }

    case csr_value_int64:
        {
            const int32_t * j = matrix->column_indices;
            const int64_t * a = (const int64_t *) matrix->values;
            int32_t i = 0;
            int64_t k = 0;
            for (i = 0; i < matrix->num_rows-1; i++) {
                for (k = matrix->row_ptr[i]; k < matrix->row_ptr[i+1]-1; k++) {
                    fprintf(f, "%*d%s%*"PRId64"%s",
                            colidx_width, j[k], value_delim,
                            value_width, a[k], nonzero_delim);
                }
                fprintf(f, "%*d%s%*"PRId64"%s",
                        colidx_width, j[k], value_delim,
                        value_width, a[k], row_delim);
            }
            for (k = matrix->row_ptr[i]; k < matrix->row_ptr[i+1]-1; k++) {
                fprintf(f, "%*d%s%*"PRId64"%s",
                        colidx_width, j[k], value_delim,
                        value_width, a[k], nonzero_delim);
            }
            fprintf(f, "%*d%s%*"PRId64,
                    colidx_width, j[k], value_delim,
                    value_width, a[k]);
            return 0;
        }

    case csr_value_f32:
        {
            const int32_t * j = matrix->column_indices;
            const float * a = (const float *) matrix->values;
            int32_t i = 0;
            int64_t k = 0;
            for (i = 0; i < matrix->num_rows-1; i++) {
                for (k = matrix->row_ptr[i]; k < matrix->row_ptr[i+1]-1; k++) {
                    fprintf(f, "%*d%s%*.*f%s",
                            colidx_width, j[k], value_delim,
                            value_width, value_precision, a[k], nonzero_delim);
                }
                fprintf(f, "%*d%s%*.*f%s",
                        colidx_width, j[k], value_delim,
                        value_width, value_precision, a[k], row_delim);
            }
            for (k = matrix->row_ptr[i]; k < matrix->row_ptr[i+1]-1; k++) {
                fprintf(f, "%*d%s%*.*f%s",
                        colidx_width, j[k], value_delim,
                        value_width, value_precision, a[k], nonzero_delim);
            }
            fprintf(f, "%*d%s%*.*f",
                    colidx_width, j[k], value_delim,
                    value_width, value_precision, a[k]);
            return 0;
        }

    case csr_value_f64:
        {
            const int32_t * j = matrix->column_indices;
            const double * a = (const double *) matrix->values;
            int32_t i = 0;
            int64_t k = 0;
            for (i = 0; i < matrix->num_rows-1; i++) {
                for (k = matrix->row_ptr[i]; k < matrix->row_ptr[i+1]-1; k++) {
                    fprintf(f, "%*d%s%*.*lf%s",
                            colidx_width, j[k], value_delim,
                            value_width, value_precision, a[k], nonzero_delim);
                }
                fprintf(f, "%*d%s%*.*lf%s",
                        colidx_width, j[k], value_delim,
                        value_width, value_precision, a[k], row_delim);
            }
            for (k = matrix->row_ptr[i]; k < matrix->row_ptr[i+1]-1; k++) {
                fprintf(f, "%*d%s%*.*lf%s",
                        colidx_width, j[k], value_delim,
                        value_width, value_precision, a[k], nonzero_delim);
            }
            fprintf(f, "%*d%s%*.*lf",
                    colidx_width, j[k], value_delim,
                    value_width, value_precision, a[k]);
            return 0;
        }

    case csr_value_complex32:
        {
            const int32_t * j = matrix->column_indices;
            const float * a = (const float *) matrix->values;
            int32_t i = 0;
            int64_t k = 0;
            for (i = 0; i < matrix->num_rows-1; i++) {
                for (k = matrix->row_ptr[i]; k < matrix->row_ptr[i+1]-1; k++) {
                    fprintf(f, "%*d%s%*.*f%s%*.*f%s",
                            colidx_width, j[k], value_delim,
                            value_width, value_precision, a[2*k+0], value_delim,
                            value_width, value_precision, a[2*k+1], nonzero_delim);
                }
                fprintf(f, "%*d%s%*.*f%s%*.*f%s",
                        colidx_width, j[k], value_delim,
                        value_width, value_precision, a[2*k+0], value_delim,
                        value_width, value_precision, a[2*k+1], row_delim);
            }
            for (k = matrix->row_ptr[i]; k < matrix->row_ptr[i+1]-1; k++) {
                fprintf(f, "%*d%s%*.*f%s%*.*f%s",
                        colidx_width, j[k], value_delim,
                        value_width, value_precision, a[2*k+0], value_delim,
                        value_width, value_precision, a[2*k+1], nonzero_delim);
            }
            fprintf(f, "%*d%s%*.*f%s%*.*f",
                    colidx_width, j[k], value_delim,
                    value_width, value_precision, a[2*k+0], value_delim,
                    value_width, value_precision, a[2*k+1]);
            return 0;
        }

    default:
        return EINVAL;
    }
}

/**
 * `csr_matrix_int32_from_matrix_market()` converts a matrix in the
 * matrix market format to a CSR matrix.
 */
int csr_matrix_int32_from_matrix_market(
    struct csr_matrix_int32 * matrix,
    const struct matrix_market * matrix_market,
    enum csr_value_format value_format)
{
    int err;
    if (matrix_market->object != matrix_market_matrix)
        return EINVAL;
    if (matrix_market->format != matrix_market_coordinate)
        return EINVAL;

    /* Determine the number of matrix nonzeros. */
    int64_t num_nonzeros;
    err = matrix_market_num_nonzeros(
        matrix_market, &num_nonzeros);
    if (err)
        return err;

    /* 1. Allocate storage for matrix. */
    int32_t num_rows = matrix_market->num_rows;
    int32_t num_columns = matrix_market->num_columns;
    err = csr_matrix_int32_alloc(
        matrix, num_rows, num_columns, num_nonzeros, value_format);
    if (err)
        return err;

    /* 2. Zero row pointers, column indices and nonzeros. */
    int64_t * row_ptr = matrix->row_ptr;
#pragma omp parallel for
    for (int32_t i = 0; i <= num_rows; i++)
        row_ptr[i] = 0;

    int32_t * column_indices = matrix->column_indices;
#pragma omp parallel for
    for (int64_t k = 0; k < num_nonzeros; k++)
        column_indices[k] = 0;

    err = csr_matrix_int32_zero(matrix);
    if (err) {
        csr_matrix_int32_free(matrix);
        return err;
    }

    /* Check if we need to convert matrix values to a different data type. */
    bool convert_values = true;
    if ((matrix_market->field == matrix_market_real && value_format == csr_value_f32) ||
        (matrix_market->field == matrix_market_double && value_format == csr_value_f64) ||
        (matrix_market->field == matrix_market_complex && value_format == csr_value_complex32) ||
        (matrix_market->field == matrix_market_integer && value_format == csr_value_int32) ||
        (matrix_market->field == matrix_market_pattern && value_format == csr_value_binary))
    {
        convert_values = false;
    }

    if (!convert_values) {
        /* 3a. Sort column indices and nonzeros by their rows. */
        err = matrix_market_sort_nonzeros(
            matrix_market, row_ptr, column_indices, matrix->values);
        if (err) {
            csr_matrix_int32_free(matrix);
            return err;
        }
    } else {
        /* 3b.1. Allocate temporary storage for the sorted matrix market values. */
        size_t value_size;
        switch (matrix_market->field) {
        case matrix_market_real: value_size = sizeof(float); break;
        case matrix_market_double: value_size = sizeof(double); break;
        case matrix_market_complex: value_size = 2*sizeof(float); break;
        case matrix_market_integer: value_size = sizeof(int); break;
        case matrix_market_pattern: value_size = 0; break;
        default:
            csr_matrix_int32_free(matrix);
            return EINVAL;
        }
        void * values = malloc(num_nonzeros * value_size);
        if (!values) {
            csr_matrix_int32_free(matrix);
            return errno;
        }

        /* 3b.2. Sort column indices and nonzeros by their rows. */
        err = matrix_market_sort_nonzeros(
            matrix_market, row_ptr, column_indices, values);
        if (err) {
            free(values);
            csr_matrix_int32_free(matrix);
            return err;
        }

        /* 3b.3. Convert values to the desired format. */
        switch (matrix_market->field) {
        case matrix_market_real:
            switch (matrix->value_format) {
            case csr_value_binary: break;
            case csr_value_int32:
                for (int64_t k = 0; k < num_nonzeros; k++)
                    ((int32_t *) matrix->values)[k] = (int32_t)((float *) values)[k];
                break;
            case csr_value_int64:
                for (int64_t k = 0; k < num_nonzeros; k++)
                    ((int64_t *) matrix->values)[k] = (int64_t)((float *) values)[k];
                break;
            case csr_value_f64:
                for (int64_t k = 0; k < num_nonzeros; k++)
                    ((double *) matrix->values)[k] = (double)((float *) values)[k];
                break;
            case csr_value_complex32:
                for (int64_t k = 0; k < num_nonzeros; k++)
                    ((float *) matrix->values)[2*k+0] = (float)((float *) values)[k];
                break;
            default:
                return EINVAL;
            }
            break;

        case matrix_market_double:
            switch (matrix->value_format) {
            case csr_value_binary: break;
            case csr_value_int32:
                for (int64_t k = 0; k < num_nonzeros; k++)
                    ((int32_t *) matrix->values)[k] = (int32_t)((double *) values)[k];
                break;
            case csr_value_int64:
                for (int64_t k = 0; k < num_nonzeros; k++)
                    ((int64_t *) matrix->values)[k] = (int64_t)((double *) values)[k];
                break;
            case csr_value_f32:
                for (int64_t k = 0; k < num_nonzeros; k++)
                    ((float *) matrix->values)[k] = (float)((double *) values)[k];
                break;
            case csr_value_complex32:
                for (int64_t k = 0; k < num_nonzeros; k++)
                    ((float *) matrix->values)[2*k+0] = (float)((double *) values)[k];
                break;
            default:
                return EINVAL;
            }
            break;

        case matrix_market_complex:
            switch (matrix->value_format) {
            case csr_value_binary: break;
            case csr_value_int32:
                for (int64_t k = 0; k < num_nonzeros; k++)
                    ((int32_t *) matrix->values)[k] = (int32_t)((float *) values)[2*k+0];
                break;
            case csr_value_int64:
                for (int64_t k = 0; k < num_nonzeros; k++)
                    ((int64_t *) matrix->values)[k] = (int64_t)((float *) values)[2*k+0];
                break;
            case csr_value_f32:
                for (int64_t k = 0; k < num_nonzeros; k++)
                    ((float *) matrix->values)[k] = (float)((float *) values)[2*k+0];
                break;
            case csr_value_f64:
                for (int64_t k = 0; k < num_nonzeros; k++)
                    ((double *) matrix->values)[k] = (double)((float *) values)[2*k+0];
                break;
            default:
                return EINVAL;
            }
            break;

        case matrix_market_integer:
            switch (matrix->value_format) {
            case csr_value_binary: break;
            case csr_value_int64:
                for (int64_t k = 0; k < num_nonzeros; k++)
                    ((int64_t *) matrix->values)[k] = (int64_t)((int *) values)[k];
                break;
            case csr_value_f32:
                for (int64_t k = 0; k < num_nonzeros; k++)
                    ((float *) matrix->values)[k] = (float)((int *) values)[k];
                break;
            case csr_value_f64:
                for (int64_t k = 0; k < num_nonzeros; k++)
                    ((double *) matrix->values)[k] = (double)((int *) values)[k];
                break;
            case csr_value_complex32:
                for (int64_t k = 0; k < num_nonzeros; k++)
                    ((float *) matrix->values)[2*k+0] = (float)((int *) values)[k];
                break;
            default:
                return EINVAL;
            }
            break;

        case matrix_market_pattern:
            switch (matrix->value_format) {
            case csr_value_binary: break;
            case csr_value_int64:
                for (int64_t k = 0; k < num_nonzeros; k++)
                    ((int64_t *) matrix->values)[k] = 1;
                break;
            case csr_value_f32:
                for (int64_t k = 0; k < num_nonzeros; k++)
                    ((float *) matrix->values)[k] = 1.0f;
                break;
            case csr_value_f64:
                for (int64_t k = 0; k < num_nonzeros; k++)
                    ((double *) matrix->values)[k] = 1.0;
                break;
            case csr_value_complex32:
                for (int64_t k = 0; k < num_nonzeros; k++)
                    ((float *) matrix->values)[2*k+0] = 1.0f;
                break;
            default:
                return EINVAL;
            }
            break;
        }
        free(values);
    }

    return 0;
}

/**
 * `csr_matrix_int32_spmv_binary_int32_int32()` multiplies a CSR matrix
 * of binary values with a 32-bit integer source vector, resulting in
 * a 32-bit floating point destination vector.
 */
static int csr_matrix_int32_spmv_binary_int32_int32(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst)
{
    if (matrix->value_format != csr_value_binary ||
        src->value_format != vector_value_int32 ||
        dst->value_format != vector_value_int32)
        return EINVAL;
    const int64_t * p = matrix->row_ptr;
    const int32_t * j = matrix->column_indices;
    const int32_t * x = (const int32_t *) src->values;
    int32_t * y = (int32_t *) dst->values;
#pragma omp for
    for (int32_t i = 0; i < matrix->num_rows; i++) {
        float z = 0.0;
        for (int64_t k = p[i]; k < p[i+1]; k++)
            z += x[j[k]];
        y[i] += z;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv_binary_int32_f32()` multiplies a CSR matrix
 * of binary values with a 32-bit integer source vector, resulting in
 * a 32-bit floating point destination vector.
 */
static int csr_matrix_int32_spmv_binary_int32_f32(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst)
{
    if (matrix->value_format != csr_value_binary ||
        src->value_format != vector_value_int32 ||
        dst->value_format != vector_value_f32)
        return EINVAL;
    const int64_t * p = matrix->row_ptr;
    const int32_t * j = matrix->column_indices;
    const int32_t * x = (const int32_t *) src->values;
    float * y = (float *) dst->values;
#pragma omp for
    for (int32_t i = 0; i < matrix->num_rows; i++) {
        float z = 0.0;
        for (int64_t k = p[i]; k < p[i+1]; k++)
            z += x[j[k]];
        y[i] += z;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv_binary_int32_f64()` multiplies a CSR matrix
 * of binary values with a 32-bit integer source vector, resulting in
 * a 64-bit floating point destination vector.
 */
static int csr_matrix_int32_spmv_binary_int32_f64(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst)
{
    if (matrix->value_format != csr_value_binary ||
        src->value_format != vector_value_int32 ||
        dst->value_format != vector_value_f64)
        return EINVAL;
    const int64_t * p = matrix->row_ptr;
    const int32_t * j = matrix->column_indices;
    const int32_t * x = (const int32_t *) src->values;
    double * y = (double *) dst->values;
#pragma omp for
    for (int32_t i = 0; i < matrix->num_rows; i++) {
        double z = 0.0;
        for (int64_t k = p[i]; k < p[i+1]; k++)
            z += x[j[k]];
        y[i] += z;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv_binary_f32_f32()` multiplies a CSR matrix of
 * binary values with a 32-bit floating point source vector, resulting
 * in a 32-bit floating point destination vector.
 */
static int csr_matrix_int32_spmv_binary_f32_f32(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst)
{
    if (matrix->value_format != csr_value_binary ||
        src->value_format != vector_value_f32 ||
        dst->value_format != vector_value_f32)
        return EINVAL;
    const int64_t * p = matrix->row_ptr;
    const int32_t * j = matrix->column_indices;
    const float * x = (const float *) src->values;
    float * y = (float *) dst->values;
#pragma omp for
    for (int32_t i = 0; i < matrix->num_rows; i++) {
        float z = 0.0;
        for (int64_t k = p[i]; k < p[i+1]; k++)
            z += x[j[k]];
        y[i] += z;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv_binary_f32_f64()` multiplies a CSR matrix of
 * binary values with a 32-bit floating point source vector, resulting
 * in a 64-bit floating point destination vector.
 */
static int csr_matrix_int32_spmv_binary_f32_f64(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst)
{
    if (matrix->value_format != csr_value_binary ||
        src->value_format != vector_value_f32 ||
        dst->value_format != vector_value_f64)
        return EINVAL;
    const int64_t * restrict p = matrix->row_ptr;
    const int32_t * restrict j = matrix->column_indices;
    const float * restrict x = (const float * restrict) src->values;
    double * restrict y = (double * restrict) dst->values;
#pragma omp for
    for (int32_t i = 0; i < matrix->num_rows; i++) {
        double z = 0.0;
        for (int64_t k = p[i]; k < p[i+1]; k++)
            z += x[j[k]];
        y[i] += z;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv_binary_f64_f32()` multiplies a CSR matrix of
 * binary values with a 64-bit floating point source vector, resulting
 * in a 32-bit floating point destination vector.
 */
static int csr_matrix_int32_spmv_binary_f64_f32(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst)
{
    if (matrix->value_format != csr_value_binary ||
        src->value_format != vector_value_f64 ||
        dst->value_format != vector_value_f32)
        return EINVAL;
    const int64_t * restrict p = matrix->row_ptr;
    const int32_t * restrict j = matrix->column_indices;
    const double * restrict x = (const double * restrict) src->values;
    float * restrict y = (float * restrict) dst->values;
#pragma omp for
    for (int32_t i = 0; i < matrix->num_rows; i++) {
        float z = 0.0;
        for (int64_t k = p[i]; k < p[i+1]; k++)
            z += x[j[k]];
        y[i] += z;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv_binary_f64_f64()` multiplies a CSR matrix of
 * binary values with a 64-bit floating point source vector, resulting
 * in a 64-bit floating point destination vector.
 */
static int csr_matrix_int32_spmv_binary_f64_f64(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst)
{
    if (matrix->value_format != csr_value_binary ||
        src->value_format != vector_value_f64 ||
        dst->value_format != vector_value_f64)
        return EINVAL;
    const int64_t * restrict p = matrix->row_ptr;
    const int32_t * restrict j = matrix->column_indices;
    const double * restrict x = (const double * restrict) src->values;
    double * restrict y = (double * restrict) dst->values;
#pragma omp for
    for (int32_t i = 0; i < matrix->num_rows; i++) {
        double z = 0.0;
        for (int64_t k = p[i]; k < p[i+1]; k++)
            z += x[j[k]];
        y[i] += z;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv_binary()` multiplies a CSR matrix of binary
 * values with a vector.
 */
static int csr_matrix_int32_spmv_binary(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst,
    int64_t * num_flops)
{
    int err;
    if (matrix->value_format != csr_value_binary)
        return EINVAL;
    if (src->value_format == vector_value_int32 && dst->value_format == vector_value_int32) {
        err = csr_matrix_int32_spmv_binary_int32_int32(matrix, src, dst);
        if (err)
            return err;
    } else if (src->value_format == vector_value_int32 && dst->value_format == vector_value_f32) {
        err = csr_matrix_int32_spmv_binary_int32_f32(matrix, src, dst);
        if (err)
            return err;
    } else if (src->value_format == vector_value_int32 && dst->value_format == vector_value_f64) {
        err = csr_matrix_int32_spmv_binary_int32_f64(matrix, src, dst);
        if (err)
            return err;
    } else if (src->value_format == vector_value_f32 && dst->value_format == vector_value_f32) {
        err = csr_matrix_int32_spmv_binary_f32_f32(matrix, src, dst);
        if (err)
            return err;
    } else if (src->value_format == vector_value_f32 && dst->value_format == vector_value_f64) {
        err = csr_matrix_int32_spmv_binary_f32_f64(matrix, src, dst);
        if (err)
            return err;
    } else if (src->value_format == vector_value_f64 && dst->value_format == vector_value_f32) {
        err = csr_matrix_int32_spmv_binary_f64_f32(matrix, src, dst);
        if (err)
            return err;
    } else if (src->value_format == vector_value_f64 && dst->value_format == vector_value_f64) {
        err = csr_matrix_int32_spmv_binary_f64_f64(matrix, src, dst);
        if (err)
            return err;
    } else {
        return EINVAL;
    }
#pragma omp master
    if (num_flops)
        *num_flops += matrix->num_nonzeros;
    return 0;
}

/**
 * `csr_matrix_int32_spmv_int32_int32_int32()` multiplies a CSR matrix
 * of 32-bit integer values with a 32-bit integer source vector,
 * resulting in a 32-bit integer destination vector.
 */
static int csr_matrix_int32_spmv_int32_int32_int32(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst)
{
    if (matrix->value_format != csr_value_int32 ||
        src->value_format != vector_value_int32 ||
        dst->value_format != vector_value_int32)
        return EINVAL;
    const int64_t * restrict p = matrix->row_ptr;
    const int32_t * restrict j = matrix->column_indices;
    const int32_t * restrict a = (const int32_t * restrict) matrix->values;
    const int32_t * restrict x = (const int32_t * restrict) src->values;
    int32_t * restrict y = (int32_t * restrict) dst->values;
#pragma omp for
    for (int32_t i = 0; i < matrix->num_rows; i++) {
        int32_t z = 0.0;
        for (int64_t k = p[i]; k < p[i+1]; k++)
            z += a[k] * x[j[k]];
        y[i] += z;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv_int32_int32_f32()` multiplies a CSR matrix
 * of 32-bit integer values with a 32-bit integer source vector,
 * resulting in a 32-bit integer destination vector.
 */
static int csr_matrix_int32_spmv_int32_int32_f32(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst)
{
    if (matrix->value_format != csr_value_int32 ||
        src->value_format != vector_value_int32 ||
        dst->value_format != vector_value_int32)
        return EINVAL;
    const int64_t * restrict p = matrix->row_ptr;
    const int32_t * restrict j = matrix->column_indices;
    const int32_t * restrict a = (const int32_t * restrict) matrix->values;
    const int32_t * restrict x = (const int32_t * restrict) src->values;
    float * restrict y = (float * restrict) dst->values;
#pragma omp for
    for (int32_t i = 0; i < matrix->num_rows; i++) {
        float z = 0.0;
        for (int64_t k = p[i]; k < p[i+1]; k++)
            z += a[k] * x[j[k]];
        y[i] += z;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv_int32_int32_f64()` multiplies a CSR matrix
 * of 32-bit integer values with a 32-bit integer source vector,
 * resulting in a 32-bit floating point destination vector.
 */
static int csr_matrix_int32_spmv_int32_int32_f64(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst)
{
    if (matrix->value_format != csr_value_int32 ||
        src->value_format != vector_value_int32 ||
        dst->value_format != vector_value_int32)
        return EINVAL;
    const int64_t * restrict p = matrix->row_ptr;
    const int32_t * restrict j = matrix->column_indices;
    const int32_t * restrict a = (const int32_t * restrict) matrix->values;
    const int32_t * restrict x = (const int32_t * restrict) src->values;
    double * restrict y = (double * restrict) dst->values;
#pragma omp for
    for (int32_t i = 0; i < matrix->num_rows; i++) {
        double z = 0.0;
        for (int64_t k = p[i]; k < p[i+1]; k++)
            z += a[k] * x[j[k]];
        y[i] += z;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv_int32_f32_f32()` multiplies a CSR matrix of
 * 32-bit integer values with a 32-bit floating point source vector,
 * resulting in a 32-bit floating point destination vector.
 */
static int csr_matrix_int32_spmv_int32_f32_f32(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst)
{
    if (matrix->value_format != csr_value_int32 ||
        src->value_format != vector_value_f32 ||
        dst->value_format != vector_value_f32)
        return EINVAL;
    const int64_t * restrict p = matrix->row_ptr;
    const int32_t * restrict j = matrix->column_indices;
    const int32_t * restrict a = (const int32_t * restrict) matrix->values;
    const float * restrict x = (const float * restrict) src->values;
    float * restrict y = (float * restrict) dst->values;
#pragma omp for
    for (int32_t i = 0; i < matrix->num_rows; i++) {
        float z = 0.0;
        for (int64_t k = p[i]; k < p[i+1]; k++)
            z += a[k] * x[j[k]];
        y[i] += z;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv_int32_f32_f64()` multiplies a CSR matrix of
 * 32-bit integer values with a 32-bit floating point source vector,
 * resulting in a 64-bit floating point destination vector.
 */
static int csr_matrix_int32_spmv_int32_f32_f64(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst)
{
    if (matrix->value_format != csr_value_int32 ||
        src->value_format != vector_value_f32 ||
        dst->value_format != vector_value_f64)
        return EINVAL;
    const int64_t * restrict p = matrix->row_ptr;
    const int32_t * restrict j = matrix->column_indices;
    const int32_t * restrict a = (const int32_t * restrict) matrix->values;
    const float * restrict x = (const float * restrict) src->values;
    double * restrict y = (double * restrict) dst->values;
#pragma omp for
    for (int32_t i = 0; i < matrix->num_rows; i++) {
        double z = 0.0;
        for (int64_t k = p[i]; k < p[i+1]; k++)
            z += a[k] * x[j[k]];
        y[i] += z;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv_int32_f64_f32()` multiplies a CSR matrix of
 * 32-bit integer values with a 64-bit floating point source vector,
 * resulting in a 32-bit floating point destination vector.
 */
static int csr_matrix_int32_spmv_int32_f64_f32(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst)
{
    if (matrix->value_format != csr_value_int32 ||
        src->value_format != vector_value_f64 ||
        dst->value_format != vector_value_f32)
        return EINVAL;
    const int64_t * restrict p = matrix->row_ptr;
    const int32_t * restrict j = matrix->column_indices;
    const int32_t * restrict a = (const int32_t * restrict) matrix->values;
    const double * restrict x = (const double * restrict) src->values;
    float * restrict y = (float * restrict) dst->values;
#pragma omp for
    for (int32_t i = 0; i < matrix->num_rows; i++) {
        float z = 0.0;
        for (int64_t k = p[i]; k < p[i+1]; k++)
            z += a[k] * x[j[k]];
        y[i] += z;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv_int32_f64_f64()` multiplies a CSR matrix of
 * 32-bit integer values with a 64-bit floating point source vector,
 * resulting in a 64-bit floating point destination vector.
 */
static int csr_matrix_int32_spmv_int32_f64_f64(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst)
{
    if (matrix->value_format != csr_value_int32 ||
        src->value_format != vector_value_f64 ||
        dst->value_format != vector_value_f64)
        return EINVAL;
    const int64_t * restrict p = matrix->row_ptr;
    const int32_t * restrict j = matrix->column_indices;
    const int32_t * restrict a = (const int32_t * restrict) matrix->values;
    const double * restrict x = (const double * restrict) src->values;
    double * restrict y = (double * restrict) dst->values;
#pragma omp for
    for (int32_t i = 0; i < matrix->num_rows; i++) {
        double z = 0.0;
        for (int64_t k = p[i]; k < p[i+1]; k++)
            z += a[k] * x[j[k]];
        y[i] += z;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv_int32()` multiplies a CSR matrix of 32-bit
 * integer values with a vector.
 */
static int csr_matrix_int32_spmv_int32(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst,
    int64_t * num_flops)
{
    int err;
    if (matrix->value_format != csr_value_int32)
        return EINVAL;

    if (src->value_format == vector_value_int32 && dst->value_format == vector_value_int32) {
        err = csr_matrix_int32_spmv_int32_int32_int32(matrix, src, dst);
        if (err)
            return err;
    } else if (src->value_format == vector_value_int32 && dst->value_format == vector_value_f32) {
        err = csr_matrix_int32_spmv_int32_int32_f32(matrix, src, dst);
        if (err)
            return err;
    } else if (src->value_format == vector_value_int32 && dst->value_format == vector_value_f64) {
        err = csr_matrix_int32_spmv_int32_int32_f64(matrix, src, dst);
        if (err)
            return err;
    } else if (src->value_format == vector_value_f32 && dst->value_format == vector_value_f32) {
        err = csr_matrix_int32_spmv_int32_f32_f32(matrix, src, dst);
        if (err)
            return err;
    } else if (src->value_format == vector_value_f32 && dst->value_format == vector_value_f64) {
        err = csr_matrix_int32_spmv_int32_f32_f64(matrix, src, dst);
        if (err)
            return err;
    } else if (src->value_format == vector_value_f64 && dst->value_format == vector_value_f32) {
        err = csr_matrix_int32_spmv_int32_f64_f32(matrix, src, dst);
        if (err)
            return err;
    } else if (src->value_format == vector_value_f64 && dst->value_format == vector_value_f64) {
        err = csr_matrix_int32_spmv_int32_f64_f64(matrix, src, dst);
        if (err)
            return err;
    } else {
        return EINVAL;
    }
#pragma omp master
    if (num_flops)
        *num_flops += 2*matrix->num_nonzeros;
    return 0;
}

/**
 * `csr_matrix_int32_spmv_int64_f32_f32()` multiplies a CSR matrix of
 * 64-bit integer values with a 32-bit floating point source vector,
 * resulting in a 32-bit floating point destination vector.
 */
static int csr_matrix_int32_spmv_int64_f32_f32(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst)
{
    if (matrix->value_format != csr_value_int64 ||
        src->value_format != vector_value_f32 ||
        dst->value_format != vector_value_f32)
        return EINVAL;
    const int64_t * restrict p = matrix->row_ptr;
    const int32_t * restrict j = matrix->column_indices;
    const int64_t * restrict a = (const int64_t * restrict) matrix->values;
    const float * restrict x = (const float * restrict) src->values;
    float * restrict y = (float * restrict) dst->values;
#pragma omp for
    for (int32_t i = 0; i < matrix->num_rows; i++) {
        float z = 0.0;
        for (int64_t k = p[i]; k < p[i+1]; k++)
            z += a[k] * x[j[k]];
        y[i] += z;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv_int64_f32_f64()` multiplies a CSR matrix of
 * 64-bit integer values with a 32-bit floating point source vector,
 * resulting in a 64-bit floating point destination vector.
 */
static int csr_matrix_int32_spmv_int64_f32_f64(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst)
{
    if (matrix->value_format != csr_value_int64 ||
        src->value_format != vector_value_f32 ||
        dst->value_format != vector_value_f64)
        return EINVAL;
    const int64_t * restrict p = matrix->row_ptr;
    const int32_t * restrict j = matrix->column_indices;
    const int64_t * restrict a = (const int64_t * restrict) matrix->values;
    const float * restrict x = (const float * restrict) src->values;
    double * restrict y = (double * restrict) dst->values;
#pragma omp for
    for (int32_t i = 0; i < matrix->num_rows; i++) {
        double z = 0.0;
        for (int64_t k = p[i]; k < p[i+1]; k++)
            z += a[k] * x[j[k]];
        y[i] += z;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv_int64_f64_f32()` multiplies a CSR matrix of
 * 64-bit integer values with a 64-bit floating point source vector,
 * resulting in a 32-bit floating point destination vector.
 */
static int csr_matrix_int32_spmv_int64_f64_f32(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst)
{
    if (matrix->value_format != csr_value_int64 ||
        src->value_format != vector_value_f64 ||
        dst->value_format != vector_value_f32)
        return EINVAL;
    const int64_t * restrict p = matrix->row_ptr;
    const int32_t * restrict j = matrix->column_indices;
    const int64_t * restrict a = (const int64_t * restrict) matrix->values;
    const double * restrict x = (const double * restrict) src->values;
    float * restrict y = (float * restrict) dst->values;
#pragma omp for
    for (int32_t i = 0; i < matrix->num_rows; i++) {
        float z = 0.0;
        for (int64_t k = p[i]; k < p[i+1]; k++)
            z += a[k] * x[j[k]];
        y[i] += z;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv_int64_f64_f64()` multiplies a CSR matrix of
 * 64-bit integer values with a 64-bit floating point source vector,
 * resulting in a 64-bit floating point destination vector.
 */
static int csr_matrix_int32_spmv_int64_f64_f64(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst)
{
    if (matrix->value_format != csr_value_int64 ||
        src->value_format != vector_value_f64 ||
        dst->value_format != vector_value_f64)
        return EINVAL;
    const int64_t * restrict p = matrix->row_ptr;
    const int32_t * restrict j = matrix->column_indices;
    const int64_t * restrict a = (const int64_t * restrict) matrix->values;
    const double * restrict x = (const double * restrict) src->values;
    double * restrict y = (double * restrict) dst->values;
#pragma omp for
    for (int32_t i = 0; i < matrix->num_rows; i++) {
        double z = 0.0;
        for (int64_t k = p[i]; k < p[i+1]; k++)
            z += a[k] * x[j[k]];
        y[i] += z;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv_int64()` multiplies a CSR matrix of 64-bit
 * integer values.
 */
static int csr_matrix_int32_spmv_int64(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst,
    int64_t * num_flops)
{
    int err;
    if (matrix->value_format != csr_value_int64)
        return EINVAL;
    if (src->value_format == vector_value_f32 && dst->value_format == vector_value_f32) {
        err = csr_matrix_int32_spmv_int64_f32_f32(matrix, src, dst);
        if (err)
            return err;
    } else if (src->value_format == vector_value_f32 && dst->value_format == vector_value_f64) {
        err = csr_matrix_int32_spmv_int64_f32_f64(matrix, src, dst);
        if (err)
            return err;
    } else if (src->value_format == vector_value_f64 && dst->value_format == vector_value_f32) {
        err = csr_matrix_int32_spmv_int64_f64_f32(matrix, src, dst);
        if (err)
            return err;
    } else if (src->value_format == vector_value_f64 && dst->value_format == vector_value_f64) {
        err = csr_matrix_int32_spmv_int64_f64_f64(matrix, src, dst);
        if (err)
            return err;
    } else {
        return EINVAL;
    }
#pragma omp master
    if (num_flops)
        *num_flops += 2*matrix->num_nonzeros;
    return 0;
}

/**
 * `csr_matrix_int32_spmv_f32_f32_f32()` multiplies a CSR matrix of
 * 32-bit floating point values with a 32-bit floating point source
 * vector, resulting in a 32-bit floating point destination vector.
 */
static int csr_matrix_int32_spmv_f32_f32_f32(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst)
{
    if (matrix->value_format != csr_value_f32 ||
        src->value_format != vector_value_f32 ||
        dst->value_format != vector_value_f32)
        return EINVAL;
    const int64_t * restrict p = matrix->row_ptr;
    const int32_t * restrict j = matrix->column_indices;
    const float * restrict a = (const float * restrict) matrix->values;
    const float * restrict x = (const float * restrict) src->values;
    float * restrict y = (float * restrict) dst->values;
#pragma omp for
    for (int32_t i = 0; i < matrix->num_rows; i++) {
        float z = 0.0;
        for (int64_t k = p[i]; k < p[i+1]; k++)
            z += a[k] * x[j[k]];
        y[i] += z;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv_f32_f32_f64()` multiplies a CSR matrix of
 * 32-bit floating point values with a 32-bit floating point source
 * vector, resulting in a 64-bit floating point destination vector.
 */
static int csr_matrix_int32_spmv_f32_f32_f64(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst)
{
    if (matrix->value_format != csr_value_f32 ||
        src->value_format != vector_value_f32 ||
        dst->value_format != vector_value_f64)
        return EINVAL;
    const int64_t * restrict p = matrix->row_ptr;
    const int32_t * restrict j = matrix->column_indices;
    const float * restrict a = (const float * restrict) matrix->values;
    const float * restrict x = (const float * restrict) src->values;
    double * restrict y = (double * restrict) dst->values;
#pragma omp for
    for (int32_t i = 0; i < matrix->num_rows; i++) {
        double z = 0.0;
        for (int64_t k = p[i]; k < p[i+1]; k++)
            z += a[k] * x[j[k]];
        y[i] += z;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv_f32_f64_f32()` multiplies a CSR matrix of
 * 32-bit floating point values with a 64-bit floating point source
 * vector, resulting in a 32-bit floating point destination vector.
 */
static int csr_matrix_int32_spmv_f32_f64_f32(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst)
{
    if (matrix->value_format != csr_value_f32 ||
        src->value_format != vector_value_f64 ||
        dst->value_format != vector_value_f32)
        return EINVAL;
    const int64_t * restrict p = matrix->row_ptr;
    const int32_t * restrict j = matrix->column_indices;
    const float * restrict a = (const float * restrict) matrix->values;
    const double * restrict x = (const double * restrict) src->values;
    float * restrict y = (float * restrict) dst->values;
#pragma omp for
    for (int32_t i = 0; i < matrix->num_rows; i++) {
        float z = 0.0;
        for (int64_t k = p[i]; k < p[i+1]; k++)
            z += a[k] * x[j[k]];
        y[i] += z;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv_f32_f64_f64()` multiplies a CSR matrix of
 * 32-bit floating point values with a 64-bit floating point source
 * vector, resulting in a 64-bit floating point destination vector.
 */
static int csr_matrix_int32_spmv_f32_f64_f64(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst)
{
    if (matrix->value_format != csr_value_f32 ||
        src->value_format != vector_value_f64 ||
        dst->value_format != vector_value_f64)
        return EINVAL;
    const int64_t * restrict p = matrix->row_ptr;
    const int32_t * restrict j = matrix->column_indices;
    const float * restrict a = (const float * restrict) matrix->values;
    const double * restrict x = (const double * restrict) src->values;
    double * restrict y = (double * restrict) dst->values;
#pragma omp for
    for (int32_t i = 0; i < matrix->num_rows; i++) {
        double z = 0.0;
        for (int64_t k = p[i]; k < p[i+1]; k++)
            z += a[k] * x[j[k]];
        y[i] += z;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv_f32()` multiplies a CSR matrix of 32-bit
 * floating point values with a vector.
 */
static int csr_matrix_int32_spmv_f32(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst,
    int64_t * num_flops)
{
    int err;
    if (matrix->value_format != csr_value_f32)
        return EINVAL;
    if (src->value_format == vector_value_f32 && dst->value_format == vector_value_f32) {
        err = csr_matrix_int32_spmv_f32_f32_f32(matrix, src, dst);
        if (err)
            return err;
    } else if (src->value_format == vector_value_f32 && dst->value_format == vector_value_f64) {
        err = csr_matrix_int32_spmv_f32_f32_f64(matrix, src, dst);
        if (err)
            return err;
    } else if (src->value_format == vector_value_f64 && dst->value_format == vector_value_f32) {
        err = csr_matrix_int32_spmv_f32_f64_f32(matrix, src, dst);
        if (err)
            return err;
    } else if (src->value_format == vector_value_f64 && dst->value_format == vector_value_f64) {
        err = csr_matrix_int32_spmv_f32_f64_f64(matrix, src, dst);
        if (err)
            return err;
    } else {
        return EINVAL;
    }
#pragma omp master
    if (num_flops)
        *num_flops += 2*matrix->num_nonzeros;
    return 0;
}

/**
 * `csr_matrix_int32_spmv_f64_f32_f32()` multiplies a CSR matrix of
 * 64-bit floating point values with a 32-bit floating point source
 * vector, resulting in a 32-bit floating point destination vector.
 */
static int csr_matrix_int32_spmv_f64_f32_f32(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst)
{
    if (matrix->value_format != csr_value_f64 ||
        src->value_format != vector_value_f32 ||
        dst->value_format != vector_value_f32)
        return EINVAL;
    const int64_t * restrict p = matrix->row_ptr;
    const int32_t * restrict j = matrix->column_indices;
    const double * restrict a = (const double * restrict) matrix->values;
    const float * restrict x = (const float * restrict) src->values;
    float * restrict y = (float * restrict) dst->values;
#pragma omp for
    for (int32_t i = 0; i < matrix->num_rows; i++) {
        float z = 0.0;
        for (int64_t k = p[i]; k < p[i+1]; k++)
            z += a[k] * x[j[k]];
        y[i] += z;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv_f64_f32_f64()` multiplies a CSR matrix of
 * 64-bit floating point values with a 32-bit floating point source
 * vector, resulting in a 64-bit floating point destination vector.
 */
static int csr_matrix_int32_spmv_f64_f32_f64(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst)
{
    if (matrix->value_format != csr_value_f64 ||
        src->value_format != vector_value_f32 ||
        dst->value_format != vector_value_f64)
        return EINVAL;
    const int64_t * restrict p = matrix->row_ptr;
    const int32_t * restrict j = matrix->column_indices;
    const double * restrict a = (const double * restrict) matrix->values;
    const float * restrict x = (const float * restrict) src->values;
    double * restrict y = (double * restrict) dst->values;
#pragma omp for
    for (int32_t i = 0; i < matrix->num_rows; i++) {
        double z = 0.0;
        for (int64_t k = p[i]; k < p[i+1]; k++)
            z += a[k] * x[j[k]];
        y[i] += z;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv_f64_f64_f32()` multiplies a CSR matrix of
 * 64-bit floating point values with a 64-bit floating point source
 * vector, resulting in a 32-bit floating point destination vector.
 */
static int csr_matrix_int32_spmv_f64_f64_f32(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst)
{
    if (matrix->value_format != csr_value_f64 ||
        src->value_format != vector_value_f64 ||
        dst->value_format != vector_value_f32)
        return EINVAL;
    const int64_t * restrict p = matrix->row_ptr;
    const int32_t * restrict j = matrix->column_indices;
    const double * restrict a = (const double * restrict) matrix->values;
    const double * restrict x = (const double * restrict) src->values;
    float * restrict y = (float * restrict) dst->values;
#pragma omp for
    for (int32_t i = 0; i < matrix->num_rows; i++) {
        float z = 0.0;
        for (int64_t k = p[i]; k < p[i+1]; k++)
            z += a[k] * x[j[k]];
        y[i] += z;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv_f64_f64_f64()` multiplies a CSR matrix of
 * 64-bit floating point values with a 64-bit floating point source
 * vector, resulting in a 64-bit floating point destination vector.
 */
static int csr_matrix_int32_spmv_f64_f64_f64(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst)
{
    if (matrix->value_format != csr_value_f64 ||
        src->value_format != vector_value_f64 ||
        dst->value_format != vector_value_f64)
        return EINVAL;
    const int64_t * restrict p = matrix->row_ptr;
    const int32_t * restrict j = matrix->column_indices;
    const double * restrict a = (const double * restrict) matrix->values;
    const double * restrict x = (const double * restrict) src->values;
    double * restrict y = (double * restrict) dst->values;
#pragma omp for
    for (int32_t i = 0; i < matrix->num_rows; i++) {
        double z = 0.0;
        for (int64_t k = p[i]; k < p[i+1]; k++)
            z += a[k] * x[j[k]];
        y[i] += z;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv_f64()` multiplies a CSR matrix of 64-bit
 * floating point values with a vector.
 */
static int csr_matrix_int32_spmv_f64(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst,
    int64_t * num_flops)
{
    int err;
    if (matrix->value_format != csr_value_f64)
        return EINVAL;
    if (src->value_format == vector_value_f32 && dst->value_format == vector_value_f32) {
        err = csr_matrix_int32_spmv_f64_f32_f32(matrix, src, dst);
        if (err)
            return err;
    } else if (src->value_format == vector_value_f32 && dst->value_format == vector_value_f64) {
        err = csr_matrix_int32_spmv_f64_f32_f64(matrix, src, dst);
        if (err)
            return err;
    } else if (src->value_format == vector_value_f64 && dst->value_format == vector_value_f32) {
        err = csr_matrix_int32_spmv_f64_f64_f32(matrix, src, dst);
        if (err)
            return err;
    } else if (src->value_format == vector_value_f64 && dst->value_format == vector_value_f64) {
        err = csr_matrix_int32_spmv_f64_f64_f64(matrix, src, dst);
        if (err)
            return err;
    } else {
        return EINVAL;
    }
#pragma omp master
    if (num_flops)
        *num_flops += 2*matrix->num_nonzeros;
    return 0;
}

/**
 * `csr_matrix_int32_spmv_complex32_f32_complex32()` multiplies a CSR
 * matrix of 32-bit floating point values with a 32-bit complex
 * floating point source vector, resulting in a 32-bit complex
 * floating point destination vector.
 */
static int csr_matrix_int32_spmv_complex32_f32_complex32(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst)
{
    if (matrix->value_format != csr_value_complex32 ||
        src->value_format != vector_value_f32 ||
        dst->value_format != vector_value_complex32)
        return EINVAL;
    const int64_t * restrict p = matrix->row_ptr;
    const int32_t * restrict j = matrix->column_indices;
    const float * restrict a = (const float * restrict) matrix->values;
    const float * restrict x = (const float * restrict) src->values;
    float * restrict y = (float * restrict) dst->values;
#pragma omp for
    for (int32_t i = 0; i < matrix->num_rows; i++) {
        float zr = 0.0;
        float zi = 0.0;
        for (int64_t k = p[i]; k < p[i+1]; k++) {
            zr += a[2*k+0] * x[j[k]];
            zi += a[2*k+1] * x[j[k]];
        }
        y[2*i+0] += zr;
        y[2*i+1] += zi;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv_complex32_f64_complex32()` multiplies a CSR
 * matrix of 32-bit floating point values with a 32-bit complex
 * floating point source vector, resulting in a 32-bit complex
 * floating point destination vector.
 */
static int csr_matrix_int32_spmv_complex32_f64_complex32(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst)
{
    if (matrix->value_format != csr_value_complex32 ||
        src->value_format != vector_value_f64 ||
        dst->value_format != vector_value_complex32)
        return EINVAL;
    const int64_t * restrict p = matrix->row_ptr;
    const int32_t * restrict j = matrix->column_indices;
    const float * restrict a = (const float * restrict) matrix->values;
    const double * restrict x = (const double * restrict) src->values;
    float * restrict y = (float * restrict) dst->values;
#pragma omp for
    for (int32_t i = 0; i < matrix->num_rows; i++) {
        float zr = 0.0;
        float zi = 0.0;
        for (int64_t k = p[i]; k < p[i+1]; k++) {
            zr += a[2*k+0] * x[j[k]];
            zi += a[2*k+1] * x[j[k]];
        }
        y[2*i+0] += zr;
        y[2*i+1] += zi;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv_complex32_complex32_complex32()` multiplies
 * a CSR matrix of 32-bit complex floating point values with a 32-bit
 * complex floating point source vector, resulting in a 32-bit complex
 * floating point destination vector.
 */
static int csr_matrix_int32_spmv_complex32_complex32_complex32(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst)
{
    if (matrix->value_format != csr_value_complex32 ||
        src->value_format != vector_value_complex32 ||
        dst->value_format != vector_value_complex32)
        return EINVAL;
    const int64_t * restrict p = matrix->row_ptr;
    const int32_t * restrict j = matrix->column_indices;
    const float * restrict a = (const float * restrict) matrix->values;
    const float * restrict x = (const float * restrict) src->values;
    float * restrict y = (float * restrict) dst->values;
#pragma omp for
    for (int32_t i = 0; i < matrix->num_rows; i++) {
        float zr = 0.0;
        float zi = 0.0;
        for (int64_t k = p[i]; k < p[i+1]; k++) {
            zr += a[2*k+0] * x[2*j[k]+0] - a[2*k+1] * x[2*j[k]+1];
            zi += a[2*k+0] * x[2*j[k]+1] + a[2*k+1] * x[2*j[k]+0];
        }
        y[2*i+0] += zr;
        y[2*i+1] += zi;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv_complex32()` multiplies a CSR matrix of
 * complex 32-bit floating point values with a vector.
 */
static int csr_matrix_int32_spmv_complex32(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst,
    int64_t * num_flops)
{
    int err;
    if (matrix->value_format != csr_value_complex32)
        return EINVAL;
    if (dst->value_format != vector_value_complex32)
        return EINVAL;

    switch (src->value_format) {
    case vector_value_f32:
        err = csr_matrix_int32_spmv_complex32_f32_complex32(matrix, src, dst);
        if (err)
            return err;
#pragma omp master
        if (num_flops)
            *num_flops += 4*matrix->num_nonzeros;
        break;
    case vector_value_f64:
        err = csr_matrix_int32_spmv_complex32_f64_complex32(matrix, src, dst);
        if (err)
            return err;
#pragma omp master
        if (num_flops)
            *num_flops += 4*matrix->num_nonzeros;
        break;
    case vector_value_complex32:
        err = csr_matrix_int32_spmv_complex32_complex32_complex32(matrix, src, dst);
        if (err)
            return err;
#pragma omp master
        if (num_flops)
            *num_flops += 8*matrix->num_nonzeros;
        break;
    default:
        return EINVAL;
    }
    return 0;
}

/**
 * `csr_matrix_int32_spmv()` multiplies a CSR matrix with a vector.
 */
int csr_matrix_int32_spmv(
    const struct csr_matrix_int32 * matrix,
    const struct vector * src,
    struct vector * dst,
    int64_t * num_flops)
{
    if (matrix->num_rows != dst->num_values)
        return EINVAL;
    if (matrix->num_columns != src->num_values)
        return EINVAL;

    switch (matrix->value_format) {
    case csr_value_binary:
        return csr_matrix_int32_spmv_binary(matrix, src, dst, num_flops);
    case csr_value_int32:
        return csr_matrix_int32_spmv_int32(matrix, src, dst, num_flops);
    case csr_value_int64:
        return csr_matrix_int32_spmv_int64(matrix, src, dst, num_flops);
    case csr_value_f32:
        return csr_matrix_int32_spmv_f32(matrix, src, dst, num_flops);
    case csr_value_f64:
        return csr_matrix_int32_spmv_f64(matrix, src, dst, num_flops);
    case csr_value_complex32:
        return csr_matrix_int32_spmv_complex32(matrix, src, dst, num_flops);
    default:
        return EINVAL;
    }
}
