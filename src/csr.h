/*
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
 * Last modified: 2020-10-17
 *
 * Sparse matrices in the compressed sparse row (CSR) format.
 */

#ifndef CSR_H
#define CSR_H

#include <stdint.h>
#include <stdio.h>

struct matrix_market;
struct vector;

/**
 * `csr_value_format` is used to enumerate different formats used for
 * nonzero values of CSR matrices.
 */
enum csr_value_format
{
    csr_value_binary,
    csr_value_int32,
    csr_value_int64,
    csr_value_f32,
    csr_value_f64,
    csr_value_complex32
};

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
    enum csr_value_format * format);

/**
 * `csr_matrix_int32` is a data structure for sparse matrices in the
 * compressed sparse row (CSR) format with 32-bit integer indices.
 */
struct csr_matrix_int32
{
    /**
     * `num_rows` is the number of matrix rows.
     */
    int32_t num_rows;

    /**
     * `num_columns` is the number of matrix columns.
     */
    int32_t num_columns;

    /**
     * `row_ptr` is an array of offsets to the first nonzero in each
     * row `i` for `i=0,1,...,num_rows-1`, and `row_ptr[num_rows]` is
     * the total number of nonzeros, which is equal to
     * `num_nonzeros`.
     *
     * The values in the `row_ptr` array may be computed as the prefix
     * sum of the number of nonzeros in each row.
     */
    int64_t * row_ptr;

    /**
     * `num_nonzeros` is the number of nonzeros in the matrix.
     */
    int64_t num_nonzeros;

    /**
     * `column_indices` is an array containing the column index of
     * each nonzero matrix entry.
     */
    int32_t * column_indices;

    /**
     * `value_format` is the type associated with nonzero values.
     */
    enum csr_value_format value_format;

    /**
     * `values` is an array containing the value associated with each
     * nonzero matrix entry.
     */
    void * values;
};

/**
 * `csr_matrix_int32_alloc()` allocates storage for a sparse matrix in
 * CSR format.
 */
int csr_matrix_int32_alloc(
    struct csr_matrix_int32 * matrix,
    int32_t num_rows,
    int32_t num_columns,
    int64_t num_nonzeros,
    enum csr_value_format value_format);

/**
 * `csr_matrix_int32_free()` frees memory and other resources
 * associated with a matrix.
 */
void csr_matrix_int32_free(
    struct csr_matrix_int32 * matrix);

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
    void * values);

/**
 * `csr_matrix_int32_zero()` zeros the values of a matrix.
 */
int csr_matrix_int32_zero(
    struct csr_matrix_int32 * matrix);

/**
 * `csr_matrix_int32_print()` prints a matrix.
 */
int csr_matrix_int32_print(
    const struct csr_matrix_int32 * matrix,
    FILE * f,
    const char * row_delimiter,
    const char * nonzero_delimiter,
    const char * column_index_and_value_separator,
    int column_index_field_width,
    int value_field_width,
    int value_precision);

/**
 * `csr_matrix_int32_from_matrix_market()` converts a matrix in the
 * Matrix Market format to a CSR matrix.
 */
int csr_matrix_int32_from_matrix_market(
    struct csr_matrix_int32 * matrix,
    const struct matrix_market * matrix_market,
    enum csr_value_format value_format);

/**
 * `csr_matrix_int32_spmv()` multiplies a CSR matrix with a vector.
 */
int csr_matrix_int32_spmv(
    const struct csr_matrix_int32 * matrix,
    const struct vector * x,
    struct vector * y,
    int64_t * num_flops);

#endif
