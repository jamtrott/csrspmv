/*
 * Benchmark program for CSR SpMV
 * Copyright (C) 2023 James D. Trotter
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
 * Last modified: 2023-02-22
 *
 * Matrix market file format.
 */

#ifndef MATRIX_MARKET_H
#define MATRIX_MARKET_H

#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>

/**
 * `matrix_market_object` is used to enumerate different kinds of
 * Matrix Market objects.
 */
enum matrix_market_object
{
    matrix_market_matrix,
    matrix_market_vector
};

/**
 * `matrix_market_object_str()` is a string representing the Matrix
 * Market object type.
 */
const char * matrix_market_object_str(
    enum matrix_market_object object);

/**
 * `matrix_market_format` is used to enumerate different kinds of
 * Matrix Market formats.
 */
enum matrix_market_format
{
    matrix_market_array,     /* array of dense matrix values */
    matrix_market_coordinate /* coordinate format of sparse matrix values */
};

/**
 * `matrix_market_format_str()` is a string representing the Matrix
 * Market format type.
 */
const char * matrix_market_format_str(
    enum matrix_market_format format);

/**
 * `matrix_market_field` is used to enumerate different kinds of
 * fields for matrix values in Matrix Market files.
 */
enum matrix_market_field
{
    matrix_market_real,
    matrix_market_double,
    matrix_market_complex,
    matrix_market_integer,
    matrix_market_pattern
};

/**
 * `matrix_market_field_str()` is a string representing the Matrix
 * Market field type.
 */
const char * matrix_market_field_str(
    enum matrix_market_field field);

/**
 * `matrix_market_symmetry` is used to enumerate different kinds of
 * symmetry for matrices in Matrix Market files.
 */
enum matrix_market_symmetry
{
    matrix_market_general,
    matrix_market_symmetric,
    matrix_market_skew_symmetric,
    matrix_market_hermitian
};

/**
 * `matrix_market_symmetry_str()` is a string representing the Matrix
 * Market symmetry type.
 */
const char * matrix_market_symmetry_str(
    enum matrix_market_symmetry symmetry);

/**
 * `matrix_market_coordinate_real` represents a nonzero matrix entry
 * in a Matrix Market file with `coordinate` format and `real` field.
 */
struct matrix_market_coordinate_real
{
    int i, j; /* row and column index */
    float a;  /* nonzero value */
};

/**
 * `matrix_market_coordinate_double` represents a nonzero matrix entry
 * in a Matrix Market file with `coordinate` format and `double` field.
 */
struct matrix_market_coordinate_double
{
    int i, j; /* row and column index */
    double a; /* nonzero value */
};

/**
 * `matrix_market_coordinate_complex` represents a nonzero matrix entry
 * in a Matrix Market file with `coordinate` format and `complex` field.
 */
struct matrix_market_coordinate_complex
{
    int i, j;     /* row and column index */
    float ar, ai; /* real and imaginary parts of nonzero value */
};

/**
 * `matrix_market_coordinate_integer` represents a nonzero matrix entry
 * in a Matrix Market file with `coordinate` format and `integer` field.
 */
struct matrix_market_coordinate_integer
{
    int i, j; /* row and column index */
    int a;    /* nonzero value */
};

/**
 * `matrix_market_coordinate_pattern` represents a nonzero matrix entry
 * in a Matrix Market file with `coordinate` format and `pattern` field.
 */
struct matrix_market_coordinate_pattern
{
    int i, j; /* row and column index */
};

/**
 * `matrix_market` is a data structure for matrices in the Matrix
 * Market format.
 */
struct matrix_market
{
    /**
     * `object` is the type of Matrix Market object: `matrix` or `vector`.
     */
    enum matrix_market_object object;

    /**
     * `format` is the matrix format: `coordinate` or `array`.
     */
    enum matrix_market_format format;

    /**
     * `field` is the matrix field: `real`, `double`, `complex`,
     * `integer` or `pattern`.
     */
    enum matrix_market_field field;

    /**
     * `symmetry` is the matrix symmetry: `general`, `symmetric`,
     * `skew-symmetric`, or `hermitian`.
     */
    enum matrix_market_symmetry symmetry;

    /**
     * `num_comment_lines` is the number of comment lines.
     */
    int num_comment_lines;

    /**
     * `comment_lines` is an array containing comment lines.
     */
    char ** comment_lines;

    /**
     * `num_rows` is the number of rows in the matrix.
     */
    int num_rows;

    /**
     * `num_columns` is the number of columns in the matrix.
     */
    int num_columns;

    /**
     * `num_nonzeros` is the number of nonzero entries in the matrix.
     *
     * The number of nonzeros depends on the matrix `symmetry`:
     *
     * - If `symmetry` is `general`, then `num_nonzeros` is the number
     *   of nonzero entries.
     *
     * - If `symmetry` is `symmetric` or `hermitian`, then
     *   `num_nonzeros` is the number of nonzero entries on or below
     *   the diagonal.
     *
     * - If `symmetry` is `skew-symmetric`, then `num_nonzeros` is the
     *   number of nonzero entries below the diagonal.
     */
    int64_t num_nonzeros;

    /**
     * `data` contains data for the nonzero matrix entries.
     *
     * The storage format of nonzero values depends on the matrix
     * `format` and `field`:
     *
     *   - If `format` is `array` and `field` is `real`, `double` or
     *     `integer`, then `data` is an array of
     *     `num_rows*num_columns` values of type `float`, `double` or
     *     `integer`, respectively.
     *
     *   - If `format` is `array` and `field` is `complex`, then
     *     `data` is an array of `2*num_rows*num_columns` values of
     *     type `float`.
     *
     *   - If `format` is `coordinate` and `field` is `real`,
     *     `double`, `complex`, `integer` or `pattern`, then `data` is
     *     an array of `num_nonzeros` values of type
     *     `matrix_market_coordinate_real`,
     *     `matrix_market_coordinate_double`,
     *     `matrix_market_coordinate_complex`,
     *     `matrix_market_coordinate_integer` or
     *     `matrix_market_coordinate_pattern`, respectively.
     */
    void * data;
};

/**
 * `matrix_market_free()` frees resources associated with a matrix in
 * Matrix Market format.
 */
void matrix_market_free(
    struct matrix_market * matrix);

/**
 * `matrix_market_error` are error codes that are used for error
 * handling when parsing files in the Matrix Market format.
 */
enum matrix_market_error
{
    MATRIX_MARKET_SUCCESS = 0,               /* no error */
    MATRIX_MARKET_ERR = -1,                  /* error code provided by errno */
    MATRIX_MARKET_ERR_EOF = -2,              /* unexpected end-of-file */
    MATRIX_MARKET_ERR_LINE_TOO_LONG = -3,    /* line exceeds maximum length */
    MATRIX_MARKET_ERR_INVALID_HEADER = -4,   /* invalid header */
    MATRIX_MARKET_ERR_INVALID_OBJECT = -5,   /* invalid object type */
    MATRIX_MARKET_ERR_INVALID_FORMAT = -6,   /* invalid format */
    MATRIX_MARKET_ERR_INVALID_FIELD = -7,    /* invalid field type */
    MATRIX_MARKET_ERR_INVALID_SYMMETRY = -8, /* invalid symmetry type */
    MATRIX_MARKET_ERR_INVALID_SIZE = -9,     /* invalid size */
    MATRIX_MARKET_ERR_INVALID_DATA = -10,    /* invalid matrix data */
};

/**
 * `matrix_market_strerror()` is a string describing an error code.
 */
const char * matrix_market_strerror(int err);

/**
 * `matrix_market_read()` reads a matrix from a stream in the Matrix
 * Market format.
 */
int matrix_market_read(
    struct matrix_market * matrix,
    FILE * f,
    int * line_number,
    int * column_number);

/**
 * `matrix_market_write()` writes a matrix to a stream in the Matrix
 * Market format.
 */
int matrix_market_write(
    const struct matrix_market * matrix,
    FILE * f,
    int field_width,
    int precision);

/**
 * `matrix_market_num_nonzeros()` is the total number of nonzeros a
 * matrix in the Matrix Market format, including strictly upper
 * triangular parts of symmetric, skew-symmetric or Hermitian
 * matrices.
 */
int matrix_market_num_nonzeros(
    const struct matrix_market * matrix,
    int64_t * nonzeros,
    bool include_symmetric_part);

/**
 * `matrix_market_num_nonzeros_diagonal()` is the number of nonzeros
 * on the main diagonal of a matrix in the Matrix Market format.
 */
int matrix_market_num_nonzeros_diagonal(
    const struct matrix_market * matrix,
    int64_t * nonzeros);

/**
 * `matrix_market_nonzeros_per_row()` counts the number of nonzeros in
 * each row of a matrix in the Matrix Market format.
 *
 * If `include_strict_upper_triangular_part` is `true` and `symmetry`
 * is `symmetric`, `skew-symmetric` or `hermitian`, then nonzeros in
 * the strict upper triangular part are also counted. Conversely, if
 * `include_strict_upper_triangular_part` is `false`, then only
 * nonzeros in the lower triangular part of the matrix are counted.
 *
 * `matrix_market_nonzeros_per_row()` returns `EINVAL` if `symmetry`
 * is `general` and `include_strict_upper_triangular_part` is `false`.
 */
int matrix_market_nonzeros_per_row(
    const struct matrix_market * matrix,
    bool include_strict_upper_triangular_part,
    int64_t * nonzeros_per_row);

/**
 * `matrix_market_sort_nonzeros()` sorts the nonzeros according to
 * their rows and columns for a matrix in the Matrix Market format.
 *
 * If `symmetry` is `symmetric`, `skew-symmetric` or `hermitian`, then
 * nonzero values of the lower and strict upper triangular parts of
 * the matrix are included.
 */
int matrix_market_sort_nonzeros(
    const struct matrix_market * matrix,
    int64_t * row_ptr,
    int32_t * column_indices,
    void * values,
    bool include_symmetric_part);

#endif
