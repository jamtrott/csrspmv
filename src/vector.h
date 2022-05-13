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
 * Dense vectors represented as arrays.
 */

#ifndef VECTOR_H
#define VECTOR_H

#include <stdint.h>
#include <stdio.h>

struct matrix_market;

/**
 * `vector_value_format` is used to enumerate different formats used
 * for vector values.
 */
enum vector_value_format
{
    vector_value_int32,
    vector_value_f32,
    vector_value_f64,
    vector_value_complex32,
};

/**
 * `vector_value_format_str()` is a string representing a given vector
 * value format.
 */
const char * vector_value_format_str(
    enum vector_value_format value_format);

/**
 * `parse_vector_value_format()` parses a string designating a format
 * for vectors.
 *
 * On success, `parse_vector_value_format()` returns `0`. If the
 * string does not correspond to a valid format, then
 * `parse_vector_value_format()` returns `EINVAL`.
 */
int parse_vector_value_format(
    const char * s,
    enum vector_value_format * format);

/**
 * `vector` is a data structure for dense vectors.
 */
struct vector
{
    /**
     * `num_values` is the number of entries in the vector.
     */
    int32_t num_values;

    /**
     * `value_format` is the type associated with each value.
     */
    enum vector_value_format value_format;

    /**
     * `values` is an array of vector values.
     */
    void * values;
};

/**
 * `vector_alloc()` allocates a dense vector.
 */
int vector_alloc(
    struct vector * vector,
    int32_t num_values,
    enum vector_value_format value_format);

/**
 * `vector_zero()` sets entries of a vector to zero.
 */
int vector_zero(
    struct vector * vector);

/**
 * `vector_ones()` sets entries of a vector to one.
 */
int vector_ones(
    struct vector * vector);

/**
 * `vector_init()` creates a dense vector from existing data.
 */
int vector_init(
    struct vector * vector,
    int32_t num_values,
    enum vector_value_format value_format,
    void * values);

/**
 * `vector_free()` frees memory and other resources associated with a
 * vector.
 */
void vector_free(
    struct vector * vector);

/**
 * `vector_print()` prints a vector.
 */
int vector_print(
    const struct vector * vector,
    FILE * f,
    const char * delimiter,
    const char * complex_delimiter,
    int field_width,
    int precision);

/**
 * `vector_from_matrix_market()` converts a vector in the Matrix
 * Market format to a dense vector.
 */
int vector_from_matrix_market(
    struct vector * vector,
    const struct matrix_market * matrix_market);

/**
 * `vector_to_matrix_market()` converts a vector to Matrix Market
 * format.
 */
int vector_to_matrix_market(
    struct vector * vector,
    struct matrix_market * matrix_market,
    int num_comment_lines,
    const char ** comment_lines);

#endif
