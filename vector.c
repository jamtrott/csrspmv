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
 * Dense vectors represented as arrays.
 */

#include "vector.h"
#include "matrix_market.h"

#include <errno.h>

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/**
 * `vector_value_format_str()` is a string representing a given vector
 * value format.
 */
const char * vector_value_format_str(
    enum vector_value_format value_format)
{
    switch (value_format) {
    case vector_value_f32:
        return "f32";
    case vector_value_f64:
        return "f64";
    case vector_value_complex32:
        return "complex32";
    default:
        return "unknown";
    }
}

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
    enum vector_value_format * format)
{
    if (strcmp("f32", s) == 0) {
        *format = vector_value_f32;
    } else if (strcmp("f64", s) == 0) {
        *format = vector_value_f64;
    } else if (strcmp("complex32", s) == 0) {
        *format = vector_value_complex32;
    } else {
        return EINVAL;
    }
    return 0;
}

/**
 * `vector_size()` computes the storage size of a vector's data.
 */
static int vector_size(
    int32_t num_values,
    enum vector_value_format value_format,
    size_t * size)
{
    switch (value_format) {
    case vector_value_f32:
        *size = num_values * sizeof(float);
        return 0;
    case vector_value_f64:
        *size = num_values * sizeof(double);
        return 0;
    case vector_value_complex32:
        *size = num_values * 2 * sizeof(float);
        return 0;
    default:
        return EINVAL;
    }
}

/**
 * `vector_alloc()` allocates a dense vector.
 */
int vector_alloc(
    struct vector * vector,
    int32_t num_values,
    enum vector_value_format value_format)
{
    int err;
    size_t size;
    err = vector_size(num_values, value_format, &size);
    if (err)
        return err;

    void * values = malloc(size);
    if (!values)
        return errno;

    vector->num_values = num_values;
    vector->value_format = value_format;
    vector->values = values;
    return 0;
}

/**
 * `vector_zero()` sets entries of a vector to zero.
 */
int vector_zero(
    struct vector * vector)
{
    switch (vector->value_format) {
    case vector_value_f32:
#pragma omp parallel for
        for (int32_t i = 0; i < vector->num_values; i++)
            ((float *) vector->values)[i] = 0.0f;
        break;

    case vector_value_f64:
#pragma omp parallel for
        for (int32_t i = 0; i < vector->num_values; i++)
            ((double *) vector->values)[i] = 0.0;
        break;

    case vector_value_complex32:
#pragma omp parallel for
        for (int32_t i = 0; i < vector->num_values; i++)
            ((float *) vector->values)[2*i+0] = ((float *) vector->values)[2*i+1] = 0.0f;
        break;

    default:
        return EINVAL;
    }
    return 0;
}

/**
 * `vector_ones()` sets entries of a vector to one.
 */
int vector_ones(
    struct vector * vector)
{
    switch (vector->value_format) {
    case vector_value_f32:
#pragma omp parallel for
        for (int32_t i = 0; i < vector->num_values; i++)
            ((float *) vector->values)[i] = 1.0f;
        break;

    case vector_value_f64:
#pragma omp parallel for
        for (int32_t i = 0; i < vector->num_values; i++)
            ((double *) vector->values)[i] = 1.0;
        break;

    case vector_value_complex32:
#pragma omp parallel for
        for (int32_t i = 0; i < vector->num_values; i++) {
            ((float *) vector->values)[2*i+0] = 1.0f;
            ((float *) vector->values)[2*i+1] = 0.0f;
        }
        break;

    default:
        return EINVAL;
    }
    return 0;
}

/**
 * `vector_init()` creates a dense vector from existing data.
 */
int vector_init(
    struct vector * vector,
    int32_t num_values,
    enum vector_value_format value_format,
    void * values)
{
    vector->num_values = num_values;
    vector->value_format = value_format;
    vector->values = values;
    return 0;
}

/**
 * `vector_free()` destroys the given vector.
 */
void vector_free(
    struct vector * vector)
{
    free(vector->values);
}

/**
 * `vector_print()` prints a vector.
 */
int vector_print(
    const struct vector * vector,
    FILE * f,
    const char * delimiter,
    const char * complex_delimiter,
    int field_width,
    int precision)
{
    int i;
    switch (vector->value_format) {
    case vector_value_f32:
        {
            const float * x = (const float *) vector->values;
            for (i = 0; i < vector->num_values-1; i++)
                fprintf(f, "%*.*f%s", field_width, precision, x[i], delimiter);
            fprintf(f, "%*.*f", field_width, precision, x[i]);
        }
        break;

    case vector_value_f64:
        {
            const double * x = (const double *) vector->values;
            for (i = 0; i < vector->num_values-1; i++)
                fprintf(f, "%*.*lf%s", field_width, precision, x[i], delimiter);
            fprintf(f, "%*.*lf", field_width, precision, x[i]);
        }
        break;

    case vector_value_complex32:
        {
            const float * x = (const float *) vector->values;
            for (i = 0; i < vector->num_values-1; i++) {
                fprintf(f, "%*.*f%s%*.*f%s",
                        field_width, precision, x[i], complex_delimiter,
                        field_width, precision, x[i], delimiter);
            }
            fprintf(f, "%*.*f%s%*.*f",
                    field_width, precision, x[i], complex_delimiter,
                    field_width, precision, x[i]);
        }
        break;

    default:
        return EINVAL;
    }

    return 0;
}

/**
 * `vector_from_matrix_market()` converts a vector in the matrix
 * market format to a dense vector.
 */
int vector_from_matrix_market(
    struct vector * vector,
    const struct matrix_market * matrix_market)
{
    int err;
    if (matrix_market->object != matrix_market_matrix &&
        matrix_market->num_columns != 1)
        return EINVAL;
    if (matrix_market->format != matrix_market_array)
        return ENOTSUP;

    switch (matrix_market->field) {
    case matrix_market_real:
        {
            /* 1. Allocate storage for vector values. */
            int32_t num_values = matrix_market->num_rows;
            enum vector_value_format value_format = vector_value_f32;
            err = vector_alloc(vector, num_values, value_format);
            if (err)
                return err;

            /* 2. Copy vector values from matrix market data. */
            float * dst = (float *) vector->values;
            const float * src = (const float *) matrix_market->data;
#pragma omp parallel for
            for (int i = 0; i < num_values; i++)
                dst[i] = src[i];
        }
        break;

    case matrix_market_double:
        {
            /* 1. Allocate storage for vector values. */
            int32_t num_values = matrix_market->num_rows;
            enum vector_value_format value_format = vector_value_f64;
            err = vector_alloc(vector, num_values, value_format);
            if (err)
                return err;

            /* 2. Copy vector values from matrix market data. */
            double * dst = (double *) vector->values;
            const double * src = (const double *) matrix_market->data;
#pragma omp parallel for
            for (int i = 0; i < num_values; i++)
                dst[i] = src[i];
        }
        break;

    case matrix_market_complex:
        {
            /* 1. Allocate storage for vector values. */
            int32_t num_values = matrix_market->num_rows;
            enum vector_value_format value_format = vector_value_complex32;
            err = vector_alloc(vector, num_values, value_format);
            if (err)
                return err;

            /* 2. Copy vector values from matrix market data. */
            float * dst = (float *) vector->values;
            const float * src = (const float *) matrix_market->data;
#pragma omp parallel for
            for (int i = 0; i < num_values; i++) {
                dst[2*i+0] = src[2*i+0];
                dst[2*i+1] = src[2*i+1];
            }
        }
        break;

    default:
        return ENOTSUP;
    }

    return 0;
}

/**
 * `vector_to_matrix_market()` converts a vector to matrix market format.
 */
int vector_to_matrix_market(
    struct vector * vector,
    struct matrix_market * matrix_market,
    int num_comment_lines,
    const char ** comment_lines)
{
    matrix_market->object = matrix_market_vector;
    matrix_market->format = matrix_market_array;
    matrix_market->symmetry = matrix_market_general;
    
    /* Copy comment lines. */
    matrix_market->num_comment_lines = num_comment_lines;
    matrix_market->comment_lines = malloc(num_comment_lines * sizeof(char *));
    if (!matrix_market->comment_lines)
        return errno;
    for (int i = 0; i < num_comment_lines; i++)
        matrix_market->comment_lines[i] = strdup(comment_lines[i]);

    matrix_market->num_rows = vector->num_values;
    matrix_market->num_columns = 1;
    matrix_market->num_nonzeros = vector->num_values;

    switch (vector->value_format) {
    case vector_value_f32:
        {
            /* 1. Allocate storage for vector values. */
            matrix_market->field = matrix_market_real;
            matrix_market->data = malloc(vector->num_values * sizeof(float));
            if (!matrix_market->data) {
                for (int i = 0; i < num_comment_lines; i++)
                    free(matrix_market->comment_lines[i]);
                free(matrix_market->comment_lines);
                return errno;
            }

            /* 2. Copy vector values from matrix market data. */
            const float * src = (const float *) vector->values;
            float * dst = (float *) matrix_market->data;
#pragma omp parallel for
            for (int i = 0; i < vector->num_values; i++)
                dst[i] = src[i];
        }
        break;

    case vector_value_f64:
        {
            /* 1. Allocate storage for vector values. */
            matrix_market->field = matrix_market_double;
            matrix_market->data = malloc(vector->num_values * sizeof(double));
            if (!matrix_market->data) {
                for (int i = 0; i < num_comment_lines; i++)
                    free(matrix_market->comment_lines[i]);
                free(matrix_market->comment_lines);
                return errno;
            }

            /* 2. Copy vector values from matrix market data. */
            const double * src = (const double *) vector->values;
            double * dst = (double *) matrix_market->data;
#pragma omp parallel for
            for (int i = 0; i < vector->num_values; i++)
                dst[i] = src[i];
        }
        break;

    case vector_value_complex32:
        {
            /* 1. Allocate storage for vector values. */
            matrix_market->field = matrix_market_real;
            matrix_market->data = malloc(vector->num_values * 2 * sizeof(float));
            if (!matrix_market->data) {
                for (int i = 0; i < num_comment_lines; i++)
                    free(matrix_market->comment_lines[i]);
                free(matrix_market->comment_lines);
                return errno;
            }

            /* 2. Copy vector values from matrix market data. */
            const float * src = (const float *) vector->values;
            float * dst = (float *) matrix_market->data;
#pragma omp parallel for
            for (int i = 0; i < vector->num_values; i++) {
                dst[2*i+0] = src[2*i+0];
                dst[2*i+1] = src[2*i+1];
            }
        }
        break;

    default:
        for (int i = 0; i < num_comment_lines; i++)
            free(matrix_market->comment_lines[i]);
        free(matrix_market->comment_lines);
        return EINVAL;
    }

    return 0;
}
