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
 * Matrix market file format.
 */

#include "matrix_market.h"
#include "parse.h"

#include <errno.h>
#include <unistd.h>

#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * `matrix_market_object_str()` is a string representing the matrix
 * market object type.
 */
const char * matrix_market_object_str(
    enum matrix_market_object object)
{
    switch (object) {
    case matrix_market_matrix: return "matrix";
    case matrix_market_vector: return "vector";
    default: return "unknown";
    }
}

/**
 * `matrix_market_format_str()` is a string representing the matrix
 * market format type.
 */
const char * matrix_market_format_str(
    enum matrix_market_format format)
{
    switch (format) {
    case matrix_market_array: return "array";
    case matrix_market_coordinate: return "coordinate";
    default: return "unknown";
    }
}

/**
 * `matrix_market_field_str()` is a string representing the matrix
 * market field type.
 */
const char * matrix_market_field_str(
    enum matrix_market_field field)
{
    switch (field) {
    case matrix_market_real: return "real";
    case matrix_market_double: return "double";
    case matrix_market_complex: return "complex";
    case matrix_market_integer: return "integer";
    case matrix_market_pattern: return "pattern";
    default: return "unknown";
    }
}

/**
 * `matrix_market_symmetry_str()` is a string representing the matrix
 * market symmetry type.
 */
const char * matrix_market_symmetry_str(
    enum matrix_market_symmetry symmetry)
{
    switch (symmetry) {
    case matrix_market_general: return "general";
    case matrix_market_symmetric: return "symmetric";
    case matrix_market_skew_symmetric: return "skew-symmetric";
    case matrix_market_hermitian: return "hermitian";
    default: return "unknown";
    }
}

/**
 * `matrix_market_free()` frees resources associated with matrix
 * market data.
 */
void matrix_market_free(
    struct matrix_market * matrix)
{
    free(matrix->data);
    for (int i = 0; i < matrix->num_comment_lines; i++)
        free(matrix->comment_lines[i]);
    free(matrix->comment_lines);
}

/**
 * `matrix_market_strerror()` is a string containing a description of
 * the given error code.
 */
const char * matrix_market_strerror(int err)
{
    switch (err) {
    case MATRIX_MARKET_SUCCESS:
        return "Success";
    case MATRIX_MARKET_ERR:
        return strerror(errno);
    case MATRIX_MARKET_ERR_EOF:
        if (errno) return strerror(errno);
        else return "Unexpected end-of-file";
    case MATRIX_MARKET_ERR_LINE_TOO_LONG:
        return "Maximum line length exceeded";
    case MATRIX_MARKET_ERR_INVALID_HEADER:
        return "Failed to parse header";
    case MATRIX_MARKET_ERR_INVALID_OBJECT:
        return "Failed to parse header - invalid object";
    case MATRIX_MARKET_ERR_INVALID_FORMAT:
        return "Failed to parse header - invalid format";
    case MATRIX_MARKET_ERR_INVALID_FIELD:
        return "Failed to parse header - invalid field";
    case MATRIX_MARKET_ERR_INVALID_SYMMETRY:
        return "Failed to parse header - invalid symmetry";
    case MATRIX_MARKET_ERR_INVALID_SIZE:
        return "Failed to parse header - invalid size";
    case MATRIX_MARKET_ERR_INVALID_DATA:
        return "Failed to parse data";
    default:
        return "Unknown error";
    }
}

/**
 * `matrix_market_read_line()` reads a single line from a stream in
 * the matrix market format.
 */
static int matrix_market_read_line(
    FILE * f,
    size_t line_max,
    char * linebuf)
{
    char * s = fgets(linebuf, line_max+1, f);
    if (!s && feof(f))
        return MATRIX_MARKET_ERR_EOF;
    else if (!s)
        return MATRIX_MARKET_ERR;
    else if (*s == '\0' || s[strlen(s)-1] != '\n')
        return MATRIX_MARKET_ERR_LINE_TOO_LONG;
    return 0;
}

/**
 * `matrix_market_read_header_line()` reads a header line from a
 * stream in the matrix market file format.
 */
static int matrix_market_read_header_line(
    enum matrix_market_object * object,
    enum matrix_market_format * format,
    enum matrix_market_field * field,
    enum matrix_market_symmetry * symmetry,
    FILE * f,
    size_t line_max,
    char * linebuf,
    int * line_number,
    int * column_number)
{
    int err;

    /* 1. Read the header line. */
    err = matrix_market_read_line(f, line_max, linebuf);
    if (err)
        return err;

    /* 2. Parse the identifier. */
    char * s = strtok(linebuf, " ");
    if (!s || strcmp("%%MatrixMarket", s) != 0)
        return MATRIX_MARKET_ERR_INVALID_HEADER;
    *column_number += strlen(s)+1;

    /* 3. Parse the object type. */
    s = strtok(NULL, " ");
    if (s && strcmp("matrix", s) == 0)
        *object = matrix_market_matrix;
    else if (s && strcmp("vector", s) == 0)
        *object = matrix_market_vector;
    else return MATRIX_MARKET_ERR_INVALID_OBJECT;
    *column_number += strlen(s)+1;

    /* 4. Parse the format. */
    s = strtok(NULL, " ");
    if (s && strcmp("array", s) == 0)
        *format = matrix_market_array;
    else if (s && strcmp("coordinate", s) == 0)
        *format = matrix_market_coordinate;
    else return MATRIX_MARKET_ERR_INVALID_FORMAT;
    *column_number += strlen(s)+1;

    /* 5. Parse the field type. */
    s = strtok(NULL, " ");
    if (s && strcmp("real", s) == 0)
        *field = matrix_market_real;
    else if (s && strcmp("double", s) == 0)
        *field = matrix_market_double;
    else if (s && strcmp("complex", s) == 0)
        *field = matrix_market_complex;
    else if (s && strcmp("integer", s) == 0)
        *field = matrix_market_integer;
    else if (s && strcmp("pattern", s) == 0)
        *field = matrix_market_pattern;
    else return MATRIX_MARKET_ERR_INVALID_FIELD;
    *column_number += strlen(s)+1;

    /* 6. Parse the symmetry type. */
    s = strtok(NULL, "\n");
    if (s && strcmp("general", s) == 0)
        *symmetry = matrix_market_general;
    else if (s && strcmp("symmetric", s) == 0)
        *symmetry = matrix_market_symmetric;
    else if (s && strcmp("skew-symmetric", s) == 0)
        *symmetry = matrix_market_skew_symmetric;
    else if (s && (strcmp("hermitian", s) == 0 || strcmp("Hermitian", s) == 0))
        *symmetry = matrix_market_hermitian;
    else return MATRIX_MARKET_ERR_INVALID_SYMMETRY;
    (*line_number)++; *column_number = 1;
    return 0;
}

/**
 * `comment_line_list` is a linked list data structure used for
 * parsing comment lines.
 */
struct comment_line_list
{
    char * comment_line;
    struct comment_line_list * prev;
    struct comment_line_list * next;
};

/**
 * `matrix_market_read_comment_lines()` reads comment lines from a
 * stream in the matrix market file format.
 */
static int matrix_market_read_comment_lines(
    int * num_comment_lines,
    char *** comment_lines,
    FILE * f,
    size_t line_max,
    char * linebuf,
    int * line_number,
    int * column_number)
{
    int err;

    /* 1. Read comment lines into a list. */
    struct comment_line_list * root = NULL;
    struct comment_line_list * node = NULL;
    *num_comment_lines = 0;
    while (true) {
        /* 1.1. Check if the line starts with '%'. */
        int c = fgetc(f);
        if (ungetc(c, f) == EOF || c != '%')
            break;

        /* Allocate a list node. */
        struct comment_line_list * next =
            malloc(sizeof(struct comment_line_list));
        if (!next) {
            while (node) {
                struct comment_line_list * prev = node->prev;
                free(node->comment_line);
                free(node);
                node = prev;
            }
            return MATRIX_MARKET_ERR;
        }
        next->comment_line = NULL;
        next->prev = next->next = NULL;

        /* 1.2. Read the next line as a comment line. */
        err = matrix_market_read_line(f, line_max, linebuf);
        if (err) {
            free(next->comment_line);
            free(next);
            while (node) {
                struct comment_line_list * prev = node->prev;
                free(node->comment_line);
                free(node);
                node = prev;
            }
            return err;
        }

        /* 1.3. Add the new node to the list. */
        next->comment_line = strdup(linebuf);
        if (!node) {
            root = node = next;
        } else {
            next->prev = node;
            node->next = next;
            node = node->next;
        }

        (*line_number)++; *column_number = 1;
        (*num_comment_lines)++;
    }

    /* 2. Allocate storage for comment lines. */
    *comment_lines = malloc(*num_comment_lines * sizeof(char *));
    if (!*comment_lines) {
        while (node) {
            struct comment_line_list * prev = node->prev;
            free(node->comment_line);
            free(node);
            node = prev;
        }
        return MATRIX_MARKET_ERR;
    }

    /* 3. Initialise the array of comment lines. */
    for (int i = 0; i < *num_comment_lines; i++) {
        (*comment_lines)[i] = root->comment_line;
        root = root->next;
    }

    /* 4. Clean up the list. */
    while (node) {
        struct comment_line_list * prev = node->prev;
        free(node);
        node = prev;
    }
    return 0;
}

/**
 * `matrix_market_read_size_line()` reads a size line from a
 * stream in the matrix market file format.
 */
static int matrix_market_read_size_line(
    enum matrix_market_object object,
    enum matrix_market_format format,
    int * num_rows,
    int * num_columns,
    int64_t * num_nonzeros,
    FILE * f,
    size_t line_max,
    char * linebuf,
    int * line_number,
    int * column_number)
{
    int err;

    /* 1. Read the size line. */
    err = matrix_market_read_line(f, line_max, linebuf);
    if (err)
        return err;

    const char * s = linebuf;
    if (format == matrix_market_array) {
        /* 2. Parse the number of rows. */
        err = parse_int32(s, " ", num_rows, &s);
        if (err == EINVAL) {
            return MATRIX_MARKET_ERR_INVALID_SIZE;
        } else if (err) {
            errno = err;
            return MATRIX_MARKET_ERR;
        }
        *column_number = s-linebuf+1;

        /* 3. Parse the number of columns. */
        err = parse_int32(s, "\n", num_columns, NULL);
        if (err == EINVAL) {
            return MATRIX_MARKET_ERR_INVALID_SIZE;
        } else if (err) {
            errno = err;
            return MATRIX_MARKET_ERR;
        }

        /* 4. Compute the number of nonzeros. */
        if (__builtin_mul_overflow(*num_rows, *num_columns, num_nonzeros)) {
            errno = EOVERFLOW;
            return MATRIX_MARKET_ERR;
        }
        (*line_number)++; *column_number = 1;

    } else if (format == matrix_market_coordinate) {
        /* 2. Parse the number of rows. */
        err = parse_int32(s, " ", num_rows, &s);
        if (err == EINVAL) {
            return MATRIX_MARKET_ERR_INVALID_SIZE;
        } else if (err) {
            errno = err;
            return MATRIX_MARKET_ERR;
        }
        *column_number = s-linebuf+1;

        /* 3. Parse the number of columns. */
        err = parse_int32(s, " ", num_columns, &s);
        if (err == EINVAL) {
            return MATRIX_MARKET_ERR_INVALID_SIZE;
        } else if (err) {
            errno = err;
            return MATRIX_MARKET_ERR;
        }
        *column_number = s-linebuf+1;

        /* 4. Parse the number of non-zeros. */
        err = parse_int64(s, "\n", num_nonzeros, NULL);
        if (err == EINVAL) {
            return MATRIX_MARKET_ERR_INVALID_SIZE;
        } else if (err) {
            errno = err;
            return MATRIX_MARKET_ERR;
        }
        (*line_number)++; *column_number = 1;
    }

    return 0;
}

/**
 * `matrix_market_parse_array_real()` parses a single nonzero for a
 * matrix whose format is `array` and field is `real`.
 */
static int matrix_market_parse_array_real(
    const char * s, float * a)
{
    int err = parse_float(s, "\n", a, NULL);
    if (err == EINVAL) {
        return MATRIX_MARKET_ERR_INVALID_DATA;
    } else if (err) {
        errno = err;
        return MATRIX_MARKET_ERR;
    }
    return 0;
}

/**
 * `matrix_market_parse_array_double()` parses a single nonzero for a
 * matrix whose format is `array` and field is `double`.
 */
static int matrix_market_parse_array_double(
    const char * s, double * a)
{
    int err = parse_double(s, "\n", a, NULL);
    if (err == EINVAL) {
        return MATRIX_MARKET_ERR_INVALID_DATA;
    } else if (err) {
        errno = err;
        return MATRIX_MARKET_ERR;
    }
    return 0;
}

/**
 * `matrix_market_parse_array_complex()` parses a single nonzero for a
 * matrix whose format is `array` and field is `complex`.
 */
static int matrix_market_parse_array_complex(
    const char * s, float * ar, float * ai)
{
    int err = parse_float(s, " ", ar, &s);
    if (err == EINVAL) {
        return MATRIX_MARKET_ERR_INVALID_DATA;
    } else if (err) {
        errno = err;
        return MATRIX_MARKET_ERR;
    }
    err = parse_float(s, "\n", ai, NULL);
    if (err == EINVAL) {
        return MATRIX_MARKET_ERR_INVALID_DATA;
    } else if (err) {
        errno = err;
        return MATRIX_MARKET_ERR;
    }
    return 0;
}

/**
 * `matrix_market_parse_array_integer()` parses a single nonzero for a
 * matrix whose format is `array` and field is `integer`.
 */
static int matrix_market_parse_array_integer(
    const char * s, int * a)
{
    int err = parse_int32(s, "\n", a, NULL);
    if (err == EINVAL) {
        return MATRIX_MARKET_ERR_INVALID_DATA;
    } else if (err) {
        errno = err;
        return MATRIX_MARKET_ERR;
    }
    return 0;
}

/**
 * `matrix_market_read_data_array()` reads lines of dense (array)
 * matrix data from a stream in the matrix market file format.
 */
static int matrix_market_read_data_array(
    enum matrix_market_field field,
    int num_rows,
    int num_columns,
    void ** out_data,
    FILE * f,
    size_t line_max,
    char * linebuf,
    int * line_number,
    int * column_number)
{
    int err;

    if (field == matrix_market_real) {
        /* 1. Allocate storage for matrix data. */
        float * data = (float *) malloc(
            num_rows * num_columns * sizeof(float));
        if (!data)
            return MATRIX_MARKET_ERR;

        /* 2. Read each line of data. */
        for (int i = 0; i < num_rows; i++) {
            for (int j = 0; j < num_columns; j++) {
                err = matrix_market_read_line(f, line_max, linebuf);
                if (err) {
                    free(data);
                    return err;
                }
                err = matrix_market_parse_array_real(
                    linebuf, &data[i*num_columns+j]);
                if (err) {
                    free(data);
                    return err;
                }
                (*line_number)++; *column_number = 1;
            }
        }
        *out_data = (void *) data;

    } else if (field == matrix_market_double) {
        /* 1. Allocate storage for matrix data. */
        double * data = (double *) malloc(
            num_rows * num_columns * sizeof(double));
        if (!data)
            return MATRIX_MARKET_ERR;

        /* 2. Read each line of data. */
        for (int i = 0; i < num_rows; i++) {
            for (int j = 0; j < num_columns; j++) {
                err = matrix_market_read_line(f, line_max, linebuf);
                if (err) {
                    free(data);
                    return err;
                }
                err = matrix_market_parse_array_double(
                    linebuf, &data[i*num_columns+j]);
                if (err) {
                    free(data);
                    return err;
                }
                (*line_number)++; *column_number = 1;
            }
        }
        *out_data = (void *) data;

    } else if (field == matrix_market_complex) {
        /* 1. Allocate storage for matrix data. */
        float * data = (float *) malloc(
            num_rows * num_columns * 2 * sizeof(float));
        if (!data)
            return MATRIX_MARKET_ERR;

        /* 2. Read each line of data. */
        for (int i = 0; i < num_rows; i++) {
            for (int j = 0; j < num_columns; j++) {
                err = matrix_market_read_line(f, line_max, linebuf);
                if (err) {
                    free(data);
                    return err;
                }
                err = matrix_market_parse_array_complex(
                    linebuf, &data[2*(i*num_columns+j)+0],
                    &data[2*(i*num_columns+j)+1]);
                if (err) {
                    free(data);
                    return err;
                }
                (*line_number)++; *column_number = 1;
            }
        }
        *out_data = (void *) data;

    } else if (field == matrix_market_integer) {
        /* 1. Allocate storage for matrix data. */
        int * data = (int *) malloc(
            num_rows * num_columns * sizeof(int));
        if (!data)
            return MATRIX_MARKET_ERR;

        /* 2. Read each line of data. */
        for (int i = 0; i < num_rows; i++) {
            for (int j = 0; j < num_columns; j++) {
                err = matrix_market_read_line(f, line_max, linebuf);
                if (err) {
                    free(data);
                    return err;
                }
                err = matrix_market_parse_array_integer(
                    linebuf, &data[i*num_columns+j]);
                if (err) {
                    free(data);
                    return err;
                }
                (*line_number)++; *column_number = 1;
            }
        }
        *out_data = (void *) data;

    } else {
        errno = EINVAL;
        return MATRIX_MARKET_ERR;
    }

    return 0;
}

/**
 * `matrix_market_parse_coordinate_real()` parses a single nonzero for
 * a matrix whose format is `coordinate` and field is `real`.
 */
static int matrix_market_parse_coordinate_real(
    const char * s,
    struct matrix_market_coordinate_real * data,
    int * line_number, int * column_number)
{
    const char * start = s;
    int err = parse_int32(s, " ", &data->i, &s);
    if (err == EINVAL) {
        return MATRIX_MARKET_ERR_INVALID_DATA;
    } else if (err) {
        errno = err;
        return MATRIX_MARKET_ERR;
    }
    *column_number += s-start; start = s;
    err = parse_int32(s, " ", &data->j, &s);
    if (err == EINVAL) {
        return MATRIX_MARKET_ERR_INVALID_DATA;
    } else if (err) {
        errno = err;
        return MATRIX_MARKET_ERR;
    }
    *column_number += s-start; start = s;
    err = parse_float(s, "\n", &data->a, NULL);
    if (err == EINVAL) {
        return MATRIX_MARKET_ERR_INVALID_DATA;
    } else if (err) {
        errno = err;
        return MATRIX_MARKET_ERR;
    }
    (*line_number)++; *column_number = 1;
    return 0;
}

/**
 * `matrix_market_parse_coordinate_double()` parses a single nonzero
 * for a matrix whose format is `coordinate` and field is `double`.
 */
static int matrix_market_parse_coordinate_double(
    const char * s,
    struct matrix_market_coordinate_double * data,
    int * line_number, int * column_number)
{
    const char * start = s;
    int err = parse_int32(s, " ", &data->i, &s);
    if (err == EINVAL) {
        return MATRIX_MARKET_ERR_INVALID_DATA;
    } else if (err) {
        errno = err;
        return MATRIX_MARKET_ERR;
    }
    *column_number += s-start; start = s;
    err = parse_int32(s, " ", &data->j, &s);
    if (err == EINVAL) {
        return MATRIX_MARKET_ERR_INVALID_DATA;
    } else if (err) {
        errno = err;
        return MATRIX_MARKET_ERR;
    }
    *column_number += s-start; start = s;
    err = parse_double(s, "\n", &data->a, NULL);
    if (err == EINVAL) {
        return MATRIX_MARKET_ERR_INVALID_DATA;
    } else if (err) {
        errno = err;
        return MATRIX_MARKET_ERR;
    }
    (*line_number)++; *column_number = 1;
    return 0;
}

/**
 * `matrix_market_parse_coordinate_complex()` parses a single nonzero
 * for a matrix whose format is `coordinate` and field is `complex`.
 */
static int matrix_market_parse_coordinate_complex(
    const char * s,
    struct matrix_market_coordinate_complex * data,
    int * line_number, int * column_number)
{
    const char * start = s;
    int err = parse_int32(s, " ", &data->i, &s);
    if (err == EINVAL) {
        return MATRIX_MARKET_ERR_INVALID_DATA;
    } else if (err) {
        errno = err;
        return MATRIX_MARKET_ERR;
    }
    *column_number += s-start; start = s;
    err = parse_int32(s, " ", &data->j, &s);
    if (err == EINVAL) {
        return MATRIX_MARKET_ERR_INVALID_DATA;
    } else if (err) {
        errno = err;
        return MATRIX_MARKET_ERR;
    }
    *column_number += s-start; start = s;
    err = parse_float(s, " ", &data->ar, &s);
    if (err == EINVAL) {
        return MATRIX_MARKET_ERR_INVALID_DATA;
    } else if (err) {
        errno = err;
        return MATRIX_MARKET_ERR;
    }
    *column_number += s-start; start = s;
    err = parse_float(s, "\n", &data->ai, NULL);
    if (err == EINVAL) {
        return MATRIX_MARKET_ERR_INVALID_DATA;
    } else if (err) {
        errno = err;
        return MATRIX_MARKET_ERR;
    }
    (*line_number)++; *column_number = 1;
    return 0;
}

/**
 * `matrix_market_parse_coordinate_integer()` parses a single nonzero
 * for a matrix whose format is `coordinate` and field is `integer`.
 */
static int matrix_market_parse_coordinate_integer(
    const char * s,
    struct matrix_market_coordinate_integer * data,
    int * line_number, int * column_number)
{
    const char * start = s;
    int err = parse_int32(s, " ", &data->i, &s);
    if (err == EINVAL) {
        return MATRIX_MARKET_ERR_INVALID_DATA;
    } else if (err) {
        errno = err;
        return MATRIX_MARKET_ERR;
    }
    *column_number += s-start; start = s;
    err = parse_int32(s, " ", &data->j, &s);
    if (err == EINVAL) {
        return MATRIX_MARKET_ERR_INVALID_DATA;
    } else if (err) {
        errno = err;
        return MATRIX_MARKET_ERR;
    }
    *column_number += s-start; start = s;
    err = parse_int32(s, "\n", &data->a, NULL);
    if (err == EINVAL) {
        return MATRIX_MARKET_ERR_INVALID_DATA;
    } else if (err) {
        errno = err;
        return MATRIX_MARKET_ERR;
    }
    (*line_number)++; *column_number = 1;
    return 0;
}

/**
 * `matrix_market_parse_coordinate_pattern()` parses a single nonzero
 * for a matrix whose format is `coordinate` and field is `pattern`.
 */
static int matrix_market_parse_coordinate_pattern(
    const char * s,
    struct matrix_market_coordinate_pattern * data,
    int * line_number, int * column_number)
{
    const char * start = s;
    int err = parse_int32(s, " ", &data->i, &s);
    if (err == EINVAL) {
        return MATRIX_MARKET_ERR_INVALID_DATA;
    } else if (err) {
        errno = err;
        return MATRIX_MARKET_ERR;
    }
    *column_number += s-start; start = s;
    err = parse_int32(s, "\n", &data->j, NULL);
    if (err == EINVAL) {
        return MATRIX_MARKET_ERR_INVALID_DATA;
    } else if (err) {
        errno = err;
        return MATRIX_MARKET_ERR;
    }
    (*line_number)++; *column_number = 1;
    return 0;
}

/**
 * `matrix_market_read_data_coordinate()` reads lines of sparse
 * (coordinate) matrix data from a stream in the matrix market file
 * format.
 */
static int matrix_market_read_data_coordinate(
    enum matrix_market_field field,
    int num_rows,
    int num_columns,
    int64_t num_nonzeros,
    void ** out_data,
    FILE * f,
    size_t line_max,
    char * linebuf,
    int * line_number,
    int * column_number)
{
    int err;

    if (field == matrix_market_real) {
        /* 1. Allocate storage for matrix data. */
        size_t nonzero_size = sizeof(struct matrix_market_coordinate_real);
        struct matrix_market_coordinate_real * data =
            (struct matrix_market_coordinate_real *) malloc(
                num_nonzeros * nonzero_size);
        if (!data)
            return MATRIX_MARKET_ERR;

        /* 2. Read each line of data. */
        for (int64_t k = 0; k < num_nonzeros; k++) {
            err = matrix_market_read_line(f, line_max, linebuf);
            if (err) {
                free(data);
                return err;
            }
            err = matrix_market_parse_coordinate_real(
                linebuf, &data[k], line_number, column_number);
            if (err) {
                free(data);
                return err;
            }
        }
        *out_data = (void *) data;

    } else if (field == matrix_market_double) {
        /* 1. Allocate storage for matrix data. */
        size_t nonzero_size = sizeof(struct matrix_market_coordinate_double);
        struct matrix_market_coordinate_double * data =
            (struct matrix_market_coordinate_double *) malloc(
                num_nonzeros * nonzero_size);
        if (!data)
            return MATRIX_MARKET_ERR;

        /* 2. Read each line of data. */
        for (int64_t k = 0; k < num_nonzeros; k++) {
            err = matrix_market_read_line(f, line_max, linebuf);
            if (err) {
                free(data);
                return err;
            }
            err = matrix_market_parse_coordinate_double(
                linebuf, &data[k], line_number, column_number);
            if (err) {
                free(data);
                return err;
            }
        }
        *out_data = (void *) data;

    } else if (field == matrix_market_complex) {
        /* 1. Allocate storage for matrix data. */
        size_t nonzero_size = sizeof(struct matrix_market_coordinate_complex);
        struct matrix_market_coordinate_complex * data =
            (struct matrix_market_coordinate_complex *) malloc(
                num_nonzeros * nonzero_size);
        if (!data)
            return MATRIX_MARKET_ERR;

        /* 2. Read each line of data. */
        for (int64_t k = 0; k < num_nonzeros; k++) {
            err = matrix_market_read_line(f, line_max, linebuf);
            if (err) {
                free(data);
                return err;
            }
            err = matrix_market_parse_coordinate_complex(
                linebuf, &data[k], line_number, column_number);
            if (err) {
                free(data);
                return err;
            }
        }
        *out_data = (void *) data;

    } else if (field == matrix_market_integer) {
        /* 1. Allocate storage for matrix data. */
        size_t nonzero_size = sizeof(struct matrix_market_coordinate_integer);
        struct matrix_market_coordinate_integer * data =
            (struct matrix_market_coordinate_integer *) malloc(
                num_nonzeros * nonzero_size);
        if (!data)
            return MATRIX_MARKET_ERR;

        /* 2. Read each line of data. */
        for (int64_t k = 0; k < num_nonzeros; k++) {
            err = matrix_market_read_line(f, line_max, linebuf);
            if (err) {
                free(data);
                return err;
            }
            err = matrix_market_parse_coordinate_integer(
                linebuf, &data[k], line_number, column_number);
            if (err) {
                free(data);
                return err;
            }
        }
        *out_data = (void *) data;

    } else if (field == matrix_market_pattern) {
        /* 1. Allocate storage for matrix data. */
        size_t nonzero_size = sizeof(struct matrix_market_coordinate_pattern);
        struct matrix_market_coordinate_pattern * data =
            (struct matrix_market_coordinate_pattern *) malloc(
                num_nonzeros * nonzero_size);
        if (!data)
            return MATRIX_MARKET_ERR;

        /* 2. Read each line of data. */
        for (int64_t k = 0; k < num_nonzeros; k++) {
            err = matrix_market_read_line(f, line_max, linebuf);
            if (err) {
                free(data);
                return err;
            }
            err = matrix_market_parse_coordinate_pattern(
                linebuf, &data[k], line_number, column_number);
            if (err) {
                free(data);
                return err;
            }
        }
        *out_data = (void *) data;

    } else {
        errno = EINVAL;
        return MATRIX_MARKET_ERR;
    }

    return 0;
}

/**
 * `matrix_market_read_data()` reads lines of matrix data from a
 * stream in the matrix market file format.
 */
static int matrix_market_read_data(
    enum matrix_market_object object,
    enum matrix_market_format format,
    enum matrix_market_field field,
    int num_rows,
    int num_columns,
    int64_t num_nonzeros,
    void ** data,
    FILE * f,
    size_t line_max,
    char * linebuf,
    int * line_number,
    int * column_number)
{
    int err;

    switch (format) {
    case matrix_market_array:
        err = matrix_market_read_data_array(
            field, num_rows, num_columns, data,
            f, line_max, linebuf, line_number, column_number);
        if (err)
            return err;
        break;

    case matrix_market_coordinate:
        err = matrix_market_read_data_coordinate(
            field, num_rows, num_columns, num_nonzeros, data,
            f, line_max, linebuf, line_number, column_number);
        if (err)
            return err;
        break;

    default:
        return EINVAL;
    }

    return 0;
}

/**
 * `matrix_market_read()` reads a matrix from a stream in the matrix
 * market file format.
 */
int matrix_market_read(
    struct matrix_market * matrix,
    FILE * f,
    int * line_number,
    int * column_number)
{
    int err;

    /* Allocate storage for reading lines from file. */
    long int line_max = sysconf(_SC_LINE_MAX);
    char * linebuf = malloc(line_max+1);
    if (!linebuf)
        return errno;

    matrix->comment_lines = NULL;
    matrix->data = NULL;

    /* 1. Parse the header line. */
    *line_number = 1;
    *column_number = 1;
    err = matrix_market_read_header_line(
        &matrix->object, &matrix->format,
        &matrix->field, &matrix->symmetry,
        f, line_max, linebuf,
        line_number, column_number);
    if (err) {
        free(linebuf);
        return err;
    }

    /* 2. Parse comment lines. */
    err = matrix_market_read_comment_lines(
        &matrix->num_comment_lines, &matrix->comment_lines,
        f, line_max, linebuf, line_number, column_number);
    if (err) {
        free(linebuf);
        return err;
    }

    /* 3. Parse the size line. */
    err = matrix_market_read_size_line(
        matrix->object, matrix->format,
        &matrix->num_rows, &matrix->num_columns, &matrix->num_nonzeros,
        f, line_max, linebuf, line_number, column_number);
    if (err) {
        free(matrix->comment_lines);
        matrix->comment_lines = NULL;
        free(linebuf);
        return err;
    }

    /* 4. Parse the data. */
    err = matrix_market_read_data(
        matrix->object, matrix->format, matrix->field,
        matrix->num_rows, matrix->num_columns, matrix->num_nonzeros,
        &matrix->data, f, line_max, linebuf,
        line_number, column_number);
    if (err) {
        free(matrix->comment_lines);
        matrix->comment_lines = NULL;
        free(linebuf);
        return err;
    }

    free(linebuf);
    return MATRIX_MARKET_SUCCESS;
}

/**
 * `matrix_market_write_matrix()` writes a matrix to a stream in the
 * matrix market format.
 */
static int matrix_market_write_matrix(
    const struct matrix_market * matrix,
    FILE * f,
    int field_width,
    int precision)
{
    if (matrix->object != matrix_market_matrix)
        return EINVAL;

    /* 1. Write the header line. */
    fprintf(f, "%%%%MatrixMarket %s %s %s %s\n",
            matrix_market_object_str(matrix->object),
            matrix_market_format_str(matrix->format),
            matrix_market_field_str(matrix->field),
            matrix_market_symmetry_str(matrix->symmetry));

    /* 2. Write comment lines. */
    for (int i = 0; i < matrix->num_comment_lines; i++)
        fputs(matrix->comment_lines[i], f);

    /* 3. Write the size line. */
    if (matrix->format == matrix_market_array)
        fprintf(f, "%d %d\n", matrix->num_rows, matrix->num_columns);
    else if (matrix->format == matrix_market_coordinate)
        fprintf(f, "%d %d %"PRId64"\n", matrix->num_rows, matrix->num_columns, matrix->num_nonzeros);
    else return EINVAL;

    /* 4. Write the data. */
    if (matrix->format == matrix_market_array) {
        if (matrix->field == matrix_market_real) {
            const float * a = (const float *) matrix->data;
            for (int i = 0; i < matrix->num_rows; i++) {
                for (int j = 0; j < matrix->num_columns; j++) {
                    fprintf(f, "%*.*f\n", field_width, precision,
                            a[i*matrix->num_columns+j]);
                }
            }
        } else if (matrix->field == matrix_market_double) {
            const double * a = (const double *) matrix->data;
            for (int i = 0; i < matrix->num_rows; i++) {
                for (int j = 0; j < matrix->num_columns; j++) {
                    fprintf(f, "%*.*lf\n", field_width, precision,
                            a[i*matrix->num_columns+j]);
                }
            }
        } else if (matrix->field == matrix_market_complex) {
            const float * a = (const float *) matrix->data;
            for (int i = 0; i < matrix->num_rows; i++) {
                for (int j = 0; j < matrix->num_columns; j++) {
                    fprintf(f, "%*.*f %*.*f\n",
                            field_width, precision, a[2*(i*matrix->num_columns+j)+0],
                            field_width, precision, a[2*(i*matrix->num_columns+j)+1]);
                }
            }
        } else if (matrix->field == matrix_market_integer) {
            const int * a = (const int *) matrix->data;
            for (int i = 0; i < matrix->num_rows; i++) {
                for (int j = 0; j < matrix->num_columns; j++) {
                    fprintf(f, "%*d\n", field_width,
                            a[i*matrix->num_columns+j]);
                }
            }
        } else {
            return EINVAL;
        }

    } else if (matrix->format == matrix_market_coordinate) {
        if (matrix->field == matrix_market_real) {
            const struct matrix_market_coordinate_real * a =
                (const struct matrix_market_coordinate_real *) matrix->data;
            for (int64_t k = 0; k < matrix->num_nonzeros; k++) {
                fprintf(f, "%d %d %*.*f\n", a[k].i, a[k].j,
                        field_width, precision, a[k].a);
            }
        } else if (matrix->field == matrix_market_double) {
            const struct matrix_market_coordinate_double * a =
                (const struct matrix_market_coordinate_double *) matrix->data;
            for (int64_t k = 0; k < matrix->num_nonzeros; k++) {
                fprintf(f, "%d %d %*.*lf\n", a[k].i, a[k].j,
                        field_width, precision, a[k].a);
            }
        } else if (matrix->field == matrix_market_complex) {
            const struct matrix_market_coordinate_complex * a =
                (const struct matrix_market_coordinate_complex *) matrix->data;
            for (int64_t k = 0; k < matrix->num_nonzeros; k++) {
                fprintf(f, "%d %d %*.*f %*.*f\n", a[k].i, a[k].j,
                        field_width, precision, a[k].ar,
                        field_width, precision, a[k].ai);
            }
        } else if (matrix->field == matrix_market_integer) {
            const struct matrix_market_coordinate_integer * a =
                (const struct matrix_market_coordinate_integer *) matrix->data;
            for (int64_t k = 0; k < matrix->num_nonzeros; k++) {
                fprintf(f, "%d %d %*d\n", a[k].i, a[k].j, field_width, a[k].a);
            }
        } else if (matrix->field == matrix_market_pattern) {
            const struct matrix_market_coordinate_pattern * a =
                (const struct matrix_market_coordinate_pattern *) matrix->data;
            for (int64_t k = 0; k < matrix->num_nonzeros; k++) {
                fprintf(f, "%d %d\n", a[k].i, a[k].j);
            }
        } else {
            return EINVAL;
        }

    } else {
        return EINVAL;
    }

    return 0;
}

/**
 * `matrix_market_write_vector()` writes a vector to a stream in the
 * matrix market format.
 */
int matrix_market_write_vector(
    const struct matrix_market * vector,
    FILE * f,
    int field_width,
    int precision)
{
    if (vector->object != matrix_market_vector)
        return EINVAL;
    if (vector->format == matrix_market_coordinate)
        return ENOTSUP;

    /* 1. Write the header line. */
    fprintf(f, "%%%%MatrixMarket %s %s %s %s\n",
            matrix_market_object_str(vector->object),
            matrix_market_format_str(vector->format),
            matrix_market_field_str(vector->field),
            matrix_market_symmetry_str(vector->symmetry));

    /* 2. Write comment lines. */
    for (int i = 0; i < vector->num_comment_lines; i++)
        fputs(vector->comment_lines[i], f);

    /* 3. Write the size line. */
    if (vector->format == matrix_market_array)
        fprintf(f, "%d\n", vector->num_rows);
    else return EINVAL;

    /* 4. Write the data. */
    if (vector->format == matrix_market_array) {
        if (vector->field == matrix_market_real) {
            const float * a = (const float *) vector->data;
            for (int i = 0; i < vector->num_rows; i++) {
                for (int j = 0; j < vector->num_columns; j++) {
                    fprintf(f, "%*.*f\n", field_width, precision,
                            a[i*vector->num_columns+j]);
                }
            }
        } else if (vector->field == matrix_market_double) {
            const double * a = (const double *) vector->data;
            for (int i = 0; i < vector->num_rows; i++) {
                for (int j = 0; j < vector->num_columns; j++) {
                    fprintf(f, "%*.*lf\n", field_width, precision,
                            a[i*vector->num_columns+j]);
                }
            }
        } else if (vector->field == matrix_market_complex) {
            const float * a = (const float *) vector->data;
            for (int i = 0; i < vector->num_rows; i++) {
                for (int j = 0; j < vector->num_columns; j++) {
                    fprintf(f, "%*.*f %*.*f\n",
                            field_width, precision, a[2*(i*vector->num_columns+j)+0],
                            field_width, precision, a[2*(i*vector->num_columns+j)+1]);
                }
            }
        } else if (vector->field == matrix_market_integer) {
            const int * a = (const int *) vector->data;
            for (int i = 0; i < vector->num_rows; i++) {
                for (int j = 0; j < vector->num_columns; j++) {
                    fprintf(f, "%*d\n", field_width,
                            a[i*vector->num_columns+j]);
                }
            }
        } else {
            return EINVAL;
        }

    } else {
        return EINVAL;
    }

    return 0;
}

/**
 * `matrix_market_write()` writes a matrix to a stream in the matrix
 * market format.
 */
int matrix_market_write(
    const struct matrix_market * matrix,
    FILE * f,
    int field_width,
    int precision)
{
    int err;
    switch (matrix->object) {
    case matrix_market_matrix:
        err = matrix_market_write_matrix(
            matrix, f, field_width, precision);
        if (err)
            return err;
        break;

    case matrix_market_vector:
        err = matrix_market_write_vector(
            matrix, f, field_width, precision);
        if (err)
            return err;
        break;

    default:
        return EINVAL;
    }
    return 0;
}

/**
 * `matrix_market_num_nonzeros()` is the total number of nonzeros a
 * matrix in the matrix market format, including strictly upper
 * triangular parts of symmetric, skew-symmetric or Hermitian
 * matrices.
 */
int matrix_market_num_nonzeros(
    const struct matrix_market * matrix,
    int64_t * nonzeros)
{
    int err;
    if (matrix->symmetry == matrix_market_general) {
        *nonzeros = matrix->num_nonzeros;
    } else if (matrix->symmetry == matrix_market_skew_symmetric) {
        *nonzeros = 2*matrix->num_nonzeros;
    } else if (matrix->symmetry == matrix_market_symmetric ||
               matrix->symmetry == matrix_market_hermitian)
    {
        int64_t num_diagonal_nonzeros;
        err = matrix_market_num_nonzeros_diagonal(matrix, &num_diagonal_nonzeros);
        if (err)
            return err;
        *nonzeros = 2*matrix->num_nonzeros - num_diagonal_nonzeros;
    } else {
        return EINVAL;
    }
    return 0;
}

/**
 * `matrix_market_num_nonzeros_diagonal()` is the number of nonzeros
 * on the main diagonal of a matrix in the matrix market format.
 */
int matrix_market_num_nonzeros_diagonal(
    const struct matrix_market * matrix,
    int64_t * nonzeros)
{
    if (matrix->format == matrix_market_array) {
        *nonzeros = matrix->num_columns < matrix->num_rows
            ? matrix->num_columns : matrix->num_rows;

    } else if (matrix->format == matrix_market_coordinate) {
        switch (matrix->field) {
        case matrix_market_real:
            {
                const struct matrix_market_coordinate_real * data =
                    (const struct matrix_market_coordinate_real *)
                    matrix->data;
                *nonzeros = 0;
                for (int64_t k = 0; k < matrix->num_nonzeros; k++) {
                    if (data[k].i == data[k].j)
                        (*nonzeros)++;
                }
            }
            break;
        case matrix_market_double:
            {
                const struct matrix_market_coordinate_double * data =
                    (const struct matrix_market_coordinate_double *)
                    matrix->data;
                *nonzeros = 0;
                for (int64_t k = 0; k < matrix->num_nonzeros; k++) {
                    if (data[k].i == data[k].j)
                        (*nonzeros)++;
                }
            }
            break;
        case matrix_market_complex:
            {
                const struct matrix_market_coordinate_complex * data =
                    (const struct matrix_market_coordinate_complex *)
                    matrix->data;
                *nonzeros = 0;
                for (int64_t k = 0; k < matrix->num_nonzeros; k++) {
                    if (data[k].i == data[k].j)
                        (*nonzeros)++;
                }
            }
            break;
        case matrix_market_integer:
            {
                const struct matrix_market_coordinate_integer * data =
                    (const struct matrix_market_coordinate_integer *)
                    matrix->data;
                *nonzeros = 0;
                for (int64_t k = 0; k < matrix->num_nonzeros; k++) {
                    if (data[k].i == data[k].j)
                        (*nonzeros)++;
                }
            }
            break;
        case matrix_market_pattern:
            {
                const struct matrix_market_coordinate_pattern * data =
                    (const struct matrix_market_coordinate_pattern *)
                    matrix->data;
                *nonzeros = 0;
                for (int64_t k = 0; k < matrix->num_nonzeros; k++) {
                    if (data[k].i == data[k].j)
                        (*nonzeros)++;
                }
            }
            break;
        default:
            return EINVAL;
        }

    } else {
        return EINVAL;
    }

    return 0;
}

/**
 * `matrix_market_nonzeros_per_row()` counts the number of nonzeros in
 * each row of a matrix in the matrix market format.
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
    int64_t * nonzeros_per_row)
{
    if (matrix->format == matrix_market_array) {
        for (int i = 0; i < matrix->num_rows; i++)
            nonzeros_per_row[i] += matrix->num_columns;

    } else if (matrix->format == matrix_market_coordinate) {
        switch (matrix->field) {
        case matrix_market_real:
            {
                const struct matrix_market_coordinate_real * data =
                    (const struct matrix_market_coordinate_real *)
                    matrix->data;
                if ((include_strict_upper_triangular_part &&
                     matrix->symmetry == matrix_market_general) ||
                    (!include_strict_upper_triangular_part &&
                     (matrix->symmetry == matrix_market_symmetric ||
                      matrix->symmetry == matrix_market_skew_symmetric ||
                      matrix->symmetry == matrix_market_hermitian)))
                {
                    for (int64_t k = 0; k < matrix->num_nonzeros; k++)
                        nonzeros_per_row[data[k].i-1]++;
                } else if (include_strict_upper_triangular_part &&
                           (matrix->symmetry == matrix_market_symmetric ||
                            matrix->symmetry == matrix_market_skew_symmetric ||
                            matrix->symmetry == matrix_market_hermitian))
                {
                    for (int64_t k = 0; k < matrix->num_nonzeros; k++) {
                        nonzeros_per_row[data[k].i-1]++;
                        if (data[k].i != data[k].j)
                            nonzeros_per_row[data[k].j-1]++;
                    }
                } else {
                    return EINVAL;
                }
            }
            break;
        case matrix_market_double:
            {
                const struct matrix_market_coordinate_double * data =
                    (const struct matrix_market_coordinate_double *)
                    matrix->data;
                for (int64_t k = 0; k < matrix->num_nonzeros; k++)
                    nonzeros_per_row[data[k].i-1]++;
            }
            break;
        case matrix_market_complex:
            {
                const struct matrix_market_coordinate_complex * data =
                    (const struct matrix_market_coordinate_complex *)
                    matrix->data;
                for (int64_t k = 0; k < matrix->num_nonzeros; k++)
                    nonzeros_per_row[data[k].i-1]++;
            }
            break;
        case matrix_market_integer:
            {
                const struct matrix_market_coordinate_integer * data =
                    (const struct matrix_market_coordinate_integer *)
                    matrix->data;
                for (int64_t k = 0; k < matrix->num_nonzeros; k++)
                    nonzeros_per_row[data[k].i-1]++;
            }
            break;
        case matrix_market_pattern:
            {
                const struct matrix_market_coordinate_pattern * data =
                    (const struct matrix_market_coordinate_pattern *)
                    matrix->data;
                for (int64_t k = 0; k < matrix->num_nonzeros; k++)
                    nonzeros_per_row[data[k].i-1]++;
            }
            break;
        default:
            return EINVAL;
        }

    } else {
        return EINVAL;
    }

    return 0;
}

/**
 * `matrix_market_sort_nonzeros_real()` sorts the nonzeros according
 * to their rows and columns for a matrix in the matrix market format
 * whose `field` is `real`.
 */
static int matrix_market_sort_nonzeros_real(
    const struct matrix_market * matrix,
    const struct matrix_market_coordinate_real * data,
    bool include_strict_upper_triangular_part,
    int64_t * row_ptr,
    int32_t * j,
    float * a)
{
    if (matrix->field != matrix_market_real)
        return EINVAL;

    if ((include_strict_upper_triangular_part &&
         matrix->symmetry == matrix_market_general) ||
        (!include_strict_upper_triangular_part &&
         (matrix->symmetry == matrix_market_symmetric ||
          matrix->symmetry == matrix_market_skew_symmetric ||
          matrix->symmetry == matrix_market_hermitian)))
    {
        for (int64_t k = 0; k < matrix->num_nonzeros; k++) {
            int i = data[k].i-1;
            int64_t l = row_ptr[i+1]-1;
            while (l >= row_ptr[i] && j[l] > data[k].j-1) {
                j[l+1] = j[l];
                a[l+1] = a[l];
                l--;
            }
            j[l+1] = data[k].j-1;
            a[l+1] = data[k].a;
            row_ptr[i+1]++;
        }

    } else if (include_strict_upper_triangular_part &&
               (matrix->symmetry == matrix_market_symmetric ||
                matrix->symmetry == matrix_market_skew_symmetric ||
                matrix->symmetry == matrix_market_hermitian))
    {
        for (int64_t k = 0; k < matrix->num_nonzeros; k++) {
            int i = data[k].i-1;
            int64_t l = row_ptr[i+1]-1;
            while (l >= row_ptr[i] && j[l] > data[k].j-1) {
                j[l+1] = j[l];
                a[l+1] = a[l];
                l--;
            }
            j[l+1] = data[k].j-1;
            a[l+1] = data[k].a;
            row_ptr[i+1]++;
            if (data[k].i == data[k].j)
                continue;
            i = data[k].j-1;
            l = row_ptr[i+1]-1;
            while (l >= row_ptr[i] && j[l] > data[k].i-1) {
                j[l+1] = j[l];
                a[l+1] = a[l];
                l--;
            }
            j[l+1] = data[k].i-1;
            if (matrix->symmetry == matrix_market_symmetric ||
                matrix->symmetry == matrix_market_hermitian)
            {
                a[l+1] = data[k].a;
            } else if (matrix->symmetry == matrix_market_skew_symmetric) {
                a[l+1] = -data[k].a;
            }
            row_ptr[i+1]++;
        }

    } else {
        return EINVAL;
    }
    return 0;
}

/**
 * `matrix_market_sort_nonzeros_double()` sorts the nonzeros according
 * to their rows and columns for a matrix in the matrix market format
 * whose `field` is `double`.
 */
static int matrix_market_sort_nonzeros_double(
    const struct matrix_market * matrix,
    const struct matrix_market_coordinate_double * data,
    bool include_strict_upper_triangular_part,
    int64_t * row_ptr,
    int32_t * j,
    double * a)
{
    if (matrix->field != matrix_market_double)
        return EINVAL;

    if ((include_strict_upper_triangular_part &&
         matrix->symmetry == matrix_market_general) ||
        (!include_strict_upper_triangular_part &&
         (matrix->symmetry == matrix_market_symmetric ||
          matrix->symmetry == matrix_market_skew_symmetric ||
          matrix->symmetry == matrix_market_hermitian)))
    {
        for (int64_t k = 0; k < matrix->num_nonzeros; k++) {
            int i = data[k].i-1;
            int64_t l = row_ptr[i+1]-1;
            while (l >= row_ptr[i] && j[l] > data[k].j-1) {
                j[l+1] = j[l];
                a[l+1] = a[l];
                l--;
            }
            j[l+1] = data[k].j-1;
            a[l+1] = data[k].a;
            row_ptr[i+1]++;
        }

    } else if (include_strict_upper_triangular_part &&
               (matrix->symmetry == matrix_market_symmetric ||
                matrix->symmetry == matrix_market_skew_symmetric ||
                matrix->symmetry == matrix_market_hermitian))
    {
        for (int64_t k = 0; k < matrix->num_nonzeros; k++) {
            int i = data[k].i-1;
            int64_t l = row_ptr[i+1]-1;
            while (l >= row_ptr[i] && j[l] > data[k].j-1) {
                j[l+1] = j[l];
                a[l+1] = a[l];
                l--;
            }
            j[l+1] = data[k].j-1;
            a[l+1] = data[k].a;
            row_ptr[i+1]++;
            if (data[k].i == data[k].j)
                continue;
            i = data[k].j-1;
            l = row_ptr[i+1]-1;
            while (l >= row_ptr[i] && j[l] > data[k].i-1) {
                j[l+1] = j[l];
                a[l+1] = a[l];
                l--;
            }
            j[l+1] = data[k].i-1;
            if (matrix->symmetry == matrix_market_symmetric ||
                matrix->symmetry == matrix_market_hermitian)
            {
                a[l+1] = data[k].a;
            } else if (matrix->symmetry == matrix_market_skew_symmetric) {
                a[l+1] = -data[k].a;
            }
            row_ptr[i+1]++;
        }

    } else {
        return EINVAL;
    }
    return 0;
}

/**
 * `matrix_market_sort_nonzeros_complex()` sorts the nonzeros according
 * to their rows and columns for a matrix in the matrix market format
 * whose `field` is `complex`.
 */
static int matrix_market_sort_nonzeros_complex(
    const struct matrix_market * matrix,
    const struct matrix_market_coordinate_complex * data,
    bool include_strict_upper_triangular_part,
    int64_t * row_ptr,
    int32_t * j,
    float * a)
{
    if (matrix->field != matrix_market_complex)
        return EINVAL;

    if ((include_strict_upper_triangular_part &&
         matrix->symmetry == matrix_market_general) ||
        (!include_strict_upper_triangular_part &&
         (matrix->symmetry == matrix_market_symmetric ||
          matrix->symmetry == matrix_market_skew_symmetric ||
          matrix->symmetry == matrix_market_hermitian)))
    {
        for (int64_t k = 0; k < matrix->num_nonzeros; k++) {
            int i = data[k].i-1;
            int64_t l = row_ptr[i+1]-1;
            while (l >= row_ptr[i] && j[l] > data[k].j-1) {
                j[l+1] = j[l];
                a[2*(l+1)+0] = a[2*l+0];
                a[2*(l+1)+1] = a[2*l+1];
                l--;
            }
            j[l+1] = data[k].j-1;
            a[2*(l+1)+0] = data[k].ar;
            a[2*(l+1)+1] = data[k].ai;
            row_ptr[i+1]++;
        }

    } else if (include_strict_upper_triangular_part &&
               (matrix->symmetry == matrix_market_symmetric ||
                matrix->symmetry == matrix_market_skew_symmetric ||
                matrix->symmetry == matrix_market_hermitian))
    {
        for (int64_t k = 0; k < matrix->num_nonzeros; k++) {
            int i = data[k].i-1;
            int64_t l = row_ptr[i+1]-1;
            while (l >= row_ptr[i] && j[l] > data[k].j-1) {
                j[l+1] = j[l];
                a[2*(l+1)+0] = a[2*l+0];
                a[2*(l+1)+1] = a[2*l+1];
                l--;
            }
            j[l+1] = data[k].j-1;
            a[2*(l+1)+0] = data[k].ar;
            a[2*(l+1)+1] = data[k].ai;
            row_ptr[i+1]++;
            if (data[k].i == data[k].j)
                continue;
            i = data[k].j-1;
            l = row_ptr[i+1]-1;
            while (l >= row_ptr[i] && j[l] > data[k].i-1) {
                j[l+1] = j[l];
                a[2*(l+1)+0] = a[2*l+0];
                a[2*(l+1)+1] = a[2*l+1];
                l--;
            }
            j[l+1] = data[k].i-1;
            if (matrix->symmetry == matrix_market_symmetric) {
                a[2*(l+1)+0] = data[k].ar;
                a[2*(l+1)+1] = data[k].ai;
            } else if (matrix->symmetry == matrix_market_skew_symmetric) {
                a[2*(l+1)+0] = -data[k].ar;
                a[2*(l+1)+1] = -data[k].ai;
            } else if (matrix->symmetry == matrix_market_hermitian) {
                a[2*(l+1)+0] = data[k].ar;
                a[2*(l+1)+1] = -data[k].ai;
            }
            row_ptr[i+1]++;
        }

    } else {
        return EINVAL;
    }
    return 0;
}

/**
 * `matrix_market_sort_nonzeros_integer()` sorts the nonzeros according
 * to their rows and columns for a matrix in the matrix market format
 * whose `field` is `integer`.
 */
static int matrix_market_sort_nonzeros_integer(
    const struct matrix_market * matrix,
    const struct matrix_market_coordinate_integer * data,
    bool include_strict_upper_triangular_part,
    int64_t * row_ptr,
    int32_t * j,
    int * a)
{
    if (matrix->field != matrix_market_integer)
        return EINVAL;

    if ((include_strict_upper_triangular_part &&
         matrix->symmetry == matrix_market_general) ||
        (!include_strict_upper_triangular_part &&
         (matrix->symmetry == matrix_market_symmetric ||
          matrix->symmetry == matrix_market_skew_symmetric ||
          matrix->symmetry == matrix_market_hermitian)))
    {
        for (int64_t k = 0; k < matrix->num_nonzeros; k++) {
            int i = data[k].i-1;
            int64_t l = row_ptr[i+1]-1;
            while (l >= row_ptr[i] && j[l] > data[k].j-1) {
                j[l+1] = j[l];
                a[l+1] = a[l];
                l--;
            }
            j[l+1] = data[k].j-1;
            a[l+1] = data[k].a;
            row_ptr[i+1]++;
        }

    } else if (include_strict_upper_triangular_part &&
               (matrix->symmetry == matrix_market_symmetric ||
                matrix->symmetry == matrix_market_skew_symmetric ||
                matrix->symmetry == matrix_market_hermitian))
    {
        for (int64_t k = 0; k < matrix->num_nonzeros; k++) {
            int i = data[k].i-1;
            int64_t l = row_ptr[i+1]-1;
            while (l >= row_ptr[i] && j[l] > data[k].j-1) {
                j[l+1] = j[l];
                a[l+1] = a[l];
                l--;
            }
            j[l+1] = data[k].j-1;
            a[l+1] = data[k].a;
            row_ptr[i+1]++;
            if (data[k].i == data[k].j)
                continue;
            i = data[k].j-1;
            l = row_ptr[i+1]-1;
            while (l >= row_ptr[i] && j[l] > data[k].i-1) {
                j[l+1] = j[l];
                a[l+1] = a[l];
                l--;
            }
            j[l+1] = data[k].i-1;
            if (matrix->symmetry == matrix_market_symmetric ||
                matrix->symmetry == matrix_market_hermitian)
            {
                a[l+1] = data[k].a;
            } else if (matrix->symmetry == matrix_market_skew_symmetric) {
                a[l+1] = -data[k].a;
            }
            row_ptr[i+1]++;
        }

    } else {
        return EINVAL;
    }
    return 0;
}

/**
 * `matrix_market_sort_nonzeros_pattern()` sorts the nonzeros according
 * to their rows and columns for a matrix in the matrix market format
 * whose `field` is `pattern`.
 */
static int matrix_market_sort_nonzeros_pattern(
    const struct matrix_market * matrix,
    const struct matrix_market_coordinate_pattern * data,
    bool include_strict_upper_triangular_part,
    int64_t * row_ptr,
    int32_t * j)
{
    if (matrix->field != matrix_market_pattern)
        return EINVAL;

    if ((include_strict_upper_triangular_part &&
         matrix->symmetry == matrix_market_general) ||
        (!include_strict_upper_triangular_part &&
         (matrix->symmetry == matrix_market_symmetric ||
          matrix->symmetry == matrix_market_skew_symmetric ||
          matrix->symmetry == matrix_market_hermitian)))
    {
        for (int64_t k = 0; k < matrix->num_nonzeros; k++) {
            int i = data[k].i-1;
            int64_t l = row_ptr[i+1]-1;
            while (l >= row_ptr[i] && j[l] > data[k].j-1) {
                j[l+1] = j[l];
                l--;
            }
            j[l+1] = data[k].j-1;
            row_ptr[i+1]++;
        }

    } else if (include_strict_upper_triangular_part &&
               (matrix->symmetry == matrix_market_symmetric ||
                matrix->symmetry == matrix_market_skew_symmetric ||
                matrix->symmetry == matrix_market_hermitian))
    {
        for (int64_t k = 0; k < matrix->num_nonzeros; k++) {
            int i = data[k].i-1;
            int64_t l = row_ptr[i+1]-1;
            while (l >= row_ptr[i] && j[l] > data[k].j-1) {
                j[l+1] = j[l];
                l--;
            }
            j[l+1] = data[k].j-1;
            row_ptr[i+1]++;
            if (data[k].i == data[k].j)
                continue;
            i = data[k].j-1;
            l = row_ptr[i+1]-1;
            while (l >= row_ptr[i] && j[l] > data[k].i-1) {
                j[l+1] = j[l];
                l--;
            }
            j[l+1] = data[k].i-1;
            row_ptr[i+1]++;
        }

    } else {
        return EINVAL;
    }
    return 0;
}

/**
 * `matrix_market_sort_nonzeros()` sorts the nonzeros according to
 * their rows and columns for a matrix in the matrix market format.
 */
int matrix_market_sort_nonzeros(
    const struct matrix_market * matrix,
    int64_t * row_ptr,
    int32_t * column_indices,
    void * values)
{
    int err;

    if (matrix->format != matrix_market_coordinate)
        return ENOTSUP;

    /* 1. Count the number of nonzeros in each row. */
    bool include_strict_upper_triangular_part = true;
    err = matrix_market_nonzeros_per_row(
        matrix, include_strict_upper_triangular_part, &row_ptr[1]);
    if (err)
        return err;

    /* 2. Adjust row pointers to point to the first entry of the
     * previous row. */
    for (int i = 1; i <= matrix->num_rows; i++)
        row_ptr[i] += row_ptr[i-1];
    for (int i = matrix->num_rows; i > 0; i--)
        row_ptr[i] = row_ptr[i-1];
    row_ptr[0] = 0;

    /* 3. Sort nonzeros by their rows and columns. */
    switch (matrix->field) {
    case matrix_market_real:
        err = matrix_market_sort_nonzeros_real(
            matrix, (const struct matrix_market_coordinate_real *) matrix->data,
            include_strict_upper_triangular_part, row_ptr, column_indices, (float *) values);
        if (err)
            return err;
        break;

    case matrix_market_double:
        err = matrix_market_sort_nonzeros_double(
            matrix, (const struct matrix_market_coordinate_double *) matrix->data,
            include_strict_upper_triangular_part, row_ptr, column_indices, (double *) values);
        if (err)
            return err;
        break;
    case matrix_market_complex:
        err = matrix_market_sort_nonzeros_complex(
            matrix, (const struct matrix_market_coordinate_complex *) matrix->data,
            include_strict_upper_triangular_part, row_ptr, column_indices, (float *) values);
        if (err)
            return err;
        break;
    case matrix_market_integer:
        err = matrix_market_sort_nonzeros_integer(
            matrix, (const struct matrix_market_coordinate_integer *) matrix->data,
            include_strict_upper_triangular_part, row_ptr, column_indices, (int *) values);
        if (err)
            return err;
        break;
    case matrix_market_pattern:
        err = matrix_market_sort_nonzeros_pattern(
            matrix, (const struct matrix_market_coordinate_pattern *) matrix->data,
            include_strict_upper_triangular_part, row_ptr, column_indices);
        if (err)
            return err;
        break;
    default:
        return EINVAL;
    }

    return 0;
}
