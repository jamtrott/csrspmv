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
 * Benchmarking program for sparse matrix-vector multiplication (SpMV)
 * with matrices in compressed sparse row (CSR) format.
 */

#include "csr.h"
#include "matrix_market.h"
#include "program_options.h"
#include "vector.h"

#include <errno.h>
#include <omp.h>

#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

const char * program_name = "csrspmv";
const char * program_version = "1.0";
const char * program_copyright =
    "Copyright (C) 2020 James D. Trotter";
const char * program_license =
    "License GPLv3+: GNU GPL version 3 or later <https://gnu.org/licenses/gpl.html>\n"
    "This is free software: you are free to change and redistribute it.\n"
    "There is NO WARRANTY, to the extent permitted by law.";

/**
 * `timespec_duration()` is the duration, in seconds, elapsed between
 * two given time points.
 */
static double timespec_duration(
    struct timespec t0,
    struct timespec t1)
{
    return (t1.tv_sec - t0.tv_sec) +
        (t1.tv_nsec - t0.tv_nsec) * 1e-9;
}

/**
 * `matrix_market_format_comment()` formats a string to provide some
 * useful information to be placed in the comments section of a matrix
 * market file.
 *
 * The string returned in `comment1` and `comment2` must be freed by
 * the caller.
 */
static int matrix_market_format_comment(
    int argc, char ** argv, char ** comment1, char ** comment2)
{
    /* Program name and version. */
    int bytes_needed;
    int bytes_written;
    *comment1 = NULL;
    const char * comment_fmt =
        "%% This file was written by %s %s\n";
    bytes_needed = snprintf(
        *comment1, 0, comment_fmt,
        program_invocation_short_name, program_version);
    if (bytes_needed == -1)
        return errno;
    *comment1 = malloc(bytes_needed+1);
    if (!comment1)
        return errno;
    bytes_written = snprintf(
        *comment1, bytes_needed+1, comment_fmt,
        program_invocation_short_name, program_version);
    if (bytes_written == -1) {
        free(*comment1);
        return errno;
    }

    /* Concatenate program options. */
    size_t len = 4;
    for (int i = 0; i < argc; i++)
        len += strlen(argv[i])+1;

    *comment2 = malloc(len);
    if (!*comment2) {
        free(*comment1);
        return errno;
    }

    char * s = *comment2;
    *s = '%'; s++;
    *s = ' '; s++;
    for (int i = 0; i < argc; i++) {
        *s = ' '; s++;
        for (int j = 0; j < strlen(argv[i]); j++) {
            *s = argv[i][j]; s++;
        }
    }
    *s = '\n'; s++;
    *s = '\0';
    return 0;
}

/*
 * Custom OpenMP reduction operator for combining errors from
 * different threads.
 */
#pragma omp declare reduction(                                          \
    err_add : int :                                                     \
    omp_out = omp_out ? omp_out : omp_in)                               \
    initializer (omp_priv=0)

/**
 * `main()`.
 */
int main(int argc, char *argv[])
{
    int err;
    struct timespec t0, t1;

    /* Parse program options. */
    struct program_options args;
    int argc_copy = argc;
    char ** argv_copy = argv;
    err = parse_program_options(&argc_copy, &argv_copy, &args);
    if (err) {
        fprintf(stderr, "%s: %s %s\n", program_invocation_short_name,
                strerror(err), argv_copy[0]);
        return EXIT_FAILURE;
    }
    if (!args.matrix_file) {
        fprintf(stderr, "%s: Please specify a matrix market file\n",
                program_invocation_short_name);
        program_options_free(&args);
        return EXIT_FAILURE;
    }

    FILE * f;
    if ((f = fopen(args.matrix_file, "r")) == NULL) {
        fprintf(stderr, "%s: %s: %s\n",
                program_invocation_short_name, strerror(errno),
            args.matrix_file);
        program_options_free(&args);
        return EXIT_FAILURE;
    }

    if (args.verbose > 0) {
        fprintf(stdout, "matrix_market_read: ");
        fflush(stdout);
        clock_gettime(CLOCK_MONOTONIC, &t0);
    }

    /* 1. Read a matrix from a file in the matrix market format. */
    struct matrix_market matrix_market;
    int line_number, column_number;
    err = matrix_market_read(
        &matrix_market, f, &line_number, &column_number);
    if (err) {
        fprintf(stderr, "%s: Error: %s:%d:%d: %s\n",
                program_invocation_name,
                args.matrix_file, line_number, column_number,
                matrix_market_strerror(err));
        fclose(f);
        program_options_free(&args);
        return EXIT_FAILURE;
    }
    fclose(f);

    if (args.verbose > 0) {
        clock_gettime(CLOCK_MONOTONIC, &t1);
        fprintf(stdout, "%.6f seconds "
                "%s object %s format %s field %s symmetry %d rows %d columns %"PRId64" nonzeros\n",
                timespec_duration(t0, t1),
                matrix_market_object_str(matrix_market.object),
                matrix_market_format_str(matrix_market.format),
                matrix_market_field_str(matrix_market.field),
                matrix_market_symmetry_str(matrix_market.symmetry),
                matrix_market.num_rows,
                matrix_market.num_columns,
                matrix_market.num_nonzeros);
        fprintf(stdout, "csr_matrix_int32_from_matrix_market: ");
        fflush(stdout);
        clock_gettime(CLOCK_MONOTONIC, &t0);
    }

    /* Determine which format to use for matrix values. */
    enum csr_value_format matrix_format;
    if (!args.matrix_format_auto) {
        matrix_format = args.matrix_format;
    } else {
        switch (matrix_market.field) {
        case matrix_market_real: matrix_format = csr_value_f32; break;
        case matrix_market_double: matrix_format = csr_value_f64; break;
        case matrix_market_complex: matrix_format = csr_value_complex32; break;
        case matrix_market_integer: matrix_format = csr_value_int32; break;
        case matrix_market_pattern: matrix_format = csr_value_binary; break;
        default:
            fprintf(stderr, "%s: %s\n", program_invocation_short_name,
                    strerror(ENOTSUP));
            matrix_market_free(&matrix_market);
            program_options_free(&args);
            return EXIT_FAILURE;
        }
    }

    /* 2. Convert to a compressed sparse row format. */
    struct csr_matrix_int32 csr_matrix;
    err = csr_matrix_int32_from_matrix_market(
        &csr_matrix, &matrix_market, matrix_format);
    if (err) {
        fprintf(stderr, "%s: %s\n", program_invocation_short_name,
                strerror(err));
        matrix_market_free(&matrix_market);
        program_options_free(&args);
        return EXIT_FAILURE;
    }
    matrix_market_free(&matrix_market);

    if (args.verbose > 0) {
        clock_gettime(CLOCK_MONOTONIC, &t1);
        fprintf(stdout, "%.6f seconds\n", timespec_duration(t0, t1));
    }

    /* 3. Load or generate a source vector. */
    struct vector x;
    if (args.source_vector_file) {
        FILE * f;
        if ((f = fopen(args.source_vector_file, "r")) == NULL) {
            fprintf(stderr, "%s: %s: %s\n",
                    program_invocation_short_name,
                    args.source_vector_file, strerror(errno));
            csr_matrix_int32_free(&csr_matrix);
            program_options_free(&args);
            return EXIT_FAILURE;
        }

        /* Read vector from file in matrix market format. */
        struct matrix_market matrix_market;
        int line_number, column_number;
        err = matrix_market_read(
            &matrix_market, f, &line_number, &column_number);
        if (err) {
            fprintf(stderr, "%s: Error: %s:%d:%d: %s\n",
                    program_invocation_name,
                    args.source_vector_file, line_number, column_number,
                    matrix_market_strerror(err));
            fclose(f);
            csr_matrix_int32_free(&csr_matrix);
            program_options_free(&args);
            return EXIT_FAILURE;
        }
        fclose(f);

        err = vector_from_matrix_market(&x, &matrix_market);
        if (err) {
            fprintf(stderr, "%s: %s: %s\n", program_invocation_name,
                    args.source_vector_file, strerror(err));
            matrix_market_free(&matrix_market);
            csr_matrix_int32_free(&csr_matrix);
            program_options_free(&args);
            return EXIT_FAILURE;
        }
        matrix_market_free(&matrix_market);

    } else {
        enum vector_value_format format;
        if (!args.source_vector_format_auto) {
            format = args.source_vector_format;
        } else if (csr_matrix.value_format == csr_value_f32) {
            format = vector_value_f32;
        } else if (csr_matrix.value_format == csr_value_f64) {
            format = vector_value_f64;
        } else if (csr_matrix.value_format == csr_value_complex32) {
            format = vector_value_complex32;
        } else {
            fprintf(stderr, "%s: %s\n", program_invocation_short_name,
                    strerror(ENOTSUP));
            csr_matrix_int32_free(&csr_matrix);
            program_options_free(&args);
            return EXIT_FAILURE;
        }

        err = vector_alloc(&x, csr_matrix.num_columns, format);
        if (err) {
            fprintf(stderr, "%s: %s\n", program_invocation_short_name,
                    strerror(err));
            csr_matrix_int32_free(&csr_matrix);
            program_options_free(&args);
            return EXIT_FAILURE;
        }

        err = vector_ones(&x);
        if (err) {
            fprintf(stderr, "%s: %s\n", program_invocation_short_name,
                    strerror(err));
            vector_free(&x);
            csr_matrix_int32_free(&csr_matrix);
            program_options_free(&args);
            return EXIT_FAILURE;
        }
    }

    /* 4. Allocate storage for a destination vector. */
    struct vector y;
    {
        enum vector_value_format format;
        if (!args.destination_vector_format_auto) {
            format = args.destination_vector_format;
        } else if (csr_matrix.value_format == csr_value_f32) {
            format = vector_value_f32;
        } else if (csr_matrix.value_format == csr_value_f64) {
            format = vector_value_f64;
        } else if (csr_matrix.value_format == csr_value_complex32) {
            format = vector_value_complex32;
        } else {
            fprintf(stderr, "%s: %s\n", program_invocation_short_name,
                    strerror(ENOTSUP));
            vector_free(&x);
            csr_matrix_int32_free(&csr_matrix);
            program_options_free(&args);
            return EXIT_FAILURE;
        }

        err = vector_alloc(&y, csr_matrix.num_rows, format);
        if (err) {
            fprintf(stderr, "%s: %s\n", program_invocation_short_name,
                    strerror(err));
            vector_free(&x);
            csr_matrix_int32_free(&csr_matrix);
            program_options_free(&args);
            return EXIT_FAILURE;
        }

        err = vector_zero(&y);
        if (err) {
            fprintf(stderr, "%s: %s\n", program_invocation_short_name,
                    strerror(err));
            vector_free(&y);
            vector_free(&x);
            csr_matrix_int32_free(&csr_matrix);
            program_options_free(&args);
            return EXIT_FAILURE;
        }
    }

    if (args.verbose > 0) {
        fprintf(stdout, "csr_matrix_int32_spmv: ");
        fflush(stdout);
        clock_gettime(CLOCK_MONOTONIC, &t0);
    }

    /* 5. Compute the sparse matrix-vector multiplication. */
#pragma omp parallel reduction(err_add:err)
    for (int i = 0; i < args.repeat; i++) {
        err = csr_matrix_int32_spmv(&csr_matrix, &x, &y);
        if (err)
            break;
    }
    if (err) {
        fprintf(stderr, "%s: %s\n", program_invocation_name,
                strerror(err));
        vector_free(&y);
        vector_free(&x);
        csr_matrix_int32_free(&csr_matrix);
        program_options_free(&args);
        return EXIT_FAILURE;
    }

    if (args.verbose > 0) {
        clock_gettime(CLOCK_MONOTONIC, &t1);
        fprintf(stdout, "%.6f seconds %d multiplications\n",
                timespec_duration(t0, t1), args.repeat);
        fflush(stdout);
    }

    /* 6. Write the destination vector to a file. */
    if (args.destination_vector_file) {

        if (args.verbose > 0) {
            fprintf(stdout, "vector_to_matrix_market: ");
            fflush(stdout);
            clock_gettime(CLOCK_MONOTONIC, &t0);
        }

        /* Format a comment. */
        int num_comment_lines = 2;
        char * comment_lines[2];
        err = matrix_market_format_comment(
            argc, argv, &comment_lines[0], &comment_lines[1]);
        if (err) {
            fprintf(stderr, "%s: %s\n", program_invocation_name,
                    strerror(err));
            vector_free(&y);
            vector_free(&x);
            csr_matrix_int32_free(&csr_matrix);
            program_options_free(&args);
            return EXIT_FAILURE;
        }

        /* Convert destination vector to matrix market format. */
        struct matrix_market matrix_market_vector;
        err = vector_to_matrix_market(
            &y, &matrix_market_vector,
            num_comment_lines, (const char **) comment_lines);
        if (err) {
            fprintf(stderr, "%s: %s\n", program_invocation_name,
                    strerror(err));
            free(comment_lines[1]);
            free(comment_lines[0]);
            vector_free(&y);
            vector_free(&x);
            csr_matrix_int32_free(&csr_matrix);
            program_options_free(&args);
            return EXIT_FAILURE;
        }
        free(comment_lines[1]);
        free(comment_lines[0]);

        if (args.verbose > 0) {
            clock_gettime(CLOCK_MONOTONIC, &t1);
            fprintf(stdout, "%.6f seconds\n", timespec_duration(t0, t1));
            fprintf(stdout, "matrix_market_write: ");
            fflush(stdout);
            clock_gettime(CLOCK_MONOTONIC, &t0);
        }

        /* Open destination vector file */
        FILE * f = fopen(args.destination_vector_file, "w");
        if (!f) {
            fprintf(stderr, "%s: %s: %s\n",
                    program_invocation_short_name, strerror(errno),
                    args.destination_vector_file);
            matrix_market_free(&matrix_market_vector);
            vector_free(&y);
            vector_free(&x);
            csr_matrix_int32_free(&csr_matrix);
            program_options_free(&args);
            return EXIT_FAILURE;
        }

        /* Write destination vector to file. */
        err = matrix_market_write(
            &matrix_market_vector, f,
            args.destination_vector_field_width,
            args.destination_vector_precision);
        if (err) {
            fprintf(stderr, "%s: %s: %s\n",
                    program_invocation_name,
                    args.destination_vector_file,
                    strerror(err));
            matrix_market_free(&matrix_market_vector);
            vector_free(&y);
            vector_free(&x);
            csr_matrix_int32_free(&csr_matrix);
            program_options_free(&args);
            return EXIT_FAILURE;
        }
        fclose(f);

        if (args.verbose > 0) {
            clock_gettime(CLOCK_MONOTONIC, &t1);
            fprintf(stdout, "%.6f seconds\n", timespec_duration(t0, t1));
            fflush(stdout);
        }

        matrix_market_free(&matrix_market_vector);
    }

    /* 7. Clean up. */
    vector_free(&y);
    vector_free(&x);
    csr_matrix_int32_free(&csr_matrix);
    program_options_free(&args);
    return EXIT_SUCCESS;
}
