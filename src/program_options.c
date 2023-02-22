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
 * Parsing of program options.
 */

#include "program_options.h"
#include "parse.h"
#include "vector.h"

#include <errno.h>

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * `program_options_init()` configures the default program options.
 */
static int program_options_init(
    struct program_options * args)
{
    args->matrix_file = NULL;
    args->source_vector_file = NULL;
    args->destination_vector_file = NULL;
    args->destination_vector_field_width = 0;
    args->destination_vector_precision = -1;
    args->matrix_format_auto = true;
    args->include_symmetric_part = false;
    args->matrix_format = csr_value_f64;
    args->source_vector_format_auto = true;
    args->source_vector_format = vector_value_f64;
    args->destination_vector_format_auto = true;
    args->destination_vector_format = vector_value_f64;
    args->repeat = 1;
    args->min_flops = 0;
    args->count_nonzeros = 0;
    args->verbose = 0;
    args->help = false;
    args->version = false;
    return 0;
}

/**
 * `program_options_free()` frees memory and other resources
 * associated with parsing program options.
 */
void program_options_free(
    struct program_options * args)
{
    if (args->matrix_file)
        free(args->matrix_file);
    if (args->source_vector_file)
        free(args->source_vector_file);
    if (args->destination_vector_file)
        free(args->destination_vector_file);
}

/**
 * `program_options_print_help()` prints a help text.
 */
void program_options_print_help(
    FILE * f)
{
    fprintf(f, "Usage: %s [OPTION..] FILE\n", program_name);
    fprintf(f, "\n");
    fprintf(f, " Multiply a sparse matrix in CSR format with a vector\n");
    fprintf(f, "\n");
    fprintf(f, " Options are:\n");
    fprintf(f, "  --matrix-format=FORMAT\tchoose one of: binary, int32, int64,\n");
    fprintf(f, "  \t\t\t\tf32, f64 or complex32.\n");
    fprintf(f, "  --include-symmetric-part\tuse both lower and upper triangular parts of a symmetric matrix\n");
    fprintf(f, "  --source-vector=FILE\t\tmatrix market file for source vector,\n");
    fprintf(f, "  \t\t\t\totherwise a vector of all ones is used.\n");
    fprintf(f, "  --source-vector-format=FORMAT\tchoose one of: int32, f32, f64 or complex32.\n");
    fprintf(f, "  --destination-vector=FILE\tmatrix market file for destination vector\n");
    fprintf(f, "  --destination-vector-width=N\tfield width for destination vector output\n");
    fprintf(f, "  --destination-vector-prec=N\tprecision for destination vector output\n");
    fprintf(f, "  --destination-vector-format=FORMAT\tchoose one of: int32, f32, f64 or complex32.\n");
    fprintf(f, "  --flops=N\t\t\tminimum number of arithmetic operations to perform\n");
    fprintf(f, "  -r, --repeat=N\t\trepeat matrix-vector multiplication\n");
    fprintf(f, "  --count-nonzeros\t\tdisplay the number of nonzeros assigned to each thread and exit\n");
    fprintf(f, "  -v, --verbose\t\t\tbe more verbose\n");
    fprintf(f, "\n");
    fprintf(f, "  -h, --help\t\t\tdisplay this help and exit\n");
    fprintf(f, "  --version\t\t\tdisplay version information and exit\n");
    fprintf(f, "\n");
    fprintf(f, "Report bugs to: <james@simula.no>\n");
}

/**
 * `program_options_print_version()` prints version information.
 */
static void program_options_print_version(
    FILE * f)
{
    fprintf(f, "%s %s\n", program_name, program_version);
    fprintf(f, "%s\n", program_copyright);
    fprintf(f, "%s\n", program_license);
}

/**
 * `parse_program_options()` parses program options.
 */
int parse_program_options(
    int * argc,
    char *** argv,
    struct program_options * args)
{
    int err;

    /* Set program invocation name. */
    program_invocation_name = (*argv)[0];
    program_invocation_short_name = (
        strrchr(program_invocation_name, '/')
        ? strrchr(program_invocation_name, '/') + 1
        : program_invocation_name);
    (*argc)--; (*argv)++;

    /* Set default program options. */
    err = program_options_init(args);
    if (err)
        return err;

    /* Parse program options. */
    int num_arguments_consumed = 0;
    while (*argc > 0) {
        *argc -= num_arguments_consumed;
        *argv += num_arguments_consumed;
        num_arguments_consumed = 0;
        if (*argc <= 0)
            break;

        /* Parse matrix format. */
        if (strcmp((*argv)[0], "--matrix-format") == 0) {
            if (*argc < 2) {
                program_options_free(args);
                return EINVAL;
            }
            err = parse_csr_value_format(
                (*argv)[1], &args->matrix_format);
            if (err) {
                program_options_free(args);
                return err;
            }
            args->matrix_format_auto = false;
            num_arguments_consumed += 2;
            continue;
        } else if (strstr((*argv)[0], "--matrix-format=") == (*argv)[0]) {
            err = parse_csr_value_format(
                (*argv)[0] + strlen("--matrix-format="),
                &args->matrix_format);
            if (err) {
                program_options_free(args);
                return err;
            }
            args->matrix_format_auto = false;
            num_arguments_consumed++;
            continue;
        }

        /* Parse source vector file name. */
        if (strcmp((*argv)[0], "--source-vector") == 0) {
            if (*argc < 2) {
                program_options_free(args);
                return EINVAL;
            }
            if (args->source_vector_file)
                free(args->source_vector_file);
            args->source_vector_file = strdup((*argv)[1]);
            if (err) {
                program_options_free(args);
                return err;
            }
            num_arguments_consumed += 2;
            continue;
        } else if (strstr((*argv)[0], "--source-vector=") == (*argv)[0]) {
            if (args->source_vector_file)
                free(args->source_vector_file);
            args->source_vector_file =
                strdup((*argv)[0] + strlen("--source-vector="));
            if (err) {
                program_options_free(args);
                return err;
            }
            num_arguments_consumed++;
            continue;
        }

        /* Parse source vector format. */
        if (strcmp((*argv)[0], "--source-vector-format") == 0) {
            if (*argc < 2) {
                program_options_free(args);
                return EINVAL;
            }
            err = parse_vector_value_format(
                (*argv)[1], &args->source_vector_format);
            if (err) {
                program_options_free(args);
                return err;
            }
            args->source_vector_format_auto = false;
            num_arguments_consumed += 2;
            continue;
        } else if (strstr((*argv)[0], "--source-vector-format=") == (*argv)[0]) {
            err = parse_vector_value_format(
                (*argv)[0] + strlen("--source-vector-format="),
                &args->source_vector_format);
            if (err) {
                program_options_free(args);
                return err;
            }
            args->source_vector_format_auto = false;
            num_arguments_consumed++;
            continue;
        }

        /* Parse destination vector file name. */
        if (strcmp((*argv)[0], "--destination-vector") == 0) {
            if (*argc < 2) {
                program_options_free(args);
                return EINVAL;
            }
            if (args->destination_vector_file)
                free(args->destination_vector_file);
            args->destination_vector_file = strdup((*argv)[1]);
            if (err) {
                program_options_free(args);
                return err;
            }
            num_arguments_consumed += 2;
            continue;
        } else if (strstr((*argv)[0], "--destination-vector=") == (*argv)[0]) {
            if (args->destination_vector_file)
                free(args->destination_vector_file);
            args->destination_vector_file =
                strdup((*argv)[0] + strlen("--destination-vector="));
            if (err) {
                program_options_free(args);
                return err;
            }
            num_arguments_consumed++;
            continue;
        }

        /* Parse destination vector format. */
        if (strcmp((*argv)[0], "--destination-vector-format") == 0) {
            if (*argc < 2) {
                program_options_free(args);
                return EINVAL;
            }
            err = parse_vector_value_format(
                (*argv)[1], &args->destination_vector_format);
            if (err) {
                program_options_free(args);
                return err;
            }
            args->destination_vector_format_auto = false;
            num_arguments_consumed += 2;
            continue;
        } else if (strstr((*argv)[0], "--destination-vector-format=") == (*argv)[0]) {
            err = parse_vector_value_format(
                (*argv)[0] + strlen("--destination-vector-format="),
                &args->destination_vector_format);
            if (err) {
                program_options_free(args);
                return err;
            }
            args->destination_vector_format_auto = false;
            num_arguments_consumed++;
            continue;
        }

        if (strcmp((*argv)[0], "--destination-vector-width") == 0) {
            if (*argc < 2) {
                program_options_free(args);
                return EINVAL;
            }
            err = parse_int32((*argv)[1], NULL, &args->destination_vector_field_width, NULL);
            if (err) {
                program_options_free(args);
                return err;
            }
            num_arguments_consumed += 2;
            continue;
        } else if (strstr((*argv)[0], "--destination-vector-width=") == (*argv)[0]) {
            err = parse_int32(
                (*argv)[0] + strlen("--destination-vector-width="), NULL,
                &args->destination_vector_field_width, NULL);
            if (err) {
                program_options_free(args);
                return err;
            }
            num_arguments_consumed++;
            continue;
        }

        if (strcmp((*argv)[0], "--destination-vector-prec") == 0) {
            if (*argc < 2) {
                program_options_free(args);
                return EINVAL;
            }
            err = parse_int32((*argv)[1], NULL, &args->destination_vector_precision, NULL);
            if (err) {
                program_options_free(args);
                return err;
            }
            num_arguments_consumed += 2;
            continue;
        } else if (strstr((*argv)[0], "--destination-vector-prec=") == (*argv)[0]) {
            err = parse_int32(
                (*argv)[0] + strlen("--destination-vector-prec="), NULL,
                &args->destination_vector_precision, NULL);
            if (err) {
                program_options_free(args);
                return err;
            }
            num_arguments_consumed++;
            continue;
        }

        if (strcmp((*argv)[0], "-r") == 0 || strcmp((*argv)[0], "--repeat") == 0) {
            if (*argc < 2) {
                program_options_free(args);
                return EINVAL;
            }
            err = parse_int32((*argv)[1], NULL, &args->repeat, NULL);
            if (err) {
                program_options_free(args);
                return err;
            }
            num_arguments_consumed += 2;
            continue;
        } else if (strstr((*argv)[0], "--repeat=") == (*argv)[0]) {
            err = parse_int32(
                (*argv)[0] + strlen("--repeat="), NULL, &args->repeat, NULL);
            if (err) {
                program_options_free(args);
                return err;
            }
            num_arguments_consumed++;
            continue;
        }

        if (strcmp((*argv)[0], "--flops") == 0) {
            if (*argc < 2) {
                program_options_free(args);
                return EINVAL;
            }
            err = parse_int64((*argv)[1], NULL, &args->min_flops, NULL);
            if (err) {
                program_options_free(args);
                return err;
            }
            num_arguments_consumed += 2;
            continue;
        } else if (strstr((*argv)[0], "--flops=") == (*argv)[0]) {
            err = parse_int64(
                (*argv)[0] + strlen("--flops="), NULL, &args->min_flops, NULL);
            if (err) {
                program_options_free(args);
                return err;
            }
            num_arguments_consumed++;
            continue;
        }

        if (strcmp((*argv)[0], "--include-symmetric-part") == 0) {
            args->include_symmetric_part = true;
            num_arguments_consumed++;
            continue;
        }

        if (strcmp((*argv)[0], "--count-nonzeros") == 0) {
            args->count_nonzeros = 1;
            num_arguments_consumed++;
            continue;
        }

        if (strcmp((*argv)[0], "-v") == 0 || strcmp((*argv)[0], "--verbose") == 0) {
            args->verbose++;
            num_arguments_consumed++;
            continue;
        }

        /* If requested, print program help text. */
        if (strcmp((*argv)[0], "-h") == 0 || strcmp((*argv)[0], "--help") == 0) {
            program_options_free(args);
            program_options_print_help(stdout);
            exit(EXIT_SUCCESS);
        }

        /* If requested, print program version information. */
        if (strcmp((*argv)[0], "--version") == 0) {
            program_options_free(args);
            program_options_print_version(stdout);
            exit(EXIT_SUCCESS);
        }

        /* Stop parsing options after '--'.  */
        if (strcmp((*argv)[0], "--") == 0) {
            (*argc)--; (*argv)++;
            break;
        }

        /* Parse matrix file name. */
        if (strlen((*argv)[0]) > 0 && (*argv)[0][0] != '-') {
            if (args->matrix_file)
                free(args->matrix_file);
            args->matrix_file = strdup((*argv)[0]);
            num_arguments_consumed++;
            continue;
        }

        /* Unrecognised option. */
        program_options_free(args);
        return EINVAL;
    }

    return 0;
}
