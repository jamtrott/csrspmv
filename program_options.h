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
 * Parsing of program options.
 */

#ifndef PROGRAM_OPTIONS_H
#define PROGRAM_OPTIONS_H

#include "csr.h"
#include "vector.h"

#include <stdbool.h>
#include <stdio.h>

extern const char * program_name;
extern const char * program_version;
extern const char * program_copyright;
extern const char * program_license;
extern const char * program_invocation_name;
extern const char * program_invocation_short_name;

/**
 * `program_options` contains data related program options.
 */
struct program_options
{
    char * matrix_file;
    char * source_vector_file;
    char * destination_vector_file;
    int destination_vector_field_width;
    int destination_vector_precision;
    bool matrix_format_auto;
    enum csr_value_format matrix_format;
    bool source_vector_format_auto;
    enum vector_value_format source_vector_format;
    bool destination_vector_format_auto;
    enum vector_value_format destination_vector_format;
    int repeat;
    int verbose;
    bool help;
    bool version;
};

/**
 * `parse_program_options()` parses program options.
 */
int parse_program_options(
    int * argc,
    char *** argv,
    struct program_options * args);

/**
 * `program_options_free()` frees any memory or resources associated
 * with parsing program options.
 */
void program_options_free(
    struct program_options * args);

#endif
