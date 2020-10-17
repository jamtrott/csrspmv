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
 *
 * String parsing.
 */

#include "parse.h"

#include <errno.h>

#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/**
 * `parse_long_long_int()` parses a string to produce a number that
 * may be represented with the type `long long int`.
 */
static int parse_long_long_int(
    const char * s,
    char ** endptr,
    int base,
    long long int * out_number)
{
    errno = 0;
    long long int number = strtoll(s, endptr, base);
    if ((errno == ERANGE && (number == LLONG_MAX || number == LLONG_MIN)) ||
        (errno != 0 && number == 0)) {
        return errno;
    }
    if (s == *endptr)
        return EINVAL;
    *out_number = number;
    return 0;
}

/**
 * `parse_int32()` parses a string to produce a number that may be
 * represented as a 32-bit integer.
 *
 * The number is parsed using `strtoll()`, following the conventions
 * documented in the man page for that function.  In addition, some
 * further error checking is performed to ensure that the number is
 * parsed correctly.  The parsed number is stored in `number`.
 *
 * `valid_delimiters` is either `NULL`, in which case it is ignored,
 * or, it may contain a string of characters that constitute valid
 * delimiters for the parsed string.  That is, after parsing a number,
 * if there are any remaining, unconsumed characters in the string,
 * `parse_int32()` checks if the next character is found in the string
 * `valid_delimiters`.  If the character is not found, then the string
 * is judged to be invalid, and `EINVAL` is returned.  Otherwise, the
 * final, delimiter character is consumed by the parser.
 *
 * If `endptr` is not `NULL`, the address stored in `endptr` points to
 * the first character beyond the characters that were consumed during
 * parsing.
 *
 * On success, `parse_int32()` returns `0`.  Otherwise, if the input
 * contained invalid characters, `parse_int32()` returns `EINVAL`.  If
 * the resulting number cannot be represented as a signed 32-bit
 * integer, `parse_int32()` returns `ERANGE`.
 */
int parse_int32(
    const char * s,
    const char * valid_delimiters,
    int32_t * out_number,
    const char ** endptr)
{
    int base = 10;
    char * s_end;
    long long int number;
    int err = parse_long_long_int(s, &s_end, base, &number);
    if (err)
        return err;

    /* Check that the number is within range. */
    if (number < INT32_MIN || number > INT32_MAX)
        return ERANGE;

    /* Check for a valid delimiter following the parsed number. */
    if (valid_delimiters && s_end && *s_end != '\0') {
        if (!strchr(valid_delimiters, *s_end)) {
            return EINVAL;
        }
        s_end++;
    }

    if (endptr)
        *endptr = s_end;
    *out_number = number;
    return 0;
}

/**
 * `parse_int64()` parses a string to produce a number that may be
 * represented as a 64-bit integer.
 *
 * The number is parsed using `strtoll()`, following the conventions
 * documented in the man page for that function.  In addition, some
 * further error checking is performed to ensure that the number is
 * parsed correctly.  The parsed number is stored in `number`.
 *
 * `valid_delimiters` is either `NULL`, in which case it is ignored,
 * or, it may contain a string of characters that constitute valid
 * delimiters for the parsed string.  That is, after parsing a number,
 * if there are any remaining, unconsumed characters in the string,
 * `parse_int64()` checks if the next character is found in the string
 * `valid_delimiters`.  If the character is not found, then the string
 * is judged to be invalid, and `EINVAL` is returned.  Otherwise, the
 * final, delimiter character is consumed by the parser.
 *
 * If `endptr` is not `NULL`, the address stored in `endptr` points to
 * the first character beyond the characters that were consumed during
 * parsing.
 *
 * On success, `parse_int64()` returns `0`.  Otherwise, if the input
 * contained invalid characters, `parse_int64()` returns `EINVAL`.  If
 * the resulting number cannot be represented as a signed 64-bit
 * integer, `parse_int64()` returns `ERANGE`.
 */
int parse_int64(
    const char * s,
    const char * valid_delimiters,
    int64_t * out_number,
    const char ** endptr)
{
    int base = 10;
    char * s_end;
    long long int number;
    int err = parse_long_long_int(s, &s_end, base, &number);
    if (err)
        return err;

    /* Check that the number is within range. */
    if (number < INT64_MIN || number > INT64_MAX)
        return ERANGE;

    /* Check for a valid delimiter following the parsed number. */
    if (valid_delimiters && s_end && *s_end != '\0') {
        if (!strchr(valid_delimiters, *s_end)) {
            return EINVAL;
        }
        s_end++;
    }

    if (endptr)
        *endptr = s_end;
    *out_number = number;
    return 0;
}

/**
 * `parse_float()` parses a string to produce a number that may be
 * represented as `float`.
 *
 * The number is parsed using `strtof()`, following the conventions
 * documented in the man page for that function.  In addition, some
 * further error checking is performed to ensure that the number is
 * parsed correctly.  The parsed number is stored in `number`.
 *
 * `valid_delimiters` is either `NULL`, in which case it is ignored,
 * or, it may contain a string of characters that constitute valid
 * delimiters for the parsed string.  That is, after parsing a number,
 * if there are any remaining, unconsumed characters in the string,
 * `parse_float()` checks if the next character is found in the string
 * `valid_delimiters`.  If the character is not found, then the string
 * is judged to be invalid, and `EINVAL` is returned.  Otherwise, the
 * final, delimiter character is consumed by the parser.
 *
 * If `endptr` is not `NULL`, the address stored in `endptr` points to
 * the first character beyond the characters that were consumed during
 * parsing.
 *
 * On success, `parse_float()` returns `0`.  Otherwise, if the input
 * contained invalid characters, `parse_float()` returns `EINVAL`.
 */
int parse_float(
    const char * s,
    const char * valid_delimiters,
    float * number,
    const char ** endptr)
{
    errno = 0;
    char * s_end;
    *number = strtof(s, &s_end);
    if ((errno == ERANGE && (*number == HUGE_VALF || *number == -HUGE_VALF)) ||
        (errno != 0 && number == 0)) {
        return errno;
    }
    if (s == s_end)
        return EINVAL;

    /* Check for a valid delimiter following the parsed number. */
    if (valid_delimiters && s_end && *s_end != '\0') {
        if (!strchr(valid_delimiters, *s_end)) {
            return EINVAL;
        }
        s_end++;
    }

    if (endptr)
        *endptr = s_end;
    return 0;
}

/**
 * `parse_double()` parses a string to produce a number that may be
 * represented as `double`.
 *
 * The number is parsed using `strtod()`, following the conventions
 * documented in the man page for that function.  In addition, some
 * further error checking is performed to ensure that the number is
 * parsed correctly.  The parsed number is stored in `number`.
 *
 * `valid_delimiters` is either `NULL`, in which case it is ignored,
 * or, it may contain a string of characters that constitute valid
 * delimiters for the parsed string.  That is, after parsing a number,
 * if there are any remaining, unconsumed characters in the string,
 * `parse_double()` checks if the next character is found in the string
 * `valid_delimiters`.  If the character is not found, then the string
 * is judged to be invalid, and `EINVAL` is returned.  Otherwise, the
 * final, delimiter character is consumed by the parser.
 *
 * If `endptr` is not `NULL`, the address stored in `endptr` points to
 * the first character beyond the characters that were consumed during
 * parsing.
 *
 * On success, `parse_double()` returns `0`.  Otherwise, if the input
 * contained invalid characters, `parse_double()` returns `EINVAL`.
 */
int parse_double(
    const char * s,
    const char * valid_delimiters,
    double * number,
    const char ** endptr)
{
    errno = 0;
    char * s_end;
    *number = strtod(s, &s_end);
    if ((errno == ERANGE && (*number == HUGE_VAL || *number == -HUGE_VAL)) ||
        (errno != 0 && number == 0)) {
        return errno;
    }
    if (s == s_end)
        return EINVAL;

    /* Check for a valid delimiter following the parsed number. */
    if (valid_delimiters && s_end && *s_end != '\0') {
        if (!strchr(valid_delimiters, *s_end)) {
            return EINVAL;
        }
        s_end++;
    }

    if (endptr)
        *endptr = s_end;
    return 0;
}
