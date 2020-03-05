/** @file sample_statistics.h
 *
 * Descriptive statistics derived from samples.
 */

#ifndef SAMPLE_STATISTICS_H
#define SAMPLE_STATISTICS_H

#include <stdbool.h>
#include <stdio.h>

/*
 * Standalone functions for computing various sample statistics.
 */

/**
 * `sample_min()` computes the sample minimum.
 */
double sample_min(
    int sample_size,
    const double * sample);

/**
 * `sample_max()` computes the sample maximum.
 */
double sample_max(
    int sample_size,
    const double * sample);

/**
 * `sample_mean()` computes the sample mean.
 */
double sample_mean(
    int sample_size,
    const double * sample);

/**
 * `sample_median()` computes the sample median.
 */
double sample_median(
    int sample_size,
    const double * sample);

/**
 * `sample_second_moment()` computes the second moment of a sample.
 */
double sample_second_moment(
    int sample_size,
    const double * sample);

/**
 * `sample_third_moment()` computes the third moment of a sample.
 */
double sample_third_moment(
    int sample_size,
    const double * sample);

/**
 * `sample_fourth_moment()` computes the fourth moment of a sample.
 */
double sample_fourth_moment(
    int sample_size,
    const double * sample);

/**
 * `sample_variance()` computes the sample variance.
 */
double sample_variance(
    int sample_size,
    const double * sample);

/**
 * `sample_standard_deviation()` computes the sample standard deviation.
 */
double sample_standard_deviation(
    int sample_size,
    const double * sample);

/**
 * `sample_skewness()` computes the sample skewness.
 */
double sample_skewness(
    int sample_size,
    const double * sample);

/**
 * `sample_kurtosis()` computes the sample kurtosis.
 */
double sample_kurtosis(
    int sample_size,
    const double * sample);

/*
 * Collections of sample statistics.
 */

/**
 * `sample_statistics` represents a sample of double-precision
 * floating point numbers drawn from a population, which may be used
 * to compute various sample statistics.
 */
struct sample_statistics
{
    /**
     * `sample_size` is the size of the sample.
     */
    int sample_size;

    /**
     * `max_sample_size` is the maximum size of the sample.
     */
    int max_sample_size;

    /**
     * `values` is an array containing the sample values drawn from
     * the population, for `i=0,1,...,max_sample_size-1`.
     */
    double * values;

    /**
     * `unit` is a string used to represent the unit associated with
     * the sample values.
     */
    char * unit;
};

/**
 * `sample_statistics_init()` creates sample statistics for a sample
 * of the given size.
 */
int sample_statistics_init(
    struct sample_statistics * sample_statistics,
    int max_sample_size,
    const char * unit);

/**
 * `sample_statistics_free()` destroys the given sample statistics.
 */
void sample_statistics_free(
    struct sample_statistics * sample_statistics);

/**
 * `sample_statistics_add_value()` adds a sample value.
 */
int sample_statistics_add_value(
    struct sample_statistics * sample_statistics,
    double value);

/**
 * `sample_statistics_reset()` resets the sample statistics, zeroing
 * the sample values.
 */
void sample_statistics_reset(
    struct sample_statistics * sample_statistics);

/**
 * `sample_statistics_print_headings()` prints a list of
 * comma-separated names, or headings, for each sample statistic.
 */
void sample_statistics_print_headings(
    const struct sample_statistics * sample_statistics,
    FILE * f,
    bool verbose);

/**
 * `sample_statistics_print()` prints sample statistics, formatted as
 * a list of comma-separated values.
 */
void sample_statistics_print(
    const struct sample_statistics * sample_statistics,
    FILE * f,
    bool verbose);

#endif
