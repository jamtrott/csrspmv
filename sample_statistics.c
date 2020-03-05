/** @file sample_statistics.c
 *
 * Descriptive statistics derived from samples.
 */

#include "sample_statistics.h"

#include <errno.h>

#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Standalone functions for computing various sample statistics.
 */

/**
 * `sample_min()` computes the sample minimum.
 */
double sample_min(int sample_size, const double * sample)
{
    if (sample_size == 0)
        return NAN;

    double min = sample[0];
    for (int i = 1; i < sample_size; i++) {
        if (min > sample[i])
            min = sample[i];
    }
    return min;
}

/**
 * `sample_max()` computes the sample maximum.
 */
double sample_max(int sample_size, const double * sample)
{
    if (sample_size == 0)
        return NAN;

    double max = sample[0];
    for (int i = 1; i < sample_size; i++) {
        if (max < sample[i])
            max = sample[i];
    }
    return max;
}

/**
 * `sample_mean()` computes the sample mean.
 */
double sample_mean(int sample_size, const double * sample)
{
    if (sample_size == 0)
        return NAN;

    double mean = 0.0;
    for (int i = 0; i < sample_size; i++) {
        mean += sample[i];
    }
    return mean / (double) sample_size;
}

/**
 * `sample_median()` computes the sample median.
 */
double sample_median(int sample_size, const double * in_sample)
{
    if (sample_size == 0)
        return NAN;

    /* Allocate temporary storage for the sample values. */
    double * sample = (double *) malloc(sample_size * sizeof(double *));
    if (!sample) {
        fprintf(stderr, "%s(): %s\n", __FUNCTION__, strerror(errno));
        return NAN;
    }
    for (int i = 0; i < sample_size; i++)
        sample[i] = in_sample[i];

    /* Sort the sample values. */
    for (int i = 0; i < sample_size; i++) {
        double x = sample[i];
        int j = i-1;
        while (j >= 0 && sample[j] > x) {
            sample[j+1] = sample[j];
            j--;
        }
        sample[j+1] = x;
    }

    /* Compute the median value. */
    double median;
    if (sample_size % 2 == 1) {
        median = sample[sample_size/2];
    } else {
        median = 0.5 * (sample[(sample_size-1)/2] + sample[sample_size/2]);
    }
    free(sample);
    return median;
}

/**
 * `sample_second_moment()` computes the second moment of a sample.
 */
double sample_second_moment(int sample_size, const double * sample)
{
    if (sample_size < 1)
        return NAN;

    double mean = sample_mean(sample_size, sample);
    double moment = 0.0;
    for (int i = 0; i < sample_size; i++)
        moment += (sample[i] - mean) * (sample[i] - mean);
    return moment / (double) sample_size;
}

/**
 * `sample_third_moment()` computes the third moment of a sample.
 */
double sample_third_moment(int sample_size, const double * sample)
{
    if (sample_size < 1)
        return NAN;

    double mean = sample_mean(sample_size, sample);
    double moment = 0.0;
    for (int i = 0; i < sample_size; i++)
        moment += (sample[i] - mean) * (sample[i] - mean) * (sample[i] - mean);
    return moment / (double) sample_size;
}

/**
 * `sample_fourth_moment()` computes the fourth moment of a sample.
 */
double sample_fourth_moment(int sample_size, const double * sample)
{
    if (sample_size < 1)
        return NAN;

    double mean = sample_mean(sample_size, sample);
    double moment = 0.0;
    for (int i = 0; i < sample_size; i++)
        moment +=
            (sample[i] - mean) * (sample[i] - mean) *
            (sample[i] - mean) * (sample[i] - mean);
    return moment / (double) sample_size;
}

/**
 * `sample_variance()` computes the sample variance.
 */
double sample_variance(int sample_size, const double * sample)
{
    if (sample_size < 2)
        return NAN;

    double mean = sample_mean(sample_size, sample);
    double variance = 0.0;
    for (int i = 0; i < sample_size; i++)
        variance += (sample[i] - mean) * (sample[i] - mean);
    return variance / (double) (sample_size - 1);
}

/**
 * `sample_standard_deviation()` computes the sample standard deviation.
 */
double sample_standard_deviation(int sample_size, const double * sample)
{
    if (sample_size < 2)
        return NAN;

    double variance = sample_variance(sample_size, sample);
    return sqrt(variance);
}

/**
 * `sample_skewness()` computes the sample skewness.
 */
double sample_skewness(int sample_size, const double * sample)
{
    if (sample_size < 2)
        return NAN;

    double second_moment = sample_second_moment(sample_size, sample);
    double third_moment = sample_third_moment(sample_size, sample);
    return third_moment / sqrt(second_moment * second_moment * second_moment);
}

/**
 * `sample_kurtosis()` computes the sample kurtosis.
 */
double sample_kurtosis(int sample_size, const double * sample)
{
    if (sample_size < 2)
        return NAN;

    double second_moment = sample_second_moment(sample_size, sample);
    double fourth_moment = sample_fourth_moment(sample_size, sample);
    return fourth_moment / (second_moment * second_moment);
}

/*
 * Collections of sample statistics.
 */

/**
 * `sample_statistics_init()` creates sample statistics for a sample
 * of the given size.
 */
int sample_statistics_init(
    struct sample_statistics * sample_statistics,
    int max_sample_size,
    const char * unit)
{
    if (max_sample_size < 0)
        return EINVAL;

    double * values = (double *) malloc(
        max_sample_size * sizeof(double));
    if (!values) {
        fprintf(stderr, "%s(): %s\n", __FUNCTION__, strerror(errno));
        return -1;
    }
    for (int i = 0; i < max_sample_size; i++)
        values[i] = 0.0;

    sample_statistics->sample_size = 0;
    sample_statistics->max_sample_size = max_sample_size;
    sample_statistics->values = values;
    sample_statistics->unit = strdup(unit);
    return 0;
}

/**
 * `sample_statistics_free()` destroys the given sample statistics.
 */
void sample_statistics_free(struct sample_statistics * sample_statistics)
{
    free(sample_statistics->unit);
    free(sample_statistics->values);
}

/**
 * `sample_statistics_add_value()` adds a sample value.
 */
int sample_statistics_add_value(
    struct sample_statistics * sample_statistics,
    double value)
{
    if (sample_statistics->sample_size >=
        sample_statistics->max_sample_size)
    {
        return EOVERFLOW;
    }

    sample_statistics->values[sample_statistics->sample_size] = value;
    sample_statistics->sample_size++;
    return 0;
}

/**
 * `sample_statistics_reset()` resets the sample statistics, zeroing
 * the sample values.
 */
void sample_statistics_reset(
    struct sample_statistics * sample_statistics)
{
    for (int i = 0; i < sample_statistics->max_sample_size; i++)
        sample_statistics->values[i] = 0.0;
    sample_statistics->sample_size = 0;
}

/**
 * `sample_statistics_print_headings()` prints a list of
 * comma-separated names, or headings, for each sample statistic.
 */
void sample_statistics_print_headings(
    const struct sample_statistics * stats,
    FILE * f,
    bool verbose)
{
    if (verbose) {
        fprintf(f, "unit,sample_size,min,max,mean,median,"
                "standard_deviation,skewness,kurtosis");
    } else {
        fprintf(f, "%s", stats->unit);
    }
}

/**
 * `sample_statistics_print()` prints sample statistics, formatted as
 * a list of comma-separated values.
 */
void sample_statistics_print(
    const struct sample_statistics * stats,
    FILE * f,
    bool verbose)
{
    if (verbose) {
        fprintf(f, "%s,%d,%g,%g,%g,%g,%g,%g,%g",
                stats->unit,
                stats->sample_size,
                sample_min(stats->sample_size, stats->values),
                sample_max(stats->sample_size, stats->values),
                sample_mean(stats->sample_size, stats->values),
                sample_median(stats->sample_size, stats->values),
                sample_standard_deviation(stats->sample_size, stats->values),
                sample_skewness(stats->sample_size, stats->values),
                sample_kurtosis(stats->sample_size, stats->values));
    } else {
        fprintf(f, "%g", sample_mean(stats->sample_size, stats->values));
    }
}
