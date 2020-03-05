/** @file perf_session.c
 *
 * Hardware performance monitoring sessions.
 */

#include "perf_session.h"
#include "perf_events.h"
#include "sample_statistics.h"

#include <errno.h>
#include <unistd.h>

#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

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
 * `perf_session_init()` creates a performance monitoring session.
 */
int perf_session_init(
    struct perf_session * perf_session,
    int num_groups,
    int * num_events_per_group,
    const char *** event_names_per_group,
    pid_t * pid_per_group,
    int * cpu_per_group,
    int * flags_per_group,
    int max_sample_size)
{
    int err;

    /* Allocate storage for sample statistics based on the duration of
     * each run. */
    struct sample_statistics * time =
        (struct sample_statistics *) malloc(
            sizeof(struct sample_statistics));
    if (!time) {
        fprintf(stderr, "%s(): %s\n", __FUNCTION__, strerror(errno));
        return -1;
    }
    err = sample_statistics_init(
        time, max_sample_size, "time (s)");
    if (err) {
        free(time);
        return -1;
    }

    /* Allocate storage for the event groups. */
    struct perf_event_group * perf_event_groups =
        (struct perf_event_group *) malloc(
            num_groups * sizeof(struct perf_event_group));
    if (!perf_event_groups) {
        fprintf(stderr, "%s(): %s\n", __FUNCTION__, strerror(errno));
        sample_statistics_free(time);
        free(time);
        return -1;
    }

    /* Create each event group. */
    for (int i = 0; i < num_groups; i++) {
        err = perf_event_group_init(
            &perf_event_groups[i],
            num_events_per_group[i],
            event_names_per_group[i],
            pid_per_group[i],
            cpu_per_group[i],
            flags_per_group[i]);
        if (err) {
            for (int ii = i-1; ii >= 0; ii--)
                perf_event_group_free(&perf_event_groups[ii]);
            free(perf_event_groups);
            sample_statistics_free(time);
            free(time);
            return -1;
        }
    }

    /* Allocate storage for sample statistics for time running ratios. */
    struct sample_statistics * time_running_ratio_per_event_group =
        (struct sample_statistics *) malloc(
            num_groups * sizeof(struct sample_statistics));
    if (!time_running_ratio_per_event_group) {
        fprintf(stderr, "%s(): %s\n", __FUNCTION__, strerror(errno));
        for (int i = num_groups-1; i >= 0; i--)
            perf_event_group_free(&perf_event_groups[i]);
        free(perf_event_groups);
        sample_statistics_free(time);
        free(time);
        return -1;
    }

    for (int i = 0; i < num_groups; i++) {
        err = sample_statistics_init(
            &time_running_ratio_per_event_group[i],
            max_sample_size, "time_running_ratio");
        if (err) {
            for (int ii = i-1; ii >= 0; ii--)
                sample_statistics_free(&time_running_ratio_per_event_group[ii]);
            free(time_running_ratio_per_event_group);
            for (int i = num_groups-1; i >= 0; i--)
                perf_event_group_free(&perf_event_groups[i]);
            free(perf_event_groups);
            sample_statistics_free(time);
            free(time);
            return err;
        }
    }

    /* Allocate storage for sample statistics for each event group. */
    struct sample_statistics ** sample_statistics_per_event_group =
        (struct sample_statistics **) malloc(
            num_groups * sizeof(struct sample_statistics *));
    if (!sample_statistics_per_event_group) {
        fprintf(stderr, "%s(): %s\n", __FUNCTION__, strerror(errno));
        for (int i = num_groups-1; i >= 0; i--)
            sample_statistics_free(&time_running_ratio_per_event_group[i]);
        free(time_running_ratio_per_event_group);
        for (int i = num_groups-1; i >= 0; i--)
            perf_event_group_free(&perf_event_groups[i]);
        free(perf_event_groups);
        sample_statistics_free(time);
        free(time);
        return -1;
    }

    /* Allocate storage for the sample statistics for each event. */
    for (int i = 0; i < num_groups; i++) {
        sample_statistics_per_event_group[i] =
            (struct sample_statistics *) malloc(
                perf_event_groups[i].num_events *
                sizeof(struct sample_statistics));
        if (!sample_statistics_per_event_group[i]) {
            fprintf(stderr, "%s(): %s\n", __FUNCTION__, strerror(errno));
            for (int ii = i-1; ii >= 0; ii--)
                free(sample_statistics_per_event_group[ii]);
            for (int ii = num_groups-1; ii >= 0; ii--)
                sample_statistics_free(&time_running_ratio_per_event_group[ii]);
            free(time_running_ratio_per_event_group);
            for (int ii = num_groups-1; ii >= 0; ii--)
                perf_event_group_free(&perf_event_groups[ii]);
            free(perf_event_groups);
            sample_statistics_free(time);
            free(time);
            return -1;
        }
    }

    /* Initialise the sample statistics for each event. */
    for (int i = 0; i < num_groups; i++) {
        for (int j = 0; j < perf_event_groups[i].num_events; j++) {
            err = sample_statistics_init(
                &sample_statistics_per_event_group[i][j],
                max_sample_size,
                perf_event_groups[i].events[j].name);
            if (err) {
                for (int jj = j-1; jj >= 0; jj--)
                    sample_statistics_free(&sample_statistics_per_event_group[i][jj]);
                for (int ii = i-1; ii >= 0; ii--) {
                    for (int jj = perf_event_groups[ii].num_events-1; jj >= 0; jj--)
                        sample_statistics_free(&sample_statistics_per_event_group[ii][jj]);
                }
                for (int ii = num_groups-1; i >= 0; i--)
                    free(sample_statistics_per_event_group[ii]);
                for (int ii = num_groups-1; ii >= 0; ii--)
                    sample_statistics_free(&time_running_ratio_per_event_group[ii]);
                free(time_running_ratio_per_event_group);
                for (int ii = num_groups-1; ii >= 0; ii--)
                    perf_event_group_free(&perf_event_groups[ii]);
                free(perf_event_groups);
                sample_statistics_free(time);
                free(time);
                return err;
            }
        }
    }

    perf_session->time = time;
    perf_session->num_perf_event_groups = num_groups;
    perf_session->perf_event_groups = perf_event_groups;
    perf_session->time_running_ratio_per_event_group =
        time_running_ratio_per_event_group;
    perf_session->sample_statistics_per_event_group =
        sample_statistics_per_event_group;
    perf_session->max_sample_size = max_sample_size;
    perf_session->sample_size = 0;
    return 0;
}

/**
 * `perf_session_free()` destroys a performance monitoring session.
 */
void perf_session_free(
    struct perf_session * perf_session)
{
    for (int i = perf_session->num_perf_event_groups-1; i >= 0; i--) {
        for (int j = perf_session->perf_event_groups[i].num_events-1; j >= 0; j--) {
            sample_statistics_free(
                &perf_session->sample_statistics_per_event_group[i][j]);
        }
    }
    for (int i = perf_session->num_perf_event_groups-1; i >= 0; i--)
        free(perf_session->sample_statistics_per_event_group[i]);
    free(perf_session->sample_statistics_per_event_group);
    for (int i = perf_session->num_perf_event_groups-1; i >= 0; i--)
        sample_statistics_free(&perf_session->time_running_ratio_per_event_group[i]);
    free(perf_session->time_running_ratio_per_event_group);
    for (int i = perf_session->num_perf_event_groups-1; i >= 0; i--)
        perf_event_group_free(&perf_session->perf_event_groups[i]);
    free(perf_session->perf_event_groups);
    sample_statistics_free(perf_session->time);
    free(perf_session->time);
}

/**
 * `perf_session_enable()` enables the event groups in the session,
 * scheduling them onto a PMU using time slicing.
 */
int perf_session_enable(
    struct perf_session * perf_session)
{
    /* Reset the event groups. */
    for (int i = 0; i < perf_session->num_perf_event_groups; i++)
        perf_event_group_reset(&perf_session->perf_event_groups[i]);

    /* Enable the event groups and start the clock. */
    for (int i = 0; i < perf_session->num_perf_event_groups; i++)
        perf_event_group_enable(&perf_session->perf_event_groups[i]);
    clock_gettime(CLOCK_MONOTONIC, &perf_session->start_time);
    return 0;
}

/**
 * `perf_session_read()` reads and updates the measurements or counts
 * associated with the event groups in the session, and updates the
 * sample statistics derived from the measurements.
 */
static int perf_session_read(
    struct perf_session * perf_session,
    struct timespec end_time)
{
    /* Update the session duration. */
    sample_statistics_add_value(
        perf_session->time,
        timespec_duration(perf_session->start_time, end_time));

    for (int i = 0; i < perf_session->num_perf_event_groups; i++) {
        /* Read the values of the hardware performance events. */
        struct perf_event_group * perf_event_group =
            &perf_session->perf_event_groups[i];
        perf_event_group_read(perf_event_group);

        /* Update the sample statistics. */
        if (perf_event_group->time_enabled > 0) {
            double time_running_ratio =
                (double) perf_event_group->time_running /
                (double) perf_event_group->time_enabled;
            sample_statistics_add_value(
                &perf_session->time_running_ratio_per_event_group[i],
                time_running_ratio);

            if (perf_event_group->time_running > 0) {
                for (int j = 0; j < perf_event_group->num_events; j++) {
                    struct sample_statistics * sample_statistics =
                        &perf_session->sample_statistics_per_event_group[i][j];
                    sample_statistics_add_value(
                        sample_statistics,
                        perf_event_group->event_values[j] / time_running_ratio);
                }
            }
        }
    }
    perf_session->sample_size++;
    return 0;
}

/**
 * `perf_session_disable()` disables the event groups in the session,
 * removing them from the PMU to stop counting or measuring events.
 *
 * In addition, the values associated with the hardware performance
 * event groups are read, and the corresponding sample statistics
 * are updated.
 */
int perf_session_disable(
    struct perf_session * perf_session)
{
    struct timespec end_time;
    clock_gettime(CLOCK_MONOTONIC, &end_time);
    for (int i = perf_session->num_perf_event_groups-1; i >= 0; i--)
        perf_event_group_disable(&perf_session->perf_event_groups[i]);

    return perf_session_read(perf_session, end_time);
}

/**
 * `perf_session_reset()` resets the event groups and statistics
 * associated with the session.
 */
int perf_session_reset(
    struct perf_session * perf_session)
{
    sample_statistics_reset(perf_session->time);

    for (int i = 0; i < perf_session->num_perf_event_groups; i++) {
        struct perf_event_group * perf_event_group =
            &perf_session->perf_event_groups[i];
        perf_event_group_reset(perf_event_group);
        for (int j = 0; j < perf_event_group->num_events; j++) {
            struct sample_statistics * sample_statistics =
                &perf_session->sample_statistics_per_event_group[i][j];
            sample_statistics_reset(sample_statistics);
        }
    }
    return 0;
}

/**
 * `perf_session_print_headings()` prints a list of comma-separated
 * names, or headings, describing each column of values that is
 * printed for the given performance monitoring session.
 */
void perf_session_print_headings(
    const struct perf_session * perf_session,
    FILE * f,
    bool verbose)
{
    sample_statistics_print_headings(perf_session->time, f, verbose);
    if (verbose)
        fprintf(f, ",num_perf_event_groups");
    for (int i = 0; i < perf_session->num_perf_event_groups; i++) {
        if (verbose) {
            fprintf(f, ",");
            sample_statistics_print_headings(
                &perf_session->time_running_ratio_per_event_group[i], f,
                verbose);
        }

        const struct perf_event_group * perf_event_group =
            &perf_session->perf_event_groups[i];
        for (int j = 0; j < perf_event_group->num_events; j++) {
            const struct sample_statistics * sample_statistics =
                &perf_session->sample_statistics_per_event_group[i][j];
            fprintf(f, ",");
            sample_statistics_print_headings(
                sample_statistics, f, verbose);
        }
    }
}

/**
 * `perf_session_print()` prints values associated with the
 * performance monitoring session, formatted as a list of
 * comma-separated values.
 */
void perf_session_print(
    const struct perf_session * perf_session,
    FILE * f,
    bool verbose)
{
    sample_statistics_print(perf_session->time, f, verbose);
    if (verbose)
        fprintf(f, ",%d", perf_session->num_perf_event_groups);
    for (int i = 0; i < perf_session->num_perf_event_groups; i++) {
        if (verbose) {
            fprintf(f, ",");
            sample_statistics_print(
                &perf_session->time_running_ratio_per_event_group[i], f,
                verbose);
        }

        const struct perf_event_group * perf_event_group =
            &perf_session->perf_event_groups[i];
        for (int j = 0; j < perf_event_group->num_events; j++) {
            const struct sample_statistics * sample_statistics =
                &perf_session->sample_statistics_per_event_group[i][j];
            fprintf(f, ",");
            sample_statistics_print(
                sample_statistics, f, verbose);
        }
    }
}
