/** @file perf_session.h
 *
 * Hardware performance monitoring sessions.
 */

#ifndef PERF_SESSION_H
#define PERF_SESSION_H

#include <sys/types.h>
#include <stdbool.h>
#include <stdio.h>
#include <time.h>

struct perf_event_group;
struct sample_statistics;

/**
 * `perf_session` is a collection of hardware performance event groups
 * together with sample statistics derived from the measurements of
 * each event group.
 *
 * For the sake of convenience, it is often useful to work with
 * multiple groups of hardware performance events, so that the groups
 * may be enabled simultaneously.  When they are enabled, a time
 * multiplexing scheme is used to schedule each event group onto a
 * hardware performance monitoring unit.
 *
 * A session may be run multiple times, and for each run, the duration
 * of the run is measured, as well as the values of the session's
 * hardware performance events.  Based on the values acquired during
 * each run, various sample statistics are provided.
 */
struct perf_session
{
    /**
     * `start_time` is the time at which the current run of the
     * session was started.
     */
    struct timespec start_time;

    /**
     * `time` is a collection of sample statistics derived from the
     * the duration of each run of the session.
     */
    struct sample_statistics * time;

    /**
     * `num_perf_event_groups` is the number of hardware performance
     * event groups in the performance monitoring session.
     */
    int num_perf_event_groups;

    /**
     * `perf_event_groups` is an array containing the performance
     * event groups of the performance monitoring session,
     * `i=0,1,...,num_perf_event_groups-1`.
     */
    struct perf_event_group * perf_event_groups;

    /**
     * `time_running_ratio_per_event_group` is an array containing
     * sample statistics representing, for each performance event
     * group, the ratio of time it is running on a PMU to the time it
     * is enabled.  This is useful for extrapolating from measurements
     * whenever event groups are time multiplexed onto PMUs.
     */
    struct sample_statistics * time_running_ratio_per_event_group;

    /**
     * `sample_statistics_per_event_group` is an array containing the
     * sample statistics derived from the hardware performance counter
     * measurements of each event group.
     *
     * For the `i`-th group, `sample_statistics_per_event_group[i][j]`
     * are the sample statistics for the `j`-th event of the group,
     * for `j=0,1,...,perf_event_groups[i].num_events-1`, and
     * `i=0,1,...,num_perf_event_groups-1`.
     */
    struct sample_statistics ** sample_statistics_per_event_group;

    /**
     * `max_sample_size` is the maximum number of times the values of
     * the hardware performance event groups may be read.
     */
    int max_sample_size;

    /**
     * `sample_size` is the number of times the values of the hardware
     * performance event groups have been read.
     */
    int sample_size;
};

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
    int max_sample_size);

/**
 * `perf_session_free()` destroys a performance monitoring session.
 */
void perf_session_free(
    struct perf_session * perf_session);

/**
 * `perf_session_enable()` enables the event groups in the session,
 * scheduling them onto a PMU using time slicing.
 */
int perf_session_enable(
    struct perf_session * session);

/**
 * `perf_session_disable()` disables the event groups in the session,
 * removing them from the PMU to stop counting or measuring events.
 *
 * In addition, the values associated with the hardware performance
 * event groups are read, and the corresponding sample statistics
 * are updated.
 */
int perf_session_disable(
    struct perf_session * session);

/**
 * `perf_session_reset()` resets the event groups and statistics
 * associated with the session.
 */
int perf_session_reset(
    struct perf_session * session);

/**
 * `perf_session_print_headings()` prints a list of comma-separated
 * names, or headings, describing each column of values that is
 * printed for the given performance monitoring session.
 */
void perf_session_print_headings(
    const struct perf_session * perf_session,
    FILE * f,
    bool verbose);

/**
 * `perf_session_print()` prints values associated with the
 * performance monitoring session, formatted as a list of
 * comma-separated values.
 */
void perf_session_print(
    const struct perf_session * perf_session,
    FILE * f,
    bool verbose);

#endif
