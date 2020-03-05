/** @file perf_events.h
 *
 * Hardware performance monitoring.
 */

#ifndef PERF_EVENTS_H
#define PERF_EVENTS_H

#include <sys/types.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>

/*
 * libpfm.
 */

/**
 * `libpfm_init()` initialises the libpfm library.
 */
int libpfm_init(void);

/**
 * `libpfm_free()` frees resources associated with the libpfm library.
 */
void libpfm_free(void);

/**
 * `libpfm_print_events()` prints a list of available events.
 */
int libpfm_print_events(FILE * f);

/*
 * Hardware performance events.
 */

/**
 * `perf_event` is a hardware performance event.
 */
struct perf_event
{
    /**
     * `name` is the name of the event.
     */
    char * name;

    /**
     * `fstr` is a string containing the fully qualified event name.
     */
    char * fstr;

    /**
     * `pid` is the process id of the event.
     */
    pid_t pid;

    /**
     * `cpu` is the processor that the event belongs to.
     */
    int cpu;

    /**
     * `groupfd` is a file descriptor corresponding to the leader of
     * the event group to which the event belongs, or `-1`, if the
     * event does not belong to an event group or is itself the leader
     * of an event group.
     */
    int groupfd;

    /**
     * `flags` are the flags that are applied to the event whenever
     * the event is opened with `perf_event_open()`.
     */
    int flags;

    /**
     * `fd` is a file descriptor used to access the event.
     */
    int fd;

    /**
     * `id` is a unique identifier used to distinguish the event from
     * other events.
     */
    uint64_t id;
};

/**
 * `perf_event_init()` creates a hardware performance event.
 */
int perf_event_init(
    struct perf_event * event,
    const char * name,
    pid_t pid,
    int cpu,
    int groupfd,
    int flags);

/**
 * `perf_event_free()` destroys an event group.
 */
void perf_event_free(
    struct perf_event * event);

/*
 * Performance event groups.
 */

/**
 * `perf_event_group` is a group of hardware performance events that
 * may be scheduled together onto a Performance Monitoring Unit (PMU).
 */
struct perf_event_group
{
    /**
     * `num_events` is the number of events in the group.
     */
    int num_events;

    /**
     * `events` is an array of events, `i=0,1,...,num_events-1`.
     */
    struct perf_event * events;

    /**
     * `enabled` is true if the event is currently enabled, and false
     * otherwise.
     */
    bool enabled;

    /**
     * `time_enabled` is the duration of time for which the event
     * group has been enabled.
     */
    uint64_t time_enabled;

    /**
     * `time_running` is the duration of time for which the event
     * group has been running, that is, the event group has been
     * scheduled onto a PMU.
     */
    uint64_t time_running;

    /**
     * `event_values` is an array containing the measurements or
     * counts associated with each event, `i=0,1,...,num_events-1`.
     */
    uint64_t * event_values;

    /**
     * `read_buffer_size` is the size, in bytes, of the read buffer.
     */
    int read_buffer_size;

    /**
     * `read_buffer` is a buffer that is used to read the measurements
     * or event counts for the group's events.
     */
    char * read_buffer;
};

/**
 * `perf_event_group_init()` creates an event group from given events.
 */
int perf_event_group_init(
    struct perf_event_group * event_group,
    int num_events,
    const char ** event_names,
    pid_t pid,
    int cpu,
    int flags);

/**
 * `perf_event_group_free()` destroys an event group.
 */
void perf_event_group_free(
    struct perf_event_group * event_group);

/**
 * `perf_event_group_enable()` enables an event group and schedules it
 * onto a PMU for its events to be counted or measured.
 */
int perf_event_group_enable(
    struct perf_event_group * event_group);

/**
 * `perf_event_group_disable()` disables an event group, removing it
 * from the PMU to stop counting or measuring its events.
 */
int perf_event_group_disable(
    struct perf_event_group * event_group);

/**
 * `perf_event_group_reset()` resets an event group, zeroing any of
 * the measurements or counts associated with its events.
 */
int perf_event_group_reset(
    struct perf_event_group * event_group);

/**
 * `perf_event_group_read()` reads and updates the measurements or
 * counts associated with a group's events.
 */
int perf_event_group_read(
    struct perf_event_group * event_group);

#endif
