/** @file perf_events.c
 *
 * Hardware performance monitoring.
 */

#include "perf_events.h"

#ifdef HAVE_LIBPFM
#include <perfmon/pfmlib.h>
#include <perfmon/pfmlib_perf_event.h>
#endif

#include <errno.h>
#include <unistd.h>

#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * libpfm.
 */

/**
 * `libpfm_init()` initialises the libpfm library.
 */
int libpfm_init(void)
{
#ifdef HAVE_LIBPFM
    int err = pfm_initialize();
    if (err != PFM_SUCCESS) {
        fprintf(stderr, "pfm_initialize(): %s\n", pfm_strerror(err));
        return -1;
    }
    return 0;
#else
    fprintf(stderr, "%s(): Please re-build with libpfm enabled\n",
            __FUNCTION__);
    return ENOTSUP;
#endif
}

/**
 * `libpfm_free()` frees resources associated with the libpfm library.
 */
void libpfm_free(void)
{
#ifdef HAVE_LIBPFM
    pfm_terminate();
#endif
}

#ifdef HAVE_LIBPFM
/**
 * `print_text_wrapped()` prints a string by wrapping the text, if it
 * is too long to be printed on a single line.
 */
static void print_text_wrapped(
    FILE * f,
    int width,
    int indent,
    int start_pos,
    const char * s)
{
    int num_lines = 0;
    if (width <= 0)
        return;

    if (width <= start_pos) {
        fprintf(f, "\n%*s", indent, "");
        num_lines++;
    } else {
        width -= start_pos;
    }

    size_t len = strlen(s);
    while (len > 0) {
        /* Locate the final whitespace before exceeding the width. */
        int line_len = len;
        bool hyphenate = false;
        bool break_line = false;
        if (line_len > width) {
            break_line = true;
            const char * prev_whitespace = s;
            const char * next_whitespace = strchr(s+1, ' ');
            while (next_whitespace && (next_whitespace - s < width)) {
                const char * t = prev_whitespace + 1;
                prev_whitespace = next_whitespace;
                next_whitespace = strchr(t, ' ');
            }
            if (prev_whitespace == NULL || prev_whitespace == s) {
                line_len = width;
                hyphenate = true;
            } else {
                line_len = prev_whitespace - s;
            }
        }

        fprintf(f, "%.*s", line_len, s);
        s += line_len;
        len -= line_len;
        if (break_line) {
            if (hyphenate) {
                fprintf(f, "-");
            } else {
                s++;
                len--;
            }
            fprintf(f, "\n%*s", indent, "");
            if (num_lines == 0)
                width += start_pos;
            num_lines++;
        }
    }
}
#endif

/**
 * `libpfm_print_events()` prints a list of available events.
 */
int libpfm_print_events(FILE * f)
{
#ifdef HAVE_LIBPFM
    pfm_err_t err;
    int first_column_width = 24;
    int width = 80;
    int tab_width = 4;
    int indent = 0;
    for (pfm_pmu_t pmu = PFM_PMU_NONE; pmu < PFM_PMU_MAX; pmu++) {
        pfm_pmu_info_t pmu_info;
        memset(&pmu_info, 0, sizeof(pmu_info));
        err = pfm_get_pmu_info(pmu, &pmu_info);
        if (err == PFM_ERR_NOTSUPP) {
            continue;
        } else if (err != PFM_SUCCESS) {
            fprintf(stderr, "pfm_get_pmu_info(): %s\n", pfm_strerror(err));
            return -1;
        }
        if (!pmu_info.is_present)
            continue;

        for (int event = pmu_info.first_event;
             event != -1;
             event = pfm_get_event_next(event))
        {
            pfm_event_info_t event_info;
            memset(&event_info, 0, sizeof(event_info));
            err = pfm_get_event_info(event, PFM_OS_NONE, &event_info);
            if (err != PFM_SUCCESS) {
                fprintf(stderr, "pfm_get_event_info(): %s\n", pfm_strerror(err));
                return -1;
            }

            indent = 0;
            fprintf(f, "%*s%s::%s\n", indent, "",
                    pmu_info.name, event_info.name);
            indent = 4;
            fprintf(f, "%*s", indent, "");
            print_text_wrapped(f, width, indent, 0, event_info.desc);
            fprintf(f, "\n");

            for (int i = 0; i < event_info.nattrs; i++) {
                pfm_event_attr_info_t event_attr;
                memset(&event_attr, 0, sizeof(event_attr));
                err = pfm_get_event_attr_info(
                    event, i, PFM_OS_NONE, &event_attr);
                if (err != PFM_SUCCESS) {
                    fprintf(stderr, "pfm_get_event_attr_info(): %s\n",
                            pfm_strerror(err));
                    return -1;
                }

                indent = 0;
                fprintf(f, "%*s", first_column_width, event_attr.name);

                size_t event_attr_name_len = strlen(event_attr.name);
                int w = first_column_width > event_attr_name_len ?
                    first_column_width : event_attr_name_len;
                indent = w + 2*tab_width - w % tab_width;
                fprintf(f, "%*s", 2*tab_width - w % tab_width, "");
                print_text_wrapped(f, width - indent, indent, 0, event_attr.desc);
                fprintf(f, "\n");
            }
            fprintf(f, "\n");
        }
        fprintf(f, "\n");
    }
    return 0;
#else
    fprintf(stderr, "%s(): Please re-build with libpfm enabled\n",
            __FUNCTION__);
    return ENOTSUP;
#endif
}

#ifdef HAVE_LIBPFM
/**
 * `libpfm_event_encoding()` creates a `perf_event_attr` object from a
 * given string that desribes a hardware performance event.
 */
static int libpfm_event_encoding(
    struct perf_event_attr * event_attr,
    const char * event_name,
    char ** fstr)
{
    int err;
    memset(event_attr, 0, sizeof(struct perf_event_attr));
    pfm_perf_encode_arg_t pfm_perf_encode_args;
    pfm_perf_encode_args.attr = event_attr;
    pfm_perf_encode_args.fstr = fstr;
    pfm_perf_encode_args.size = sizeof(pfm_perf_encode_arg_t);
    pfm_perf_encode_args.idx = 0;
    pfm_perf_encode_args.cpu = 0;
    pfm_perf_encode_args.flags = 0;
    err = pfm_get_os_event_encoding(
        event_name, PFM_PLM3, PFM_OS_PERF_EVENT,
        &pfm_perf_encode_args);
    if (err != PFM_SUCCESS) {
        fprintf(stderr,
                "pfm_get_os_event_encoding(): %s: %s\n",
                event_name, pfm_strerror(err));
        return -1;
    }
    return 0;
}
#endif

/*
 * Hardware performance events.
 */

/**
 * `perf_event_init()` creates a hardware performance event.
 */
int perf_event_init(
    struct perf_event * event,
    const char * event_name,
    pid_t pid,
    int cpu,
    int groupfd,
    int flags)
{
#ifdef HAVE_LIBPFM
    int err;

    /* Use libpfm to obtain an encoding of the event from its name. */
    struct perf_event_attr event_attr;
    char * fstr = NULL;
    err = libpfm_event_encoding(
        &event_attr, event_name, &fstr);
    if (err) {
        if (fstr)
            free(fstr);
        return err;
    }

    /* Use `perf_event_open()` to obtain a file descriptor for the
     * hardware performance event. */
    event_attr.read_format =
        PERF_FORMAT_GROUP
        | PERF_FORMAT_ID
        | PERF_FORMAT_TOTAL_TIME_ENABLED
        | PERF_FORMAT_TOTAL_TIME_RUNNING;
    event_attr.disabled = 1;
    int fd = perf_event_open(
        &event_attr, pid, cpu, groupfd, flags);
    if (fd < 0) {
        fprintf(stderr, "perf_event_open(): %s: %s\n",
                event_name, strerror(errno));
        if (fstr)
            free(fstr);
        return -1;
    }

    /* Obtain a unique identifier for the event. */
    uint64_t id;
    err = ioctl(fd, PERF_EVENT_IOC_ID, &id);
    if (err) {
        fprintf(stderr, "ioctl(): %s: %s\n",
                event_name, strerror(errno));
        close(fd);
        if (fstr)
            free(fstr);
        return -1;
    }

    event->name = strdup(event_name);
    event->fstr = fstr;
    event->pid = pid;
    event->cpu = cpu;
    event->groupfd = groupfd;
    event->flags = flags;
    event->fd = fd;
    event->id = id;
    return 0;
#else
    fprintf(stderr, "%s(): Please re-build with libpfm enabled\n",
            __FUNCTION__);
    return ENOTSUP;
#endif
}

/**
 * `perf_event_free()` destroys an event group.
 */
void perf_event_free(
    struct perf_event * event)
{
    close(event->fd);
    free(event->fstr);
    free(event->name);
}

/*
 * Performance event groups.
 */

/**
 * `perf_event_group_read_format` is used to read event values from a
 * file descriptor corresponding to an event group leader.
 */
struct perf_event_group_read_format
{
    uint64_t num_events;    /* The number of events */
    uint64_t time_enabled;  /* if PERF_FORMAT_TOTAL_TIME_ENABLED */
    uint64_t time_running;  /* if PERF_FORMAT_TOTAL_TIME_RUNNING */
    struct {
        uint64_t value;     /* The value of the event */
        uint64_t id;        /* if PERF_FORMAT_ID */
    } values[];
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
    int flags)
{
    int err;

    /* Allocate storage for events. */
    struct perf_event * events = (struct perf_event *)
        malloc(num_events * sizeof(struct perf_event));
    if (!events) {
        fprintf(stderr, "%s(): %s\n", __FUNCTION__, strerror(errno));
        return -1;
    }

    /* Create the first event, which will become the group leader. */
    int groupfd = -1;
    if (num_events > 0) {
        err = perf_event_init(
            &events[0], event_names[0],
            pid, cpu, groupfd, flags);
        if (err) {
            free(events);
            return err;
        }
        groupfd = events[0].fd;
    }

    /* Create the remaining events. */
    for (int i = 1; i < num_events; i++) {
        err = perf_event_init(
            &events[i], event_names[i],
            pid, cpu, groupfd, flags);
        if (err) {
            for (int j = i-1; j >= 0; j--)
                perf_event_free(&events[j]);
            free(events);
            return err;
        }
    }

#ifdef HAVE_LIBPFM
    /* Reset the event group. */
    if (num_events > 0) {
        err = ioctl(
            groupfd,
            PERF_EVENT_IOC_RESET,
            PERF_IOC_FLAG_GROUP);
        if (err < 0) {
            fprintf(stderr, "ioctl(%d, %s): %s: %s\n",
                    groupfd, "PERF_EVENT_IOC_RESET",
                    event_names[0], strerror(errno));
            for (int i = 0; i < num_events; i++)
                perf_event_free(&events[i]);
            free(events);
            return -1;
        }
    }
#endif

    /* Allocate storage for event values. */
    uint64_t * event_values = (uint64_t *) malloc(
        num_events * sizeof(uint64_t));
    if (!event_values) {
        fprintf(stderr, "%s(): %s\n", __FUNCTION__, strerror(errno));
        for (int i = 0; i < num_events; i++)
            perf_event_free(&events[i]);
        free(events);
        return -1;
    }
    for (int i = 0; i < num_events; i++)
        event_values[i] = 0;

    /* Allocate storage for the buffer used to read event counts. */
    int read_buffer_size =
        sizeof(struct perf_event_group_read_format) +
        2 * num_events * sizeof(uint64_t);
    char * read_buffer = (char *) malloc(read_buffer_size);
    if (!read_buffer) {
        fprintf(stderr, "%s(): %s\n", __FUNCTION__, strerror(errno));
        free(event_values);
        for (int i = 0; i < num_events; i++)
            perf_event_free(&events[i]);
        free(events);
        return -1;
    }

    event_group->num_events = num_events;
    event_group->events = events;
    event_group->enabled = false;
    event_group->time_enabled = 0;
    event_group->time_running = 0;
    event_group->event_values = event_values;
    event_group->read_buffer_size = read_buffer_size;
    event_group->read_buffer = read_buffer;
    return 0;
}

/**
 * `perf_event_group_free()` destroys an event group.
 */
void perf_event_group_free(
    struct perf_event_group * event_group)
{
    free(event_group->read_buffer);
    free(event_group->event_values);
    for (int i = 0; i < event_group->num_events; i++)
        perf_event_free(&event_group->events[i]);
    free(event_group->events);
}

/**
 * `perf_event_group_enable()` enables an event group and schedules it
 * into a PMU for its events to be counted or measured.
 */
int perf_event_group_enable(
    struct perf_event_group * event_group)
{
    if (event_group->enabled)
        return 0;
    if (event_group->num_events <= 0) {
        event_group->enabled = true;
        return 0;
    }

#ifdef HAVE_LIBPFM
    int err = ioctl(
        event_group->events[0].fd,
        PERF_EVENT_IOC_ENABLE,
        PERF_IOC_FLAG_GROUP);
    if (err < 0) {
        fprintf(stderr, "ioctl(%d, %s): %s: %s\n",
                event_group->events[0].fd, "PERF_EVENT_IOC_ENABLE",
                event_group->events[0].name, strerror(errno));
        return -1;
    }
    event_group->enabled = true;
    return 0;
#else
    fprintf(stderr, "%s(): Please re-build with libpfm enabled\n",
            __FUNCTION__);
    return ENOTSUP;
#endif
}

/**
 * `perf_event_group_disable()` disables an event group, removing it
 * from the PMU to stop counting or measuring its events.
 */
int perf_event_group_disable(
    struct perf_event_group * event_group)
{
    if (!event_group->enabled)
        return 0;
    if (event_group->num_events <= 0) {
        event_group->enabled = false;
        return 0;
    }

#ifdef HAVE_LIBPFM
    int err = ioctl(
        event_group->events[0].fd,
        PERF_EVENT_IOC_DISABLE,
        PERF_IOC_FLAG_GROUP);
    if (err < 0) {
        fprintf(stderr, "ioctl(%d, %s): %s: %s\n",
                event_group->events[0].fd, "PERF_EVENT_IOC_DISABLE",
                event_group->events[0].name, strerror(errno));
        return -1;
    }
    event_group->enabled = false;
    return 0;
#else
    fprintf(stderr, "%s(): Please re-build with libpfm enabled\n",
            __FUNCTION__);
    return ENOTSUP;
#endif
}

/**
 * `perf_event_group_reset()` resets an event group, zeroing any of
 * the measurements or counts associated with its events.
 */
int perf_event_group_reset(
    struct perf_event_group * event_group)
{
    if (event_group->num_events <= 0)
        return 0;

#ifdef HAVE_LIBPFM
    int err = ioctl(
        event_group->events[0].fd,
        PERF_EVENT_IOC_RESET,
        PERF_IOC_FLAG_GROUP);
    if (err < 0) {
        fprintf(stderr, "ioctl(%d, %s): %s: %s\n",
                event_group->events[0].fd, "PERF_EVENT_IOC_RESET",
                event_group->events[0].name, strerror(errno));
        return -1;
    }
    return 0;
#else
    fprintf(stderr, "%s(): Please re-build with libpfm enabled\n",
            __FUNCTION__);
    return ENOTSUP;
#endif
}

/**
 * `perf_event_group_read()` reads and updates the measurements or
 * counts associated with a group's events.
 */
int perf_event_group_read(
    struct perf_event_group * event_group)
{
    if (event_group->num_events <= 0)
        return 0;

    /* Read from the file descriptor of the event group leader. */
    ssize_t bytes_read = read(
        event_group->events[0].fd,
        event_group->read_buffer,
        event_group->read_buffer_size);
    if (bytes_read == -1) {
        fprintf(stderr, "%s(): %s: read(): %s\n",
                __FUNCTION__, event_group->events[0].name,
                strerror(errno));
        return -1;
    }

    const struct perf_event_group_read_format * read_format =
        (const struct perf_event_group_read_format *)
        event_group->read_buffer;
    if (event_group->num_events > read_format->num_events) {
        fprintf(stderr, "%s(): %s: Expected at most %d events\n",
                __FUNCTION__, event_group->events[0].name,
                event_group->num_events);
        return -1;
    }

    /* Update the event counts. */
    for (int i = 0; i < read_format->num_events; i++) {
        uint64_t value = read_format->values[i].value;
        uint64_t id = read_format->values[i].id;

        for (int j = 0; j < event_group->num_events; j++) {
            const struct perf_event * event = &event_group->events[j];
            if (event->id == id) {
                event_group->event_values[j] = value;
            }
        }
    }
    event_group->time_enabled = read_format->time_enabled;
    event_group->time_running = read_format->time_running;
    return 0;
}
