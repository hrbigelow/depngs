#ifndef _TIMER_H
#define _TIMER_H

#include <time.h>
#include <stdio.h>

#define ELAPSED_MS \
    ((((end_time).tv_sec * 1000000000 + (end_time).tv_nsec) -           \
      ((beg_time).tv_sec * 1000000000 + (beg_time).tv_nsec)) / 1000000)

#ifdef _TIMER
#define PROGRESS_DECL() struct timespec beg_time, end_time
#define PROGRESS_START() clock_gettime(CLOCK_REALTIME, &beg_time)
#define PROGRESS_MSG(msg)                                   \
    do {                                                    \
        clock_gettime(CLOCK_REALTIME, &end_time);           \
        fprintf(stderr, "%li\t%s\n",                        \
                ELAPSED_MS,                                 \
                (msg));                                     \
        fflush(stderr);                                     \
    } while (0)
#else
#define PROGRESS_DECL()
#define PROGRESS_START()
#define PROGRESS_MSG(msg)
#endif


/* call at start of program to initialize the start time */
void
timer_init();

/* return a statically allocated string showing current time and
   elapsed time. */
const char *
timer_progress();

#endif /* _TIMER_H */
