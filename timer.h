#include <time.h>

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
