#include "timer.h"
#include <string.h>

static struct timespec program_start_time;

static char progress_string[200];

void
timer_init()
{
    clock_gettime(CLOCK_REALTIME, &program_start_time);
}


const char *
timer_progress()
{
    struct timespec now;
    clock_gettime(CLOCK_REALTIME, &now);
    unsigned elapsed = now.tv_sec - program_start_time.tv_sec;
    time_t cal = time(NULL);
    
    strcpy(progress_string, ctime(&cal));
    progress_string[strlen(progress_string) - 1] = '\0';
    sprintf(progress_string + strlen(progress_string), 
            " (%02i:%02i:%02i elapsed)", 
            elapsed / 3600,
            (elapsed % 3600) / 60,
            elapsed % 60);
    return progress_string;
}
