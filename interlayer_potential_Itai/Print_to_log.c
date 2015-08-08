#include "declarations.h"

void Print_to_log_begin(char *string, FILE* logfile, time_t &t0, clock_t &c0, int print)
{
  if(print){
    t0 = time(NULL);
    c0 = clock();
    
    fprintf(logfile,"***************************************\n");
    fprintf(logfile,"Begining the calculation of %s\n",string);
    fprintf(logfile,"***************************************\n\n");
    fprintf(logfile,"Begining wall time:              %ld\n", (long) t0);
    fprintf(logfile,"Begining CPU time:               %d\n\n", (int) c0);
    fflush(logfile);
  }
}

void Print_to_log_end(char *string, FILE* logfile, time_t &t0, clock_t &c0, int print)
{ 
  int sec, min, hour, day;
  if(print){
    time_t  t1 = time(NULL);
    clock_t c1 = clock();
    
    fprintf (logfile,"Ending wall:                    %ld\n", (long) t1);
    fprintf (logfile,"Ending CPU:                     %d\n\n", (int) c1);
    
    day  = int((t1-t0)/86400.0);
    hour = int(((t1-t0)-day*86400.0)/3600.0);
    min  = int(((t1-t0)-day*86400.0-hour*3600.0)/60.0);
    sec  = int((t1-t0)-day*86400.0-hour*3600.0-min*60.0);
    
    fprintf (logfile,"Elapsed wall clock time:        %ld sec.\n", (t1-t0));
    fprintf (logfile,"Elapsed wall clock time:        %i:%i:%i:%i (days:hours:minutes:seconds).\n", day, hour, min, sec);
    fprintf (logfile,"Elapsed CPU time:               %f sec.\n\n", (float) (c1 - c0)/CLOCKS_PER_SEC);
    fprintf(logfile,"***************************************\n");
    fprintf(logfile,"Ended the calculation of %s\n",string);
    fprintf(logfile,"***************************************\n\n");
    fflush(logfile);
  }
}
