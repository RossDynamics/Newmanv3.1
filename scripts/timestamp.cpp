/* Utility to insert date strings into a Tecplot FTLE output file.
The output from newman2D consists of zones that are labelled with a floating point time value.  
It is often convenient when plotting to replace the time label with a date and time stamp.  This utility allows you to do this. 
When the program runs, you will need to provide the date and time of your time origin.*/

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>

#define PST (-8)


int main() {

char  *datestr;
time_t timezerosec;
time_t timenowsec;
double  Time_Scale;
struct tm timezero;
struct tm timenow;
int   year0, month0 ,day0, hour0, minute0, second0;
int   year1, month1 ,day1, hour1, minute1, second1;
double  zonehours;
double timediff;

 
  /* prompt user for origin date */
  printf ("\nEnter date and time  of time origin (YYYY-MM-DD HH-MM-SS): "); 
  scanf ("%d-%d-%d %d-%d-%d",&year0, &month0, &day0, &hour0, &minute0, &second0 );

  timezero.tm_year = year0-1900;
  timezero.tm_mon = month0-1;
  timezero.tm_mday = day0;
  timezero.tm_hour = hour0;
  timezero.tm_min = minute0;
  timezero.tm_sec = second0;
  timezero.tm_isdst = 0;


  timezerosec=mktime(&timezero);
  datestr=asctime(&timezero);
  printf("Time origin: %s",datestr);
  
   /* prompt user for current date */
  printf ("\nEnter current date and time (YYYY-MM-DD HH-MM-SS): "); 
  scanf ("%d-%d-%d %d-%d-%d",&year1, &month1, &day1, &hour1, &minute1, &second1 );

  timenow.tm_year = year1-1900;
  timenow.tm_mon = month1-1;
  timenow.tm_mday = day1;
  timenow.tm_hour = hour1;
  timenow.tm_min = minute1;
  timenow.tm_sec = second1;
  timenow.tm_isdst = 0;


  timenowsec=mktime(&timenow);
  datestr=asctime(&timenow);
  printf("Time now: %s",datestr);

  printf("\nAll time data will be converted to seconds. \n Enter conversion factor (e.g. Hours = 3600): ");
  scanf("%lf",&Time_Scale);

  timediff=difftime(timenowsec,timezerosec);
  zonehours=timediff/Time_Scale;
  
  printf("Time since origin: %12f\n",zonehours);

printf("DONE\a\n");

return(0);
}