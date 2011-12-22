#include <time.h>

#define SECS_PER_CLOCK (1./CLOCKS_PER_SEC)

//void seconds_( double *usecs, double *ssecs ) 
void seconds_( double *usecs ) 
{
/*   struct tms t;
   clock_t utime;
   clock_t stime;

   times(&t);
   utime = t.tms_utime;
   stime = t.tms_stime;

   *usecs = utime * SECS_PER_CLOCK;*/
   // *ssecs = stime * SECS_PER_CLOCK;
	*usecs = 0.;
}


