#include <sys/time.h>
//#include <time.h>
                            
int gettime(double *sec)
{
  /*  extern int gettimeofday(); */
  struct timeval t;
  struct timezone tz;
                              
  int stat;
  stat = gettimeofday(&t, &tz);
  *sec = (double)(t.tv_sec + t.tv_usec/1000000.0);
  return(stat);
}
