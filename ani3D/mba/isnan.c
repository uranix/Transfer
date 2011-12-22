//to call floating-point check from the Fortran routine
#include <math.h>

int isnan(double x) { return x != x; }
int isinf(double x) { return !isnan(x) && isnan(x - x); }

void fpcheck_( double *a, int *flag ) 
{ 
   *flag = isnan(*a) && !isinf(*a); 
}

