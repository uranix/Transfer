#ifndef H_USER3_MESH3D
#define H_USER3_MESH3D


#include"struct3.h"


/* exported  function  function */
int  bounSurf( int i, double u, double v, double *x, double *y, double *z );
void  V_U( int i, double u, double *v );
double userSizeFace( double x, double y, double z );
double  periodic(int dir);

int setsizefunction(double (*f) (double, double, double));
int setbounsurffunction(int (*f) (int, double, double, double*, double*, double*));
int setvufunction(void (*f) (int, double, double*));
int setcgmfunction(int (*f) (int, double, double, double*, double*, double*));
int setperiodicfunction(double (*f) (int));

#endif








