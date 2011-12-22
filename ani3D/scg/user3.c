#include"struct3.h"
#include"tree3.h"
#include"region3.h"
#include"user3.h"
#include"error3.h"

extern StrucRegion3  reg3;
double S0 = 0.1;

static double (*usersizefunction) (double, double, double) = NULL;
static int (*userbounsurffunction) (int, double, double, double*, double*, double*) = NULL;
static void (*uservufunction) (int, double, double*) = NULL;
static int (*cgmfunction) (int, double, double, double*, double*, double*) = NULL;
static double (*periodicfunction) (int) = NULL;

int bounSurf(int i, double u, double v, double *x, double *y, double *z) {
	if (i==-1) {
		if (cgmfunction) return cgmfunction(i, u, v, x, y, z);
		else return -1;
	}
	if (i==0) return bounSurf0(u, v, x, y, z);
	if (userbounsurffunction) return userbounsurffunction(i, u, v, x, y, z);
	errorExit3(3,"Use set_param_functions() to define parametrization functions");
   return  -1;
} /*bounSurf*/
double  periodic(int dir) {
	if (reg3.iSurf == 0)  return 0.0;
	if (periodicfunction)  return periodicfunction(dir);
	else  return 0.0;
}

void V_U(int i, double u, double *v) {
	if (i==0) {
		V_U0(u, v);
		return ;
	}
	if (uservufunction) {
		uservufunction(i, u, v);
		return ;
	}
	errorExit3(3,"Use set_param_functions() to define parametrization functions");
} /*V_U*/

double userSizeFace(double x, double y, double z) {
	if (usersizefunction) return usersizefunction(x, y, z);
	return S0;
} /*userSizeFace*/


/* Register user functions */
int setsizefunction(double (*f) (double, double, double)) {
	usersizefunction = f;
	return 0;
}
int setbounsurffunction(int (*f) (int, double, double, double*, double*, double*)) {
	userbounsurffunction = f;
	return 0;
}
int setvufunction(void (*f) (int, double, double*)) {
	uservufunction = f;
	return 0;
}
int setcgmfunction(int (*f) (int, double, double, double*, double*, double*)) {
	cgmfunction = f;
	return 0;
}
int setperiodicfunction(double (*f) (int)) {
	periodicfunction = f;
	return 0;
}
