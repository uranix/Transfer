#include "mba_fixshape.h"

int mbaFixShape(
	int *np, int maxp, int *nf, int maxf, int *ne, int maxe, 
	double *xyp, int *ipf, int *ipe, int *lbf, int *lbe, 
	int npv, int nfv, int nev, int *ipv, int *ifv, int *iev, 
	int flagauto, int status, int maxskipe, int maxqitr, 
	int	(*metricfunction)(double *, double *, double *, double *),
	double quality, double *rquality, 
	int maxwr, int maxwi, double *rw, int *iw, int iprint) 
{
	int ierr;
	mbafixshape_(np, &maxp, nf, &maxf, ne, &maxe, 
				xyp, ipf, ipe, lbf, lbe, 
				&npv, &nfv, &nev, ipv, ifv, iev, 
				&flagauto, &status, &maxskipe, &maxqitr, 
				(I_fp)metricfunction, 
				&quality, rquality,
				&maxwr, &maxwi, rw, iw, 
				&iprint, &ierr);
	return ierr;
}