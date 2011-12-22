#ifndef __LIBMBA_H__
#define __LIBMBA_H__

#ifdef __cplusplus
extern "C" {
#endif

int mbaFixShape(
	int *np, int maxp, int *nf, int maxf, int *ne, int maxe, 
	double *xyp, int *ipf, int *ipe, int *lbf, int *lbe, 
	int npv, int nfv, int nev, int *ipv, int *ifv, int *iev, 
	int flagauto, int status, int maxskipe, int maxqitr, 
	int	(*metricfunction)(double *, double *, double *, double*),
	double quality, double *rquality, 
	int maxwr, int maxwi, double *rw, int *iw, int iprint);

enum : int {
	ANI_METRIC_SCALAR = 0,
	ANI_METRIC_TENSOR = 1
};

#ifdef __cplusplus
}
#endif

#endif