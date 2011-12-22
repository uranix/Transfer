/* f2h postprocessed file */
#ifndef __MBA_ANALYTIC_H__
#define __MBA_ANALYTIC_H__

#include "f2c.h"

/* @f2h@ */ /* Subroutine */ int mbaanalytic_(integer *np, integer *maxp, 
	integer *nf, integer *maxf, integer *ne, integer *maxe, doublereal *
	xyp, integer *ipf, integer *ipe, integer *lbf, integer *lbe, integer *
	nestar, integer *npv, integer *nfv, integer *nev, integer *ipv, 
	integer *ifv, integer *iev, logical *flagauto, integer *status, 
	integer *maxskipe, integer *maxqitr, I_fp metricfunction, doublereal *
	quality, doublereal *rquality, integer *maxwr, integer *maxwi, 
	doublereal *rw, integer *iw, integer *iprint, integer *ierr)
;

/* @f2h@ */ /* Subroutine */ int mbaanalyticshort_(integer *np, integer *maxp,
	 integer *nf, integer *maxf, integer *ne, integer *maxe, doublereal *
	xyp, integer *ipf, integer *ipe, integer *lbf, integer *nestar, 
	integer *status, integer *maxwr, integer *maxwi, doublereal *rw, 
	integer *iw, integer *iprint, integer *ierr)
;

#endif
