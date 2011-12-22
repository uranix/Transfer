/* f2h postprocessed file */
#ifndef __MOVEP_H__
#define __MOVEP_H__

#include "f2c.h"

/* @f2h@ */ /* Subroutine */ int movep_(integer *iwp, integer *iwe, integer *
	ne, doublereal *xyp, integer *ipf, integer *ipe, doublereal *hstar, 
	integer *icp, integer *ife, integer *l1e, integer *l2e, integer *nl2, 
	integer *nstep, integer *status, doublereal *hesp, doublereal *
	rquality, doublereal *detg, doublereal *qe, I_fp metricfunction, 
	logical *flaganalytic, integer *lfu, integer *leu, integer *ifu, 
	integer *ieu, integer *ipfs, integer *ipes, doublereal *qeu, integer *
	npw, integer *new__, doublereal *xypw, doublereal *hespw, integer *
	ipew, integer *milintrp, integer *mrlintrp, integer *ise, doublereal *
	rse, integer *icontrol, doublereal *rmove, logical *flag__)
;

/* @f2h@ */ doublereal distsr_(integer *ipo, doublereal *nx, doublereal *ny, 
	doublereal *nz, integer *lf, integer *ifs, integer *ipfs, doublereal *
	xyp)
;

/* @f2h@ */ doublereal distsf_(integer *ipo, doublereal *nx, doublereal *ny, 
	doublereal *nz, integer *le, integer *ies, integer *ipes, doublereal *
	xyp)
;

#endif
