/* f2h postprocessed file */
#ifndef __CLPSR_H__
#define __CLPSR_H__

#include "f2c.h"

/* @f2h@ */ /* Subroutine */ int clpsr_(integer *iwr, integer *iwe, integer *
	np, integer *ne, doublereal *xyp, integer *ipe, doublereal *hstar, 
	integer *icp, integer *iep, integer *ife, integer *iee, integer *l1e, 
	integer *l2e, integer *nl2, integer *nstep, integer *iholp, integer *
	ihole, integer *status, doublereal *hesp, doublereal *rquality, 
	doublereal *detg, doublereal *qe, I_fp metricfunction, logical *
	flaganalytic, integer *lfu, integer *leu, integer *ifu, integer *ieu, 
	integer *ipfu, integer *ipeu, doublereal *qeu, integer *npw, integer *
	new__, doublereal *xypw, doublereal *hespw, integer *ipew, integer *
	milintrp, integer *mrlintrp, integer *ise, doublereal *rse, integer *
	icontrol, logical *flag__)
;

/* @f2h@ */ /* Subroutine */ int chksotet_(integer *ip1, integer *ip2, 
	integer *ipa, integer *ipb, doublereal *xyps, doublereal *xyp, 
	integer *ipe, doublereal *v1, logical *flag__)
;

#endif
