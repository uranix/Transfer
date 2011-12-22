/* f2h postprocessed file */
#ifndef __AUXSR_H__
#define __AUXSR_H__

#include "f2c.h"

/* @f2h@ */ /* Subroutine */ int maksr_(integer *ipa, integer *ipb, integer *
	le, integer *ies, integer *ipes, integer *lr, integer *irs, logical *
	flag__)
;

/* @f2h@ */ /* Subroutine */ int clrsr_(integer *ipa, integer *ipb, integer *
	icp, integer *ipf, integer *ife, integer *lf, integer *ifs, integer *
	le, integer *ies, integer *icrab)
;

/* @f2h@ */ doublereal angle2edges_(doublereal *xya, doublereal *xyb, 
	doublereal *xyc)
;

/* @f2h@ */ doublereal caledge_(doublereal *xy1, doublereal *xy2)
;

/* @f2h@ */ doublereal sqredge_(doublereal *xy1, doublereal *xy2)
;

#endif
