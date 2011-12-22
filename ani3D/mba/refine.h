/* f2h postprocessed file */
#ifndef __REFINE_H__
#define __REFINE_H__

#include "f2c.h"

/* @f2h@ */ /* Subroutine */ int initializerefinement_(integer *np, integer *
	ne, doublereal *xyp, integer *ipe, doublereal *mapmtr, integer *
	ref2mapmtr)
;

/* @f2h@ */ /* Subroutine */ int uniformrefinement_(integer *np, integer *
	maxp, integer *nf, integer *maxf, integer *ne, integer *maxe, 
	doublereal *xyp, integer *ipf, integer *ipe, integer *lbf, integer *
	lbe, doublereal *mapmtr, integer *ref2mapmtr, integer *iw, integer *
	maxwi)
;

/* @f2h@ */ /* Subroutine */ int localrefinement_(integer *np, integer *maxp, 
	integer *nf, integer *maxf, integer *ne, integer *maxe, doublereal *
	xyp, integer *ipf, integer *ipe, integer *lbf, integer *lbe, 
	doublereal *mapmtr, integer *ref2mapmtr, logical *splitflag, integer *
	iw, integer *maxwi)
;

/* @f2h@ */ /* Subroutine */ int calreg_(doublereal *xy1, doublereal *xy2, 
	doublereal *xy3, doublereal *xy4, doublereal *qe)
;

/* @f2h@ */ /* Subroutine */ int fixedmap_(doublereal *xyp, doublereal *
	mapmatrix, doublereal *xypt)
;

#endif
