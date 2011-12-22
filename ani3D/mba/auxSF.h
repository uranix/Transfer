/* f2h postprocessed file */
#ifndef __AUXSF_H__
#define __AUXSF_H__

#include "f2c.h"

/* @f2h@ */ /* Subroutine */ int maksf_(integer *ip, integer *lf, integer *
	ifs, integer *le, integer *ies, integer *ipes, integer *ipf, integer *
	ife, integer *maxs, integer *lp1, integer *ip1s, integer *lp2, 
	integer *ip2s, integer *icp2s, integer *ls, integer *ipss, integer *
	iess, integer *ldf, integer *idfs, integer *lde, integer *ides, 
	logical *flagface, logical *flagedge)
;

/* @f2h@ */ /* Subroutine */ int shutf_(doublereal *zz1, doublereal *zz2, 
	doublereal *xyp1, doublereal *xyp2, doublereal *xyp3, logical *flag__)
;

/* @f2h@ */ doublereal projectf_(doublereal *xyp, doublereal *xy1, doublereal 
	*xy2, doublereal *xy3, doublereal *xyt)
;

/* @f2h@ */ /* Subroutine */ int baricoords_(doublereal *xyp, doublereal *xy1,
	 doublereal *xy2, doublereal *xy3, doublereal *b1, doublereal *b2, 
	doublereal *b3)
;

/* @f2h@ */ doublereal angle2faces_(doublereal *xya, doublereal *xyb, 
	doublereal *xyc, doublereal *xyd)
;

/* @f2h@ */ /* Subroutine */ int trinormal_(doublereal *xy1, doublereal *xy2, 
	doublereal *xy3, doublereal *normal)
;

/* @f2h@ */ /* Subroutine */ int vecmul_(doublereal *a, doublereal *b, 
	doublereal *c__)
;

/* @f2h@ */ doublereal dotmul_(doublereal *a, doublereal *b)
;

#endif
