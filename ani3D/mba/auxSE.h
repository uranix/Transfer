/* f2h postprocessed file */
#ifndef __AUXSE_H__
#define __AUXSE_H__

#include "f2c.h"

/* @f2h@ */ /* Subroutine */ int makse_(integer *ie, integer *icp, integer *
	iep, integer *ipf, integer *ipe, integer *ife, integer *iee, 
	doublereal *qe, integer *maxs, integer *lf, integer *le, integer *ifs,
	 integer *ies, integer *ipfs, integer *ipes, doublereal *qes)
;

/* @f2h@ */ /* Subroutine */ int calso_(doublereal *xyp, integer *ipe, 
	integer *le, integer *ies, integer *ios)
;

/* @f2h@ */ /* Subroutine */ int chkso_(integer *ip, doublereal *xyps, 
	doublereal *xyp, integer *ipe, integer *le, integer *ies, integer *
	ios, logical *flag__)
;

/* @f2h@ */ /* Subroutine */ int findse_(integer *le, integer *ies, integer *
	ie, integer *nes)
;

/* @f2h@ */ /* Subroutine */ int copyse_(integer *lfu, integer *leu, integer *
	ifu, integer *ieu, integer *ipfu, integer *ipeu, doublereal *qeu, 
	integer *lf, integer *le, integer *ifs, integer *ies, integer *ipfs, 
	integer *ipes, doublereal *qes)
;

/* @f2h@ */ /* Subroutine */ int copysemove_(integer *lfu, integer *leu, 
	integer *ifu, integer *ieu, doublereal *qeu, integer *lf, integer *le,
	 integer *ifs, integer *ies, doublereal *qes)
;

/* @f2h@ */ /* Subroutine */ int copysq_(integer *le, doublereal *qeu, 
	doublereal *xypu, doublereal *hespu, doublereal *detgu, doublereal *
	qes, doublereal *xyps, doublereal *hesps, doublereal *detgs)
;

/* @f2h@ */ logical check1j_(integer *i__, integer *j)
;

/* @f2h@ */ logical check2j_(integer *i1, integer *i2, integer *j)
;

/* @f2h@ */ logical check22_(integer *i1, integer *i2, integer *j1, integer *
	j2)
;

/* @f2h@ */ logical check3j_(integer *i1, integer *i2, integer *i3, integer *
	j)
;

/* @f2h@ */ logical check33_(integer *i1, integer *i2, integer *i3, integer *
	j1, integer *j2, integer *j3)
;

/* @f2h@ */ logical check13_(integer *i1, integer *j1, integer *j2, integer *
	j3)
;

/* @f2h@ */ logical check14_(integer *i1, integer *j1, integer *j2, integer *
	j3, integer *j4)
;

/* @f2h@ */ /* Subroutine */ int swapdd_(doublereal *d1, doublereal *d2)
;

/* @f2h@ */ /* Subroutine */ int swapii_(integer *i1, integer *i2)
;

#endif
