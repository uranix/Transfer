/* f2h postprocessed file */
#ifndef __LIST_H__
#define __LIST_H__

#include "f2c.h"

/* @f2h@ */ /* Subroutine */ int lstmak_(integer *nestar, integer *ne, 
	integer *l1e, integer *l2e, integer *nl2, integer *icntl, integer *
	lhol)
;

/* @f2h@ */ /* Subroutine */ int lstupd_(integer *ne, integer *l1e, integer *
	nl2, integer *l2e, integer *icntl, doublereal *qe, integer *ie, 
	doublereal *qie)
;

/* @f2h@ */ /* Subroutine */ int lstdel_(integer *ne, integer *l1e, integer *
	nl2, integer *l2e, integer *icntl, integer *lhol, doublereal *qe, 
	integer *ie)
;

/* @f2h@ */ /* Subroutine */ int lstadd_(integer *ne, integer *l1e, integer *
	nl2, integer *l2e, integer *icntl, integer *lhol, doublereal *qe, 
	doublereal *qie, integer *ie)
;

/* @f2h@ */ integer prevl2_(integer *nl2, integer *l2e, doublereal *qe, 
	doublereal *qie)
;

/* @f2h@ */ integer prevl2ie_(integer *nl2, integer *l2e, doublereal *qe, 
	integer *ie, integer *l1e)
;

/* @f2h@ */ integer prevl1_(integer *l1e, integer *l2e, integer *iprevl2, 
	doublereal *qe, doublereal *qie)
;

#endif
