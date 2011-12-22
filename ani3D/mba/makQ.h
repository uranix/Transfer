/* f2h postprocessed file */
#ifndef __MAKQ_H__
#define __MAKQ_H__

#include "f2c.h"

/* @f2h@ */ /* Subroutine */ int makq_(integer *nloop, integer *np, integer *
	ne, doublereal *xyp, integer *ipe, integer *nev, integer *iev, 
	integer *nestar, doublereal *hstar, doublereal *hesp, doublereal *
	detg, doublereal *qe)
;

/* @f2h@ */ /* Subroutine */ int updqa_(integer *n, doublereal *xyp, integer *
	ipe, integer *iee, doublereal *qe)
;

/* @f2h@ */ /* Subroutine */ int updqb_(integer *nes, integer *le, integer *
	ies, doublereal *xyp, integer *ipes, doublereal *qes)
;

/* @f2h@ */ doublereal calvol_(doublereal *xy1, doublereal *xy2, doublereal *
	xy3, doublereal *xy4)
;

/* @f2h@ */ doublereal mutualorientation_(doublereal *xy1, doublereal *xy2, 
	doublereal *xy3, doublereal *xy4, doublereal *xy5)
;

/* @f2h@ */ doublereal calsqr_(doublereal *xy1, doublereal *xy2, doublereal *
	xy3)
;

/* @f2h@ */ /* Subroutine */ int caldet_(doublereal *hesp, doublereal *detg)
;

/* @f2h@ */ /* Subroutine */ int spectralmodule_(doublereal *hesp, doublereal 
	*detg)
;

/* @f2h@ */ /* Subroutine */ int calqe_(doublereal *hes1, doublereal *xy1, 
	doublereal *hes2, doublereal *xy2, doublereal *hes3, doublereal *xy3, 
	doublereal *hes4, doublereal *xy4, doublereal *hstar, doublereal *qe, 
	doublereal *volume)
;

/* @f2h@ */ /* Subroutine */ int calqf_(doublereal *hes1, doublereal *xy1, 
	doublereal *hes2, doublereal *xy2, doublereal *hes3, doublereal *xy3, 
	doublereal *hes4, doublereal *xy4, doublereal *hstar, integer *if__, 
	integer *ir)
;

/* @f2h@ */ /* Subroutine */ int hesbnd_(integer *lp, integer *ips, integer *
	icp, doublereal *hesp, doublereal *hesps)
;

/* @f2h@ */ /* Subroutine */ int iniq_analytic__(integer *np, doublereal *xyp,
	 I_fp metricfunction, doublereal *hesp)
;

/* @f2h@ */ doublereal avgq_(integer *ne, doublereal *qe, integer *l1e, 
	integer *l2e)
;

#endif
