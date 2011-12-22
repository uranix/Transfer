/* f2h postprocessed file */
#ifndef __MAKM_H__
#define __MAKM_H__

#include "f2c.h"

/* @f2h@ */ /* Subroutine */ int makm_(integer *np, integer *nf, integer *
	maxf, integer *ne, integer *maxe, doublereal *xyp, integer *ipf, 
	integer *ipe, integer *icp, integer *ipp, integer *iep, integer *ife, 
	integer *iee, integer *iholp, integer *iholf, integer *ihole, integer 
	*iepw, integer *nepw, integer *status, integer *npv, integer *nfv, 
	integer *nev, integer *ipv, integer *ifv, integer *iev, integer *ierr)
;

/* @f2h@ */ /* Subroutine */ int updm_(integer *np, integer *nf, integer *ne, 
	doublereal *xyp, integer *ipf, integer *ipe, integer *icp, integer *
	ipp, integer *ife, integer *iee, integer *iholp, integer *iholf, 
	integer *ihole, integer *status, integer *npv, integer *nfv, integer *
	nev, integer *ipv, integer *ifv, integer *iev, doublereal *hesp, 
	doublereal *qe, integer *ipw)
;

/* @f2h@ */ logical cmpe_(integer *i1, integer *i2, integer *i3, integer *iep,
	 integer *nep, integer *ie1, integer *ie2)
;

/* @f2h@ */ logical cmpf_(integer *i1, integer *i2, integer *ifp, integer *
	nfp, integer *if1, integer *if2)
;

/* @f2h@ */ logical cmpp_(integer *ip, integer *ifp, integer *nfp, integer *
	ipf)
;

/* @f2h@ */ logical cmpr_(integer *ip, integer *ifp, integer *nfp, doublereal 
	*xyp, integer *ipf)
;

/* @f2h@ */ logical crosspoint_(doublereal *xyp, integer *nf, integer *ifp, 
	integer *ip, integer *ipf)
;

/* @f2h@ */ integer countcolors_(integer *n, integer *ice)
;

/* @f2h@ */ /* Subroutine */ int scale2cube_(integer *np, doublereal *xyp, 
	logical *flag__)
;

/* @f2h@ */ /* Subroutine */ int scaleback_(doublereal *xypi, doublereal *
	xypo)
;

/* @f2h@ */ /* Subroutine */ int randr_(doublereal *xy1, doublereal *xy2, 
	doublereal *xy3, doublereal *xy4, doublereal *rout, doublereal *rin)
;

/* @f2h@ */ /* Subroutine */ int copymeshdata_(integer *np, integer *ne, 
	doublereal *xyp, doublereal *hesp, integer *ipe, doublereal *xypw, 
	doublereal *hespw, integer *ipew)
;

/* @f2h@ */ doublereal surfacearea_(integer *nf, doublereal *xyp, integer *
	ipf, integer *ic)
;

/* @f2h@ */ doublereal fixedarea_(integer *nfv, doublereal *xyp, integer *ipf,
	 integer *ifv)
;

/* @f2h@ */ doublereal domainvolume_(integer *ne, doublereal *xyp, integer *
	ipe, integer *ic)
;

/* @f2h@ */ doublereal fixedvolume_(integer *nev, doublereal *xyp, integer *
	ipe, integer *iev)
;

#endif
