/* f2h postprocessed file */
#ifndef __UTILS_H__
#define __UTILS_H__

#include "f2c.h"

/* @f2h@ */ /* Subroutine */ int listp2p_(integer *np, integer *ne, integer *
	maxlist, integer *ipe, integer *npp, integer *ipp, integer *iw)
;

/* @f2h@ */ /* Subroutine */ int listr2p_(integer *np, integer *nr, integer *
	ne, integer *maxr, integer *ipe, integer *ipr, integer *iw)
;

/* @f2h@ */ /* Subroutine */ int listr2r_(integer *np, integer *nr, integer *
	ne, integer *maxl, integer *ipe, integer *nrr, integer *irr, integer *
	iw, integer *ierr)
;

/* @f2h@ */ /* Subroutine */ int liste2r_(integer *np, integer *nr, integer *
	ne, integer *ipe, integer *ire, integer *nep, integer *iep)
;

/* @f2h@ */ /* Subroutine */ int liste2f_(integer *np, integer *nf, integer *
	ne, integer *ipe, integer *ife, integer *nep, integer *iep)
;

/* @f2h@ */ /* Subroutine */ int liste2fb_(integer *np, integer *nfb, integer 
	*ne, integer *ipf, integer *ipe, integer *ife, integer *nep, integer *
	iep)
;

/* @f2h@ */ /* Subroutine */ int listconv_(integer *np, integer *nr, integer *
	ne, integer *nep, integer *iep, integer *l, integer *ire, integer *nx,
	 integer *maxx, integer *nrp, integer *irp, integer *iw, integer *
	ierr)
;

/* @f2h@ */ /* Subroutine */ int reversemap_(integer *np, integer *nr, 
	integer *nrp, integer *irp, integer *npr, integer *ipr)
;

/* @f2h@ */ /* Subroutine */ int backreferences_(integer *np, integer *ne, 
	integer *l, integer *m, integer *ipe, integer *nep, integer *iep)
;

/* @f2h@ */ /* Subroutine */ int dualnormals_(integer *np, integer *nr, 
	integer *ne, integer *ipe, integer *ipr, doublereal *xyp, doublereal *
	nrm, integer *iw)
;

/* @f2h@ */ /* Subroutine */ int addboundaryfaces_(integer *np, integer *nf, 
	integer *maxf, integer *ne, doublereal *xyp, integer *ipf, integer *
	ipe, integer *lbf, integer *lbe, integer *iw)
;

/* @f2h@ */ /* Subroutine */ int addmaterialfaces_(integer *np, integer *nf, 
	integer *maxf, integer *ne, doublereal *xyp, integer *ipf, integer *
	ipe, integer *lbf, integer *lbe, integer *iw)
;

/* @f2h@ */ /* Subroutine */ int global2local_(integer *myid, integer *ice, 
	integer *np, integer *nf, integer *maxf, integer *ne, doublereal *xyp,
	 integer *ipf, integer *ipe, integer *lbf, integer *lbe, integer *npv,
	 integer *nfv, integer *nev, integer *ipv, integer *ifv, integer *iev,
	 integer *npl, integer *nfl, integer *nel, doublereal *xypl, integer *
	ipfl, integer *ipel, integer *lbfl, integer *lbel, integer *npvl, 
	integer *nfvl, integer *nevl, integer *ipvl, integer *ifvl, integer *
	ievl, integer *nfvi, integer *ippl, integer *iffl, integer *ipw, 
	integer *ifw, integer *iew)
;

/* @f2h@ */ /* Subroutine */ int local2global_(integer *nmeshes, integer *
	maxp, integer *maxf, integer *maxe, integer *np, integer *nf, integer 
	*ne, doublereal *xyp, integer *ipf, integer *ipe, integer *lbf, 
	integer *lbe, integer *npv, integer *nfv, integer *nev, integer *ipv, 
	integer *ifv, integer *iev, integer *npl, integer *nfl, integer *nel, 
	doublereal *xypl, integer *ipfl, integer *ipel, integer *lbfl, 
	integer *lbel, integer *npvl, integer *nfvl, integer *nevl, integer *
	ipvl, integer *ifvl, integer *ievl, integer *ippl, integer *iffl, 
	integer *ips, integer *ipw, integer *ifw, integer *iew)
;

/* @f2h@ */ /* Subroutine */ int maktnode_(integer *np, integer *maxp, 
	integer *ne, integer *ipp, integer *ipe, integer *ice)
;

#endif
