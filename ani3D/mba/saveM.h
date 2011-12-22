/* f2h postprocessed file */
#ifndef __SAVEM_H__
#define __SAVEM_H__

#include "f2c.h"

/* @f2h@ */ /* Subroutine */ int savemani_(integer *np, integer *nf, integer *
	ne, doublereal *xyp, integer *ipf, integer *ipe, integer *lbf, 
	integer *lbe, integer *npv, integer *nfv, integer *nev, integer *ipv, 
	integer *ifv, integer *iev, logical *flagi, integer *ipp, integer *
	iff, char *fname, ftnlen fname_len)
;

/* @f2h@ */ /* Subroutine */ int saves_(integer *np, doublereal *sol, char *
	fname, ftnlen fname_len)
;

/* @f2h@ */ /* Subroutine */ int savemgmv_(integer *np, integer *nf, integer *
	ne, doublereal *xyp, integer *ipf, integer *ipe, integer *lbf, 
	integer *lbe, char *fname, integer *iw, ftnlen fname_len)
;

/* @f2h@ */ /* Subroutine */ int savemgeo_(integer *np, integer *maxp, 
	integer *nf, integer *maxf, integer *ne, integer *maxe, doublereal *
	xyp, integer *ipf, integer *ipe, integer *lbp, integer *lbf, integer *
	lbe, integer *npv, integer *nfv, integer *nev, integer *ipv, integer *
	ifv, integer *iev, char *fname, ftnlen fname_len)
;

/* @f2h@ */ /* Subroutine */ int savemaft_(integer *np, integer *nf, integer *
	ne, doublereal *xyp, integer *ipf, integer *ipe, integer *lbf, 
	integer *lbe, char *fname, ftnlen fname_len)
;

#endif
