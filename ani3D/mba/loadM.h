/* f2h postprocessed file */
#ifndef __LOADM_H__
#define __LOADM_H__

#include "f2c.h"

/* @f2h@ */ /* Subroutine */ int loadmani_(integer *maxp, integer *maxf, 
	integer *maxe, integer *np, integer *nf, integer *ne, doublereal *xyp,
	 integer *ipf, integer *ipe, integer *lbf, integer *lbe, integer *npv,
	 integer *nfv, integer *nev, integer *ipv, integer *ifv, integer *iev,
	 integer *ipp, integer *iff, char *fname, ftnlen fname_len)
;

/* @f2h@ */ /* Subroutine */ int loadmani_header__(integer *np, integer *nf, 
	integer *ne, integer *npv, integer *nfv, integer *nev, char *fname, 
	ftnlen fname_len)
;

/* @f2h@ */ /* Subroutine */ int loads_(integer *np, doublereal *sol, char *
	fname, ftnlen fname_len)
;

/* @f2h@ */ /* Subroutine */ int loadmaft_(integer *maxp, integer *maxf, 
	integer *maxe, integer *np, integer *nf, integer *ne, doublereal *xyp,
	 integer *ipf, integer *ipe, integer *lbf, integer *lbe, integer *npv,
	 integer *nfv, integer *nev, integer *ipv, integer *ifv, integer *iev,
	 integer *ipp, integer *iff, char *fname, ftnlen fname_len)
;

/* @f2h@ */ /* Subroutine */ int loadmgmv_(integer *maxp, integer *maxf, 
	integer *maxe, integer *np, integer *nf, integer *ne, doublereal *xyp,
	 integer *ipf, integer *ipe, integer *lbf, integer *lbe, integer *npv,
	 integer *nfv, integer *nev, integer *ipv, integer *ifv, integer *iev,
	 integer *ipp, integer *iff, char *fname, ftnlen fname_len)
;

#endif
