/* f2h postprocessed file */
#ifndef __UPDATE_H__
#define __UPDATE_H__

#include "f2c.h"

/* @f2h@ */ /* Subroutine */ int pntadd_(integer *ip, integer *np, integer *
	maxp, integer *icp, doublereal *xyp, doublereal *hesp, doublereal *
	detg, integer *iholp, integer *icps, doublereal *xyps, doublereal *
	hesps, doublereal *detgs)
;

/* @f2h@ */ /* Subroutine */ int pntupd_(integer *ip, integer *icp, 
	doublereal *xyp, doublereal *hesp, doublereal *detg, integer *icps, 
	doublereal *xyps, doublereal *hesps, doublereal *detgs)
;

/* @f2h@ */ /* Subroutine */ int pntdel_(integer *ip, integer *np, integer *
	icp, integer *iholp)
;

/* @f2h@ */ /* Subroutine */ int facadd_(integer *if__, integer *nf, integer *
	maxf, integer *iholf)
;

/* @f2h@ */ /* Subroutine */ int facupd_(integer *nfs, integer *ipf, integer *
	ifs, integer *ipfs)
;

/* @f2h@ */ /* Subroutine */ int facdel_(integer *if__, integer *nf, integer *
	ipf, integer *iholf)
;

/* @f2h@ */ /* Subroutine */ int eleadd_(integer *ne, integer *maxe, integer *
	ihole)
;

/* @f2h@ */ /* Subroutine */ int eleupd_(integer *nes, integer *iep, integer *
	ipe, integer *ife, integer *iee, integer *lf, integer *le, integer *
	ifs, integer *ies, integer *ipfs, integer *ipes)
;

/* @f2h@ */ /* Subroutine */ int eledel_(integer *ie, integer *ipe, integer *
	iee)
;

#endif
