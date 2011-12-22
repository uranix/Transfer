/* nlnfnc.f -- translated by f2c (version 20090411).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

#include "ani2.h"
#include "utils.h"
#include "error.h"
#include "auxSE.h"
#include "auxSF.h"
#include "auxSP.h"
#include "auxSR.h"
#include "update.h"

#include "lintrp3d.h"
#include "makQ.h"
#include "makM.h"

/* Table of constant values */

static integer c__64 = 64;
static integer c__1 = 1;

/* ================================================================ */
/* @f2h@ */ doublereal nlnfnc_(integer *nu, doublereal *u, doublereal *xyp, 
	doublereal *hesp, doublereal *detg, doublereal *hstar, integer *ips, 
	integer *le, integer *ies, doublereal *xyps, integer *icps, integer *
	ipes, doublereal *hesps, doublereal *detgs, doublereal *qes, integer *
	npw, integer *new__, doublereal *xypw, doublereal *hespw, integer *
	ipew, I_fp metricfunction, logical *flaganalytic, integer *milintrp, 
	integer *mrlintrp, integer *ise, doublereal *rse, integer *icontrol)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2;

    /* Local variables */
    static integer i__, n;
    static doublereal v;
    static integer i1, i2, i3, i4, ie, ip[5], ldh, ipa, ipb, ipc, nxy;
    static doublereal xypo[3];

/* ================================================================ */
/* group (F) */
/* group (ANI) */
/* ================================================================ */
/* ================================================================ */
/* group (F) */
/* ===== POINTS === POINTS === POINTS === POINTS ====================== */
/*  Old definition of the colors */

/*  Colors:   T - tough point: can not be changed anyway, */
/*                             boundary edges with T end points */
/*                                can not be changed */
/*                Remark: the T-color appears only in a */
/*                        mixture with another color. */

/*            V - fixed point:       can not be changed anyway */
/*           VB - edge point:        has to live on the edge */
/*            B - surface point:     has to live on the surface */
/*          VBI - internal VB-point: similar to VB */
/*           BI - internal B-point:  similar to  B */
/*            I - internal point:    free point inside the domain */

/*  Working colors which will be replaced by one from the above list: */
/*           TV - tough  V-point */
/*          TVB - tough VB-point */
/*           TB - tough  B-point */

/* ========================================================== */
/*     Integer   iTVnode, iTVBnode, iTBnode */
/*     Parameter(iTVnode = 1, iTVBnode = 2, iTBnode = 3) */

/*     Integer   iVnode, iVBnode, iBnode, iVBInode, iBInode, iInode */
/*     Parameter(iVnode   = 4, iVBnode = 5, iBnode = 6) */
/*     Parameter(iVBInode = 7, iBInode = 8, iInode = 9) */


/* ==================================================================== */
/*  New binary definition of colors. Each bit stands for specific */
/*  location of a point: boundary, surface, vertex, temporary */
/*  frosen point, etc. Before mesh regeneration, we try to recover */
/*  as many point characteristics as possible. For example, a fix */
/*  point on the domain boundary is vertex, boundary point and surface */
/*  point simulataneously. It could be also the temporary frosen point. */

/*  The difefrence between V abd T nodes is as follows. We can put */
/*  a new point on V-V or V-T edge, but we can not split T-T edge. */

/*  The characteristics are recovered with binary operations IAND */
/*  and IOR. To distinguish between old and new colors, we shall use */
/*  prefix "j" instead of "i". */

/*   bit 01  -  vertex                (V) */
/*   bit 02  -  edge node             (R) */
/*   bit 03  -  boundary node         (B) only outer boundary */
/*   bit 04  -  surface node          (S) */
/*   bit 05  -  reserved */
/*   bit 06  -  reserved */
/*   bit 07  -  inner node            (I) */
/*   bit 08  -  temporary frosen node (T) */

/* REMARK 1. All colors should be powers of 2. */

/* REMARK 2. The contradictory colors are jBnode and jInode. */

/* ==================================================================== */
/* ========================================================== */
/* Colors for faces: */
/*     F - fictitious faces: used for dummy face in IFE */
/*                           supposed to be negative */
/*     V - fixed faces:      face can not be modified */
/*     M - material faces:   interface between different materials */
/*                           interval of colors [iMface; MaxS] */
/*                           is reserved for 100 different interfaces */
/* ========================================================== */
/*     Integer   iVface, iMface */
/*     Parameter(iVface = MaxS - 100 , iMface = MaxS - 99) */
/* group (ANI) */
/* group (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --rse;
    --ise;
    ipew -= 5;
    hespw -= 7;
    xypw -= 4;
    --qes;
    --hesps;
    ipes -= 6;
    --xyps;
    --ies;
    --detg;
    hesp -= 7;
    xyp -= 4;
    --u;

    /* Function Body */
    ip[0] = 1;
    ip[1] = 2;
    ip[2] = 3;
    ip[3] = 4;
    ip[4] = 1;
    for (i__ = 1; i__ <= 3; ++i__) {
	xyps[i__] = u[i__];
    }
    if (ifxnode_(icps, &c__64)) {
	for (i__ = 1; i__ <= 6; ++i__) {
	    hesps[i__] = hesp[i__ + *ips * 6];
	}
	ldh = 6;
	nxy = 1;
	if (! (*flaganalytic)) {
	    lintrp3d_(new__, &ipew[5], npw, &xypw[4], &ldh, &hespw[7], &nxy, &
		    xyps[1], &hesps[1], &ise[1], milintrp, &rse[1], mrlintrp, 
		    icontrol);
	} else {
	    scaleback_(&xyps[1], xypo);
	    iniq_analytic__(&c__1, xypo, (I_fp)metricfunction, &hesps[1]);
	}
	caldet_(&hesps[1], detgs);
    } else {
	for (i__ = 1; i__ <= 6; ++i__) {
	    hesps[i__] = hesp[i__ + *ips * 6];
	}
	*detgs = detg[*ips];
    }
    ret_val = 1.;
    i__1 = *le;
    for (n = 1; n <= i__1; ++n) {
	ie = ies[n];
	if (ie <= 0) {
	    goto L10;
	}
	for (i1 = 1; i1 <= 4; ++i1) {
	    if (ipes[i1 + n * 5] == *ips) {
		i2 = ip[i1];
		i3 = ip[i2];
		i4 = ip[i3];
		ipa = ipes[i2 + n * 5];
		ipb = ipes[i3 + n * 5];
		ipc = ipes[i4 + n * 5];
		calqe_(&hesps[1], &xyps[1], &hesp[ipa * 6 + 1], &xyp[ipa * 3 
			+ 1], &hesp[ipb * 6 + 1], &xyp[ipb * 3 + 1], &hesp[
			ipc * 6 + 1], &xyp[ipc * 3 + 1], hstar, &qes[n], &v);
		goto L5;
	    }
	}
L5:
/* Computing MIN */
	d__1 = ret_val, d__2 = qes[n];
	ret_val = min(d__1,d__2);
L10:
	;
    }
    ret_val = 1. - ret_val;
    return ret_val;
} /* nlnfnc_ */

