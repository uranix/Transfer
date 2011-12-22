/* splitE.f -- translated by f2c (version 20090411).
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

#include "utils.h"
#include "error.h"
#include "auxSE.h"
#include "auxSF.h"
#include "auxSP.h"
#include "auxSR.h"
#include "makM.h"
#include "makQ.h"
#include "update.h"
#include "list.h"

/* Table of constant values */

static integer c__1007 = 1007;

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int splite_(integer *iwe, integer *np, integer *
	maxp, integer *ne, integer *maxe, doublereal *xyp, integer *ipe, 
	doublereal *hstar, integer *icp, integer *iep, integer *ife, integer *
	iee, integer *iholp, integer *ihole, doublereal *hesp, doublereal *
	detg, doublereal *qe, integer *lf, integer *le, integer *ifs, integer 
	*ies, integer *ipfs, integer *ipes, doublereal *qes, logical *flag__)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, n;
    static doublereal v;
    static integer i1, i2, i3, ke, ip[5], ie1, ne1, ip1, ip2, ip3, ip4, ips, 
	    icps;
    static doublereal xyps[3];
    static integer leold;
    static doublereal detgs, hesps[6];

/* ================================================================ */
/* group (M) */
/* group (Q) */
/* group (S) */
/* group (W) */
/* ================================================================ */
/* ====================================================================== */
/* MaxS evaluates the number of elements in a superelement. */
/* Additionaly, it bounds the number of different boundary identificators. */
/* ====================================================================== */
/* ====================================================================== */
/* Colors: iVface - a fixed face */
/*         iMface - a material face */
/* ====================================================================== */
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
/* ================================================================ */
/* The operation does not change the quality list, since the last */
/* does not work properly when the quality is descreased. */

/* The metric in the new point is interpolated, since 8-tree is not */
/* yet available. The latter comes from an attempt to save memory. */
/* ================================================================ */
/* group (M) */
/* ========================================================== */
/* The local operations with grid: */
/*   F2E   : face is replaced by an edge. Two elements are */
/*           replaced by three elements. */

/*   E2F   : reverse to F2E operation. It works if 3 elements */
/*           have the common edge. */

/*   SWAP  : a generalization of E2F. The edge is destroyed */
/*           and a hole is filled with new tetrahedra. */

/*   DELET : delete a point. */

/*   MOVE  : movement of a point. */

/*   SPLTE : the elements is splitted into 4 elements by an */
/*           internal point. The operation is used only to */
/*           make grid to satisfy some FE restrictions. */
/* ========================================================== */
/* ...  available operations and their order (from 1 to 7 ONLY!) */
/* ...  number of operations */
/* group (Q) */
/* group (S) */
/* group (Flag) */
/* ================================================================ */
/* group (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --qes;
    ipes -= 6;
    ipfs -= 5;
    --ies;
    --ifs;
    --qe;
    --detg;
    hesp -= 7;
    --ihole;
    --iholp;
    iee -= 5;
    ife -= 5;
    --iep;
    --icp;
    ipe -= 6;
    xyp -= 4;

    /* Function Body */
    *flag__ = TRUE_;
    ip[0] = 1;
    ip[1] = 2;
    ip[2] = 3;
    ip[3] = 4;
    ip[4] = 1;
    ie1 = *iwe;
    findse_(le, &ies[1], &ie1, &ne1);
    ip1 = ipes[ne1 * 5 + 1];
    ip2 = ipes[ne1 * 5 + 2];
    ip3 = ipes[ne1 * 5 + 3];
    ip4 = ipes[ne1 * 5 + 4];
/* ... creating an inner point */
    icps = 64;
    for (i__ = 1; i__ <= 3; ++i__) {
	xyps[i__ - 1] = (xyp[i__ + ip1 * 3] + xyp[i__ + ip2 * 3] + xyp[i__ + 
		ip3 * 3] + xyp[i__ + ip4 * 3]) / 4;
    }
    for (i__ = 1; i__ <= 6; ++i__) {
	hesps[i__ - 1] = (hesp[i__ + ip1 * 6] + hesp[i__ + ip2 * 6] + hesp[
		i__ + ip3 * 6] + hesp[i__ + ip4 * 6]) / 4;
    }
/* !!!  LDH = 6 */
/* !!!  nXY = 1 */
/* !!!  Call LINTRP3D(nEw, IPEw, nPw, XYPw, LDH, HesPw, nXY, XYPs, */
/* !!! &              HesPs, iSE, miLINTRP, rSE, mrLINTRP, iControl) */
    caldet_(hesps, &detgs);
    pntadd_(&ips, np, maxp, &icp[1], &xyp[4], &hesp[7], &detg[1], &iholp[1], &
	    icps, xyps, hesps, &detgs);
/* ... creating 4 elements */
    leold = *le;
    for (i1 = 1; i1 <= 4; ++i1) {
	i2 = ip[i1];
	i3 = ip[i2];
	ip1 = ipe[i1 + *iwe * 5];
	ip2 = ipe[i2 + *iwe * 5];
	ip3 = ipe[i3 + *iwe * 5];
	if (i1 == 1) {
	    ke = ne1;
	} else {
	    ++(*le);
	    if (*le > 1000) {
		goto L9000;
	    }
	    ke = *le;
	}
	calqe_(&hesp[ip1 * 6 + 1], &xyp[ip1 * 3 + 1], &hesp[ip2 * 6 + 1], &
		xyp[ip2 * 3 + 1], &hesp[ip3 * 6 + 1], &xyp[ip3 * 3 + 1], 
		hesps, xyps, hstar, &qes[ke], &v);
	ipes[ke * 5 + 1] = ip1;
	ipes[ke * 5 + 2] = ip2;
	ipes[ke * 5 + 3] = ip3;
	ipes[ke * 5 + 4] = ips;
	ipes[ke * 5 + 5] = ipe[ie1 * 5 + 5];
    }
/* ... updating the grid */
/* !!!  Call lstUpd(nE, L1E, nL2, L2E, nStep, qE, iEs(nE1), qEs(nE1)) */
/* !!!  next line simulates lstAdd */
    qe[ies[ne1]] = qes[ne1];
    eledel_(&ies[ne1], &ipe[6], &iee[5]);
    i__1 = *le;
    for (n = leold + 1; n <= i__1; ++n) {
	eleadd_(ne, maxe, &ihole[1]);
/* !!!     Call lstAdd(nE, L1E, nL2, L2E, nStep, IHolE, */
/* !!! &               qE, qEs(n), iEs(n)) */
/* !!!     3 next lines simulate lstAdd */
	++(*ne);
	ies[n] = *ne;
	qe[ies[n]] = qes[n];
	eledel_(&ies[n], &ipe[6], &iee[5]);
    }
    eleupd_(&ne1, &iep[1], &ipe[6], &ife[5], &iee[5], lf, le, &ifs[1], &ies[1]
	    , &ipfs[5], &ipes[6]);
    i__1 = *le;
    for (n = leold + 1; n <= i__1; ++n) {
	eleupd_(&n, &iep[1], &ipe[6], &ife[5], &iee[5], lf, le, &ifs[1], &ies[
		1], &ipfs[5], &ipes[6]);
    }
    return 0;
L9000:
    errmes_(&c__1007, "splitE", "local parameter MaxS is small", (ftnlen)6, (
	    ftnlen)29);
    return 0;
} /* splite_ */

