/* auxSR.f -- translated by f2c (version 20090411).
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

#include "error.h"
#include "makQ.h"
#include "auxSP.h"
#include "auxSE.h"
#include "auxSF.h"
#include "auxSR.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;
static integer c__4 = 4;

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int maksr_(integer *ipa, integer *ipb, integer *
	le, integer *ies, integer *ipes, integer *lr, integer *irs, logical *
	flag__)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, m, n, m1, ip1, ip2, ip3, ip4, iwk, ipt, itr, icnt;

/* ================================================================ */
/* The ordered sequence of edges "orthogonal" to the edge {iPa, iPb} */
/* is computed. */

/* Remark: The superlement is not modified. */
/* ================================================================ */
/* group (S) */
/* group (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    irs -= 4;
    ipes -= 6;
    --ies;

    /* Function Body */
    *flag__ = TRUE_;
    *lr = 0;
    i__1 = *le;
    for (n = 1; n <= i__1; ++n) {
	ip1 = ipes[n * 5 + 1];
	ip2 = ipes[n * 5 + 2];
	ip3 = ipes[n * 5 + 3];
	ip4 = ipes[n * 5 + 4];
	if (check14_(ipa, &ip1, &ip2, &ip3, &ip4) && check14_(ipb, &ip1, &ip2,
		 &ip3, &ip4)) {
	    ++(*lr);
	    irs[*lr * 3 + 1] = n;
	    icnt = 1;
	    for (i__ = 1; i__ <= 4; ++i__) {
		ipt = ipes[i__ + n * 5];
		if (ipt != *ipa && ipt != *ipb) {
		    ++icnt;
		    irs[icnt + *lr * 3] = ipt;
		}
	    }
	}
    }
/* ... reodering of the edges (2 iterations are allowed) */
    itr = 0;
L10:
    i__1 = *lr - 1;
    for (n = 1; n <= i__1; ++n) {
	ip1 = irs[n * 3 + 2];
	ip2 = irs[n * 3 + 3];
	i__2 = *lr;
	for (m = n + 1; m <= i__2; ++m) {
	    ip3 = irs[m * 3 + 2];
	    ip4 = irs[m * 3 + 3];
	    m1 = m;
	    if (ip2 == ip3) {
		goto L20;
	    } else if (ip2 == ip4) {
		irs[m * 3 + 2] = ip4;
		irs[m * 3 + 3] = ip3;
		goto L20;
	    }
	}
	for (i__ = 1; i__ <= 3; ++i__) {
	    iwk = irs[i__ + 3];
	    irs[i__ + 3] = irs[i__ + n * 3];
	    irs[i__ + n * 3] = iwk;
	}
	irs[5] = ip2;
	irs[6] = ip1;
	++itr;
	if (itr == 3) {
/*  ...  there is no connected path */
	    *flag__ = FALSE_;
	    goto L1000;
	}
	goto L10;
L20:
	for (i__ = 1; i__ <= 3; ++i__) {
	    iwk = irs[i__ + (n + 1) * 3];
	    irs[i__ + (n + 1) * 3] = irs[i__ + m1 * 3];
	    irs[i__ + m1 * 3] = iwk;
	}
    }
L1000:
    return 0;
} /* maksr_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int clrsr_(integer *ipa, integer *ipb, integer *
	icp, integer *ipf, integer *ife, integer *lf, integer *ifs, integer *
	le, integer *ies, integer *icrab)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k, n, i1, i2, ip[4], ic1, ip1, ip2, ifp, ift, iet;
    static logical flag__;
    static integer icnt[1000], iclrf;

/* ================================================================ */
/* The color of edge {iPa, iPb} is computed. */

/* Remark: The superlement is not modified. */
/* ================================================================ */
/* ====================================================================== */
/* MaxS evaluates the number of elements in a superelement. */
/* Additionaly, it bounds the number of different boundary identificators. */
/* ====================================================================== */
/* ====================================================================== */
/* Colors: iVface - a fixed face */
/*         iMface - a material face */
/* ====================================================================== */
/* ================================================================ */
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
/* group (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --ies;
    --ifs;
    ife -= 5;
    ipf -= 5;
    --icp;

    /* Function Body */
    ip[0] = 1;
    ip[1] = 2;
    ip[2] = 3;
    ip[3] = 1;
    flag__ = FALSE_;
    i__1 = *lf;
    for (n = 1; n <= i__1; ++n) {
	icnt[n - 1] = -1;
    }
    *icrab = 0;
    i__1 = *lf;
    for (n = 1; n <= i__1; ++n) {
	ift = ifs[n];
	if (ift <= 0) {
	    goto L10;
	}
	for (i1 = 1; i1 <= 3; ++i1) {
	    i2 = ip[i1];
	    ip1 = ipf[i1 + (ift << 2)];
	    ip2 = ipf[i2 + (ift << 2)];
	    if (*ipa == ip1 && *ipb == ip2 || *ipb == ip1 && *ipa == ip2) {
		iclrf = ipf[(ift << 2) + 4];
		icnt[n - 1] = 0;
		if (*icrab == 0) {
		    ifp = ift;
		    *icrab = iclrf;
		} else if (*icrab != iclrf) {
		    *icrab = minclr_(&icp[*ipa], &icp[*ipb]);
/*  ...  a temporary fix for V-V edges */
		    delxnode_(icrab, &c__1);
		    addxnode_(icrab, &c__2);
		    goto L1000;
		}
	    }
	}
L10:
	;
    }
    if (*icrab != 0) {
	*icrab = 8;
/*  ...  checking for edges on boundary surfaces */
	ic1 = 0;
	i__1 = *le;
	for (k = 1; k <= i__1; ++k) {
	    iet = ies[k];
	    if (iet > 0) {
		for (i__ = 1; i__ <= 4; ++i__) {
		    if (ife[i__ + (iet << 2)] == ifp) {
			++ic1;
		    }
		    if (ic1 == 2) {
			goto L1000;
		    }
		}
	    }
	}
	addxnode_(icrab, &c__4);
    } else {
	*icrab = 64;
    }
L1000:
    return 0;
} /* clrsr_ */

/* ================================================================ */
/* @f2h@ */ doublereal angle2edges_(doublereal *xya, doublereal *xyb, 
	doublereal *xyc)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3;

    /* Local variables */
    static integer i__;
    static doublereal d1, d2, xym[3], xyn[3];

/* ================================================================ */
/* Cosine of angle between two edges {A,B} and {A,C}. */
/* ================================================================ */
    /* Parameter adjustments */
    --xyc;
    --xyb;
    --xya;

    /* Function Body */
    for (i__ = 1; i__ <= 3; ++i__) {
	xyn[i__ - 1] = xyb[i__] - xya[i__];
	xym[i__ - 1] = xyc[i__] - xya[i__];
    }
/* Computing 2nd power */
    d__1 = xyn[0];
/* Computing 2nd power */
    d__2 = xyn[1];
/* Computing 2nd power */
    d__3 = xyn[2];
    d1 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
/* Computing 2nd power */
    d__1 = xym[0];
/* Computing 2nd power */
    d__2 = xym[1];
/* Computing 2nd power */
    d__3 = xym[2];
    d2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
    ret_val = (xyn[0] * xym[0] + xyn[1] * xym[1] + xyn[2] * xym[2]) / sqrt(d1 
	    * d2);
    return ret_val;
} /* angle2edges_ */

/* ================================================================ */
/* @f2h@ */ doublereal caledge_(doublereal *xy1, doublereal *xy2)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3;

/* ================================================================ */
    /* Parameter adjustments */
    --xy2;
    --xy1;

    /* Function Body */
/* Computing 2nd power */
    d__1 = xy1[1] - xy2[1];
/* Computing 2nd power */
    d__2 = xy1[2] - xy2[2];
/* Computing 2nd power */
    d__3 = xy1[3] - xy2[3];
    ret_val = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
    return ret_val;
} /* caledge_ */

/* ================================================================ */
/* @f2h@ */ doublereal sqredge_(doublereal *xy1, doublereal *xy2)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3;

/* ================================================================ */
    /* Parameter adjustments */
    --xy2;
    --xy1;

    /* Function Body */
/* Computing 2nd power */
    d__1 = xy1[1] - xy2[1];
/* Computing 2nd power */
    d__2 = xy1[2] - xy2[2];
/* Computing 2nd power */
    d__3 = xy1[3] - xy2[3];
    ret_val = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
    return ret_val;
} /* sqredge_ */

