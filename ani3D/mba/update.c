/* update.f -- translated by f2c (version 20090411).
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
#include "makM.h"
#include "update.h"

/* Table of constant values */

static integer c__1003 = 1003;
static integer c__1004 = 1004;
static integer c__1006 = 1006;
static integer c__1000 = 1000;

/* ========================================================== */
/* PONITS ... POINT ... POINTS ... POINT */
/* ========================================================== */
/* @f2h@ */ /* Subroutine */ int pntadd_(integer *ip, integer *np, integer *
	maxp, integer *icp, doublereal *xyp, doublereal *hesp, doublereal *
	detg, integer *iholp, integer *icps, doublereal *xyps, doublereal *
	hesps, doublereal *detgs)
{
    static integer nholp;

/* ========================================================== */
    /* Parameter adjustments */
    hesps -= 7;
    xyps -= 4;
    --iholp;
    --detg;
    hesp -= 7;
    xyp -= 4;
    --icp;

    /* Function Body */
    nholp = iholp[1];
    ++(*np);
    if (*np > *maxp) {
	errmes_(&c__1003, "pntAdd", "local parameter MaxP is small", (ftnlen)
		6, (ftnlen)29);
    }
    if (nholp == 0) {
	*ip = *np;
    } else {
	*ip = iholp[nholp + 1];
	--nholp;
	iholp[1] = nholp;
    }
    pntupd_(ip, &icp[1], &xyp[4], &hesp[7], &detg[1], icps, &xyps[4], &hesps[
	    7], detgs);
    return 0;
} /* pntadd_ */

/* ========================================================== */
/* @f2h@ */ /* Subroutine */ int pntupd_(integer *ip, integer *icp, 
	doublereal *xyp, doublereal *hesp, doublereal *detg, integer *icps, 
	doublereal *xyps, doublereal *hesps, doublereal *detgs)
{
    static integer i__;

/* ========================================================== */
    /* Parameter adjustments */
    --hesps;
    --xyps;
    --detg;
    hesp -= 7;
    xyp -= 4;
    --icp;

    /* Function Body */
    icp[*ip] = *icps;
    for (i__ = 1; i__ <= 3; ++i__) {
	xyp[i__ + *ip * 3] = xyps[i__];
    }
    for (i__ = 1; i__ <= 6; ++i__) {
	hesp[i__ + *ip * 6] = hesps[i__];
    }
    detg[*ip] = *detgs;
    return 0;
} /* pntupd_ */

/* ========================================================== */
/* @f2h@ */ /* Subroutine */ int pntdel_(integer *ip, integer *np, integer *
	icp, integer *iholp)
{
    static integer nholp;

/* ========================================================== */
    /* Parameter adjustments */
    --iholp;
    --icp;

    /* Function Body */
    nholp = iholp[1];
    ++nholp;
    iholp[nholp + 1] = *ip;
    icp[*ip] = -icp[*ip];
    iholp[1] = nholp;
    --(*np);
    return 0;
} /* pntdel_ */

/* ========================================================== */
/* FACES ... FACE ... FACES ... FACE */
/* ========================================================== */
/* @f2h@ */ /* Subroutine */ int facadd_(integer *if__, integer *nf, integer *
	maxf, integer *iholf)
{
    static integer nholf;

/* ========================================================== */
    /* Parameter adjustments */
    --iholf;

    /* Function Body */
    nholf = iholf[1];
    ++(*nf);
    if (*nf > *maxf) {
	errmes_(&c__1004, "facAdd", "local parameter MaxF is small", (ftnlen)
		6, (ftnlen)29);
    }
    if (nholf == 0) {
	*if__ = *nf;
    } else {
	*if__ = iholf[nholf + 1];
	--nholf;
	iholf[1] = nholf;
    }
    return 0;
} /* facadd_ */

/* ========================================================== */
/* @f2h@ */ /* Subroutine */ int facupd_(integer *nfs, integer *ipf, integer *
	ifs, integer *ipfs)
{
    static integer i__, if__;

/* ========================================================== */
    /* Parameter adjustments */
    ipfs -= 5;
    --ifs;
    ipf -= 5;

    /* Function Body */
    if__ = ifs[*nfs];
    for (i__ = 1; i__ <= 4; ++i__) {
	ipf[i__ + (if__ << 2)] = ipfs[i__ + (*nfs << 2)];
    }
    return 0;
} /* facupd_ */

/* ========================================================== */
/* @f2h@ */ /* Subroutine */ int facdel_(integer *if__, integer *nf, integer *
	ipf, integer *iholf)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, nholf;

/* ========================================================== */
    /* Parameter adjustments */
    --iholf;
    ipf -= 5;

    /* Function Body */
    nholf = iholf[1];
    ++nholf;
    iholf[nholf + 1] = *if__;
    for (i__ = 1; i__ <= 4; ++i__) {
	ipf[i__ + (*if__ << 2)] = -(i__1 = ipf[i__ + (*if__ << 2)], dabs(i__1))
		;
    }
    iholf[1] = nholf;
    --(*nf);
    return 0;
} /* facdel_ */

/* ========================================================== */
/* ELEMENTS ... ELEMENT ... ELEMENTS ... ELEMENT */
/* ========================================================== */
/* @f2h@ */ /* Subroutine */ int eleadd_(integer *ne, integer *maxe, integer *
	ihole)
{
/* ========================================================== */
    /* Parameter adjustments */
    --ihole;

    /* Function Body */
    if (*ne >= *maxe && ihole[1] == 0) {
	errmes_(&c__1006, "eleAdd", "local parameter MaxE is small", (ftnlen)
		6, (ftnlen)29);
    }
    return 0;
} /* eleadd_ */

/* ========================================================== */
/* @f2h@ */ /* Subroutine */ int eleupd_(integer *nes, integer *iep, integer *
	ipe, integer *ife, integer *iee, integer *lf, integer *le, integer *
	ifs, integer *ies, integer *ipfs, integer *ipes)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, k, n, i1, i2, i3, j1, j2, j3, ie, ke, ip1, ip2, ip3, 
	    jp1, jp2, jp3, ieu[1000], iet, iref[5];

/* ========================================================== */
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
/* group (Local variables) */
/* ========================================================== */
    /* Parameter adjustments */
    ipes -= 6;
    ipfs -= 5;
    --ies;
    --ifs;
    iee -= 5;
    ife -= 5;
    ipe -= 6;
    --iep;

    /* Function Body */
    iref[0] = 1;
    iref[1] = 2;
    iref[2] = 3;
    iref[3] = 4;
    iref[4] = 1;
    ie = ies[*nes];
    for (i__ = 1; i__ <= 4; ++i__) {
	ipe[i__ + ie * 5] = ipes[i__ + *nes * 5];
	iep[ipe[i__ + ie * 5]] = ie;
	ife[i__ + (ie << 2)] = 0;
    }
    ipe[ie * 5 + 5] = ipes[*nes * 5 + 5];
    for (i1 = 1; i1 <= 4; ++i1) {
	i2 = iref[i1];
	i3 = iref[i2];
	ip1 = ipes[i1 + *nes * 5];
	ip2 = ipes[i2 + *nes * 5];
	ip3 = ipes[i3 + *nes * 5];
	i__1 = *lf;
	for (n = 1; n <= i__1; ++n) {
	    if (ifs[n] > 0) {
		jp1 = ipfs[(n << 2) + 1];
		jp2 = ipfs[(n << 2) + 2];
		jp3 = ipfs[(n << 2) + 3];
		if (check33_(&ip1, &ip2, &ip3, &jp1, &jp2, &jp3)) {
		    ife[i1 + (ie << 2)] = ifs[n];
		    goto L10;
		}
	    }
	}
L10:
	if (ife[i1 + (ie << 2)] != 0) {
/*  ...  analyzing the detailed structure of this face */
/*  ...  the face may be interior and with boundary points */
	    maksp_(&ip1, &iep[1], &ipe[6], &iee[5], &c__1000, &ke, ieu);
	    i__1 = ke;
	    for (k = 1; k <= i__1; ++k) {
		iet = ieu[k - 1];
		i__2 = *le;
		for (n = 1; n <= i__2; ++n) {
		    if (iet == (i__3 = ies[n], dabs(i__3))) {
			goto L15;
		    }
		}
		for (j1 = 1; j1 <= 4; ++j1) {
		    j2 = iref[j1];
		    j3 = iref[j2];
		    jp1 = ipe[j1 + iet * 5];
		    jp2 = ipe[j2 + iet * 5];
		    jp3 = ipe[j3 + iet * 5];
		    if (check33_(&ip1, &ip2, &ip3, &jp1, &jp2, &jp3)) {
			goto L30;
		    }
		}
L15:
		;
	    }
	    iee[i1 + (ie << 2)] = 0;
	}
	i__1 = *le;
	for (k = 1; k <= i__1; ++k) {
	    ke = ies[k];
	    if (ke <= 0) {
		goto L20;
	    }
	    if (! check3j_(&ip1, &ip2, &ip3, &ipes[k * 5 + 1])) {
		goto L20;
	    }
	    if (ke == ie) {
		goto L20;
	    }
	    for (j1 = 1; j1 <= 4; ++j1) {
		j2 = iref[j1];
		j3 = iref[j2];
		jp1 = ipes[j1 + k * 5];
		jp2 = ipes[j2 + k * 5];
		jp3 = ipes[j3 + k * 5];
		if (check33_(&ip1, &ip2, &ip3, &jp1, &jp2, &jp3)) {
		    iee[i1 + (ie << 2)] = ke;
		    iee[j1 + (ke << 2)] = ie;
		    ife[j1 + (ke << 2)] = ife[i1 + (ie << 2)];
		    goto L30;
		}
	    }
L20:
	    ;
	}
L30:
	;
    }
    return 0;
} /* eleupd_ */

/* ========================================================== */
/* @f2h@ */ /* Subroutine */ int eledel_(integer *ie, integer *ipe, integer *
	iee)
{
    static integer i__, j, iet;

/* ========================================================== */
    /* Parameter adjustments */
    iee -= 5;
    ipe -= 6;

    /* Function Body */
    ipe[*ie * 5 + 1] = 0;
    for (i__ = 1; i__ <= 4; ++i__) {
	iet = iee[i__ + (*ie << 2)];
	if (iet != 0) {
	    for (j = 1; j <= 4; ++j) {
		if (iee[j + (iet << 2)] == *ie) {
		    iee[j + (iet << 2)] = 0;
		}
	    }
	    iee[i__ + (*ie << 2)] = 0;
	}
    }
    return 0;
} /* eledel_ */

