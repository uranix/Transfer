/* ZZ.f -- translated by f2c (version 20090411).
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
#include "zz.h"
#include "lapack.h"

/* Table of constant values */

static integer c__1001 = 1001;
static integer c__1002 = 1002;
static integer c__4 = 4;
static integer c__1007 = 1007;
static integer c__2011 = 2011;
static integer c__1 = 1;
static integer c__40 = 40;
static integer c__3011 = 3011;

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int p02p1_(integer *np, integer *nf, integer *ne,
	 doublereal *xyp, integer *ipf, integer *ipe, doublereal *fp0, 
	doublereal *fp1, integer *maxwr, integer *maxwi, doublereal *rw, 
	integer *iw)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, n;
    static doublereal s;
    static integer i1, i2, ie, ip;
    static doublereal xy[3000]	/* was [3][1000] */;
    static integer ip1, ip2, ip3, ip4, iiw;
    static doublereal vol;
    static integer irw, ibnd, iiep, kbuf[1000], iarm, lbuf, inep, imass;
    static doublereal values[1000], weights[1000];

/* ================================================================ */
/* ================================================================ */
/* ================================================================ */
/*  The routine maps a discontinuous piece-wise constant function */
/*  with d.o.f in elements onto a continuous piece-wise linear */
/*  function with d.o.f. at vertices. We use the ZZ method for that. */

/*  Limitation to the mesh: Each boundary node can be connected with */
/*                          an interior note by at most 2 mesh edges. */

/*  *** Remarks */
/*         1. Extrapolation at boundary nodes adds some diffusion */
/*            for functions with sharp gradients near boundary. */
/* ================================================================ */
/*  *** Input */
/*         nP  - the number of mesh nodes */
/*         nF  - the number of boundary edges */
/*         nE  - the number of triangles */

/*         XYP - coordinates of the mesh nodes */
/*         IPF - map: boundary edge -> vertices */
/*         IPE - map: triangle -> vertices */

/*         fP0 - discontinuous piece-wise constant function */

/*  *** Output: */
/*         fP1 - continuous piece-wice linear function */

/*  *** Work memory: */
/*          rW - real*8  array of size MaxWr */
/*          iW - integer array of size MaxWi */
/* ================================================================ */
/* ====================================================================== */
/* MaxS evaluates the number of elements in a superelement. */
/* Additionaly, it bounds the number of different boundary identificators. */
/* ====================================================================== */
/* ====================================================================== */
/* Colors: iVface - a fixed face */
/*         iMface - a material face */
/* ====================================================================== */
/* ... external functions */
/* ... local variables */
/* ================================================================ */
/* ... distribute memory */
    /* Parameter adjustments */
    --iw;
    --rw;
    --fp1;
    --fp0;
    ipe -= 5;
    ipf -= 4;
    xyp -= 4;

    /* Function Body */
    ibnd = 0;
    iarm = ibnd + *np;
    inep = iarm + *np;
    iiep = inep + *np;
    iiw = iiep + (*ne << 2);
    imass = 0;
    irw = imass + *np;
    if (*maxwi < iiw) {
	errmes_(&c__1001, "P02P1", "Increase size of MaxWi", (ftnlen)5, (
		ftnlen)22);
    }
    if (*maxwr < irw) {
	errmes_(&c__1002, "P02P1", "Increase size of MaxWr", (ftnlen)5, (
		ftnlen)22);
    }
/* ... create maps: vertix -> triangles */
    backreferences_(np, ne, &c__4, &c__4, &ipe[5], &iw[inep + 1], &iw[iiep + 
	    1]);
/* ... mark boundary nodes */
    i__1 = *np << 1;
    for (n = 1; n <= i__1; ++n) {
	iw[n] = 0;
    }
    i__1 = *nf;
    for (n = 1; n <= i__1; ++n) {
	for (i__ = 1; i__ <= 3; ++i__) {
	    ip = ipf[i__ + n * 3];
	    iw[ibnd + ip] = 1;
	}
    }
/* ... initialize Mass and fP1 */
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	rw[imass + n] = 0.;
	fp1[n] = 0.;
    }
/* ... interpolate in interior points with least squares */
    i2 = 0;
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	i1 = i2 + 1;
	i2 = iw[inep + n];
	if (i2 - i1 + 1 > 1000) {
	    errmes_(&c__1007, "DG2P1", "Local parameter MaxS is small", (
		    ftnlen)5, (ftnlen)29);
	}
	i__2 = i2;
	for (i__ = i1; i__ <= i__2; ++i__) {
	    ie = iw[iiep + i__];
	    ip1 = ipe[(ie << 2) + 1];
	    ip2 = ipe[(ie << 2) + 2];
	    ip3 = ipe[(ie << 2) + 3];
	    ip4 = ipe[(ie << 2) + 4];
	    vol = calvol_(&xyp[ip1 * 3 + 1], &xyp[ip2 * 3 + 1], &xyp[ip3 * 3 
		    + 1], &xyp[ip4 * 3 + 1]);
	    vol = dabs(vol);
	    if (iw[ibnd + n] > 0) {
		for (j = 1; j <= 4; ++j) {
		    ip = ipe[j + (ie << 2)];
		    if (iw[ibnd + ip] == 0) {
			rw[imass + n] += vol / 20;
		    }
		}
	    } else {
		rw[imass + n] += vol / 4;
		k = i__ - i1 + 1;
		values[k - 1] = fp0[ie];
		for (j = 1; j <= 3; ++j) {
		    xy[j + k * 3 - 4] = (xyp[j + ip1 * 3] + xyp[j + ip2 * 3] 
			    + xyp[j + ip3 * 3] + xyp[j + ip4 * 3]) / 4;
		}
	    }
	}
	if (iw[ibnd + n] == 0) {
	    k = i2 - i1 + 1;
	    fp1[n] = lsvalue_(&k, xy, values, &xyp[n * 3 + 1]);
	}
    }
/* ... extrapolate fP1 at boundary nodes from nearest interior nodes */
    i2 = 0;
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	i1 = i2 + 1;
	i2 = iw[inep + n];
	if (iw[ibnd + n] == 1) {
	    lbuf = 0;
	    i__2 = i2;
	    for (i__ = i1; i__ <= i__2; ++i__) {
		ie = iw[iiep + i__];
		ip1 = ipe[(ie << 2) + 1];
		ip2 = ipe[(ie << 2) + 2];
		ip3 = ipe[(ie << 2) + 3];
		ip4 = ipe[(ie << 2) + 4];
		vol = calvol_(&xyp[ip1 * 3 + 1], &xyp[ip2 * 3 + 1], &xyp[ip3 *
			 3 + 1], &xyp[i__ + ip4 * 3]);
		vol = dabs(vol);
		for (j = 1; j <= 4; ++j) {
		    ip = ipe[j + (ie << 2)];
		    findse_(&lbuf, kbuf, &ip, &k);
		    if (k > 0) {
			weights[k - 1] += vol / 20;
		    } else if (iw[ibnd + ip] == 0) {
			++lbuf;
			kbuf[lbuf - 1] = ip;
			weights[lbuf - 1] = vol / 20;
		    }
		}
	    }
	    i__2 = lbuf;
	    for (k = 1; k <= i__2; ++k) {
		ip = kbuf[k - 1];
		fp1[n] += fp1[ip] * weights[k - 1] / rw[imass + n];
	    }
	    if (lbuf > 0) {
		iw[iarm + n] = 1;
	    }
	}
    }
/* ... remove boundary node from the list */
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	if (iw[iarm + n] == 1) {
	    iw[ibnd + n] = 0;
	}
    }
/* ... the second level of extrapolation */
    i2 = 0;
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	i1 = i2 + 1;
	i2 = iw[inep + n];
	if (iw[ibnd + n] == 1) {
	    k = 0;
	    s = 0.;
	    i__2 = i2;
	    for (i__ = i1; i__ <= i__2; ++i__) {
		ie = iw[iiep + i__];
		for (j = 1; j <= 4; ++j) {
		    ip = ipe[j + (ie << 2)];
		    if (iw[ibnd + ip] == 0) {
			++k;
			s += fp1[ip];
		    }
		}
	    }
	    if (k == 0) {
		errmes_(&c__2011, "P02P1", "Mesh violates 2-arm rule", (
			ftnlen)5, (ftnlen)24);
	    }
	    fp1[n] = s / k;
	}
    }
    return 0;
} /* p02p1_ */

/* ================================================================ */
/* @f2h@ */ doublereal lsvalue_(integer *k, doublereal *xy, doublereal *
	values, doublereal *xy0)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static doublereal a[16]	/* was [4][4] */;
    static integer i__, j;
    static doublereal s[4];
    static integer info, ipiv[4];
    static doublereal work[40];

/* ================================================================ */
/*  This routine uses least square linear fit to points xy(2,*) and */
/*  evaluates the value of the linear function at point xy0. */
/* ================================================================ */
/* (local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --xy0;
    --values;
    xy -= 4;

    /* Function Body */
    for (i__ = 1; i__ <= 4; ++i__) {
	for (j = 1; j <= 4; ++j) {
	    a[i__ + (j << 2) - 5] = 0.;
	}
	s[i__ - 1] = 0.;
    }
/* ... generate the least squares matrix */
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a[0] += xy[i__ * 3 + 1] * xy[i__ * 3 + 1];
	a[4] += xy[i__ * 3 + 1] * xy[i__ * 3 + 2];
	a[8] += xy[i__ * 3 + 1] * xy[i__ * 3 + 3];
	a[12] += xy[i__ * 3 + 1];
	a[5] += xy[i__ * 3 + 2] * xy[i__ * 3 + 2];
	a[9] += xy[i__ * 3 + 2] * xy[i__ * 3 + 3];
	a[13] += xy[i__ * 3 + 2];
	a[10] += xy[i__ * 3 + 3] * xy[i__ * 3 + 3];
	a[14] += xy[i__ * 3 + 3];
    }
    a[15] = (doublereal) (*k);
/* ... generate the RHS */
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s[0] += xy[i__ * 3 + 1] * values[i__];
	s[1] += xy[i__ * 3 + 2] * values[i__];
	s[2] += xy[i__ * 3 + 3] * values[i__];
	s[3] += values[i__];
    }
    dsysv_("U", &c__4, &c__1, a, &c__4, ipiv, s, &c__4, work, &c__40, &info, (ftnlen)1);
    if (info != 0) {
	errmes_(&c__3011, "LSvalue", "Error in Lapack routine dsysv", (ftnlen)7, (ftnlen)29);
    }
/* ... evaluate the linear function */
    ret_val = s[0] * xy0[1] + s[1] * xy0[2] + s[2] * xy0[3] + s[3];
    return ret_val;
} /* lsvalue_ */

