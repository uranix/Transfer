/* lintrp3D.f -- translated by f2c (version 20090411).
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
#include "makQ.h"
#include "makM.h"
#include "auxSP.h"
#include "auxSE.h"
#include "auxSF.h"
#include "auxSR.h"
#include "lintrp3d.h"
#include "clpSR.h"
#include "update.h"
#include "list.h"

/* Table of constant values */

static integer c__6101 = 6101;
static integer c__6102 = 6102;
static integer c__1001 = 1001;
static integer c__1002 = 1002;
static integer c__1009 = 1009;
static integer c__6103 = 6103;
static integer c__6104 = 6104;
static integer c__1 = 1;
static doublereal c_b60 = 1e-6;
static integer c__6105 = 6105;

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int lintrp3d_(integer *nt, integer *tet, integer 
	*nv, doublereal *vrt, integer *ldf, doublereal *f, integer *nxyz, 
	doublereal *xyz, doublereal *g, integer *imem, integer *nimem, 
	doublereal *dmem, integer *ndmem, integer *icontrol)
{
    /* System generated locals */
    integer f_dim1, f_offset, g_dim1, g_offset, i__1, i__2;

    /* Local variables */
    static doublereal h__[32];
    static integer i__, qt, ref, nqt, itet;
    static logical flags[2];
    static integer nunit;

/* ================================================================ */
/* ================================================================ */
/*  Subrouine  LINTRP3D  interpolates the peicewise linear  function */
/*  determined at the nodes of 3D triangulation onto a given set of points. */

/*  Parameters of subroutine in order to locate (without typing): */
/*       NT         - number of tetrahedra */
/*       tet(4, NT) - list of tetrahedra ( 4 vertices ) */
/*       NV         - number of nodes (vetices) of  triangulation */
/*       vrt(3, NV) - coords of vertices */

/*       LDF          - leading dimension of vector function F and G */
/*       F(LDF, NV)   - the values of vector-function in the nodes of */
/*                      triangulation */
/*       NXY          - number of points where the value of function to be */
/*                      determined */
/*       xyz(3, NXYZ) - coords of points */
/*       G(LDF, NXYZ) - the returned values of vector-function in given */
/*                      points */

/*       imem(Nimem) - auxiliary array of integers of length Nimem */
/*       dmem(Ndmem) - auxiliary array of double precision numbers of */
/*                     length Ndmem */

/*       iContol - 4-digit integer representing 3 control parameters */
/*            1st digit = 1 means "Initialisation is needed" */
/*                      = 2 means "Initialisation is not needed" */
/*                      otherwise the input is erroneous */
/*            2nd digit  = 1 means  "Nearby a true curved boundary patch */
/*                                   the domain is convex" */
/*                         2 means  "Nearby a true curved boundary patch */
/*                                   the domain is not convex" */
/*                        otherwise  the input is erroneous */
/*            3,4 digits = 00 means "No output information is provided" */
/*                  ab   > 00 means "Important warnings are written to */
/*                                   the channel (NUNIT=ab )" */
/* ================================================================ */
/*  REMARK 1. */
/*    If the domain is polyhedral, and for a given point there is no an */
/*    element containing it within tolerance PREC=10^{-6} , then the ERROR */
/*    is assumed. The code could not find an element within PREC tolerance */
/*    for MaxTrials. */

/*  REMARK 2. */
/*    If the domain has a curved boundary, it is possible that a given */
/*    point is in the domain but out of the mesh. In this case the tolerance */
/*    PREC will be relaxed until an element will be found. Yet, the warning */
/*    about the point will be issued if NUNIT > 0. That is why for the first */
/*    usages of the code NUNIT > 0 are recommended. */

/* ================================================================ */
/* ================================================================ */
/*      Parameter(RelativeTrials = 1D0) */
/* ================================================================ */
/* ... decoding control parameters */
    /* Parameter adjustments */
    tet -= 5;
    vrt -= 4;
    f_dim1 = *ldf;
    f_offset = 1 + f_dim1;
    f -= f_offset;
    g_dim1 = *ldf;
    g_offset = 1 + g_dim1;
    g -= g_offset;
    xyz -= 4;
    --imem;
    --dmem;

    /* Function Body */
    i__ = *icontrol / 1000;
    if (i__ == 1) {
	flags[0] = TRUE_;
    } else if (i__ == 2) {
	flags[0] = FALSE_;
    } else {
	errmes_(&c__6101, "lintrp3D", "Wrong 1st digit in iControl", (ftnlen)
		8, (ftnlen)27);
    }
    i__ = (*icontrol - *icontrol / 1000 * 1000) / 100;
    if (i__ == 1) {
	flags[1] = TRUE_;
    } else if (i__ == 2) {
	flags[1] = FALSE_;
    } else {
	errmes_(&c__6102, "lintrp3D", "Wrong 2nd digit in iControl", (ftnlen)
		8, (ftnlen)27);
    }
    nunit = *icontrol % 100;
/* ... distributing working memory */
    if (flags[0]) {
	qt = 5;
	initqt_(&nqt, &imem[qt], h__, &dmem[33]);
	i__1 = *nv;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *nimem - 4;
	    drops_(&nqt, &imem[qt], &i__2, h__, &dmem[33], &vrt[4], &i__);
	}
	ref = qt + (nqt << 3);
	itet = ref + *nv + 1;
    } else {
	qt = imem[1];
	nqt = imem[2];
	ref = imem[3];
	itet = imem[4];
	for (i__ = 1; i__ <= 32; ++i__) {
	    h__[i__ - 1] = dmem[i__];
	}
    }
    if (itet + (*nt << 2) > *nimem) {
	errmes_(&c__1001, "lintrp3D", "not enough memory for Integer arrays", 
		(ftnlen)8, (ftnlen)36);
    }
    if (nqt * 3 + 32 > *ndmem) {
	errmes_(&c__1002, "lintrp3D", "not enough memory for Real*8 arrays", (
		ftnlen)8, (ftnlen)35);
    }
/* ... calling the main module */
    i__1 = *nimem - itet - (*nt << 2);
    restore_(nt, &tet[5], nv, &vrt[4], ldf, &f[f_offset], nxyz, &xyz[4], &g[
	    g_offset], &imem[qt], &dmem[33], &imem[ref], &imem[itet], h__, 
	    flags, &imem[itet + (*nt << 2)], &i__1, &nunit);
/* ... saving memory distribution */
    if (flags[0]) {
	imem[1] = qt;
	imem[2] = nqt;
	imem[3] = ref;
	imem[4] = itet;
	for (i__ = 1; i__ <= 32; ++i__) {
	    dmem[i__] = h__[i__ - 1];
	}
    }
    return 0;
} /* lintrp3d_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int initqt_(integer *nqt, integer *qt, 
	doublereal *h__, doublereal *xyzc)
{
    static integer k;

/* ================================================================ */
/*     Initializing of octtree sructure */
/* ================================================================ */
/* ================================================================ */
/*      Parameter(RelativeTrials = 1D0) */
/* ================================================================ */
    /* Parameter adjustments */
    xyzc -= 4;
    --h__;
    qt -= 9;

    /* Function Body */
    h__[1] = .5f;
    for (k = 2; k <= 32; ++k) {
	h__[k] = h__[k - 1] / 2;
    }
    *nqt = 0;
    *nqt = newqt_(nqt, &qt[9]);
    xyzc[*nqt * 3 + 1] = .5f;
    xyzc[*nqt * 3 + 2] = .5f;
    xyzc[*nqt * 3 + 3] = .5f;
    return 0;
} /* initqt_ */

/* ================================================================ */
/* @f2h@ */ integer newqt_(integer *nqt, integer *qt)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer i__;

/* ================================================================ */
/*  The function allocates a new octtree */
/* ================================================================ */
/* ================================================================ */
/* ================================================================ */
    /* Parameter adjustments */
    qt -= 9;

    /* Function Body */
    ++(*nqt);
    for (i__ = 1; i__ <= 8; ++i__) {
	qt[i__ + (*nqt << 3)] = 0;
    }
    ret_val = *nqt;
    return ret_val;
} /* newqt_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int setij_(doublereal *xyzc, doublereal *xyz, 
	integer *i__, integer *j, integer *k)
{
/* ================================================================ */
/*  Definition of octant of point XYZ with respect to the centre of */
/*  quadtree XYZc. */
/* ================================================================ */
/* ================================================================ */
    /* Parameter adjustments */
    --xyz;
    --xyzc;

    /* Function Body */
    if (xyz[1] > xyzc[1]) {
	*i__ = 2;
    } else {
	*i__ = 1;
    }
    if (xyz[2] > xyzc[2]) {
	*j = 2;
    } else {
	*j = 1;
    }
    if (xyz[3] > xyzc[3]) {
	*k = 2;
    } else {
	*k = 1;
    }
    return 0;
} /* setij_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int drops_(integer *nqt, integer *qt, integer *
	nqtav, doublereal *h__, doublereal *xyzc, doublereal *xyz, integer *
	idx)
{
    static integer i__, j, k, l, ip, dir[2], new__, ptr;

/* ================================================================ */
/*     The arrangement of a new point in octtree */
/* ================================================================ */
/* ================================================================ */
/*      Parameter(RelativeTrials = 1D0) */
/* ================================================================ */
    /* Parameter adjustments */
    xyz -= 4;
    xyzc -= 4;
    --h__;
    qt -= 15;

    /* Function Body */
    dir[0] = -1;
    dir[1] = 1;
    ip = 1;
    l = 1;
    while(TRUE_) {
	++l;
	if (l > 32) {
	    goto L1000;
	}
	setij_(&xyzc[ip * 3 + 1], &xyz[*idx * 3 + 1], &i__, &j, &k);
	ptr = qt[i__ + (j + (k + (ip << 1) << 1) << 1)];
	if (ptr == 0) {
	    qt[i__ + (j + (k + (ip << 1) << 1) << 1)] = *idx;
	    return 0;
	} else if (ptr < 0) {
	    ip = -ptr;
	} else if (ptr > 0) {
	    new__ = newqt_(nqt, &qt[15]);
	    if (*nqt > *nqtav) {
		errmes_(&c__1001, "drops", "not enough memory for Integer ar"
			"rays", (ftnlen)5, (ftnlen)36);
	    }
	    qt[i__ + (j + (k + (ip << 1) << 1) << 1)] = -new__;
	    xyzc[new__ * 3 + 1] = xyzc[ip * 3 + 1] + dir[i__ - 1] * h__[l];
	    xyzc[new__ * 3 + 2] = xyzc[ip * 3 + 2] + dir[j - 1] * h__[l];
	    xyzc[new__ * 3 + 3] = xyzc[ip * 3 + 3] + dir[k - 1] * h__[l];
	    setij_(&xyzc[new__ * 3 + 1], &xyz[ptr * 3 + 1], &i__, &j, &k);
	    qt[i__ + (j + (k + (new__ << 1) << 1) << 1)] = ptr;
	    ip = new__;
	}
    }
    return 0;
L1000:
    errmes_(&c__1009, "drops", "local parameter MaxH is small", (ftnlen)5, (
	    ftnlen)29);
    return 0;
} /* drops_ */

/* ================================================================ */
/* @f2h@ */ doublereal sqrdst_(doublereal *a, doublereal *b)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3;

/* ================================================================ */
/*  Function returns the square of distance between points of a space */
/* ================================================================ */
    /* Parameter adjustments */
    --b;
    --a;

    /* Function Body */
/* Computing 2nd power */
    d__1 = a[1] - b[1];
/* Computing 2nd power */
    d__2 = a[2] - b[2];
/* Computing 2nd power */
    d__3 = a[3] - b[3];
    ret_val = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
    return ret_val;
} /* sqrdst_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int order2_(doublereal *xyzc, doublereal *h__, 
	doublereal *xyz, integer *ord, doublereal *sqrd)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k;
    static doublereal sqr[6]	/* was [3][2] */, dist[3];

/* ================================================================ */
/*  Ordering of subocts of octtree QT by their distanse from */
/*  the given point XYZ */
/* ================================================================ */
/* ================================================================ */
/* ================================================================ */
/* ================================================================ */
    /* Parameter adjustments */
    --sqrd;
    --ord;
    --xyz;
    --xyzc;

    /* Function Body */
    for (i__ = 1; i__ <= 3; ++i__) {
	dist[i__ - 1] = (d__1 = xyz[i__] - xyzc[i__], dabs(d__1));
/* Computing 2nd power */
	d__1 = dist[i__ - 1];
	sqr[i__ - 1] = d__1 * d__1;
/* Computing 2nd power */
	d__1 = d_dim(&dist[i__ - 1], h__);
	sqr[i__ + 2] = d__1 * d__1;
    }
    setij_(&xyzc[1], &xyz[1], &i__, &j, &k);
    ord[1] = i__ + 2 * (j - 1) + 4 * (k - 1);
    sqrd[1] = sqr[0] + sqr[1] + sqr[2];
    if (dist[0] <= dist[1] && dist[1] <= dist[2]) {
	i__1 = 3 - i__;
	ord[2] = i__1 + 2 * (j - 1) + 4 * (k - 1);
	i__1 = 3 - j;
	ord[3] = i__ + 2 * (i__1 - 1) + 4 * (k - 1);
	i__1 = 3 - i__;
	i__2 = 3 - j;
	ord[4] = i__1 + 2 * (i__2 - 1) + 4 * (k - 1);
	i__1 = 3 - k;
	ord[5] = i__ + 2 * (j - 1) + 4 * (i__1 - 1);
	i__1 = 3 - i__;
	i__2 = 3 - k;
	ord[6] = i__1 + 2 * (j - 1) + 4 * (i__2 - 1);
	i__1 = 3 - j;
	i__2 = 3 - k;
	ord[7] = i__ + 2 * (i__1 - 1) + 4 * (i__2 - 1);
	sqrd[2] = sqr[0] + sqr[4] + sqr[5];
	sqrd[3] = sqr[3] + sqr[1] + sqr[5];
	sqrd[4] = sqr[0] + sqr[1] + sqr[5];
	sqrd[5] = sqr[3] + sqr[4] + sqr[2];
	sqrd[6] = sqr[0] + sqr[4] + sqr[2];
	sqrd[7] = sqr[3] + sqr[1] + sqr[2];
	goto L1;
    }
    if (dist[1] <= dist[0] && dist[0] <= dist[2]) {
	i__1 = 3 - j;
	ord[2] = i__ + 2 * (i__1 - 1) + 4 * (k - 1);
	i__1 = 3 - i__;
	ord[3] = i__1 + 2 * (j - 1) + 4 * (k - 1);
	i__1 = 3 - i__;
	i__2 = 3 - j;
	ord[4] = i__1 + 2 * (i__2 - 1) + 4 * (k - 1);
	i__1 = 3 - k;
	ord[5] = i__ + 2 * (j - 1) + 4 * (i__1 - 1);
	i__1 = 3 - j;
	i__2 = 3 - k;
	ord[6] = i__ + 2 * (i__1 - 1) + 4 * (i__2 - 1);
	i__1 = 3 - i__;
	i__2 = 3 - k;
	ord[7] = i__1 + 2 * (j - 1) + 4 * (i__2 - 1);
	sqrd[2] = sqr[3] + sqr[1] + sqr[5];
	sqrd[3] = sqr[0] + sqr[4] + sqr[5];
	sqrd[4] = sqr[0] + sqr[1] + sqr[5];
	sqrd[5] = sqr[3] + sqr[4] + sqr[2];
	sqrd[6] = sqr[3] + sqr[1] + sqr[2];
	sqrd[7] = sqr[0] + sqr[4] + sqr[2];
	goto L1;
    }
    if (dist[0] <= dist[2] && dist[2] <= dist[1]) {
	i__1 = 3 - i__;
	ord[2] = i__1 + 2 * (j - 1) + 4 * (k - 1);
	i__1 = 3 - k;
	ord[3] = i__ + 2 * (j - 1) + 4 * (i__1 - 1);
	i__1 = 3 - i__;
	i__2 = 3 - k;
	ord[4] = i__1 + 2 * (j - 1) + 4 * (i__2 - 1);
	i__1 = 3 - j;
	ord[5] = i__ + 2 * (i__1 - 1) + 4 * (k - 1);
	i__1 = 3 - i__;
	i__2 = 3 - j;
	ord[6] = i__1 + 2 * (i__2 - 1) + 4 * (k - 1);
	i__1 = 3 - j;
	i__2 = 3 - k;
	ord[7] = i__ + 2 * (i__1 - 1) + 4 * (i__2 - 1);
	sqrd[2] = sqr[0] + sqr[4] + sqr[5];
	sqrd[3] = sqr[3] + sqr[4] + sqr[2];
	sqrd[4] = sqr[0] + sqr[4] + sqr[2];
	sqrd[5] = sqr[3] + sqr[1] + sqr[5];
	sqrd[6] = sqr[0] + sqr[1] + sqr[5];
	sqrd[7] = sqr[3] + sqr[1] + sqr[2];
	goto L1;
    }
    if (dist[1] <= dist[2] && dist[2] <= dist[0]) {
	i__1 = 3 - j;
	ord[2] = i__ + 2 * (i__1 - 1) + 4 * (k - 1);
	i__1 = 3 - k;
	ord[3] = i__ + 2 * (j - 1) + 4 * (i__1 - 1);
	i__1 = 3 - j;
	i__2 = 3 - k;
	ord[4] = i__ + 2 * (i__1 - 1) + 4 * (i__2 - 1);
	i__1 = 3 - i__;
	ord[5] = i__1 + 2 * (j - 1) + 4 * (k - 1);
	i__1 = 3 - i__;
	i__2 = 3 - j;
	ord[6] = i__1 + 2 * (i__2 - 1) + 4 * (k - 1);
	i__1 = 3 - i__;
	i__2 = 3 - k;
	ord[7] = i__1 + 2 * (j - 1) + 4 * (i__2 - 1);
	sqrd[2] = sqr[3] + sqr[1] + sqr[5];
	sqrd[3] = sqr[3] + sqr[4] + sqr[2];
	sqrd[4] = sqr[3] + sqr[1] + sqr[2];
	sqrd[5] = sqr[0] + sqr[4] + sqr[5];
	sqrd[6] = sqr[0] + sqr[1] + sqr[5];
	sqrd[7] = sqr[0] + sqr[4] + sqr[2];
	goto L1;
    }
    if (dist[2] <= dist[1] && dist[1] <= dist[0]) {
	i__1 = 3 - k;
	ord[2] = i__ + 2 * (j - 1) + 4 * (i__1 - 1);
	i__1 = 3 - j;
	ord[3] = i__ + 2 * (i__1 - 1) + 4 * (k - 1);
	i__1 = 3 - j;
	i__2 = 3 - k;
	ord[4] = i__ + 2 * (i__1 - 1) + 4 * (i__2 - 1);
	i__1 = 3 - i__;
	ord[5] = i__1 + 2 * (j - 1) + 4 * (k - 1);
	i__1 = 3 - i__;
	i__2 = 3 - k;
	ord[6] = i__1 + 2 * (j - 1) + 4 * (i__2 - 1);
	i__1 = 3 - i__;
	i__2 = 3 - j;
	ord[7] = i__1 + 2 * (i__2 - 1) + 4 * (k - 1);
	sqrd[2] = sqr[3] + sqr[4] + sqr[2];
	sqrd[3] = sqr[3] + sqr[1] + sqr[5];
	sqrd[4] = sqr[3] + sqr[1] + sqr[2];
	sqrd[5] = sqr[0] + sqr[4] + sqr[5];
	sqrd[6] = sqr[0] + sqr[4] + sqr[2];
	sqrd[7] = sqr[0] + sqr[1] + sqr[5];
	goto L1;
    }
    if (dist[2] <= dist[0] && dist[0] <= dist[1]) {
	i__1 = 3 - k;
	ord[2] = i__ + 2 * (j - 1) + 4 * (i__1 - 1);
	i__1 = 3 - i__;
	ord[3] = i__1 + 2 * (j - 1) + 4 * (k - 1);
	i__1 = 3 - i__;
	i__2 = 3 - k;
	ord[4] = i__1 + 2 * (j - 1) + 4 * (i__2 - 1);
	i__1 = 3 - j;
	ord[5] = i__ + 2 * (i__1 - 1) + 4 * (k - 1);
	i__1 = 3 - j;
	i__2 = 3 - k;
	ord[6] = i__ + 2 * (i__1 - 1) + 4 * (i__2 - 1);
	i__1 = 3 - i__;
	i__2 = 3 - j;
	ord[7] = i__1 + 2 * (i__2 - 1) + 4 * (k - 1);
	sqrd[2] = sqr[3] + sqr[4] + sqr[2];
	sqrd[3] = sqr[0] + sqr[4] + sqr[5];
	sqrd[4] = sqr[0] + sqr[4] + sqr[2];
	sqrd[5] = sqr[3] + sqr[1] + sqr[5];
	sqrd[6] = sqr[3] + sqr[1] + sqr[2];
	sqrd[7] = sqr[0] + sqr[1] + sqr[5];
	goto L1;
    }
L1:
    i__1 = 3 - i__;
    i__2 = 3 - j;
    i__3 = 3 - k;
    ord[8] = i__1 + 2 * (i__2 - 1) + 4 * (i__3 - 1);
    sqrd[8] = sqr[3] + sqr[4] + sqr[5];
    return 0;
} /* order2_ */

/* ================================================================ */
/* @f2h@ */ integer nearst_(integer *qt, doublereal *xyzc, doublereal *xyz, 
	doublereal *point, doublereal *h__)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer k[32], l, ip[32];
    static doublereal min__;
    static integer ord[256]	/* was [8][32] */, ptr;
    static doublereal sqrd[256]	/* was [8][32] */, sqdist;

/* ================================================================ */
/*     This function returns the index of the octtree cell which contains */
/*     the nearest to given point  grid node */
/* ================================================================ */
/* ================================================================ */
/*      Parameter(RelativeTrials = 1D0) */
/* ================================================================ */
    /* Parameter adjustments */
    --h__;
    --point;
    xyz -= 4;
    xyzc -= 4;
    qt -= 9;

    /* Function Body */
    l = 1;
    ptr = -1;
    min__ = 2.f;
    ret_val = 0;
    while(TRUE_) {
	if (ptr < 0) {
	    ++l;
	    if (l > 32) {
		goto L1000;
	    }
	    ip[l - 1] = -ptr;
	    order2_(&xyzc[-ptr * 3 + 1], &h__[l], &point[1], &ord[(l << 3) - 
		    8], &sqrd[(l << 3) - 8]);
	    k[l - 1] = 1;
	} else {
	    if (ptr > 0) {
		sqdist = sqrdst_(&point[1], &xyz[ptr * 3 + 1]);
		if (min__ > sqdist) {
		    min__ = sqdist;
		    ret_val = ord[k[l - 1] + (l << 3) - 9] + (ip[l - 1] - 1 <<
			     3);
		}
	    }
	    while(TRUE_) {
		++k[l - 1];
		if (k[l - 1] <= 8) {
		    if (sqrd[k[l - 1] + (l << 3) - 9] < min__) {
			goto L101;
		    }
		}
		--l;
		if (l == 1) {
		    goto L102;
		}
	    }
	}
L101:
	ptr = qt[ord[k[l - 1] + (l << 3) - 9] + (ip[l - 1] << 3)];
    }
L102:
    return ret_val;
L1000:
    errmes_(&c__1009, "nearest", "local parameter MaxH is small", (ftnlen)7, (
	    ftnlen)29);
    return ret_val;
} /* nearst_ */

/* ================================================================ */
/* @f2h@ */ integer basetet_(integer *qt, doublereal *xyzc, doublereal *vrt, 
	integer *nt, integer *tet, integer *ref, integer *itet, doublereal *
	xyz, doublereal *h__, integer *buf, integer *nbuf, doublereal *prec)
{
    /* System generated locals */
    integer ret_val, i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__, k, n, maxtrials, ip, idx;

/* ================================================================ */
/*  The function  determines the underlying tetrahedron */
/*  for a point XYZ inside of triangulation */
/* ================================================================ */
/* ================================================================ */
/*      Parameter(RelativeTrials = 1D0) */
/* ================================================================ */
    /* Parameter adjustments */
    buf -= 3;
    --h__;
    --xyz;
    --itet;
    --ref;
    tet -= 5;
    vrt -= 4;
    xyzc -= 4;
    --qt;

    /* Function Body */
    ret_val = 0;
    k = 0;
/* Computing MAX */
    d__1 = *nt * .1;
    maxtrials = (integer) max(d__1,20.);
    while(k < maxtrials) {
	ip = nearst_(&qt[1], &xyzc[4], &vrt[4], &xyz[1], &h__[1]);
	if (ip == 0) {
	    goto L101;
	}
	idx = qt[ip];
	i__1 = ref[idx + 1] - 1;
	for (i__ = ref[idx]; i__ <= i__1; ++i__) {
	    if (enclose_(&xyz[1], &vrt[4], &tet[(itet[i__] << 2) + 1], prec)) 
		    {
		ret_val = itet[i__];
		goto L101;
	    }
	}
/* Yuri Vassilevski: instead of searching in a neighborhood, take */
/*                   another cell of the octtree */
/* ...  checking for neighbooring tetrahedra */
/*         do i = ref(idx), ref(idx+1)-1 */
/*           iT = itet(i) */
/*           Do j = 1, 4 */
/*              idx2 = tet(j, iT) */
/*              Do m = ref(idx2), ref(idx2+1)-1 */
/*                 if ( ENCLOSE(XYZ,vrt,tet(1,itet(m)),PREC) ) then */
/*                    BASETET = itet(m) */
/*                    goto 101 */
/*                 endif */
/*              End do */
/*           End do */
/*         enddo */
/* ...  end of */
	++k;
	if (k << 1 > *nbuf) {
	    goto L1000;
	}
	buf[(k << 1) + 1] = ip;
	buf[(k << 1) + 2] = idx;
	qt[ip] = 0;
    }
L101:
    for (n = k; n >= 1; --n) {
	qt[buf[(n << 1) + 1]] = buf[(n << 1) + 2];
    }
    return ret_val;
L1000:
    errmes_(&c__1001, "basetet", "not enough memory for Integer arrays", (
	    ftnlen)7, (ftnlen)36);
    return ret_val;
} /* basetet_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int restore_(integer *nt, integer *tet, integer *
	nv, doublereal *vrt, integer *ldf, doublereal *f, integer *nxyz, 
	doublereal *xyz, doublereal *g, integer *qt, doublereal *xyzc, 
	integer *ref, integer *itet, doublereal *h__, logical *flags, integer 
	*buf, integer *nbuf, integer *nunit)
{
    /* Format strings */
    static char fmt_21[] = "(\002Caution! PREC1 was relaxed: \002,f11.7,\002"
	    " since\002)";
    static char fmt_22[] = "(\002the point \002,3f9.5,\002 is probably out o"
	    "f the mesh\002)";
    static char fmt_23[] = "(\002After MaxTrials failed to find a host withi"
	    "n \002,f8.6,\002 tolerance\002)";
    static char fmt_24[] = "(\002NEED TO REINSTALL constants in lintrp.fd"
	    "\002)";

    /* System generated locals */
    integer f_dim1, f_offset, g_dim1, g_offset, i__1, i__2;

    /* Local variables */
    static doublereal a, b, c__, d__;
    static integer i__, j, ip[4];
    static doublereal det;
    static integer idx, irel;
    static doublereal prec1;

    /* Fortran I/O blocks */
    static cilist io___43 = { 0, 0, 0, fmt_21, 0 };
    static cilist io___44 = { 0, 0, 0, fmt_22, 0 };
    static cilist io___45 = { 0, 0, 0, fmt_23, 0 };
    static cilist io___46 = { 0, 0, 0, fmt_24, 0 };


/* ================================================================ */
/*  Code builds the transposed list of tetrahedra grid and restores */
/*  the values of vector-function F in the points of interest with the aid */
/*  of the octtree structure */
/* ================================================================ */
/* ================================================================ */
/* Predefined tolerance for checking whether the point belongs to an element */
/*      Parameter(RelativeTrials = 1D0) */
/* Predefined number of possible relaxations of the above tolerance (6 orders of 10, here) */
/* ================================================================ */
/* ================================================================ */
/* (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    tet -= 5;
    vrt -= 4;
    f_dim1 = *ldf;
    f_offset = 1 + f_dim1;
    f -= f_offset;
    g_dim1 = *ldf;
    g_offset = 1 + g_dim1;
    g -= g_offset;
    xyz -= 4;
    qt -= 15;
    xyzc -= 4;
    --ref;
    --itet;
    --h__;
    --flags;
    --buf;

    /* Function Body */
    if (flags[1]) {
	i__1 = *nv + 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ref[i__] = 0;
	}
	i__1 = *nt;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (j = 1; j <= 4; ++j) {
		idx = tet[j + (i__ << 2)];
		++ref[idx];
	    }
	}
	++ref[1];
	i__1 = *nv + 1;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    ref[i__] = ref[i__ - 1] + ref[i__];
	}
	i__1 = *nt;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (j = 1; j <= 4; ++j) {
		idx = tet[j + (i__ << 2)];
		--ref[idx];
		itet[ref[idx]] = i__;
	    }
	}
    }
    i__1 = *nxyz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	prec1 = 1e-6;
	irel = 0;
L11:
	idx = basetet_(&qt[15], &xyzc[4], &vrt[4], nt, &tet[5], &ref[1], &
		itet[1], &xyz[i__ * 3 + 1], &h__[1], &buf[1], nbuf, &prec1);
	if (idx <= 0) {
/* if Curvelinear boundaries exist, then probably the node is out of the mesh */
/* Relaxation of PREC is applied */
	    if (flags[2]) {
		if (prec1 < .05) {
		    prec1 *= 10;
		} else {
		    prec1 *= 3.16f;
		}
		if (irel < 7) {
		    ++irel;
		    goto L11;
		} else {
		    errmes_(&c__6103, "lintrp3D", "Max number of relaxations"
			    " for PREC is reached", (ftnlen)8, (ftnlen)45);
		}
/*  if no Curvelinear boundaries, then error is assumed or inconsistency */
/* of parameters in lintrp.fd. REINSTALL constants in lintrp.fd! */
	    } else {
		errmes_(&c__6104, "lintrp3D", "Failed to find element within"
			" PREC tolerance", (ftnlen)8, (ftnlen)44);
	    }
	}
	if (prec1 != 1e-6) {
	    if (*nunit > 0) {
		io___43.ciunit = *nunit;
		s_wsfe(&io___43);
		do_fio(&c__1, (char *)&prec1, (ftnlen)sizeof(doublereal));
		e_wsfe();
		io___44.ciunit = *nunit;
		s_wsfe(&io___44);
		do_fio(&c__1, (char *)&xyz[i__ * 3 + 1], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&xyz[i__ * 3 + 2], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&xyz[i__ * 3 + 3], (ftnlen)sizeof(
			doublereal));
		e_wsfe();
		io___45.ciunit = *nunit;
		s_wsfe(&io___45);
		do_fio(&c__1, (char *)&c_b60, (ftnlen)sizeof(doublereal));
		e_wsfe();
		io___46.ciunit = *nunit;
		s_wsfe(&io___46);
		e_wsfe();
	    }
	}
	for (j = 1; j <= 4; ++j) {
	    ip[j - 1] = tet[j + (idx << 2)];
	}
	a = calvol_(&xyz[i__ * 3 + 1], &vrt[ip[1] * 3 + 1], &vrt[ip[2] * 3 + 
		1], &vrt[ip[3] * 3 + 1]);
	b = calvol_(&xyz[i__ * 3 + 1], &vrt[ip[0] * 3 + 1], &vrt[ip[3] * 3 + 
		1], &vrt[ip[2] * 3 + 1]);
	c__ = calvol_(&xyz[i__ * 3 + 1], &vrt[ip[0] * 3 + 1], &vrt[ip[1] * 3 
		+ 1], &vrt[ip[3] * 3 + 1]);
	d__ = calvol_(&xyz[i__ * 3 + 1], &vrt[ip[0] * 3 + 1], &vrt[ip[2] * 3 
		+ 1], &vrt[ip[1] * 3 + 1]);
	det = a + b + c__ + d__;
	i__2 = *ldf;
	for (j = 1; j <= i__2; ++j) {
	    g[j + i__ * g_dim1] = (a * f[j + ip[0] * f_dim1] + b * f[j + ip[1]
		     * f_dim1] + c__ * f[j + ip[2] * f_dim1] + d__ * f[j + ip[
		    3] * f_dim1]) / det;
	}
/* L10: */
    }
    return 0;
} /* restore_ */

/* ================================================================ */
/* @f2h@ */ logical enclose_(doublereal *xyz, doublereal *vrt, integer *tet, 
	doublereal *prec)
{
    /* Initialized data */

    static integer face[12]	/* was [3][4] */ = { 2,3,4,1,4,3,1,2,4,1,3,2 }
	    ;

    /* System generated locals */
    logical ret_val;

    /* Local variables */
    static integer i__;
    static doublereal vlm, frac;

/* ================================================================ */
/*  Determines does the point XYZ belong to the tetrahedra */
/* ================================================================ */
/* ================================================================ */
    /* Parameter adjustments */
    --tet;
    vrt -= 4;
    --xyz;

    /* Function Body */
/* ================================================================ */
    vlm = calvol_(&vrt[(0 + (0 + (1 + tet[1] * 3 << 3))) / 8], &vrt[(0 + (0 + 
	    (1 + tet[(0 + (0 + (face[0] << 2))) / 4] * 3 << 3))) / 8], &vrt[(
	    0 + (0 + (1 + tet[(0 + (0 + (face[1] << 2))) / 4] * 3 << 3))) / 8]
	    , &vrt[(0 + (0 + (1 + tet[(0 + (0 + (face[2] << 2))) / 4] * 3 << 
	    3))) / 8]);
    for (i__ = 1; i__ <= 4; ++i__) {
	frac = calvol_(&xyz[1], &vrt[tet[face[i__ * 3 - 3]] * 3 + 1], &vrt[
		tet[face[i__ * 3 - 2]] * 3 + 1], &vrt[tet[face[i__ * 3 - 1]] *
		 3 + 1]) / vlm;
	if (frac <= -(*prec)) {
	    ret_val = FALSE_;
	    return ret_val;
	}
    }
    ret_val = TRUE_;
    return ret_val;
} /* enclose_ */

/* ================================================================ */
/* @f2h@ */ doublereal sizeqt_(doublereal *point, integer *imem, doublereal *
	dmem)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static integer iqt;

/* ================================================================ */
/* ================================================================ */
/*  The function returns the size of the octtree cell which contains */
/*  the given point with coords point(3). */

/*  PARAMETERS: */
/*     point(3) - user given point inside the unit cube (0,1)^3 */

/*     dmem(*)  - real*8  working memory of LINTRP3D */
/*     imem(*)  - integer working memory of LINTRP3D */

/*  REMARK 1. This function is the wrapping for function SizeHost. */

/* ================================================================ */
/*      Parameter(RelativeTrials = 1D0) */
    /* Parameter adjustments */
    --dmem;
    --imem;
    --point;

    /* Function Body */
    iqt = imem[1];
    ret_val = sizehost_(&imem[iqt], &dmem[33], &point[1], &dmem[1]);
    return ret_val;
} /* sizeqt_ */

/* ================================================================ */
/* @f2h@ */ doublereal sizehost_(integer *qt, doublereal *xyzc, doublereal *
	point, doublereal *h__)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static integer l, ip[32], ord[256]	/* was [8][32] */, ptr;
    static doublereal sqrd[256]	/* was [8][32] */;

/* ================================================================ */
/* ================================================================ */
/*      Parameter(RelativeTrials = 1D0) */
/* ================================================================ */
    /* Parameter adjustments */
    --h__;
    --point;
    xyzc -= 4;
    qt -= 9;

    /* Function Body */
    l = 1;
    ptr = -1;
    ret_val = 1.;
    while(TRUE_) {
	if (ptr < 0) {
	    ++l;
	    if (l > 32) {
		goto L1000;
	    }
	    ip[l - 1] = -ptr;
	    order2_(&xyzc[-ptr * 3 + 1], &h__[l], &point[1], &ord[(l << 3) - 
		    8], &sqrd[(l << 3) - 8]);
	} else {
	    ret_val = h__[l];
	    return ret_val;
	}
	ptr = qt[ord[(l << 3) - 8] + (ip[l - 1] << 3)];
    }
    errmes_(&c__6105, "SizeHost", "host cell was not found", (ftnlen)8, (
	    ftnlen)23);
    return ret_val;
L1000:
    errmes_(&c__1009, "SizeHost", "local parameter MaxH is small", (ftnlen)8, 
	    (ftnlen)29);
    return ret_val;
} /* sizehost_ */

