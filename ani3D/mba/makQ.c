/* makQ.f -- translated by f2c (version 20090411).
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
#include "makM.h"
#include "update.h"

#include "makQ.h"
#include "lapack.h"

/* Table of constant values */

static doublereal c_b2 = .3333;
static integer c__3 = 3;
static integer c__15 = 15;
static integer c__3011 = 3011;
static integer c__64 = 64;
static integer c__6001 = 6001;

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int makq_(integer *nloop, integer *np, integer *
	ne, doublereal *xyp, integer *ipe, integer *nev, integer *iev, 
	integer *nestar, doublereal *hstar, doublereal *hesp, doublereal *
	detg, doublereal *qe)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__, n, i1, i2, i3, i4;
    static doublereal vn, vstar, detavg, hesavg[6];

/* ================================================================ */
/* group (M) */
/* group (Q) */
/* ================================================================ */
/* ============================================================ */
/* Magic numbers are used to adjust the code performance: */

/* magicNumber      - corection to the final number of tetrahedrons */

/* MaxItrSwap       - the maximal permissible number of local */
/*                    iterations in swapR.f */

/* cMax2DAngle      - the cosine of the minimal angle between neighboring */
/*                    triangles living on a curvilinear surface */

/* ***PREC          - under reconstruction at the moment */

/* MinMPIElements   - the minimal number of tetrahedrons per processor */

/* LocalItrFactor   - the relative number of local sweeps with respect */
/*                    to the number of tetrahedrons */

/* MaxEdgeRefine    - refinement of curvilinear boundary edges */

/* MaxCutSearch     - the maximal number of cuts of a search direction */

/* MaxInertiaTensor - the maximal number of inertial tensors */

/* MaxBaskets       - the maximal number of baskets with bad tetrahedra */
/*                    (for the whole grid) */

/* AniRatio         - the maximal ratio of Hessian eigenvalues */
/* AniEigenvalue    - the eigenvalue used when H is zero matrix */

/* QualityFactor    - increase in the quality (set to 0.1%) */
/* ============================================================ */
/* ================================================================ */
/* Quality computation for mesh elements. */

/* Pre-conditions:  1. connectivity structure {IPE(4, *), XYP(3, *)} */
/*                  2. tensor metric field HesP(6, *) */

/* Post-conditions: 1. determinant of the tensor metric, detG(*) */
/*                  2. quality of mesh elements, qE(*) */

/* Remark: nLoop is used in parallel codes. */
/*         It should be equal to 1 in serial nodal/analytic codes. */
/*         It should be equal to 0 in serial fixquality code. */
/* ================================================================ */
/* group (M) */
/* 2nd order quadrature with 4 points */
/* group (Q) */
/* group (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --qe;
    --detg;
    hesp -= 7;
    --iev;
    ipe -= 5;
    xyp -= 4;

    /* Function Body */
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	caldet_(&hesp[n * 6 + 1], &detg[n]);
    }
    if (*nloop == 1) {
	vstar = 0.;
	i__1 = *ne;
	for (n = 1; n <= i__1; ++n) {
	    i1 = ipe[(n << 2) + 1];
	    i2 = ipe[(n << 2) + 2];
	    i3 = ipe[(n << 2) + 3];
	    i4 = ipe[(n << 2) + 4];
/*  ...  1-point quadrature */
	    for (i__ = 1; i__ <= 6; ++i__) {
		hesavg[i__ - 1] = (hesp[i__ + i1 * 6] + hesp[i__ + i2 * 6] + 
			hesp[i__ + i3 * 6] + hesp[i__ + i4 * 6]) / 4;
	    }
	    caldet_(hesavg, &detavg);
	    vn = calvol_(&xyp[i1 * 3 + 1], &xyp[i2 * 3 + 1], &xyp[i3 * 3 + 1],
		     &xyp[i4 * 3 + 1]);
	    vstar += dabs(vn) * sqrt((dabs(detavg)));
/*  ...  4-point quadrature */
/*           d1 = detG(i1) */
/*           d2 = detG(i2) */
/*           d3 = detG(i3) */
/*           d4 = detG(i4) */

/*           dsum = dsqrt(d1 * T2B + (d2 + d3 + d4) * T2A) */
/*    &           + dsqrt(d2 * T2B + (d3 + d4 + d1) * T2A) */
/*    &           + dsqrt(d3 * T2B + (d4 + d1 + d2) * T2A) */
/*    &           + dsqrt(d4 * T2B + (d1 + d2 + d3) * T2A) */
/*           Vn = calVol(XYP(1, i1), XYP(1, i2), XYP(1, i3), XYP(1, i4)) */
/*           VStar = VStar + dabs(Vn) * dsum / 4 */
	}
	d__1 = vstar / *nestar * 1. * 12. / sqrt(2.);
	*hstar = pow_dd(&d__1, &c_b2);
    }
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	i1 = ipe[(n << 2) + 1];
	i2 = ipe[(n << 2) + 2];
	i3 = ipe[(n << 2) + 3];
	i4 = ipe[(n << 2) + 4];
	calqe_(&hesp[i1 * 6 + 1], &xyp[i1 * 3 + 1], &hesp[i2 * 6 + 1], &xyp[
		i2 * 3 + 1], &hesp[i3 * 6 + 1], &xyp[i3 * 3 + 1], &hesp[i4 * 
		6 + 1], &xyp[i4 * 3 + 1], hstar, &qe[n], &vn);
    }
/* ... set quality of fixed elements to 1 */
    i__1 = *nev;
    for (n = 1; n <= i__1; ++n) {
	qe[iev[n]] = 1.;
    }
    return 0;
} /* makq_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int updqa_(integer *n, doublereal *xyp, integer *
	ipe, integer *iee, doublereal *qe)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static integer i1, i2, i3, j1, j2, j3, i4, j4;
    static doublereal v1, v2;
    static integer ie, ip[5], ip1, ip2, ip3, jp1, jp2, jp3, ip4, jp4;

/* ================================================================ */
/* Initial quality modification for tangled elements and */
/* their closest (face-) neighboors. */
/* ================================================================ */
/* (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --qe;
    iee -= 5;
    ipe -= 6;
    xyp -= 4;

    /* Function Body */
    ip[0] = 1;
    ip[1] = 2;
    ip[2] = 3;
    ip[3] = 4;
    ip[4] = 1;
    for (i1 = 1; i1 <= 4; ++i1) {
	ie = iee[i1 + (*n << 2)];
	if (ie <= 0) {
	    goto L100;
	}
	i2 = ip[i1];
	i3 = ip[i2];
	ip1 = ipe[i1 + *n * 5];
	ip2 = ipe[i2 + *n * 5];
	ip3 = ipe[i3 + *n * 5];
	for (j1 = 1; j1 <= 4; ++j1) {
	    j2 = ip[j1];
	    j3 = ip[j2];
	    jp1 = ipe[j1 + ie * 5];
	    jp2 = ipe[j2 + ie * 5];
	    jp3 = ipe[j3 + ie * 5];
	    if (check33_(&ip1, &ip2, &ip3, &jp1, &jp2, &jp3)) {
		i4 = ip[i3];
		ip4 = ipe[i4 + *n * 5];
		j4 = ip[j3];
		jp4 = ipe[j4 + ie * 5];
		v1 = calvol_(&xyp[ip1 * 3 + 1], &xyp[ip2 * 3 + 1], &xyp[ip3 * 
			3 + 1], &xyp[ip4 * 3 + 1]);
		v2 = calvol_(&xyp[ip1 * 3 + 1], &xyp[ip2 * 3 + 1], &xyp[ip3 * 
			3 + 1], &xyp[jp4 * 3 + 1]);
		if (v1 * v2 >= 0.) {
		    qe[*n] = -(d__1 = qe[*n], dabs(d__1));
		    qe[ie] = -(d__1 = qe[ie], dabs(d__1));
		}
		goto L100;
	    }
	}
L100:
	;
    }
    return 0;
} /* updqa_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int updqb_(integer *nes, integer *le, integer *
	ies, doublereal *xyp, integer *ipes, doublereal *qes)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer k, i1, i2, i3, j1, j2, j3, i4, j4;
    static doublereal v1, v2;
    static integer ip[5], ip1, ip2, ip3, jp1, jp2, jp3, ip4, jp4;

/* ================================================================ */
/* Dynamic quality modification for tangled elements inside */
/* a super-element. */

/* Remark: non-efficient, time-consuming, but robust algorithm. */
/* ================================================================ */
/* group (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --qes;
    ipes -= 6;
    xyp -= 4;
    --ies;

    /* Function Body */
    ip[0] = 1;
    ip[1] = 2;
    ip[2] = 3;
    ip[3] = 4;
    ip[4] = 1;
    for (i1 = 1; i1 <= 4; ++i1) {
	i2 = ip[i1];
	i3 = ip[i2];
	ip1 = ipes[i1 + *nes * 5];
	ip2 = ipes[i2 + *nes * 5];
	ip3 = ipes[i3 + *nes * 5];
	i__1 = *le;
	for (k = 1; k <= i__1; ++k) {
	    if (ies[k] < 0) {
		goto L20;
	    }
	    if (k == *nes) {
		goto L20;
	    }
	    for (j1 = 1; j1 <= 4; ++j1) {
		j2 = ip[j1];
		j3 = ip[j2];
		jp1 = ipes[j1 + k * 5];
		jp2 = ipes[j2 + k * 5];
		jp3 = ipes[j3 + k * 5];
		if (check33_(&ip1, &ip2, &ip3, &jp1, &jp2, &jp3)) {
		    i4 = ip[i3];
		    ip4 = ipes[i4 + *nes * 5];
		    j4 = ip[j3];
		    jp4 = ipes[j4 + k * 5];
		    v1 = calvol_(&xyp[ip1 * 3 + 1], &xyp[ip2 * 3 + 1], &xyp[
			    ip3 * 3 + 1], &xyp[ip4 * 3 + 1]);
		    v2 = calvol_(&xyp[ip1 * 3 + 1], &xyp[ip2 * 3 + 1], &xyp[
			    ip3 * 3 + 1], &xyp[jp4 * 3 + 1]);
		    if (v1 * v2 >= 0.) {
			qes[*nes] = -(d__1 = qes[*nes], dabs(d__1));
			qes[k] = -(d__1 = qes[k], dabs(d__1));
		    }
		    goto L100;
		}
	    }
L20:
	    ;
	}
L100:
	;
    }
    return 0;
} /* updqb_ */

/* ================================================================ */
/* @f2h@ */ doublereal calvol_(doublereal *xy1, doublereal *xy2, doublereal *
	xy3, doublereal *xy4)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static integer i__;
    static doublereal v1[3], v2[3], v3[3];

/* ================================================================ */
/* The oriented volume of the tetrahedron given by 4 vertices */
/* ================================================================ */
/* group (Local variables) */
    /* Parameter adjustments */
    --xy4;
    --xy3;
    --xy2;
    --xy1;

    /* Function Body */
    for (i__ = 1; i__ <= 3; ++i__) {
	v1[i__ - 1] = xy1[i__] - xy4[i__];
	v2[i__ - 1] = xy2[i__] - xy4[i__];
	v3[i__ - 1] = xy3[i__] - xy4[i__];
    }
    ret_val = (v1[0] * (v2[1] * v3[2] - v2[2] * v3[1]) + v1[1] * (v2[2] * v3[
	    0] - v2[0] * v3[2]) + v1[2] * (v2[0] * v3[1] - v2[1] * v3[0])) / 
	    6.;
    return ret_val;
} /* calvol_ */

/* ================================================================ */
/* @f2h@ */ doublereal mutualorientation_(doublereal *xy1, doublereal *xy2, 
	doublereal *xy3, doublereal *xy4, doublereal *xy5)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static integer i__;
    static doublereal v2[3], v3[3], v4[3], v5[3], ax, ay, az, vol1, vol2;

/* ================================================================ */
/* Mutual orientation of points p4 an p5 with respect to plane p1-p2-p3. */
/* ================================================================ */
/* group (Local variables) */
    /* Parameter adjustments */
    --xy5;
    --xy4;
    --xy3;
    --xy2;
    --xy1;

    /* Function Body */
    for (i__ = 1; i__ <= 3; ++i__) {
	v2[i__ - 1] = xy2[i__] - xy1[i__];
	v3[i__ - 1] = xy3[i__] - xy1[i__];
	v4[i__ - 1] = xy4[i__] - xy1[i__];
	v5[i__ - 1] = xy5[i__] - xy1[i__];
    }
    ax = v2[1] * v3[2] - v2[2] * v3[1];
    ay = v2[2] * v3[0] - v2[0] * v3[2];
    az = v2[0] * v3[1] - v2[1] * v3[0];
    vol1 = ax * v4[0] + ay * v4[1] + az * v4[2];
    vol2 = ax * v5[0] + ay * v5[1] + az * v5[2];
    ret_val = vol1 * vol2;
    return ret_val;
} /* mutualorientation_ */

/* ================================================================ */
/* @f2h@ */ doublereal calsqr_(doublereal *xy1, doublereal *xy2, doublereal *
	xy3)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3;

    /* Local variables */
    static doublereal ax, ay, az, bx, by, bz;

/* ================================================================ */
/*  Routine calculates area of a triangle formed by three points */
/* ================================================================ */
/* Local variables */
    /* Parameter adjustments */
    --xy3;
    --xy2;
    --xy1;

    /* Function Body */
    ax = xy1[1] - xy3[1];
    ay = xy1[2] - xy3[2];
    az = xy1[3] - xy3[3];
    bx = xy2[1] - xy3[1];
    by = xy2[2] - xy3[2];
    bz = xy2[3] - xy3[3];
/* Computing 2nd power */
    d__1 = ay * bz - az * by;
/* Computing 2nd power */
    d__2 = ax * bz - az * bx;
/* Computing 2nd power */
    d__3 = ax * by - ay * bx;
    ret_val = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3) * .5;
    return ret_val;
} /* calsqr_ */

/* ====================================================================== */
/* @f2h@ */ /* Subroutine */ int caldet_(doublereal *hesp, doublereal *detg)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */

/* ====================================================================== */
/* Routine computes the determinant of the metric HesP. */

/* Remark: robustness of the overall code was increased by */
/*         replacing errMes() with the computation of |H|. */
/* ====================================================================== */
/* ====================================================================== */
    /* Parameter adjustments */
    --hesp;

    /* Function Body */
/* Computing 2nd power */
    d__1 = hesp[5];
    *detg = hesp[1] * (hesp[2] * hesp[3] - d__1 * d__1) - hesp[4] * (hesp[4] *
	     hesp[3] - hesp[5] * hesp[6]) + hesp[6] * (hesp[4] * hesp[5] - 
	    hesp[6] * hesp[2]);
    if (*detg <= 0.) {
	spectralmodule_(&hesp[1], detg);
    }
    return 0;
} /* caldet_ */

/* ====================================================================== */
/* @f2h@ */ /* Subroutine */ int spectralmodule_(doublereal *hesp, doublereal 
	*detg)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal a[9]	/* was [3][3] */, e[3];
    static integer i__;
    static doublereal rw[15];
    static integer info;

/* ====================================================================== */
/* ====================================================================== */
/* Routines sets the minimal eigenvalue of HesP to some constant */
/* ====================================================================== */
/* ============================================================ */
/* Magic numbers are used to adjust the code performance: */

/* magicNumber      - corection to the final number of tetrahedrons */

/* MaxItrSwap       - the maximal permissible number of local */
/*                    iterations in swapR.f */

/* cMax2DAngle      - the cosine of the minimal angle between neighboring */
/*                    triangles living on a curvilinear surface */

/* ***PREC          - under reconstruction at the moment */

/* MinMPIElements   - the minimal number of tetrahedrons per processor */

/* LocalItrFactor   - the relative number of local sweeps with respect */
/*                    to the number of tetrahedrons */

/* MaxEdgeRefine    - refinement of curvilinear boundary edges */

/* MaxCutSearch     - the maximal number of cuts of a search direction */

/* MaxInertiaTensor - the maximal number of inertial tensors */

/* MaxBaskets       - the maximal number of baskets with bad tetrahedra */
/*                    (for the whole grid) */

/* AniRatio         - the maximal ratio of Hessian eigenvalues */
/* AniEigenvalue    - the eigenvalue used when H is zero matrix */

/* QualityFactor    - increase in the quality (set to 0.1%) */
/* ============================================================ */
/* Local arrays for LAPACK */
/* ====================================================================== */
    /* Parameter adjustments */
    --hesp;

    /* Function Body */
    a[0] = hesp[1];
    a[4] = hesp[2];
    a[8] = hesp[3];
    a[3] = hesp[4];
    a[7] = hesp[5];
    a[6] = hesp[6];
    dsyev_("V", "U", &c__3, a, &c__3, e, rw, &c__15, &info, (ftnlen)1, (
	    ftnlen)1);
    if (info != 0) {
	errmes_(&c__3011, "calDet", "Error in LAPACK routine dsyev", (ftnlen)
		6, (ftnlen)29);
    }
    e[0] = dabs(e[0]);
    e[1] = dabs(e[1]);
    e[2] = dabs(e[2]);
/* Computing MAX */
    d__1 = e[0], d__2 = e[2] * 1e-8;
    e[0] = max(d__1,d__2);
    if (e[0] == 0.) {
	for (i__ = 1; i__ <= 3; ++i__) {
	    e[i__ - 1] = 1e-4;
	}
    }
    for (i__ = 1; i__ <= 6; ++i__) {
	hesp[i__] = 0.;
    }
    for (i__ = 1; i__ <= 3; ++i__) {
/* Computing 2nd power */
	d__1 = a[i__ * 3 - 3];
	hesp[1] += e[i__ - 1] * (d__1 * d__1);
/* Computing 2nd power */
	d__1 = a[i__ * 3 - 2];
	hesp[2] += e[i__ - 1] * (d__1 * d__1);
/* Computing 2nd power */
	d__1 = a[i__ * 3 - 1];
	hesp[3] += e[i__ - 1] * (d__1 * d__1);
	hesp[4] += e[i__ - 1] * a[i__ * 3 - 3] * a[i__ * 3 - 2];
	hesp[5] += e[i__ - 1] * a[i__ * 3 - 2] * a[i__ * 3 - 1];
	hesp[6] += e[i__ - 1] * a[i__ * 3 - 3] * a[i__ * 3 - 1];
    }
    *detg = e[0] * e[1] * e[2];
    return 0;
} /* spectralmodule_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int calqe_(doublereal *hes1, doublereal *xy1, 
	doublereal *hes2, doublereal *xy2, doublereal *hes3, doublereal *xy3, 
	doublereal *hes4, doublereal *xy4, doublereal *hstar, doublereal *qe, 
	doublereal *volume)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal f;
    static integer i__, n;
    static doublereal x1, y1, z1, pk, vk, detavg, hesavg[6];

/* ================================================================ */
/* Computing quality of tetrahedron {xy1, ..., xy4} assuming that */
/* the metric field is linear. */
/* ================================================================ */
/* ================================================================ */
/* 2nd order quadrature with 4 points */
/* ================================================================ */
/*  local variables */
/* ================================================================ */
    /* Parameter adjustments */
    --xy4;
    --hes4;
    --xy3;
    --hes3;
    --xy2;
    --hes2;
    --xy1;
    --hes1;

    /* Function Body */
    pk = 0.;
    for (n = 1; n <= 6; ++n) {
	if (n == 1) {
	    x1 = xy1[1] - xy4[1];
	    y1 = xy1[2] - xy4[2];
	    z1 = xy1[3] - xy4[3];
	    for (i__ = 1; i__ <= 6; ++i__) {
		hesavg[i__ - 1] = (hes1[i__] + hes4[i__]) / 2;
	    }
	} else if (n == 2) {
	    x1 = xy2[1] - xy4[1];
	    y1 = xy2[2] - xy4[2];
	    z1 = xy2[3] - xy4[3];
	    for (i__ = 1; i__ <= 6; ++i__) {
		hesavg[i__ - 1] = (hes2[i__] + hes4[i__]) / 2;
	    }
	} else if (n == 3) {
	    x1 = xy3[1] - xy4[1];
	    y1 = xy3[2] - xy4[2];
	    z1 = xy3[3] - xy4[3];
	    for (i__ = 1; i__ <= 6; ++i__) {
		hesavg[i__ - 1] = (hes3[i__] + hes4[i__]) / 2;
	    }
	} else if (n == 4) {
	    x1 = xy1[1] - xy3[1];
	    y1 = xy1[2] - xy3[2];
	    z1 = xy1[3] - xy3[3];
	    for (i__ = 1; i__ <= 6; ++i__) {
		hesavg[i__ - 1] = (hes1[i__] + hes3[i__]) / 2;
	    }
	} else if (n == 5) {
	    x1 = xy2[1] - xy3[1];
	    y1 = xy2[2] - xy3[2];
	    z1 = xy2[3] - xy3[3];
	    for (i__ = 1; i__ <= 6; ++i__) {
		hesavg[i__ - 1] = (hes2[i__] + hes3[i__]) / 2;
	    }
	} else if (n == 6) {
	    x1 = xy1[1] - xy2[1];
	    y1 = xy1[2] - xy2[2];
	    z1 = xy1[3] - xy2[3];
	    for (i__ = 1; i__ <= 6; ++i__) {
		hesavg[i__ - 1] = (hes1[i__] + hes2[i__]) / 2;
	    }
	}
	pk += sqrt(hesavg[0] * x1 * x1 + hesavg[1] * y1 * y1 + hesavg[2] * z1 
		* z1 + hesavg[3] * 2 * x1 * y1 + hesavg[4] * 2 * y1 * z1 + 
		hesavg[5] * 2 * x1 * z1);
    }
/* ... 1-point quadrature rule */
    for (i__ = 1; i__ <= 6; ++i__) {
	hesavg[i__ - 1] = (hes1[i__] + hes2[i__] + hes3[i__] + hes4[i__]) / 4;
    }
    caldet_(hesavg, &detavg);
    *volume = calvol_(&xy1[1], &xy2[1], &xy3[1], &xy4[1]);
    vk = *volume * sqrt((dabs(detavg)));
/* ... 4-point quadrature rule */
/*     d1 = det1 */
/*     d2 = det2 */
/*     d3 = det3 */
/*     d4 = det4 */

/*     dsum = dsqrt(d1 * T2B + (d2 + d3 + d4) * T2A) */
/*    &     + dsqrt(d2 * T2B + (d3 + d4 + d1) * T2A) */
/*    &     + dsqrt(d3 * T2B + (d4 + d1 + d2) * T2A) */
/*    &     + dsqrt(d4 * T2B + (d1 + d2 + d3) * T2A) */
/*     volume = calVol(xy1, xy2, xy3, xy4) */
/*     Vk = volume * dsum / 4 */
/* Computing 3rd power */
    d__1 = pk;
    *qe = dabs(vk) * 1832.8208 / (d__1 * (d__1 * d__1));
    if (*hstar > 0.) {
	x1 = pk / (*hstar * 6);
/* Computing MIN */
	d__1 = x1, d__2 = 1. / x1;
	x1 = min(d__1,d__2);
/* Computing 5th power */
	d__1 = x1 * (2. - x1), d__2 = d__1, d__1 *= d__1;
	f = d__2 * (d__1 * d__1);
	*qe *= f;
    }
    return 0;
} /* calqe_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int calqf_(doublereal *hes1, doublereal *xy1, 
	doublereal *hes2, doublereal *xy2, doublereal *hes3, doublereal *xy3, 
	doublereal *hes4, doublereal *xy4, doublereal *hstar, integer *if__, 
	integer *ir)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k, n;
    static doublereal x1, y1, z1, qf[4], qr[6];
    static integer imin;
    static doublereal hesavg[6];

/* ================================================================ */
/* group (Local variables) */
/* group (Function) */
/* ================================================================ */
    /* Parameter adjustments */
    --ir;
    --if__;
    --xy4;
    --hes4;
    --xy3;
    --hes3;
    --xy2;
    --hes2;
    --xy1;
    --hes1;

    /* Function Body */
    for (n = 1; n <= 6; ++n) {
	if (n == 1) {
	    x1 = xy1[1] - xy2[1];
	    y1 = xy1[2] - xy2[2];
	    z1 = xy1[3] - xy2[3];
	    for (i__ = 1; i__ <= 6; ++i__) {
		hesavg[i__ - 1] = (hes1[i__] + hes2[i__]) / 2;
	    }
	} else if (n == 2) {
	    x1 = xy1[1] - xy3[1];
	    y1 = xy1[2] - xy3[2];
	    z1 = xy1[3] - xy3[3];
	    for (i__ = 1; i__ <= 6; ++i__) {
		hesavg[i__ - 1] = (hes1[i__] + hes3[i__]) / 2;
	    }
	} else if (n == 3) {
	    x1 = xy1[1] - xy4[1];
	    y1 = xy1[2] - xy4[2];
	    z1 = xy1[3] - xy4[3];
	    for (i__ = 1; i__ <= 6; ++i__) {
		hesavg[i__ - 1] = (hes1[i__] + hes4[i__]) / 2;
	    }
	} else if (n == 4) {
	    x1 = xy2[1] - xy3[1];
	    y1 = xy2[2] - xy3[2];
	    z1 = xy2[3] - xy3[3];
	    for (i__ = 1; i__ <= 6; ++i__) {
		hesavg[i__ - 1] = (hes2[i__] + hes3[i__]) / 2;
	    }
	} else if (n == 5) {
	    x1 = xy2[1] - xy4[1];
	    y1 = xy2[2] - xy4[2];
	    z1 = xy2[3] - xy4[3];
	    for (i__ = 1; i__ <= 6; ++i__) {
		hesavg[i__ - 1] = (hes2[i__] + hes4[i__]) / 2;
	    }
	} else if (n == 6) {
	    x1 = xy3[1] - xy4[1];
	    y1 = xy3[2] - xy4[2];
	    z1 = xy3[3] - xy4[3];
	    for (i__ = 1; i__ <= 6; ++i__) {
		hesavg[i__ - 1] = (hes3[i__] + hes4[i__]) / 2;
	    }
	}
	x1 = sqrt(hesavg[0] * x1 * x1 + hesavg[1] * y1 * y1 + hesavg[2] * z1 *
		 z1 + hesavg[3] * 2 * x1 * y1 + hesavg[4] * 2 * y1 * z1 + 
		hesavg[5] * 2 * x1 * z1) / *hstar;
/* Computing MIN */
	d__1 = x1, d__2 = 1. / x1;
	x1 = min(d__1,d__2);
	ir[n] = n;
/* Computing 5th power */
	d__1 = x1 * (2. - x1), d__2 = d__1, d__1 *= d__1;
	qr[n - 1] = d__2 * (d__1 * d__1);
    }
    qf[0] = min(qr[0],qr[1]);
    qf[0] = min(qf[0],qr[3]);
    qf[1] = min(qr[3],qr[4]);
    qf[1] = min(qf[1],qr[5]);
    qf[2] = min(qr[1],qr[2]);
    qf[2] = min(qf[2],qr[5]);
    qf[3] = min(qr[0],qr[2]);
    qf[3] = min(qf[3],qr[4]);
    for (i__ = 1; i__ <= 4; ++i__) {
	if__[i__] = i__;
    }
    for (i__ = 1; i__ <= 3; ++i__) {
	imin = i__;
	for (j = i__ + 1; j <= 4; ++j) {
	    if (qf[j - 1] < qf[imin - 1]) {
		imin = j;
	    }
	}
	x1 = qf[i__ - 1];
	qf[i__ - 1] = qf[imin - 1];
	qf[imin - 1] = x1;
	k = if__[i__];
	if__[i__] = if__[imin];
	if__[imin] = k;
    }
    for (i__ = 1; i__ <= 5; ++i__) {
	imin = i__;
	for (j = i__ + 1; j <= 6; ++j) {
	    if (qr[j - 1] < qr[imin - 1]) {
		imin = j;
	    }
	}
	x1 = qr[i__ - 1];
	qr[i__ - 1] = qr[imin - 1];
	qr[imin - 1] = x1;
	k = ir[i__];
	ir[i__] = ir[imin];
	ir[imin] = k;
    }
    return 0;
} /* calqf_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int hesbnd_(integer *lp, integer *ips, integer *
	icp, doublereal *hesp, doublereal *hesps)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, n, nb, ni, ip;
    static doublereal hesb[6], hesi[6];

/* ================================================================ */
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
/* (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --hesps;
    hesp -= 7;
    --icp;
    --ips;

    /* Function Body */
    ni = 0;
    nb = 0;
    for (i__ = 1; i__ <= 6; ++i__) {
	hesi[i__ - 1] = 0.;
	hesb[i__ - 1] = 0.;
    }
    i__1 = *lp;
    for (n = 1; n <= i__1; ++n) {
	ip = ips[n];
	if (ifxnode_(&icp[ip], &c__64)) {
	    ++ni;
	    for (i__ = 1; i__ <= 6; ++i__) {
		hesi[i__ - 1] += hesp[i__ + ip * 6];
	    }
	} else {
	    ++nb;
	    for (i__ = 1; i__ <= 6; ++i__) {
		hesb[i__ - 1] += hesp[i__ + ip * 6];
	    }
	}
    }
    if (ni > 0) {
	for (i__ = 1; i__ <= 6; ++i__) {
	    hesps[i__] = hesi[i__ - 1] / ni;
	}
    } else if (nb > 0) {
	for (i__ = 1; i__ <= 6; ++i__) {
	    hesps[i__] = hesb[i__ - 1] / nb;
	}
    } else {
	errmes_(&c__6001, "HesBnd", "system error", (ftnlen)6, (ftnlen)12);
    }
    return 0;
} /* hesbnd_ */

/* ====================================================================== */
/* @f2h@ */ /* Subroutine */ int iniq_analytic__(integer *np, doublereal *xyp,
	 I_fp metricfunction, doublereal *hesp)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, n;
    static doublereal x, y, z__, metric[9]	/* was [3][3] */;

/* ====================================================================== */
/*  Three Fortran routines below create a metric field which */
/*  is 3x3 variable positive definite diagonal tensor HesP, */

/*             M11   M12   M13 */
/*      HesP = M12   M22   M23 */
/*             M13   M23   M33 */

/*  where Mij = Mij(x, y, z). */

/*  The tensor element are enumerated in the following order: */
/*  HesP_{11}, HesP_{22}, HesP_{33}, HesP_{12}, HesP_{23}, HesP_{13} */

/* ====================================================================== */
/* type of the metric: isotropic/scalar or full tensor */
/* ===================================================================== */
    /* Parameter adjustments */
    hesp -= 7;
    xyp -= 4;

    /* Function Body */
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	x = xyp[n * 3 + 1];
	y = xyp[n * 3 + 2];
	z__ = xyp[n * 3 + 3];
	i__ = (*metricfunction)(&x, &y, &z__, metric);
	hesp[n * 6 + 1] = metric[0];
	hesp[n * 6 + 2] = metric[4];
	hesp[n * 6 + 3] = metric[8];
	hesp[n * 6 + 4] = metric[3];
	hesp[n * 6 + 5] = metric[7];
	hesp[n * 6 + 6] = metric[6];
    }
    return 0;
} /* iniq_analytic__ */

/* ================================================================ */
/* @f2h@ */ doublereal avgq_(integer *ne, doublereal *qe, integer *l1e, 
	integer *l2e)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer n, ie;

/* ================================================================ */
    /* Parameter adjustments */
    --l2e;
    l1e -= 3;
    --qe;

    /* Function Body */
    ret_val = 0.;
    ie = l2e[1];
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	ret_val += qe[ie];
	ie = l1e[(ie << 1) + 2];
    }
    ret_val /= *ne;
    return ret_val;
} /* avgq_ */

