/* auxSF.f -- translated by f2c (version 20090411).
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

static integer c__6001 = 6001;
static integer c__1007 = 1007;

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int maksf_(integer *ip, integer *lf, integer *
	ifs, integer *le, integer *ies, integer *ipes, integer *ipf, integer *
	ife, integer *maxs, integer *lp1, integer *ip1s, integer *lp2, 
	integer *ip2s, integer *icp2s, integer *ls, integer *ipss, integer *
	iess, integer *ldf, integer *idfs, integer *lde, integer *ides, 
	logical *flagface, logical *flagedge)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, k, m, n, i1, i2, i3, ift, nft, ipt[3], mpt, npt;
    static logical flag__;
    static integer iref[5], iclrs;

/* ================================================================ */
/* ================================================================ */
/* Compute: */
/*   (a) list of faces IPSs(3, lS) around point iP is created. */

/*   (b) list of surface points connected with point iP */

/*   (c) edge case : list of points on an edge passing through iP */
/*       face case : list of points on the face passing through iP */
/* ================================================================ */
/* group (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --ides;
    --idfs;
    --iess;
    ipss -= 4;
    --icp2s;
    --ip2s;
    --ip1s;
    ife -= 5;
    ipf -= 5;
    ipes -= 6;
    --ies;
    --ifs;

    /* Function Body */
    iref[0] = 1;
    iref[1] = 2;
    iref[2] = 3;
    iref[3] = 4;
    iref[4] = 1;
    *lp1 = 0;
    *lp2 = 0;
    *ls = 0;
    *ldf = 0;
    *lde = 0;
    *flagface = FALSE_;
    *flagedge = FALSE_;
    i__1 = *le;
    for (n = 1; n <= i__1; ++n) {
	flag__ = FALSE_;
	for (i1 = 1; i1 <= 4; ++i1) {
	    if (ipes[i1 + n * 5] == *ip) {
		flag__ = TRUE_;
		++(*lde);
		ides[*lde] = n;
		i2 = i1;
		for (i__ = 1; i__ <= 3; ++i__) {
		    i2 = iref[i2];
		    ipt[i__ - 1] = ipes[i2 + n * 5];
		    findse_(lp2, &ip2s[1], &ipt[i__ - 1], &npt);
		    if (npt == 0) {
			++(*lp2);
			if (*lp2 > *maxs) {
			    goto L1000;
			}
			ip2s[*lp2] = ipt[i__ - 1];
			icp2s[*lp2] = 0;
		    }
		}
		i__2 = *ls;
		for (m = 1; m <= i__2; ++m) {
		    if (check33_(ipt, &ipt[1], &ipt[2], &ipss[m * 3 + 1], &
			    ipss[m * 3 + 2], &ipss[m * 3 + 3])) {
			goto L10;
		    }
		}
		++(*ls);
		if (*ls > *maxs) {
		    goto L1000;
		}
		for (i__ = 1; i__ <= 3; ++i__) {
		    ipss[i__ + *ls * 3] = ipt[i__ - 1];
		}
		iess[*ls] = n;
		goto L15;
	    }
L10:
	    ;
	}
L15:
	if (flag__) {
	    for (i1 = 1; i1 <= 4; ++i1) {
		ift = ife[i1 + (ies[n] << 2)];
		if (ift != 0) {
		    i2 = iref[i1];
		    i3 = iref[i2];
		    ipt[0] = ipes[i1 + n * 5];
		    ipt[1] = ipes[i2 + n * 5];
		    ipt[2] = ipes[i3 + n * 5];
		    if (check13_(ip, ipt, &ipt[1], &ipt[2])) {
			iclrs = ipf[(ift << 2) + 4];
			if (iclrs < 0) {
			    errmes_(&c__6001, "auxSF", "system error", (
				    ftnlen)5, (ftnlen)12);
			}
			for (i__ = 1; i__ <= 3; ++i__) {
			    findse_(lp2, &ip2s[1], &ipt[i__ - 1], &npt);
			    if (npt > 0) {
				if (icp2s[npt] == 0) {
				    *flagface = TRUE_;
				    icp2s[npt] = iclrs;
				} else if (icp2s[npt] != iclrs) {
				    *flagedge = TRUE_;
				    findse_(lp1, &ip1s[1], &ipt[i__ - 1], &
					    mpt);
				    if (mpt == 0) {
					++(*lp1);
					if (*lp1 > *maxs) {
					    goto L1000;
					}
					ip1s[*lp1] = ipt[i__ - 1];
				    }
				    icp2s[npt] += iclrs;
				}
			    }
			}
			findse_(lf, &ifs[1], &ift, &nft);
			findse_(ldf, &idfs[1], &nft, &k);
			if (k <= 0) {
			    ++(*ldf);
			    if (*ldf > *maxs) {
				goto L1000;
			    }
			    idfs[*ldf] = nft;
			}
		    }
		}
/* L20: */
	    }
	}
    }
    if (*flagedge) {
	goto L9000;
    } else if (*flagface) {
	i__1 = *lp2;
	for (n = 1; n <= i__1; ++n) {
	    if (icp2s[n] != 0) {
		++(*lp1);
		ip1s[*lp1] = ip2s[n];
	    }
	}
    } else {
	*lp1 = *lp2;
	i__1 = *lp1;
	for (n = 1; n <= i__1; ++n) {
	    ip1s[n] = ip2s[n];
	}
    }
L9000:
    return 0;
L1000:
    errmes_(&c__1007, "makSF", "local variable MaxS is small", (ftnlen)5, (
	    ftnlen)28);
    return 0;
} /* maksf_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int shutf_(doublereal *zz1, doublereal *zz2, 
	doublereal *xyp1, doublereal *xyp2, doublereal *xyp3, logical *flag__)
{
    static logical flagopen;
    static doublereal v1, v2;

/* ================================================================ */
/* Routine returns .TRUE. if the ray [ZZ1, ZZ2) intersects */
/* triangle {xyz1, xyz2, xyz3}. If flag = .TRUE. the triangle is */
/* assumed to be the open domain. */
/* ================================================================ */
/* group (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --xyp3;
    --xyp2;
    --xyp1;
    --zz2;
    --zz1;

    /* Function Body */
    flagopen = *flag__;
    *flag__ = FALSE_;
    v1 = calvol_(&zz1[1], &zz2[1], &xyp3[1], &xyp1[1]);
    v2 = calvol_(&zz1[1], &zz2[1], &xyp3[1], &xyp2[1]);
    if (flagopen) {
	if (v1 * v2 >= 0.) {
	    goto L1000;
	}
    } else {
	if (v1 * v2 > 0.) {
	    goto L1000;
	}
    }
    v1 = calvol_(&zz1[1], &zz2[1], &xyp2[1], &xyp1[1]);
    v2 = calvol_(&zz1[1], &zz2[1], &xyp2[1], &xyp3[1]);
    if (flagopen) {
	if (v1 * v2 >= 0.) {
	    goto L1000;
	}
    } else {
	if (v1 * v2 > 0.) {
	    goto L1000;
	}
    }
    v1 = calvol_(&zz1[1], &zz2[1], &xyp1[1], &xyp2[1]);
    v2 = calvol_(&zz1[1], &zz2[1], &xyp1[1], &xyp3[1]);
    if (flagopen) {
	if (v1 * v2 >= 0.) {
	    goto L1000;
	}
    } else {
	if (v1 * v2 > 0.) {
	    goto L1000;
	}
    }
    *flag__ = TRUE_;
L1000:
    return 0;
} /* shutf_ */

/* ================================================================ */
/* @f2h@ */ doublereal projectf_(doublereal *xyp, doublereal *xy1, doublereal 
	*xy2, doublereal *xy3, doublereal *xyt)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Local variables */
    static doublereal a, b;
    static integer i__;
    static doublereal f1, f2, v1[3], v2[3], m11, m12, m22, vp[3], det;

/* ================================================================ */
/* Routine finds a point xyt lying inside triangular face */
/* {xy1, xy2, xy3} which is the closest point to xyp. It returns */
/* the distance to this point. */

/* *** Remarks: */
/*        1. We minimize the distance ||a * x21 + b * x31 - p||^2 */
/*           w.r.t. a and b. */
/* ================================================================ */
/* group (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --xyt;
    --xy3;
    --xy2;
    --xy1;
    --xyp;

    /* Function Body */
    for (i__ = 1; i__ <= 3; ++i__) {
	v1[i__ - 1] = xy2[i__] - xy1[i__];
	v2[i__ - 1] = xy3[i__] - xy1[i__];
	vp[i__ - 1] = xyp[i__] - xy1[i__];
    }
    m11 = dotmul_(v1, v1);
    m12 = dotmul_(v1, v2);
    m22 = dotmul_(v2, v2);
    f1 = dotmul_(vp, v1);
    f2 = dotmul_(vp, v2);
/* ... solve a linear system (Wronskian is always non-zero) */
    det = m22 * m11 - m12 * m12;
    a = (m22 * f1 - m12 * f2) / det;
    b = (m11 * f2 - m12 * f1) / det;
/* ... analyze few cases (seach along b, a and a + b = 1) */
    if (a < 0.) {
	a = 0.;
/* Computing MAX */
	d__1 = 0., d__2 = f2 / m22;
	b = max(d__1,d__2);
	b = min(1.,b);
    } else if (b < 0.) {
	b = 0.;
/* Computing MAX */
	d__1 = 0., d__2 = f1 / m11;
	a = max(d__1,d__2);
	a = min(1.,a);
    } else if (a + b > 1.) {
/* Computing MAX */
	d__1 = 0., d__2 = (f1 - m12) / (m11 - m12);
	a = max(d__1,d__2);
	a = min(1.,a);
	b = 1. - a;
    }
/* ... compute point the closest point */
    for (i__ = 1; i__ <= 3; ++i__) {
	xyt[i__] = xy1[i__] + a * v1[i__ - 1] + b * v2[i__ - 1];
    }
    ret_val = caledge_(&xyp[1], &xyt[1]);
    return ret_val;
} /* projectf_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int baricoords_(doublereal *xyp, doublereal *xy1,
	 doublereal *xy2, doublereal *xy3, doublereal *b1, doublereal *b2, 
	doublereal *b3)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static doublereal a[3], b[3], c__[3];
    static integer i__;
    static doublereal ar1, art;

/* ================================================================ */
/* Routine computes barycentric coordinates of point xyp lying */
/* inside triangle {xy1, xy2, xy3}. They coordinates are always */
/* positive and sum to 1 upto round-off errors. */
/* ================================================================ */
/* group (Local variables) */
/* ================================================================ */
/* ... square of area of the triangle */
    /* Parameter adjustments */
    --xy3;
    --xy2;
    --xy1;
    --xyp;

    /* Function Body */
    for (i__ = 1; i__ <= 3; ++i__) {
	a[i__ - 1] = xy2[i__] - xy1[i__];
	b[i__ - 1] = xy3[i__] - xy1[i__];
    }
    vecmul_(a, b, c__);
/* Computing 2nd power */
    d__1 = c__[0];
/* Computing 2nd power */
    d__2 = c__[1];
/* Computing 2nd power */
    d__3 = c__[2];
    art = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
/* ... first baricentric coordinate */
    for (i__ = 1; i__ <= 3; ++i__) {
	a[i__ - 1] = xy2[i__] - xyp[i__];
	b[i__ - 1] = xy3[i__] - xyp[i__];
    }
    vecmul_(a, b, c__);
/* Computing 2nd power */
    d__1 = c__[0];
/* Computing 2nd power */
    d__2 = c__[1];
/* Computing 2nd power */
    d__3 = c__[2];
    ar1 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
    *b1 = sqrt(ar1 / art);
/* ... second baricentric coordinate */
    for (i__ = 1; i__ <= 3; ++i__) {
	a[i__ - 1] = xy1[i__] - xyp[i__];
    }
    vecmul_(a, b, c__);
/* Computing 2nd power */
    d__1 = c__[0];
/* Computing 2nd power */
    d__2 = c__[1];
/* Computing 2nd power */
    d__3 = c__[2];
    ar1 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
    *b2 = sqrt(ar1 / art);
/* ... third baricentric coordinate */
    *b3 = 1 - *b2 - *b1;
    if (*b1 + *b2 > 1.) {
	*b1 = min(1.,*b1);
/* Computing MAX */
	d__1 = 0., d__2 = 1 - *b1;
	*b2 = max(d__1,d__2);
	*b3 = 0.;
    }
    return 0;
} /* baricoords_ */

/* ================================================================ */
/* @f2h@ */ doublereal angle2faces_(doublereal *xya, doublereal *xyb, 
	doublereal *xyc, doublereal *xyd)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3;

    /* Local variables */
    static integer i__;
    static doublereal d1, d2, v1[3], v2[3], v3[3], xym[3], xyn[3];

/* ================================================================ */
/* Cosine of angle between two faces {A,B,C} and {A,B,D}. */

/* Remark: the angle is equal to -1 for flat faces. */
/* ================================================================ */
    /* Parameter adjustments */
    --xyd;
    --xyc;
    --xyb;
    --xya;

    /* Function Body */
    for (i__ = 1; i__ <= 3; ++i__) {
	v1[i__ - 1] = xyb[i__] - xya[i__];
	v2[i__ - 1] = xyc[i__] - xyb[i__];
	v3[i__ - 1] = xyd[i__] - xyb[i__];
    }
    vecmul_(v1, v2, xyn);
    vecmul_(v1, v3, xym);
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
} /* angle2faces_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int trinormal_(doublereal *xy1, doublereal *xy2, 
	doublereal *xy3, doublereal *normal)
{
    static doublereal a[3], b[3];
    static integer i__;

/* ================================================================ */
    /* Parameter adjustments */
    --normal;
    --xy3;
    --xy2;
    --xy1;

    /* Function Body */
    for (i__ = 1; i__ <= 3; ++i__) {
	a[i__ - 1] = xy2[i__] - xy1[i__];
	b[i__ - 1] = xy3[i__] - xy1[i__];
    }
    vecmul_(a, b, &normal[1]);
    return 0;
} /* trinormal_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int vecmul_(doublereal *a, doublereal *b, 
	doublereal *c__)
{
/* ================================================================ */
    /* Parameter adjustments */
    --c__;
    --b;
    --a;

    /* Function Body */
    c__[1] = a[2] * b[3] - b[2] * a[3];
    c__[2] = b[1] * a[3] - a[1] * b[3];
    c__[3] = a[1] * b[2] - b[1] * a[2];
    return 0;
} /* vecmul_ */

/* ================================================================ */
/* @f2h@ */ doublereal dotmul_(doublereal *a, doublereal *b)
{
    /* System generated locals */
    doublereal ret_val;

/* ================================================================ */
    /* Parameter adjustments */
    --b;
    --a;

    /* Function Body */
    ret_val = a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
    return ret_val;
} /* dotmul_ */

