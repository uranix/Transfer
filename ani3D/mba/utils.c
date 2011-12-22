/* utils.f -- translated by f2c (version 20090411).
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

/* Table of constant values */

static integer c__4 = 4;
static integer c__2011 = 2011;
static integer c__6 = 6;
static integer c__3 = 3;
static integer c__0 = 0;
static integer c__1011 = 1011;
static integer c__1004 = 1004;
static integer c__1010 = 1010;
static integer c__1003 = 1003;
static integer c__1006 = 1006;

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int listp2p_(integer *np, integer *ne, integer *
	maxlist, integer *ipe, integer *npp, integer *ipp, integer *iw)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j, m, n, i1, i2, ie, nl, nlo, ipt, iiep, inep;

/* ================================================================ */
/*  The routine creates connectivity lists P->P for mesh points. */

/*  *** Remarks: */
/*         1. iW(*) - working memory of size 4 * nE + nP */
/* ================================================================ */
/* ================================================================ */
    /* Parameter adjustments */
    --iw;
    --ipp;
    --npp;
    ipe -= 5;

    /* Function Body */
    inep = 0;
    iiep = inep + *np;
    backreferences_(np, ne, &c__4, &c__4, &ipe[5], &iw[inep + 1], &iw[iiep + 
	    1]);
/* ... main algorithm: array nEP is overloaded inside */
    nl = 0;
    i2 = 0;
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	nlo = nl;
	i1 = i2 + 1;
	i2 = iw[inep + n];
	i__2 = i2;
	for (m = i1; m <= i__2; ++m) {
	    ie = iw[iiep + m];
	    for (j = 1; j <= 4; ++j) {
		ipt = ipe[j + (ie << 2)];
		if (iw[inep + ipt] > 0) {
		    ++nl;
		    if (nl > *maxlist) {
			errmes_(&c__2011, "listP2P", "user parameter MaxList"
				" is small", (ftnlen)7, (ftnlen)31);
		    }
		    ipp[nl] = ipt;
		    iw[inep + ipt] = -iw[inep + ipt];
		}
	    }
	}
	npp[n] = nl;
/*  ...  recovering values of array nEP */
	i__2 = nl;
	for (m = nlo + 1; m <= i__2; ++m) {
	    ipt = ipp[m];
	    iw[inep + ipt] = -iw[inep + ipt];
	}
    }
    return 0;
} /* listp2p_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int listr2p_(integer *np, integer *nr, integer *
	ne, integer *maxr, integer *ipe, integer *ipr, integer *iw)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j, m, n, i1, i2, ie, ipt, nro, iiep, inep;

/* ================================================================ */
/*  The routine reates connectivity list R -> P for mesh edges. */
/*  The algorithm has a linear arithmetical complexity. */

/*  Parameters: */
/*      nP, nR, nE - the number of points, edges and elements */
/*      MaxR       - the maximal number of edges */

/*      IPE(4, nE) - connectivity list of tetrahedra: E -> P */
/*      IPR(2, nR) - connectivity list of edges:      R -> P */

/*      iW(*) - working memory of size 4 * nE + nP */
/* ================================================================ */
/* ================================================================ */
    /* Parameter adjustments */
    --iw;
    ipr -= 3;
    ipe -= 5;

    /* Function Body */
    inep = 0;
    iiep = inep + *np;
    backreferences_(np, ne, &c__4, &c__4, &ipe[5], &iw[inep + 1], &iw[iiep + 
	    1]);
/* ... main algorithm: array nEP is overloaded inside */
    *nr = 0;
    i2 = 0;
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	nro = *nr;
	i1 = i2 + 1;
	i2 = iw[inep + n];
	i__2 = i2;
	for (m = i1; m <= i__2; ++m) {
	    ie = iw[iiep + m];
	    for (j = 1; j <= 4; ++j) {
		ipt = ipe[j + (ie << 2)];
		if (ipt > n && iw[inep + ipt] > 0) {
		    ++(*nr);
		    if (*nr > *maxr) {
			errmes_(&c__2011, "listR2P", "user parameter MaxR is"
				" small", (ftnlen)7, (ftnlen)28);
		    }
		    ipr[(*nr << 1) + 1] = n;
		    ipr[(*nr << 1) + 2] = ipt;
		    iw[inep + ipt] = -iw[inep + ipt];
		}
	    }
	}
/*  ...  recovering values of array nEP */
	i__2 = *nr;
	for (m = nro + 1; m <= i__2; ++m) {
	    ipt = ipr[(m << 1) + 2];
	    iw[inep + ipt] = -iw[inep + ipt];
	}
    }
    return 0;
} /* listr2p_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int listr2r_(integer *np, integer *nr, integer *
	ne, integer *maxl, integer *ipe, integer *nrr, integer *irr, integer *
	iw, integer *ierr)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer j, k, m, n, i1, i2, ie, nl, nlo, irt, iend, iiep, iire, 
	    iier, inep, iner;

/* ================================================================ */
/*  The routine creates connectivity lists R->R for mesh edges. */
/*  Routine returns 0 upon successful completion. */

/*  *** Remarks: */
/*         1. iW(*) - working memory of size 12 * nE + nR */
/* ================================================================ */
/* ================================================================ */
    /* Parameter adjustments */
    --iw;
    --irr;
    --nrr;
    ipe -= 5;

    /* Function Body */
    *ierr = 0;
    iire = 1;
    inep = iire + *ne * 6;
    iiep = inep + *np;
    iend = iiep + (*ne << 2);
    liste2r_(np, nr, ne, &ipe[5], &iw[iire], &iw[inep], &iw[iiep]);
    iner = inep;
    iier = iner + *nr;
    iend = iier + *ne * 6;
    backreferences_(nr, ne, &c__6, &c__6, &iw[iire], &iw[iner], &iw[iier]);
    nl = 0;
    i2 = 0;
    i__1 = *nr;
    for (n = 1; n <= i__1; ++n) {
	nlo = nl;
	i1 = i2 + 1;
	i2 = iw[iner + n - 1];
	i__2 = i2;
	for (m = i1; m <= i__2; ++m) {
	    ie = iw[iier + m - 1];
	    for (j = 1; j <= 6; ++j) {
		irt = iw[iire + (ie - 1) * 6 + j - 1];
		i__3 = nl;
		for (k = nlo + 1; k <= i__3; ++k) {
		    if (irt == irr[k]) {
			goto L100;
		    }
		}
		++nl;
		if (nl > *maxl) {
		    *ierr = n;
		    goto L9000;
		}
		irr[nl] = irt;
L100:
		;
	    }
	}
	nrr[n] = nl;
    }
L9000:
    return 0;
} /* listr2r_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int liste2r_(integer *np, integer *nr, integer *
	ne, integer *ipe, integer *ire, integer *nep, integer *iep)
{
    /* Initialized data */

    static integer ipr[12]	/* was [2][6] */ = { 1,2,1,3,1,4,2,3,2,4,3,4 }
	    ;

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, m, n, m1, ie, ip1, ip2, ip3, jp1, jp2;

/* ================================================================ */
/*  The routine computes connectivity lists for mesh edges */

/*  nEP(*) - working array of size nP */
/*  IEP(*) - working array of size 4*nE */
/* ================================================================ */
/* ================================================================ */
/* group (Local variables) */
    /* Parameter adjustments */
    --iep;
    --nep;
    ire -= 7;
    ipe -= 5;

    /* Function Body */
/* ================================================================ */
    backreferences_(np, ne, &c__4, &c__4, &ipe[5], &nep[1], &iep[1]);
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	for (i__ = 1; i__ <= 6; ++i__) {
	    ire[i__ + n * 6] = 0;
	}
    }
    *nr = 0;
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	for (i__ = 1; i__ <= 6; ++i__) {
	    if (ire[i__ + n * 6] == 0) {
		++(*nr);
		ire[i__ + n * 6] = *nr;
		ip1 = ipe[ipr[(i__ << 1) - 2] + (n << 2)];
		ip2 = ipe[ipr[(i__ << 1) - 1] + (n << 2)];
		ip3 = max(ip1,ip2);
		m1 = 1;
		if (ip3 > 1) {
		    m1 = nep[ip3 - 1] + 1;
		}
		i__2 = nep[ip3];
		for (m = m1; m <= i__2; ++m) {
		    ie = iep[m];
		    for (j = 1; j <= 6; ++j) {
			jp1 = ipe[ipr[(j << 1) - 2] + (ie << 2)];
			jp2 = ipe[ipr[(j << 1) - 1] + (ie << 2)];
			if (check22_(&ip1, &ip2, &jp1, &jp2)) {
			    if (ire[j + ie * 6] == 0) {
				ire[j + ie * 6] = *nr;
			    }
			}
		    }
		}
	    }
/* L20: */
	}
    }
    return 0;
} /* liste2r_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int liste2f_(integer *np, integer *nf, integer *
	ne, integer *ipe, integer *ife, integer *nep, integer *iep)
{
    /* Initialized data */

    static integer iref[5] = { 1,2,3,4,1 };

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, n, i1, i2, i3, j1, j2, j3, ie2, ip1, ip2, ip3, jp1, 
	    jp2, jp3;

/* ================================================================ */
/*  The routine creates connectivity list E->F for mesh elements */

/*  *** Remarks: */
/*         1. Working memory is nEP(nP), IEP(4 * nE) */
/* ================================================================ */
/* ================================================================ */
/* group (Local variables) */
    /* Parameter adjustments */
    --iep;
    --nep;
    ife -= 5;
    ipe -= 5;

    /* Function Body */
/* ================================================================ */
    backreferences_(np, ne, &c__4, &c__4, &ipe[5], &nep[1], &iep[1]);
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	for (i__ = 1; i__ <= 4; ++i__) {
	    ife[i__ + (n << 2)] = 0;
	}
    }
    *nf = 0;
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	for (i1 = 1; i1 <= 4; ++i1) {
	    if (ife[i1 + (n << 2)] == 0) {
		++(*nf);
		ife[i1 + (n << 2)] = *nf;
		i2 = iref[i1];
		i3 = iref[i2];
		ip1 = ipe[i1 + (n << 2)];
		ip2 = ipe[i2 + (n << 2)];
		ip3 = ipe[i3 + (n << 2)];
		if (cmpe_(&ip1, &ip2, &ip3, &iep[1], &nep[1], &n, &ie2)) {
		    for (j1 = 1; j1 <= 4; ++j1) {
			j2 = iref[j1];
			j3 = iref[j2];
			jp1 = ipe[j1 + (ie2 << 2)];
			jp2 = ipe[j2 + (ie2 << 2)];
			jp3 = ipe[j3 + (ie2 << 2)];
			if (check33_(&ip1, &ip2, &ip3, &jp1, &jp2, &jp3)) {
			    ife[j1 + (ie2 << 2)] = *nf;
			    goto L10;
			}
		    }
		}
	    }
L10:
	    ;
	}
    }
    return 0;
} /* liste2f_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int liste2fb_(integer *np, integer *nfb, integer 
	*ne, integer *ipf, integer *ipe, integer *ife, integer *nep, integer *
	iep)
{
    /* Initialized data */

    static integer iref[5] = { 1,2,3,4,1 };

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer n, i1, i2, i3, if__, ip1, ip2, ip3;

/* ================================================================ */
/*  The routine computes connectivity list E->F for BOUNDARY faces */

/*  *** Remarks: */
/*         1. Working memory is nEP(nP), IEP(3 * nFb) */
/* ================================================================ */
/* ================================================================ */
/* group (Local variables) */
    /* Parameter adjustments */
    --iep;
    --nep;
    ife -= 5;
    ipe -= 5;
    ipf -= 4;

    /* Function Body */
/* ================================================================ */
    backreferences_(np, nfb, &c__3, &c__3, &ipf[4], &nep[1], &iep[1]);
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	for (i1 = 1; i1 <= 4; ++i1) {
	    ife[i1 + (n << 2)] = 0;
	    i2 = iref[i1];
	    i3 = iref[i2];
	    ip1 = ipe[i1 + (n << 2)];
	    ip2 = ipe[i2 + (n << 2)];
	    ip3 = ipe[i3 + (n << 2)];
	    if (cmpe_(&ip1, &ip2, &ip3, &iep[1], &nep[1], &c__0, &if__)) {
		ife[i1 + (n << 2)] = if__;
	    }
	}
    }
    return 0;
} /* liste2fb_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int listconv_(integer *np, integer *nr, integer *
	ne, integer *nep, integer *iep, integer *l, integer *ire, integer *nx,
	 integer *maxx, integer *nrp, integer *irp, integer *iw, integer *
	ierr)
{
    /* System generated locals */
    integer ire_dim1, ire_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, n, i1, i2, ie, ir, nx0;

/* ================================================================ */
/*  The routine convolutes unstructured map X->Y and structured map */
/*  Y->Z to get the map X->Z. For examples, if X means points (P), */
/*  Y means elements (E), and Z means edges (R), we get the map from */
/*  a point to all edges in the elements having this point. */

/*  Routine returns 0 upon successful completion. */

/*  *** Remarks: */
/*         1. Working memory is iW(nR) */
/* ================================================================ */
/* ================================================================ */
    /* Parameter adjustments */
    --nep;
    --iep;
    ire_dim1 = *l;
    ire_offset = 1 + ire_dim1;
    ire -= ire_offset;
    --nrp;
    --irp;
    --iw;

    /* Function Body */
    *ierr = 0;
    i__1 = *nr;
    for (n = 1; n <= i__1; ++n) {
	iw[n] = 1;
    }
    *nx = 0;
    i2 = 0;
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	nx0 = *nx;
	i1 = i2 + 1;
	i2 = nep[n];
	i__2 = i2;
	for (i__ = i1; i__ <= i__2; ++i__) {
	    ie = iep[i__];
	    i__3 = *l;
	    for (j = 1; j <= i__3; ++j) {
		ir = ire[j + ie * ire_dim1];
		if (iw[ir] > 0) {
		    ++(*nx);
		    if (*nx > *maxx) {
			*ierr = n;
			goto L9000;
		    }
		    irp[*nx] = ir;
		    iw[ir] = -iw[ir];
		}
	    }
	}
	nrp[n] = *nx;
/* ...    restore array iW */
	i__2 = *nx;
	for (i__ = nx0 + 1; i__ <= i__2; ++i__) {
	    ir = irp[i__];
	    iw[ir] = -iw[ir];
	}
    }
L9000:
    return 0;
} /* listconv_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int reversemap_(integer *np, integer *nr, 
	integer *nrp, integer *irp, integer *npr, integer *ipr)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, n, i1, i2, ir, iref;

/* ================================================================ */
/* Routine creates map R->P reverse to the map P->R. */
/* ================================================================ */
    /* Parameter adjustments */
    --ipr;
    --npr;
    --irp;
    --nrp;

    /* Function Body */
    i__1 = *nr;
    for (n = 1; n <= i__1; ++n) {
	npr[n] = 0;
    }
    i2 = 0;
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	i1 = i2 + 1;
	i2 = nrp[n];
	i__2 = i2;
	for (i__ = i1; i__ <= i__2; ++i__) {
	    ir = irp[i__];
	    ++npr[ir];
	}
    }
    i__1 = *nr;
    for (n = 2; n <= i__1; ++n) {
	npr[n] += npr[n - 1];
    }
    i2 = 0;
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	i1 = i2 + 1;
	i2 = nrp[n];
	i__2 = i2;
	for (i__ = i1; i__ <= i__2; ++i__) {
	    ir = irp[i__];
	    iref = npr[ir];
	    ipr[iref] = n;
	    npr[ir] = iref - 1;
	}
    }
    i2 = 0;
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	i1 = i2 + 1;
	i2 = nrp[n];
	i__2 = i2;
	for (i__ = i1; i__ <= i__2; ++i__) {
	    ir = irp[i__];
	    ++npr[ir];
	}
    }
    return 0;
} /* reversemap_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int backreferences_(integer *np, integer *ne, 
	integer *l, integer *m, integer *ipe, integer *nep, integer *iep)
{
    /* System generated locals */
    integer ipe_dim1, ipe_offset, i__1, i__2;

    /* Local variables */
    static integer i__, n, i1, iref;

/* ================================================================ */
/* Routine creates map P->E reverse to the map E->P. */
/*     nEP(P) - nEP(P-1) = number of elements having common */
/*                         point P. */
/*     IPE([nEP(P-1) + 1 : nEP(P)]) = list of elements having */
/*                         common point P. */
/* ================================================================ */
    /* Parameter adjustments */
    ipe_dim1 = *m;
    ipe_offset = 1 + ipe_dim1;
    ipe -= ipe_offset;
    --nep;
    --iep;

    /* Function Body */
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	nep[n] = 0;
    }
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	i__2 = *l;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i1 = ipe[i__ + n * ipe_dim1];
	    ++nep[i1];
	}
    }
    i__1 = *np;
    for (n = 2; n <= i__1; ++n) {
	nep[n] += nep[n - 1];
    }
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	i__2 = *l;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i1 = ipe[i__ + n * ipe_dim1];
	    iref = nep[i1];
	    iep[iref] = n;
	    nep[i1] = iref - 1;
	}
    }
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	i__2 = *l;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i1 = ipe[i__ + n * ipe_dim1];
	    ++nep[i1];
	}
    }
    return 0;
} /* backreferences_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int dualnormals_(integer *np, integer *nr, 
	integer *ne, integer *ipe, integer *ipr, doublereal *xyp, doublereal *
	nrm, integer *iw)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, m, n;
    static doublereal s;
    static integer i1, i2, ie, ip1, ip2, ip3, ip4, ipt;
    static doublereal xya[3], xyb[3], xyc[3], xyd[3], xye[3], xym[3], xyn[3], 
	    xyr[3];
    static integer iiep, inep, icnt;

/* ================================================================ */
/*  The routine computes normals to surfaces of a dual mesh */
/*  which separate points of a primary mesh. */

/*  Parameters: */
/*      nP, nR, nE - the number of points, edges and elements in the */
/*                   primary mesh */

/*      IPE(4, nE) - connectivity list of tetrahedra: E -> P */
/*      IPR(2, nR) - connectivity list of edges:      R -> P */

/*      XYP(3, nP) - Cartesian coordinates of mesh points */
/*      NRM(3, nR) - array of unit vectors normal to surfaces of a */
/*                   dual mesh separating point of the primary mesh */

/*      iW(*) - working memory of size 4 * nE + nP */
/* ================================================================ */
/* ================================================================ */
/* (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --iw;
    nrm -= 4;
    xyp -= 4;
    ipr -= 3;
    ipe -= 5;

    /* Function Body */
    inep = 0;
    iiep = inep + *np;
    backreferences_(np, ne, &c__4, &c__4, &ipe[5], &iw[inep + 1], &iw[iiep + 
	    1]);
    i__1 = *nr;
    for (n = 1; n <= i__1; ++n) {
	ip1 = ipr[(n << 1) + 1];
	ip2 = ipr[(n << 1) + 2];
	for (i__ = 1; i__ <= 3; ++i__) {
	    nrm[i__ + n * 3] = 0.;
	    xye[i__ - 1] = (xyp[i__ + ip1 * 3] + xyp[i__ + ip2 * 3]) / 2;
	    xyr[i__ - 1] = xyp[i__ + ip2 * 3] - xyp[i__ + ip1 * 3];
	}
	i1 = 1;
	if (ip1 > 1) {
	    i1 = iw[inep + ip1 - 1] + 1;
	}
	i2 = iw[inep + ip1];
	s = 0.;
	i__2 = i2;
	for (m = i1; m <= i__2; ++m) {
	    ie = iw[iiep + m];
	    icnt = 0;
	    for (j = 1; j <= 4; ++j) {
		ipt = ipe[j + (ie << 2)];
		if (ipt != ip1 && ipt != ip2) {
		    ++icnt;
		    if (icnt == 1) {
			ip3 = ipt;
		    }
		    if (icnt == 2) {
			ip4 = ipt;
		    }
		}
	    }
	    if (icnt != 2) {
		goto L10;
	    }
/*   ...   computing vertices of tet differ from iP1 & iP2 */
	    for (i__ = 1; i__ <= 3; ++i__) {
		xyc[i__ - 1] = (xyp[i__ + ip1 * 3] + xyp[i__ + ip2 * 3] + xyp[
			i__ + ip3 * 3] + xyp[i__ + ip4 * 3]) / 4;
		xya[i__ - 1] = (xyp[i__ + ip1 * 3] + xyp[i__ + ip2 * 3] + xyp[
			i__ + ip3 * 3]) / 3 - xyc[i__ - 1];
		xyb[i__ - 1] = (xyp[i__ + ip1 * 3] + xyp[i__ + ip2 * 3] + xyp[
			i__ + ip4 * 3]) / 3 - xyc[i__ - 1];
		xyd[i__ - 1] = xye[i__ - 1] - xyc[i__ - 1];
	    }
	    vecmul_(xyd, xya, xym);
	    vecmul_(xyd, xyb, xyn);
	    s = s + calnorm_(xym) + calnorm_(xyn);
/*    ...    orienting the normal vector from iP1 to iP2 */
	    if (dotmul_(xyr, xym) < 0.) {
		for (i__ = 1; i__ <= 3; ++i__) {
		    nrm[i__ + n * 3] -= xym[i__ - 1];
		}
	    } else {
		for (i__ = 1; i__ <= 3; ++i__) {
		    nrm[i__ + n * 3] += xym[i__ - 1];
		}
	    }
	    if (dotmul_(xyr, xyn) < 0.) {
		for (i__ = 1; i__ <= 3; ++i__) {
		    nrm[i__ + n * 3] -= xyn[i__ - 1];
		}
	    } else {
		for (i__ = 1; i__ <= 3; ++i__) {
		    nrm[i__ + n * 3] += xyn[i__ - 1];
		}
	    }
L10:
	    ;
	}
/*  ...  rescaling the area of the interface separating iP1 & iP2 */
/*        s = s / calNorm(NRM(1, n)) */
/*        Do i = 1, 3 */
/*           NRM(i, n) = NRM(i, n) * s */
/*        End do */
    }
    return 0;
} /* dualnormals_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int addboundaryfaces_(integer *np, integer *nf, 
	integer *maxf, integer *ne, doublereal *xyp, integer *ipf, integer *
	ipe, integer *lbf, integer *lbe, integer *iw)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k, n, i1, i2, i3, ip[5], ip1, ip2, ip3, iet, ift, iiep, 
	    iifp, infp, inep, kmax, ierr;
    static logical flage, flagf;

/* ================================================================ */
/* ================================================================ */
/*     iW(*) - working memory of size 2 * nP + 3 * nF + 4 * nE */
/* ================================================================ */
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
/* group (Functions) */
/* group (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --iw;
    --lbe;
    --lbf;
    ipe -= 5;
    ipf -= 4;
    xyp -= 4;

    /* Function Body */
    ierr = 0;
    ip[0] = 1;
    ip[1] = 2;
    ip[2] = 3;
    ip[3] = 4;
    ip[4] = 1;
    iifp = 1;
    iiep = iifp + *nf * 3;
    infp = iiep + (*ne << 2);
    inep = infp + *np;
/* ... checking that [iVface, MaxS] is clear */
    i__1 = *nf;
    for (n = 1; n <= i__1; ++n) {
	if (lbf[n] >= 900) {
	    errmes_(&c__1011, "addBoundaryFaces", "reserved boundary identif"
		    "icator is used", (ftnlen)16, (ftnlen)39);
	}
    }
/* ... creating an auxiliary structure */
    backreferences_(np, nf, &c__3, &c__3, &ipf[4], &iw[infp], &iw[iifp]);
    backreferences_(np, ne, &c__4, &c__4, &ipe[5], &iw[inep], &iw[iiep]);
/* ... creating material and missing boundaries */
    k = 0;
    kmax = 99;
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	for (i1 = 1; i1 <= 4; ++i1) {
	    i2 = ip[i1];
	    i3 = ip[i2];
	    ip1 = ipe[i1 + (n << 2)];
	    ip2 = ipe[i2 + (n << 2)];
	    ip3 = ipe[i3 + (n << 2)];
	    flagf = cmpe_(&ip1, &ip2, &ip3, &iw[iifp], &iw[infp], &c__0, &ift)
		    ;
	    flage = cmpe_(&ip1, &ip2, &ip3, &iw[iiep], &iw[inep], &n, &iet);
	    if (! flagf && ! flage) {
		++(*nf);
		if (*nf > *maxf) {
		    errmes_(&c__1004, "addBoundaryFaces", "local parameter M"
			    "axF is small", (ftnlen)16, (ftnlen)29);
		}
		ipf[*nf * 3 + 1] = ipe[i1 + (n << 2)];
		ipf[*nf * 3 + 2] = ipe[i2 + (n << 2)];
		ipf[*nf * 3 + 3] = ipe[i3 + (n << 2)];
		lbf[*nf] = 900;
	    }
	}
    }
    return 0;
} /* addboundaryfaces_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int addmaterialfaces_(integer *np, integer *nf, 
	integer *maxf, integer *ne, doublereal *xyp, integer *ipf, integer *
	ipe, integer *lbf, integer *lbe, integer *iw)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, k, n, i1, i2, i3, ip[5], ip1, ip2, ip3, iet, ift, 
	    mat1, mat2, iiep, iifp, infp, inep, kmax, ierr;
    static logical flage, flagf;
    static integer mlist[200]	/* was [2][100] */, icface;

/* ================================================================ */
/* ================================================================ */
/*     iW(*) - working memory of size 2 * nP + 3 * nF + 4 * nE */
/* ================================================================ */
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
/* group (Functions) */
/* group (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --iw;
    --lbe;
    --lbf;
    ipe -= 5;
    ipf -= 4;
    xyp -= 4;

    /* Function Body */
    ierr = 0;
    ip[0] = 1;
    ip[1] = 2;
    ip[2] = 3;
    ip[3] = 4;
    ip[4] = 1;
    iifp = 1;
    iiep = iifp + *nf * 3;
    infp = iiep + (*ne << 2);
    inep = infp + *np;
/* ... checking that [iVface, MaxS] is clear */
    i__1 = *nf;
    for (n = 1; n <= i__1; ++n) {
	if (lbf[n] >= 900) {
	    errmes_(&c__1011, "addMaterialFaces", "reserved boundary identif"
		    "icator is used", (ftnlen)16, (ftnlen)39);
	}
    }
/* ... creating an auxiliary structure */
    backreferences_(np, nf, &c__3, &c__3, &ipf[4], &iw[infp], &iw[iifp]);
    backreferences_(np, ne, &c__4, &c__4, &ipe[5], &iw[inep], &iw[iiep]);
/* ... creating material and missing boundaries */
    k = 0;
    kmax = 99;
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	for (i1 = 1; i1 <= 4; ++i1) {
	    i2 = ip[i1];
	    i3 = ip[i2];
	    ip1 = ipe[i1 + (n << 2)];
	    ip2 = ipe[i2 + (n << 2)];
	    ip3 = ipe[i3 + (n << 2)];
	    flagf = cmpe_(&ip1, &ip2, &ip3, &iw[iifp], &iw[infp], &c__0, &ift)
		    ;
	    flage = cmpe_(&ip1, &ip2, &ip3, &iw[iiep], &iw[inep], &n, &iet);
	    if (! flagf && flage && iet > n) {
		mat1 = lbe[n];
		mat2 = lbe[iet];
		if (mat1 != mat2) {
/*  ...  searching for this pair in the list of material interfaces */
		    i__2 = k;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			if (mlist[(i__ << 1) - 2] == mat1 && mlist[(i__ << 1) 
				- 1] == mat2 || mlist[(i__ << 1) - 2] == mat2 
				&& mlist[(i__ << 1) - 1] == mat1) {
			    icface = i__ + 900;
			    goto L1;
			}
		    }
/*  ...  making the new material interface */
		    ++k;
		    if (k > kmax) {
			errmes_(&c__1010, "addMaterialFaces", "not enough me"
				"mory for material faces", (ftnlen)16, (ftnlen)
				36);
		    }
		    mlist[(k << 1) - 2] = mat1;
		    mlist[(k << 1) - 1] = mat2;
		    icface = k + 900;
L1:
		    ++(*nf);
		    if (*nf > *maxf) {
			errmes_(&c__1004, "addMaterialFaces", "local paramet"
				"er MaxF is small", (ftnlen)16, (ftnlen)29);
		    }
		    ipf[*nf * 3 + 1] = ipe[i1 + (n << 2)];
		    ipf[*nf * 3 + 2] = ipe[i2 + (n << 2)];
		    ipf[*nf * 3 + 3] = ipe[i3 + (n << 2)];
		    lbf[*nf] = icface;
		}
	    }
	}
    }
    return 0;
} /* addmaterialfaces_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int global2local_(integer *myid, integer *ice, 
	integer *np, integer *nf, integer *maxf, integer *ne, doublereal *xyp,
	 integer *ipf, integer *ipe, integer *lbf, integer *lbe, integer *npv,
	 integer *nfv, integer *nev, integer *ipv, integer *ifv, integer *iev,
	 integer *npl, integer *nfl, integer *nel, doublereal *xypl, integer *
	ipfl, integer *ipel, integer *lbfl, integer *lbel, integer *npvl, 
	integer *nfvl, integer *nevl, integer *ipvl, integer *ifvl, integer *
	ievl, integer *nfvi, integer *ippl, integer *iffl, integer *ipw, 
	integer *ifw, integer *iew)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, m, n, ic, ie, if__, ip, ip1, ip2, ipt;
    static logical flag__;
    static integer nflo, mfvi, nfvt, icfree;

/* ================================================================ */
/* group (Mg) */
/* group (Ml) */
/* group (I) */
/* group (W) */
/* ================================================================ */
/*  The subroutine extracts a submesh from the global mesh using */
/*  tetrahedra colored by myID color in array ICE. */

/*  The fixed triangles are placed at the beginning of corresponding */
/*  list. The interfaces triangles are created with utility */
/*  addBoundaryFaces and are placed right after the fixed surface */
/*  triangles. The number of these triangles equals to nFvi. Note */
/*  that these triangles are not added to the list of fixed triangles. */

/*  The intereface points are also placed at the beginning of list. */
/*  The user don't need to include them in the list of fixed trianges. */
/*  It will be done automatically when the interface triangles will */
/*  be added to the list of fixed triangles. */

/*  IPFl(MaxF) - since we need to add interface triangles, nFl may be */
/*               bigger than nF. */

/*  IPPl(nPl)  - references to a global enumeration necessary to */
/*               glue submeshes in one global mesh: */
/*                  0 - point is located outside the DD interfaces; */
/*                 >0 - splitted points from different subdomains */
/*                      have the same value of IPPl. */

/* IFFl(nFl)   - references to global enumeration of surface triangles: */
/*                  0 - the auxiliary interface triangle */
/*                 >0 - splitted triangles from different meshes have */
/*                      the same value of IFFl. */

/*  In order to use the information inside IPPl, the interface */
/*  triangles should not be modified and be always at the beginning */
/*  of the corresponding list. The user has to add the interface */
/*  triangles (nFvi triangles) to the list of fixed triangles. */

/*  Working memory: IPw(2*nP), IFw(nF), IEw(2*nP + 3*nF + 4*nE) */

/*  Remarks: The number of interface triangles is correct if and */
/*  only if the global mesh is complete, i.e. there are no missing */
/*  surface triangles. If necessary, the user has to use utility */
/*  addBoundaryFaces before splitting the global mesh. */

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
/* group (Mg) */
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
/* group (Ml) */
/* group (I) */
/* group (W) */
/* group (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --iew;
    --ifw;
    --ipw;
    --iffl;
    --ippl;
    --ievl;
    --ifvl;
    --ipvl;
    --lbel;
    --lbfl;
    ipel -= 5;
    ipfl -= 4;
    xypl -= 4;
    --iev;
    --ifv;
    --ipv;
    --lbe;
    --lbf;
    ipe -= 5;
    ipf -= 4;
    xyp -= 4;
    --ice;

    /* Function Body */
    ip1 = 1;
    ip2 = ip1 + *np;
    maktnode_(np, np, ne, &ipw[ip2], &ipe[5], &ice[1]);
/* ... marking points of the subdomain with myID */
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	ipw[n] = 0;
    }
    *nel = 0;
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	if (ice[n] == *myid) {
	    ++(*nel);
	    iew[n] = *nel;
	    for (i__ = 1; i__ <= 4; ++i__) {
		ip = ipe[i__ + (n << 2)];
		ipw[ip] = ip;
		ipel[i__ + (*nel << 2)] = ip;
	    }
	    lbel[*nel] = lbe[n];
	}
    }
    *nevl = 0;
    i__1 = *nev;
    for (n = 1; n <= i__1; ++n) {
	ie = iev[n];
	if (ice[ie] == *myid) {
	    ++(*nevl);
	    ievl[*nevl] = iew[ie];
	}
    }
/* ... coping points */
    *npl = 0;
    for (m = 1; m <= 2; ++m) {
	i__1 = *np;
	for (n = 1; n <= i__1; ++n) {
	    ipt = ip2 + n - 1;
	    flag__ = m == 1 && ipw[ipt] > 0 || m == 2 && ipw[ipt] == 0;
	    if (ipw[n] != 0 && flag__) {
		++(*npl);
		ipw[n] = *npl;
		for (i__ = 1; i__ <= 3; ++i__) {
		    xypl[i__ + *npl * 3] = xyp[i__ + n * 3];
		}
		ippl[*npl] = ipw[ip2 + n - 1];
	    }
	}
    }
    *npvl = 0;
    i__1 = *npv;
    for (n = 1; n <= i__1; ++n) {
	if (ipw[ipv[n]] != 0) {
	    ++(*npvl);
	    ipvl[*npvl] = ipw[ipv[n]];
	}
    }
/* ... coping faces (fixed faces) */
    i__1 = *maxf;
    for (n = 1; n <= i__1; ++n) {
	iffl[n] = 0;
    }
    i__1 = *nf;
    for (n = 1; n <= i__1; ++n) {
	ifw[n] = 0;
    }
    *nfl = 0;
    i__1 = *nfv;
    for (n = 1; n <= i__1; ++n) {
	if__ = ifv[n];
	for (i__ = 1; i__ <= 3; ++i__) {
	    if (ipw[ipf[i__ + if__ * 3]] == 0) {
		goto L100;
	    }
	}
	++(*nfl);
	ifw[if__] = *nfl;
	for (i__ = 1; i__ <= 3; ++i__) {
	    ipfl[i__ + *nfl * 3] = ipw[ipf[i__ + if__ * 3]];
	}
	lbfl[*nfl] = lbf[if__];
	iffl[*nfl] = if__;
L100:
	;
    }
/* ... converting some surface triangles to interface triangles */
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	ipw[ip2 + n] = ipw[n];
    }
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	if (ice[n] != *myid) {
	    for (i__ = 1; i__ <= 4; ++i__) {
		ip1 = ipe[i__ + (n << 2)];
		ipw[ip2 + ip1] = -ipw[ip2 + ip1];
	    }
	}
    }
    mfvi = 0;
    i__1 = *nf;
    for (n = 1; n <= i__1; ++n) {
	if (ifw[n] > 0) {
	    goto L200;
	}
	for (i__ = 1; i__ <= 3; ++i__) {
	    if (ipw[ip2 + ipf[i__ + n * 3]] >= 0) {
		goto L200;
	    }
	}
	++mfvi;
	++(*nfl);
	ifw[n] = *nfl;
	for (i__ = 1; i__ <= 3; ++i__) {
	    ipfl[i__ + *nfl * 3] = ipw[ipf[i__ + n * 3]];
	}
	lbfl[*nfl] = lbf[n];
	iffl[*nfl] = n;
L200:
	;
    }
    nfvt = *nfl;
/* ... coping the rest of faces */
    i__1 = *nf;
    for (n = 1; n <= i__1; ++n) {
	if (ifw[n] > 0) {
	    goto L300;
	}
	for (i__ = 1; i__ <= 3; ++i__) {
	    if (ipw[ipf[i__ + n * 3]] == 0) {
		goto L300;
	    }
	}
	++(*nfl);
	ifw[n] = *nfl;
	for (i__ = 1; i__ <= 3; ++i__) {
	    ipfl[i__ + *nfl * 3] = ipw[ipf[i__ + n * 3]];
	}
	lbfl[*nfl] = lbf[n];
	iffl[*nfl] = n;
L300:
	;
    }
    *nfvl = 0;
    i__1 = *nfv;
    for (n = 1; n <= i__1; ++n) {
	if__ = ifv[n];
	if (ifw[if__] != 0) {
	    ++(*nfvl);
	    ifvl[*nfvl] = ifw[if__];
	}
    }
/* ... computing local coordinates for elements */
    i__1 = *nel;
    for (n = 1; n <= i__1; ++n) {
	for (i__ = 1; i__ <= 4; ++i__) {
	    ipel[i__ + (n << 2)] = ipw[ipel[i__ + (n << 2)]];
	}
    }
/* ... computing the interface triangles */
    nflo = *nfl;
    addboundaryfaces_(npl, nfl, maxf, nel, &xypl[4], &ipfl[4], &ipel[5], &
	    lbfl[1], &lbel[1], &iew[1]);
    *nfvi = *nfl - nflo;
/* ... moving the fixed and interface triangles */
    m = *nfvl + mfvi + 1;
    i__1 = nflo + 1;
    for (n = *nfl; n >= i__1; --n) {
	for (i__ = 1; i__ <= 3; ++i__) {
	    swapii_(&ipfl[i__ + n * 3], &ipfl[i__ + m * 3]);
	}
	swapii_(&lbfl[n], &lbfl[m]);
	swapii_(&iffl[n], &iffl[m]);
	++m;
	if (m > nflo) {
	    goto L500;
	}
    }
/* ... changing color of interface triangles (if possible) */
L500:
    *nfvi += mfvi;
    icfree = 900;
    for (ic = 1; ic <= 899; ++ic) {
	i__1 = *nfl;
	for (n = 1; n <= i__1; ++n) {
	    if (lbfl[n] == ic) {
		goto L600;
	    }
	}
	icfree = ic;
	goto L700;
L600:
	;
    }
L700:
    i__1 = *nfl;
    for (n = 1; n <= i__1; ++n) {
	if (lbfl[n] == 900) {
	    lbfl[n] = icfree;
	}
    }
    return 0;
} /* global2local_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int local2global_(integer *nmeshes, integer *
	maxp, integer *maxf, integer *maxe, integer *np, integer *nf, integer 
	*ne, doublereal *xyp, integer *ipf, integer *ipe, integer *lbf, 
	integer *lbe, integer *npv, integer *nfv, integer *nev, integer *ipv, 
	integer *ifv, integer *iev, integer *npl, integer *nfl, integer *nel, 
	doublereal *xypl, integer *ipfl, integer *ipel, integer *lbfl, 
	integer *lbel, integer *npvl, integer *nfvl, integer *nevl, integer *
	ipvl, integer *ifvl, integer *ievl, integer *ippl, integer *iffl, 
	integer *ips, integer *ipw, integer *ifw, integer *iew)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, k, n, if__, me, mf, kp, mp, mv, ip1, ip2, ip3;

/* ================================================================ */
/* group (Mg) */
/* group (Ml) */
/* group (I) */
/* group (W) */
/* ================================================================ */
/*  The subroutine gathers submeshes into a global mesh using */
/*  references to global enumeration of mesh points kept in array */
/*  IPPl(*) and the global enumeration of faces in array IFFl(*). */

/*  The meshes are kept as continuous lists of data. For instance, */
/*  array XYP(1:3, 1:nP1) keeps coordinates of mesh points in the */
/*  first grid. Continuation of the array, XYP(1:3, nP1 : nP1 + nP2), */
/*  keeps coordinates of the mesh points of the second grid. And */
/*  so on. The array having such a structure are: */

/*    XYPl(3, *) - mesh coordinates */
/*    IPFl(3, *) - connectivity list for triangles */
/*    IPEl(4, *) - connectivity list for tetrahedra */

/*    IPVl(*) - list of fixed points */
/*    IFVl(*) - list of fixed triangles */
/*    IEVl(*) - list of fixed tetrahedra */

/*    IPPl(*) - references to global enumeration of mesh point */
/*    IFFl(*) - references to global enumeration of mesh triangles */


/*  The necessary information to extract submesh data are kept in */
/*  the following arrays: */

/*    nPl(nMeshes)  - numbers of points */
/*    nFl(nMeshes)  - numbers of faces */
/*    nEl(nMeshes)  - numbers of tetrahedra */

/*    nPvl(nMeshes) - numbers of fixed points */
/*    nFvl(nMeshes) - numbers of fixed triangles */
/*    nEvl(nMeshes) - numbers of fixed tetrahedra */

/*  Working memory: IPs(MaxP), IPw(MaxP), IFw(MaxF), IEw(4*MaxE) */

/*  Remark: Array IPFl is destroyed in the algorithm. */
/* ================================================================ */
/* group (Mg) */
/* group (Ml) */
/* group (I) */
/* group (W) */
/* group (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --iew;
    --ifw;
    --ipw;
    --ips;
    --iffl;
    --ippl;
    --ievl;
    --ifvl;
    --ipvl;
    --nevl;
    --nfvl;
    --npvl;
    --lbel;
    --lbfl;
    ipel -= 5;
    ipfl -= 4;
    xypl -= 4;
    --nel;
    --nfl;
    --npl;
    --iev;
    --ifv;
    --ipv;
    --lbe;
    --lbf;
    ipe -= 5;
    ipf -= 4;
    xyp -= 4;

    /* Function Body */
    mp = 0;
    mf = 0;
    me = 0;
    i__1 = *nmeshes;
    for (n = 1; n <= i__1; ++n) {
	mp += npl[n];
	mf += nfl[n];
	me += nel[n];
    }
    if (mp > *maxp) {
	errmes_(&c__1003, "local2global", "local parameter MaxP is small", (
		ftnlen)12, (ftnlen)29);
    }
    if (mf > *maxf) {
	errmes_(&c__1004, "local2global", "local parameter MaxF is small", (
		ftnlen)12, (ftnlen)29);
    }
    if (me > *maxe) {
	errmes_(&c__1006, "local2global", "local parameter MaxE is small", (
		ftnlen)12, (ftnlen)29);
    }
    i__1 = mp;
    for (n = 1; n <= i__1; ++n) {
	ips[n] = 0;
	ipw[n] = 0;
    }
/* ... counting points on interfeices */
    kp = 0;
    i__1 = mp;
    for (n = 1; n <= i__1; ++n) {
	if (ippl[n] != 0) {
	    if (ips[ippl[n]] == 0) {
		++kp;
		ips[ippl[n]] = kp;
	    }
	}
    }
/* ... counting the rest of points */
    i__1 = mp;
    for (n = 1; n <= i__1; ++n) {
	if (ippl[n] == 0) {
	    ++kp;
	    ipw[n] = kp;
	} else {
	    ipw[n] = ips[ippl[n]];
	}
    }
/* ... making points */
    *np = kp;
    i__1 = mp;
    for (n = 1; n <= i__1; ++n) {
	for (i__ = 1; i__ <= 3; ++i__) {
	    xyp[i__ + ipw[n] * 3] = xypl[i__ + n * 3];
	}
    }
/* ... making vertices */
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	ips[n] = 0;
    }
    mp = 0;
    mv = 0;
    i__1 = *nmeshes;
    for (k = 1; k <= i__1; ++k) {
	i__2 = npvl[k];
	for (n = 1; n <= i__2; ++n) {
	    ++mv;
	    ips[ipw[mp + ipvl[mv]]] = 1;
	}
	mp += npl[k];
    }
    *npv = 0;
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	if (ips[n] != 0) {
	    ++(*npv);
	    ipv[*npv] = n;
	}
    }
/* ... making faces */
    mp = 0;
    mf = 0;
    i__1 = *nmeshes;
    for (k = 1; k <= i__1; ++k) {
	i__2 = nfl[k];
	for (n = 1; n <= i__2; ++n) {
	    ++mf;
	    ifw[mf] = 0;
	    for (i__ = 1; i__ <= 3; ++i__) {
		ipfl[i__ + mf * 3] = ipw[mp + ipfl[i__ + mf * 3]];
	    }
	}
	mp += npl[k];
    }
    backreferences_(np, &mf, &c__3, &c__3, &ipfl[4], &ips[1], &iew[1]);
    i__1 = mf;
    for (n = 1; n <= i__1; ++n) {
	if (ifw[n] == 0) {
	    ip1 = ipfl[n * 3 + 1];
	    ip2 = ipfl[n * 3 + 2];
	    ip3 = ipfl[n * 3 + 3];
	    if (cmpe_(&ip1, &ip2, &ip3, &iew[1], &ips[1], &n, &if__)) {
		if (if__ > n) {
		    ifw[n] = if__;
		}
	    }
	}
    }
    *nf = 0;
    i__1 = mf;
    for (n = 1; n <= i__1; ++n) {
	if (ifw[n] == 0 && iffl[n] > 0) {
	    ++(*nf);
	    ifw[n] = *nf;
	    for (i__ = 1; i__ <= 3; ++i__) {
		ipf[i__ + *nf * 3] = ipfl[i__ + n * 3];
	    }
	    lbf[*nf] = lbfl[n];
	} else {
	    ifw[n] = 0;
	}
    }
/* ... making fixed faces */
    i__1 = *nf;
    for (n = 1; n <= i__1; ++n) {
	iew[n] = 0;
    }
    mf = 0;
    mv = 0;
    i__1 = *nmeshes;
    for (k = 1; k <= i__1; ++k) {
	i__2 = nfvl[k];
	for (n = 1; n <= i__2; ++n) {
	    ++mv;
	    iew[ifw[mf + ifvl[mv]]] = 1;
	}
	mf += nfl[k];
    }
    *nfv = 0;
    i__1 = *nf;
    for (n = 1; n <= i__1; ++n) {
	if (iew[n] != 0) {
	    ++(*nfv);
	    ifv[*nfv] = n;
	}
    }
/* ... making elements */
    mp = 0;
    *ne = 0;
    i__1 = *nmeshes;
    for (k = 1; k <= i__1; ++k) {
	i__2 = nel[k];
	for (n = 1; n <= i__2; ++n) {
	    ++(*ne);
	    for (i__ = 1; i__ <= 4; ++i__) {
		ipe[i__ + (*ne << 2)] = ipw[mp + ipel[i__ + (*ne << 2)]];
	    }
	    lbe[*ne] = lbel[*ne];
	}
	mp += npl[k];
    }
/* ... making fixed elements */
    *nev = 0;
    me = 0;
    mv = 0;
    i__1 = *nmeshes;
    for (k = 1; k <= i__1; ++k) {
	i__2 = nevl[k];
	for (n = 1; n <= i__2; ++n) {
	    ++(*nev);
	    iev[*nev] = me + ievl[mv + n];
	}
	me += nel[k];
	mv += nevl[k];
    }
    return 0;
} /* local2global_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int maktnode_(integer *np, integer *maxp, 
	integer *ne, integer *ipp, integer *ipe, integer *ice)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, n, ip, icn, ict;

/* ================================================================ */
/*   The T-nodes are marked. The original connectivity */
/*   list IPE is used. */
/* ================================================================ */
    /* Parameter adjustments */
    --ice;
    ipe -= 5;
    --ipp;

    /* Function Body */
    i__1 = *maxp;
    for (n = 1; n <= i__1; ++n) {
	ipp[n] = 0;
    }
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	ict = ice[n] + 1;
	for (i__ = 1; i__ <= 4; ++i__) {
	    ip = ipe[i__ + (n << 2)];
	    icn = ipp[ip];
	    if (icn == 0) {
		ipp[ip] = ict;
	    } else if (icn != ict) {
		ipp[ip] = -1;
	    }
	}
    }
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	if (ipp[n] == -1) {
	    ipp[n] = n;
	} else {
	    ipp[n] = 0;
	}
    }
    return 0;
} /* maktnode_ */

