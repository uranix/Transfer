/* refine.f -- translated by f2c (version 20090411).
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
#include "refine.h"
#include "lapack.h"

/* Common Block Declarations */

struct {
    doublereal refxyp[3], scaxyp[3];
} anicrv_;

#define anicrv_1 anicrv_

/* Table of constant values */

static integer c__3 = 3;
static integer c__1001 = 1001;
static integer c__1003 = 1003;
static integer c__1004 = 1004;
static integer c__1006 = 1006;
static integer c__4 = 4;
static integer c__0 = 0;
static integer c__1002 = 1002;
static integer c__1005 = 1005;
static integer c__1007 = 1007;
static integer c__1008 = 1008;

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int initializerefinement_(integer *np, integer *
	ne, doublereal *xyp, integer *ipe, doublereal *mapmtr, integer *
	ref2mapmtr)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal a[9]	/* was [3][3] */, b[9]	/* was [3][3] */;
    static integer i__, j, k, n;
    static doublereal aa[9]	/* was [3][3] */;
    static integer info, ipiv[3];
    static doublereal xypet[12]	/* was [3][4] */;

/* ================================================================ */
/* ================================================================ */
/*   Fill mapping matrix for each tetrahedron to make it equilateral */
/*   Fill reference array to MapMtr */
/* ================================================================ */
/* Local */
/* Equilateral tet */
    /* Parameter adjustments */
    --ref2mapmtr;
    mapmtr -= 13;
    ipe -= 5;
    xyp -= 4;

    /* Function Body */
    xypet[0] = 0.;
    xypet[1] = 0.;
    xypet[2] = 0.;
    xypet[3] = sqrt(2.) * 2 / 3;
    xypet[4] = 0.;
    xypet[5] = -1.3333333333333333;
    xypet[6] = -sqrt(2.) / 3;
    xypet[7] = sqrt(6.) / 3;
    xypet[8] = -1.3333333333333333;
    xypet[9] = -sqrt(2.) / 3;
    xypet[10] = -sqrt(6.) / 3;
    xypet[11] = -1.3333333333333333;
/* Find map Ortotet to Equitet */
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    aa[i__ + j * 3 - 4] = xypet[i__ + (j + 1) * 3 - 4] - xypet[i__ - 
		    1];
	}
    }
/* Loop over tets */
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
/* ...  find map  Ortotet to Tet */
	for (i__ = 1; i__ <= 3; ++i__) {
	    for (j = 1; j <= 3; ++j) {
		a[i__ + j * 3 - 4] = xyp[i__ + ipe[j + 1 + (n << 2)] * 3] - 
			xyp[i__ + ipe[(n << 2) + 1] * 3];
	    }
	}
/* ...  find map Tet to Ortotet */
	for (i__ = 1; i__ <= 3; ++i__) {
	    for (j = 1; j <= 3; ++j) {
		b[i__ + j * 3 - 4] = 0.;
	    }
	    b[i__ + i__ * 3 - 4] = 1.;
	}
	dgesv_(&c__3, &c__3, a, &c__3, ipiv, b, &c__3, &info);
	if (info != 0) {
	    errmes_(&c__1001, "InitializeRefinement", "error in dgesv", (
		    ftnlen)20, (ftnlen)14);
	}
/* ...  find map Tet to Equitet */
	for (i__ = 1; i__ <= 3; ++i__) {
	    for (j = 1; j <= 3; ++j) {
		mapmtr[i__ + (j + n * 3) * 3] = 0.;
		for (k = 1; k <= 3; ++k) {
		    mapmtr[i__ + (j + n * 3) * 3] += aa[i__ + k * 3 - 4] * b[
			    k + j * 3 - 4];
		}
	    }
	}
	ref2mapmtr[n] = n;
    }
    return 0;
} /* initializerefinement_ */

/* ================================================================ */
/* This routine may be called ONLY after initializeRefinement */
/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int uniformrefinement_(integer *np, integer *
	maxp, integer *nf, integer *maxf, integer *ne, integer *maxe, 
	doublereal *xyp, integer *ipf, integer *ipe, integer *lbf, integer *
	lbe, doublereal *mapmtr, integer *ref2mapmtr, integer *iw, integer *
	maxwi)
{
    /* Initialized data */

    static integer refoct[48]	/* was [4][4][3] */ = { 4,5,6,1,2,3,6,1,1,3,5,
	    6,1,2,4,6,1,3,2,5,3,6,2,5,6,4,2,5,4,1,2,5,1,2,3,4,2,6,3,4,6,5,3,4,
	    5,1,3,4 };
    static integer refcrnr[16]	/* was [4][4] */ = { 1,2,3,1,1,4,5,2,2,4,6,3,
	    3,5,6,4 };
    static integer loop[4] = { 1,2,3,1 };

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k, n, i1, i2, j1, j2, i3, ie, if__, jf, ii, jj, ir,
	     kr, nr;
    static doublereal qt;
    static integer ip1, ip2, ip3, jp1, jp2, neo, nfo, npo, iend, iiep, iipe, 
	    iire, inep, ipes[4];
    static doublereal xypt[12]	/* was [3][4] */;
    static logical flage;
    static doublereal minqt[3];

/* ================================================================ */
/* ================================================================ */
/*    MapMtr - data array of length 9*nEcoarse */
/*    Ref2MapMtr - data array (lengths on input and output are 9*nE) */
/*    iW(*) - working array of size MaxWi which is at least */
/*            nP + 14 * nE */
/* ================================================================ */
/* (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --iw;
    --ref2mapmtr;
    mapmtr -= 10;
    --lbe;
    --lbf;
    ipe -= 5;
    ipf -= 4;
    xyp -= 4;

    /* Function Body */
/* ================================================================ */
    inep = 1;
    iiep = inep + *np;
    iire = iiep + (*ne << 2);
    iipe = iire + *ne * 6;
    iend = iipe + (*ne << 2);
    liste2r_(np, &nr, ne, &ipe[5], &iw[iire], &iw[inep], &iw[iiep]);
    npo = *np;
    nfo = *nf;
    neo = *ne;
    *np = npo + nr;
    *nf = nfo << 2;
    *ne = neo << 3;
    if (*np > *maxp) {
	errmes_(&c__1003, "uniformRefinement", "local parameter MaxP is small"
		, (ftnlen)17, (ftnlen)29);
    }
    if (*nf > *maxf) {
	errmes_(&c__1004, "uniformRefinement", "local parameter MaxF is small"
		, (ftnlen)17, (ftnlen)29);
    }
    if (*ne > *maxe) {
	errmes_(&c__1006, "uniformRefinement", "local parameter MaxE is small"
		, (ftnlen)17, (ftnlen)29);
    }
/* ... mapping P -> E */
    backreferences_(&npo, &neo, &c__4, &c__4, &ipe[5], &iw[inep], &iw[iiep]);
    k = 0;
    i__1 = neo;
    for (n = 1; n <= i__1; ++n) {
	for (i__ = 1; i__ <= 4; ++i__) {
	    ++k;
	    iw[iipe + k - 1] = ipe[i__ + (n << 2)];
	}
    }
/* ... splitting elements */
    kr = 0;
    i__1 = neo;
    for (n = 1; n <= i__1; ++n) {
	for (i1 = 1; i1 <= 3; ++i1) {
	    for (i2 = i1 + 1; i2 <= 4; ++i2) {
		ip1 = ipe[i1 + (n << 2)];
		ip2 = ipe[i2 + (n << 2)];
		++kr;
		ir = iw[iire + kr - 1];
		for (i__ = 1; i__ <= 3; ++i__) {
		    xyp[i__ + (npo + ir) * 3] = (xyp[i__ + ip1 * 3] + xyp[i__ 
			    + ip2 * 3]) / 2;
		}
	    }
	}
	for (i__ = 1; i__ <= 4; ++i__) {
	    ipes[i__ - 1] = ipe[i__ + (n << 2)];
	}
/* ... find optimal splitting of inner octahedron (ii=1,2,3) */
	for (ii = 1; ii <= 3; ++ii) {
	    minqt[ii - 1] = 1.;
	    for (i__ = 1; i__ <= 4; ++i__) {
		for (k = 1; k <= 4; ++k) {
		    j = npo + iw[iire + n * 6 + refoct[k + (i__ + (ii << 2) <<
			     2) - 21] - 7];
		    jj = ref2mapmtr[n];
		    if (jj < 0) {
			jj = -jj;
		    }
		    fixedmap_(&xyp[j * 3 + 1], &mapmtr[jj * 9 + 1], &xypt[k * 
			    3 - 3]);
		}
		calreg_(xypt, &xypt[3], &xypt[6], &xypt[9], &qt);
/* Computing MIN */
		d__1 = qt, d__2 = minqt[ii - 1];
		minqt[ii - 1] = min(d__1,d__2);
	    }
	}
	if (minqt[0] >= minqt[1] && minqt[0] >= minqt[2]) {
	    ii = 1;
	} else if (minqt[1] >= minqt[0] && minqt[1] >= minqt[2]) {
	    ii = 2;
	} else if (minqt[2] >= minqt[0] && minqt[2] >= minqt[1]) {
	    ii = 3;
	}
/* ... apply optimal splitting */
	for (i__ = 1; i__ <= 4; ++i__) {
	    ie = (i__ + 3) * neo + n;
	    for (k = 1; k <= 4; ++k) {
		ipe[k + (ie << 2)] = npo + iw[iire + n * 6 + refoct[k + (i__ 
			+ (ii << 2) << 2) - 21] - 7];
	    }
	    lbe[ie] = lbe[n];
	    ref2mapmtr[ie] = ref2mapmtr[n];
	    ie = (4 - i__) * neo + n;
	    for (k = 1; k <= 3; ++k) {
		ipe[k + (ie << 2)] = npo + iw[iire + n * 6 + refcrnr[k + (i__ 
			<< 2) - 5] - 7];
	    }
	    ipe[(ie << 2) + 4] = ipes[i__ - 1];
	    lbe[ie] = lbe[n];
	    ref2mapmtr[ie] = ref2mapmtr[n];
	}
    }
/* ... splitting faces */
    i__1 = nfo;
    for (n = 1; n <= i__1; ++n) {
	ip1 = ipf[n * 3 + 1];
	ip2 = ipf[n * 3 + 2];
	ip3 = ipf[n * 3 + 3];
	flage = cmpe_(&ip1, &ip2, &ip3, &iw[iiep], &iw[inep], &c__0, &ie);
	jf = nfo * 3 + n;
	k = 0;
	for (i1 = 1; i1 <= 3; ++i1) {
	    for (i2 = i1 + 1; i2 <= 4; ++i2) {
		++k;
		ip1 = iw[iipe + (ie << 2) + i1 - 5];
		ip2 = iw[iipe + (ie << 2) + i2 - 5];
		for (j1 = 1; j1 <= 3; ++j1) {
		    j2 = loop[j1];
		    jp1 = ipf[j1 + n * 3];
		    jp2 = ipf[j2 + n * 3];
		    if (check22_(&ip1, &ip2, &jp1, &jp2)) {
			ipf[j1 + jf * 3] = npo + iw[iire + ie * 6 + k - 7];
			goto L10;
		    }
		}
L10:
		;
	    }
	}
	lbf[jf] = lbf[n];
	for (i1 = 1; i1 <= 3; ++i1) {
	    if__ = (3 - i1) * nfo + n;
	    i2 = loop[i1];
	    i3 = loop[i2];
	    ipf[if__ * 3 + 1] = ipf[i1 + n * 3];
	    ipf[if__ * 3 + 2] = ipf[i1 + jf * 3];
	    ipf[if__ * 3 + 3] = ipf[i3 + jf * 3];
	    lbf[if__] = lbf[n];
	}
    }
    return 0;
} /* uniformrefinement_ */

/* ================================================================ */
/* This routine may be called ONLY after initializeRefinement */
/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int localrefinement_(integer *np, integer *maxp, 
	integer *nf, integer *maxf, integer *ne, integer *maxe, doublereal *
	xyp, integer *ipf, integer *ipe, integer *lbf, integer *lbe, 
	doublereal *mapmtr, integer *ref2mapmtr, logical *splitflag, integer *
	iw, integer *maxwi)
{
    /* Initialized data */

    static integer refoct[48]	/* was [4][4][3] */ = { 4,5,6,1,2,3,6,1,1,3,5,
	    6,1,2,4,6,1,3,2,5,3,6,2,5,6,4,2,5,4,1,2,5,1,2,3,4,2,6,3,4,6,5,3,4,
	    5,1,3,4 };
    static integer refcrnr[16]	/* was [4][4] */ = { 1,2,3,1,1,4,5,2,2,4,6,3,
	    3,5,6,4 };
    static integer loop[4] = { 1,2,3,1 };
    static integer facedge[12]	/* was [3][4] */ = { 1,2,4,2,3,6,1,3,5,4,5,6 }
	    ;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer opposite[3], i__, j, k, n, i1, i2, j1, j2, k1, k2, ie, ii, 
	    jj, nr, mr;
    static doublereal qt;
    static integer ip0, ip1, ip2, ip3, ip4, jp1, jp2, ipa, ipb, ipc, irc, neo,
	     nfo, npo, ipt, irt, iend, iipe, iiep, iire, ilbr, inep, iipf, 
	    ipes[6], newp;
    static doublereal xypt[12]	/* was [3][4] */;
    static logical flage;
    static integer medge[6];
    static doublereal minqt[3];
    static integer iendsr, facount[4];

/* ================================================================ */
/* ================================================================ */
/*  Routine refines uniformly elements marked by SplitFlag. */
/*  SplitFlag=.true. is ignored in case of negative sign of Ref2MapMtr. */
/*  Ref2MapMtr is assigned  negative sign at an element if it was split non-uniformly. */

/*  Remarks: */
/*    MapMtr - data array of length 9*nEcoarse */
/*    Ref2MapMtr - data array (lengths on input and output are 9*nE) */
/*    iW(*) - working array of size MaxWi which is at least */
/*            nP + 14 * nE + 3*nR + 4*min(MaxF,4*nF) */
/* ================================================================ */
/* (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --iw;
    --splitflag;
    --ref2mapmtr;
    mapmtr -= 10;
    --lbe;
    --lbf;
    ipe -= 5;
    ipf -= 4;
    xyp -= 4;

    /* Function Body */
/* ================================================================ */
    iire = 1;
    iipe = iire + *ne * 6;
    inep = iipe + (*ne << 2);
    iiep = inep + *np;
    iend = iiep + (*ne << 2);
    liste2r_(np, &nr, ne, &ipe[5], &iw[iire], &iw[inep], &iw[iiep]);
    ilbr = iend;
    iendsr = ilbr + nr;
    iipf = iendsr + (nr << 1);
/* Computing MIN */
    i__1 = *maxf, i__2 = *nf << 2;
    iend = iipf + (min(i__1,i__2) << 2);
    if (iend > *maxwi) {
	errmes_(&c__1001, "localRefinement", "not enough working memory", (
		ftnlen)15, (ftnlen)25);
    }
    npo = *np;
    nfo = *nf;
    neo = *ne;
/* ... mark new points on edges */
    i__1 = nr;
    for (n = 1; n <= i__1; ++n) {
	iw[ilbr + n - 1] = 0;
    }
/* ... mark mother edges */
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	if (splitflag[n] && ref2mapmtr[n] > 0) {
	    for (i__ = 1; i__ <= 6; ++i__) {
		irt = iw[iire + (n - 1) * 6 + i__ - 1];
		iw[ilbr + irt - 1] = 1;
	    }
	}
    }
/* ... count mother edges */
    mr = 0;
    i__1 = nr;
    for (irt = 1; irt <= i__1; ++irt) {
	if (iw[ilbr + irt - 1] > 0) {
	    ++(*np);
	    ++mr;
	    iw[ilbr + irt - 1] = mr;
	}
    }
/* ... assign father nodes to mother edges */
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	k = 0;
	for (i1 = 1; i1 <= 3; ++i1) {
	    for (i2 = i1 + 1; i2 <= 4; ++i2) {
		ip1 = ipe[i1 + (n << 2)];
		ip2 = ipe[i2 + (n << 2)];
		++k;
		irt = iw[iire + (n - 1) * 6 + k - 1];
		irc = iw[ilbr + irt - 1];
		if (irc > 0) {
		    iw[iendsr + (irc - 1 << 1)] = ip1;
		    iw[iendsr + (irc - 1 << 1) + 1] = ip2;
		}
	    }
	}
    }
/* ... add more points so that each face have 0, 1, or 3 new points */
    i__1 = neo;
    for (n = 1; n <= i__1; ++n) {
	k = 0;
/*        count new points for faces of each tet */
	for (j1 = 1; j1 <= 4; ++j1) {
	    facount[j1 - 1] = 0;
	}
	for (i1 = 1; i1 <= 3; ++i1) {
	    for (i2 = i1 + 1; i2 <= 4; ++i2) {
/* loop over edges */
		++k;
		irt = iw[iire + (n - 1) * 6 + k - 1];
		irc = iw[ilbr + irt - 1];
		if (irc > 0) {
/*                 Loop over faces (j1) having k-th edge */
		    for (j1 = 1; j1 <= 4; ++j1) {
			for (j2 = 1; j2 <= 3; ++j2) {
			    if (facedge[j2 + j1 * 3 - 4] == k) {
				++facount[j1 - 1];
			    }
			}
		    }
		}
	    }
	}
/*        Add new point to faces with 2 new points */
	k = 0;
	for (i1 = 1; i1 <= 3; ++i1) {
	    for (i2 = i1 + 1; i2 <= 4; ++i2) {
/* loop over edges */
		++k;
		irt = iw[iire + (n - 1) * 6 + k - 1];
		irc = iw[ilbr + irt - 1];
		if (irc == 0) {
/*                 Loop over faces (j1) having k-th edge */
		    for (j1 = 1; j1 <= 4; ++j1) {
			for (j2 = 1; j2 <= 3; ++j2) {
			    if (facedge[j2 + j1 * 3 - 4] == k && facount[j1 - 
				    1] == 2) {
/*                        add extra mother edges */
				++(*np);
				++mr;
				iw[ilbr + irt - 1] = mr;
				if (k == 1) {
				    ip1 = ipe[(n << 2) + 1];
				    ip2 = ipe[(n << 2) + 2];
				} else if (k == 2) {
				    ip1 = ipe[(n << 2) + 1];
				    ip2 = ipe[(n << 2) + 3];
				} else if (k == 3) {
				    ip1 = ipe[(n << 2) + 1];
				    ip2 = ipe[(n << 2) + 4];
				} else if (k == 4) {
				    ip1 = ipe[(n << 2) + 2];
				    ip2 = ipe[(n << 2) + 3];
				} else if (k == 5) {
				    ip1 = ipe[(n << 2) + 2];
				    ip2 = ipe[(n << 2) + 4];
				} else if (k == 6) {
				    ip1 = ipe[(n << 2) + 3];
				    ip2 = ipe[(n << 2) + 4];
				}
				iw[iendsr + (mr - 1 << 1)] = ip1;
				iw[iendsr + (mr - 1 << 1) + 1] = ip2;
				facount[j1 - 1] = 3;
				for (k1 = 1; k1 <= 4; ++k1) {
				    for (k2 = 1; k2 <= 3; ++k2) {
					if (facedge[k2 + k1 * 3 - 4] == k && 
						k1 != j1) {
					    ++facount[k1 - 1];
					}
				    }
				}
				goto L1;
			    }
			}
		    }
		}
L1:
		;
	    }
	}
    }
/* ... check for faces with 2 new points */
    i__1 = neo;
    for (n = 1; n <= i__1; ++n) {
	k = 0;
/*        count new points for faces of each tet */
	for (j1 = 1; j1 <= 4; ++j1) {
	    facount[j1 - 1] = 0;
	}
	for (i1 = 1; i1 <= 3; ++i1) {
	    for (i2 = i1 + 1; i2 <= 4; ++i2) {
		++k;
		irt = iw[iire + (n - 1) * 6 + k - 1];
		irc = iw[ilbr + irt - 1];
		if (irc > 0) {
/*                 Loop over faces (j1) having k-th edge */
		    for (j1 = 1; j1 <= 4; ++j1) {
			for (j2 = 1; j2 <= 3; ++j2) {
			    if (facedge[j2 + j1 * 3 - 4] == k) {
				++facount[j1 - 1];
			    }
			}
		    }
		}
	    }
	}
	for (j1 = 1; j1 <= 4; ++j1) {
	    if (facount[j1 - 1] == 2) {
		errmes_(&c__1002, "localRefinement", "facount=2", (ftnlen)15, 
			(ftnlen)9);
	    }
/* If this error happen, it implies that algorithm for removing situation */
/* with 2 new points at a face has to be refined! */
	}
    }
    if (*np > *maxp) {
	errmes_(&c__1003, "localRefinement", "local parameter MaxP is small", 
		(ftnlen)15, (ftnlen)29);
    }
/* ... create new points, the point number is nPo plus the edge number */
    i__1 = nr;
    for (irt = 1; irt <= i__1; ++irt) {
	irc = iw[ilbr + irt - 1];
	if (irc > 0) {
	    ipt = npo + irc;
	    ip1 = iw[iendsr + (irc - 1 << 1)];
	    ip2 = iw[iendsr + (irc - 1 << 1) + 1];
	    for (i__ = 1; i__ <= 3; ++i__) {
		xyp[i__ + ipt * 3] = (xyp[i__ + ip1 * 3] + xyp[i__ + ip2 * 3])
			 / 2;
	    }
	}
    }
/* ... mapping P -> E */
    backreferences_(&npo, &neo, &c__4, &c__4, &ipe[5], &iw[inep], &iw[iiep]);
    k = 0;
    i__1 = neo;
    for (n = 1; n <= i__1; ++n) {
	for (i__ = 1; i__ <= 4; ++i__) {
	    ++k;
	    iw[iipe + k - 1] = ipe[i__ + (n << 2)];
	}
    }
/* ... splitting elements */
    i__1 = neo;
    for (n = 1; n <= i__1; ++n) {
/*        count new points */
	newp = 0;
	for (i__ = 1; i__ <= 6; ++i__) {
	    irt = iw[iire + (n - 1) * 6 + i__ - 1];
	    irc = iw[ilbr + irt - 1];
	    if (irc > 0) {
		++newp;
		medge[newp - 1] = irc;
/* index of mother edge */
	    }
	}
	for (i__ = 1; i__ <= 4; ++i__) {
	    ipes[i__ - 1] = ipe[i__ + (n << 2)];
	}
	ipes[4] = lbe[n];
	ipes[5] = ref2mapmtr[n];
	if (newp == 1) {
/*           tet is split into 2 tets,n and nE+1. iP1,iP2 are endpoints of mother edge */
	    ++(*ne);
	    if (*ne > *maxe) {
		errmes_(&c__1006, "localRefinement", "local parameter MaxE i"
			"s small", (ftnlen)15, (ftnlen)29);
	    }
	    for (k = 1; k <= 4; ++k) {
		ipe[k + (*ne << 2)] = ipes[k - 1];
	    }
	    lbe[n] = ipes[4];
	    lbe[*ne] = ipes[4];
	    ref2mapmtr[n] = -dabs(ipes[5]);
	    ref2mapmtr[*ne] = -dabs(ipes[5]);
	    ip1 = iw[iendsr + (medge[0] - 1 << 1)];
	    ip2 = iw[iendsr + (medge[0] - 1 << 1) + 1];
	    for (i__ = 1; i__ <= 4; ++i__) {
		if (ip2 == ipe[i__ + (*ne << 2)]) {
		    ipe[i__ + (*ne << 2)] = npo + medge[0];
		}
		if (ip1 == ipe[i__ + (n << 2)]) {
		    ipe[i__ + (n << 2)] = npo + medge[0];
		}
	    }
	} else if (newp == 2) {
/*           tet is split into 4 tets, n, nE+1,nE+2,nE+3. */
	    *ne += 3;
	    if (*ne > *maxe) {
		errmes_(&c__1006, "localRefinement", "local parameter MaxE i"
			"s small", (ftnlen)15, (ftnlen)29);
	    }
	    for (k = 1; k <= 4; ++k) {
		ipe[k + (*ne - 2 << 2)] = ipes[k - 1];
		ipe[k + (*ne - 1 << 2)] = ipes[k - 1];
		ipe[k + (*ne << 2)] = ipes[k - 1];
	    }
	    lbe[n] = ipes[4];
	    lbe[*ne - 2] = ipes[4];
	    lbe[*ne - 1] = ipes[4];
	    lbe[*ne] = ipes[4];
	    ref2mapmtr[n] = -dabs(ipes[5]);
	    ref2mapmtr[*ne - 2] = -dabs(ipes[5]);
	    ref2mapmtr[*ne - 1] = -dabs(ipes[5]);
	    ref2mapmtr[*ne] = -dabs(ipes[5]);
/*           iP1,iP2,IP3,iP4 are endpoints of mother edge */
	    ip1 = iw[iendsr + (medge[0] - 1 << 1)];
	    ip2 = iw[iendsr + (medge[0] - 1 << 1) + 1];
	    ip3 = iw[iendsr + (medge[1] - 1 << 1)];
	    ip4 = iw[iendsr + (medge[1] - 1 << 1) + 1];
	    for (i__ = 1; i__ <= 4; ++i__) {
		if (ip2 == ipe[i__ + (n << 2)]) {
		    ipe[i__ + (n << 2)] = npo + medge[0];
		}
		if (ip4 == ipe[i__ + (n << 2)]) {
		    ipe[i__ + (n << 2)] = npo + medge[1];
		}
		if (ip2 == ipe[i__ + (*ne - 2 << 2)]) {
		    ipe[i__ + (*ne - 2 << 2)] = npo + medge[0];
		}
		if (ip3 == ipe[i__ + (*ne - 2 << 2)]) {
		    ipe[i__ + (*ne - 2 << 2)] = npo + medge[1];
		}
		if (ip1 == ipe[i__ + (*ne - 1 << 2)]) {
		    ipe[i__ + (*ne - 1 << 2)] = npo + medge[0];
		}
		if (ip4 == ipe[i__ + (*ne - 1 << 2)]) {
		    ipe[i__ + (*ne - 1 << 2)] = npo + medge[1];
		}
		if (ip1 == ipe[i__ + (*ne << 2)]) {
		    ipe[i__ + (*ne << 2)] = npo + medge[0];
		}
		if (ip3 == ipe[i__ + (*ne << 2)]) {
		    ipe[i__ + (*ne << 2)] = npo + medge[1];
		}
	    }
	} else if (newp == 3) {
/*           tet is split into 4 tets, n, nE+1,nE+2,nE+3. */
	    *ne += 3;
	    if (*ne > *maxe) {
		errmes_(&c__1006, "localRefinement", "local parameter MaxE i"
			"s small", (ftnlen)15, (ftnlen)29);
	    }
	    for (i__ = 1; i__ <= 3; ++i__) {
		ip1 = iw[iendsr + (medge[i__ - 1] - 1 << 1)];
		ip2 = iw[iendsr + (medge[i__ - 1] - 1 << 1) + 1];
		for (k = 1; k <= 4; ++k) {
		    if (ip1 == ipes[k - 1] || ip2 == ipes[k - 1]) {
			ipes[k - 1] = 0;
		    }
		}
		jp1 = iw[iendsr + (medge[loop[i__] - 1] - 1 << 1)];
		jp2 = iw[iendsr + (medge[loop[i__] - 1] - 1 << 1) + 1];
		if (ip1 == jp1 || ip1 == jp2) {
		    opposite[i__ - 1] = ip1;
		}
		if (ip2 == jp1 || ip2 == jp2) {
		    opposite[i__ - 1] = ip2;
		}
	    }
	    for (k = 1; k <= 4; ++k) {
		if (ipes[k - 1] != 0) {
		    ip0 = ipes[k - 1];
		}
	    }
	    ipe[(n << 2) + 1] = ip0;
	    ipe[(n << 2) + 2] = npo + medge[0];
	    ipe[(n << 2) + 3] = npo + medge[1];
	    ipe[(n << 2) + 4] = npo + medge[2];
	    lbe[n] = ipes[4];
	    ref2mapmtr[n] = -dabs(ipes[5]);
	    for (k = 1; k <= 3; ++k) {
		ipe[(*ne - k + 1 << 2) + 1] = ip0;
		ipe[(*ne - k + 1 << 2) + 2] = npo + medge[k - 1];
		ipe[(*ne - k + 1 << 2) + 3] = npo + medge[loop[k] - 1];
		ipe[(*ne - k + 1 << 2) + 4] = opposite[k - 1];
		lbe[*ne - k + 1] = ipes[4];
		ref2mapmtr[*ne - k + 1] = -dabs(ipes[5]);
	    }
	} else if (newp == 4) {
	    errmes_(&c__1004, "localRefinement", "four new points in tet", (
		    ftnlen)15, (ftnlen)22);
	} else if (newp == 5) {
	    errmes_(&c__1005, "localRefinement", "five new points in tet", (
		    ftnlen)15, (ftnlen)22);
	} else if (newp == 6) {
/* ... find optimal splitting of inner octahedron (ii=1,2,3) */
	    for (ii = 1; ii <= 3; ++ii) {
		minqt[ii - 1] = 1.;
		for (i__ = 1; i__ <= 4; ++i__) {
		    for (k = 1; k <= 4; ++k) {
			irt = iw[iire + n * 6 + refoct[k + (i__ + (ii << 2) <<
				 2) - 21] - 7];
			irc = iw[ilbr + irt - 1];
			j = npo + irc;
			jj = ref2mapmtr[n];
			if (jj < 0) {
			    jj = -jj;
			}
			fixedmap_(&xyp[j * 3 + 1], &mapmtr[jj * 9 + 1], &xypt[
				k * 3 - 3]);
		    }
		    calreg_(xypt, &xypt[3], &xypt[6], &xypt[9], &qt);
/* Computing MIN */
		    d__1 = qt, d__2 = minqt[ii - 1];
		    minqt[ii - 1] = min(d__1,d__2);
		}
	    }
	    if (minqt[0] >= minqt[1] && minqt[0] >= minqt[2]) {
		ii = 1;
	    } else if (minqt[1] >= minqt[0] && minqt[1] >= minqt[2]) {
		ii = 2;
	    } else if (minqt[2] >= minqt[0] && minqt[2] >= minqt[1]) {
		ii = 3;
	    }
/* ... apply optimal splitting */
	    for (i__ = 1; i__ <= 4; ++i__) {
		++(*ne);
/*  ...  interior tet */
		if (*ne > *maxe) {
		    errmes_(&c__1006, "localRefinement", "local parameter Ma"
			    "xE is small", (ftnlen)15, (ftnlen)29);
		}
		for (k = 1; k <= 4; ++k) {
		    irt = iw[iire + n * 6 + refoct[k + (i__ + (ii << 2) << 2) 
			    - 21] - 7];
		    irc = iw[ilbr + irt - 1];
		    ipe[k + (*ne << 2)] = npo + irc;
		}
		lbe[*ne] = ipes[4];
		ref2mapmtr[*ne] = ipes[5];
/*  ...  vertex-based tet (i<4) or original tet (i=4) */
		if (i__ == 4) {
		    ie = n;
		} else {
		    ++(*ne);
		    if (*ne > *maxe) {
			errmes_(&c__1006, "localRefinement", "local paramete"
				"r MaxE is small", (ftnlen)15, (ftnlen)29);
		    }
		    ie = *ne;
		}
		for (k = 1; k <= 3; ++k) {
		    irt = iw[iire + n * 6 + refcrnr[k + (i__ << 2) - 5] - 7];
		    irc = iw[ilbr + irt - 1];
		    ipe[k + (ie << 2)] = npo + irc;
		}
		ipe[(ie << 2) + 4] = ipes[i__ - 1];
		lbe[ie] = ipes[4];
		ref2mapmtr[ie] = ipes[5];
	    }
	}
    }
/* ... splitting faces */
    *nf = 0;
    i__1 = nfo;
    for (n = 1; n <= i__1; ++n) {
	ip1 = ipf[n * 3 + 1];
	ip2 = ipf[n * 3 + 2];
	ip3 = ipf[n * 3 + 3];
	flage = cmpe_(&ip1, &ip2, &ip3, &iw[iiep], &iw[inep], &c__0, &ie);
/*        count new points at face and fill medge */
	newp = 0;
	for (i__ = 1; i__ <= 6; ++i__) {
	    irt = iw[iire + (ie - 1) * 6 + i__ - 1];
	    irc = iw[ilbr + irt - 1];
	    if (irc > 0) {
		jp1 = iw[iendsr + (irc - 1 << 1)];
		jp2 = iw[iendsr + (irc - 1 << 1) + 1];
		if (check13_(&jp1, &ip1, &ip2, &ip3) && check13_(&jp2, &ip1, &
			ip2, &ip3)) {
		    ++newp;
		    medge[newp - 1] = irc;
		}
	    }
	}
	if (newp == 0) {
	    ++(*nf);
	    if (*nf > *maxf) {
		errmes_(&c__1007, "localRefinement", "local parameter MaxF i"
			"s small", (ftnlen)15, (ftnlen)29);
	    }
	    iw[iipf + (*nf - 1 << 2)] = ip1;
	    iw[iipf + (*nf - 1 << 2) + 1] = ip2;
	    iw[iipf + (*nf - 1 << 2) + 2] = ip3;
	    iw[iipf + (*nf - 1 << 2) + 3] = lbf[n];
	} else if (newp == 1) {
	    *nf += 2;
	    if (*nf > *maxf) {
		errmes_(&c__1007, "localRefinement", "local parameter MaxF i"
			"s small", (ftnlen)15, (ftnlen)29);
	    }
	    jp1 = iw[iendsr + (medge[0] - 1 << 1)];
	    jp2 = iw[iendsr + (medge[0] - 1 << 1) + 1];
	    if (check22_(&ip1, &ip2, &jp1, &jp2)) {
		ipa = ip3;
		ipb = ip1;
		ipc = ip2;
	    }
	    if (check22_(&ip1, &ip3, &jp1, &jp2)) {
		ipa = ip2;
		ipb = ip1;
		ipc = ip3;
	    }
	    if (check22_(&ip2, &ip3, &jp1, &jp2)) {
		ipa = ip1;
		ipb = ip2;
		ipc = ip3;
	    }
	    iw[iipf + (*nf - 2 << 2)] = ipa;
	    iw[iipf + (*nf - 2 << 2) + 1] = ipb;
	    iw[iipf + (*nf - 2 << 2) + 2] = npo + medge[0];
	    iw[iipf + (*nf - 2 << 2) + 3] = lbf[n];
	    iw[iipf + (*nf - 1 << 2)] = ipa;
	    iw[iipf + (*nf - 1 << 2) + 1] = ipc;
	    iw[iipf + (*nf - 1 << 2) + 2] = npo + medge[0];
	    iw[iipf + (*nf - 1 << 2) + 3] = lbf[n];
	} else if (newp == 3) {
/* ... define opposite */
	    for (i__ = 1; i__ <= 3; ++i__) {
		ip1 = iw[iendsr + (medge[i__ - 1] - 1 << 1)];
		ip2 = iw[iendsr + (medge[i__ - 1] - 1 << 1) + 1];
		jp1 = iw[iendsr + (medge[loop[i__] - 1] - 1 << 1)];
		jp2 = iw[iendsr + (medge[loop[i__] - 1] - 1 << 1) + 1];
		if (ip1 == jp1 || ip1 == jp2) {
		    opposite[i__ - 1] = ip1;
		}
		if (ip2 == jp1 || ip2 == jp2) {
		    opposite[i__ - 1] = ip2;
		}
	    }
/* ... */
	    *nf += 4;
	    if (*nf > *maxf) {
		errmes_(&c__1007, "localRefinement", "local parameter MaxF i"
			"s small", (ftnlen)15, (ftnlen)29);
	    }
/*           interior triangle */
	    ipa = npo + medge[0];
	    ipb = npo + medge[1];
	    ipc = npo + medge[2];
	    iw[iipf + (*nf - 4 << 2)] = ipa;
	    iw[iipf + (*nf - 4 << 2) + 1] = ipb;
	    iw[iipf + (*nf - 4 << 2) + 2] = ipc;
	    iw[iipf + (*nf - 4 << 2) + 3] = lbf[n];
	    for (k = 1; k <= 3; ++k) {
/* corner triangles */
		iw[iipf + (*nf - k << 2)] = opposite[k - 1];
		iw[iipf + (*nf - k << 2) + 1] = npo + medge[k - 1];
		iw[iipf + (*nf - k << 2) + 2] = npo + medge[loop[k] - 1];
		iw[iipf + (*nf - k << 2) + 3] = lbf[n];
	    }
	} else {
	    errmes_(&c__1008, "localRefinement", "newP for face  is wrong", (
		    ftnlen)15, (ftnlen)23);
	}
    }
    i__1 = *nf;
    for (n = 1; n <= i__1; ++n) {
	ipf[n * 3 + 1] = iw[iipf + (n - 1 << 2)];
	ipf[n * 3 + 2] = iw[iipf + (n - 1 << 2) + 1];
	ipf[n * 3 + 3] = iw[iipf + (n - 1 << 2) + 2];
	lbf[n] = iw[iipf + (n - 1 << 2) + 3];
    }
    return 0;
} /* localrefinement_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int calreg_(doublereal *xy1, doublereal *xy2, 
	doublereal *xy3, doublereal *xy4, doublereal *qe)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static integer n;
    static doublereal x1, y1, z1, pk, vk, volume;

/* ================================================================ */
/* Computing regularity quality of tetrahedron {xy1, ..., xy4} */
/* ================================================================ */
/* group (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --xy4;
    --xy3;
    --xy2;
    --xy1;

    /* Function Body */
    pk = 0.;
    for (n = 1; n <= 6; ++n) {
	if (n == 1) {
	    x1 = xy1[1] - xy4[1];
	    y1 = xy1[2] - xy4[2];
	    z1 = xy1[3] - xy4[3];
	} else if (n == 2) {
	    x1 = xy2[1] - xy4[1];
	    y1 = xy2[2] - xy4[2];
	    z1 = xy2[3] - xy4[3];
	} else if (n == 3) {
	    x1 = xy3[1] - xy4[1];
	    y1 = xy3[2] - xy4[2];
	    z1 = xy3[3] - xy4[3];
	} else if (n == 4) {
	    x1 = xy1[1] - xy3[1];
	    y1 = xy1[2] - xy3[2];
	    z1 = xy1[3] - xy3[3];
	} else if (n == 5) {
	    x1 = xy2[1] - xy3[1];
	    y1 = xy2[2] - xy3[2];
	    z1 = xy2[3] - xy3[3];
	} else if (n == 6) {
	    x1 = xy1[1] - xy2[1];
	    y1 = xy1[2] - xy2[2];
	    z1 = xy1[3] - xy2[3];
	}
	pk += sqrt(x1 * x1 + y1 * y1 + z1 * z1);
    }
    volume = calvol_(&xy1[1], &xy2[1], &xy3[1], &xy4[1]);
    vk = dabs(volume);
/* Computing 3rd power */
    d__1 = pk;
    *qe = vk * 1832.8208 / (d__1 * (d__1 * d__1));
    return 0;
} /* calreg_ */

/* @f2h@ */ /* Subroutine */ int fixedmap_(doublereal *xyp, doublereal *
	mapmatrix, doublereal *xypt)
{
    static integer i__, j;

    /* Parameter adjustments */
    --xypt;
    mapmatrix -= 4;
    --xyp;

    /* Function Body */
    for (i__ = 1; i__ <= 3; ++i__) {
	xypt[i__] = 0.;
	for (j = 1; j <= 3; ++j) {
	    xypt[i__] += mapmatrix[i__ + j * 3] * xyp[j];
	}
    }
    return 0;
} /* fixedmap_ */

