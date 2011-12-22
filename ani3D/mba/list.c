/* list.f -- translated by f2c (version 20090411).
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
#include "isnan.h"

/* Table of constant values */

static integer c__5101 = 5101;
static integer c__5111 = 5111;
static integer c__5102 = 5102;
static integer c__5103 = 5103;
static integer c__5104 = 5104;
static integer c__5113 = 5113;
static integer c__5005 = 5005;

/* ================================================================ */
/* The routines below work with an ordered list of reals. */
/* The available operations are: */
/*     (a) intialize the list */
/*     (b) add new element in the list */
/*     (c) update an element value */
/*     (d) exclude an element from the list */
/*     (e) careful check the list data */

/* *** Remarsk: */
/*         1. This is the new version of list routines. It requires */
/*            double memory for L2E (compared to the old version), */
/*            and use iCntl(4) instead of nstep. Besides, it uses */
/*            iCntl as the input data to be defined outside: */

/*            iCntl(1) is typical interval length, we recommEnd sqrt(nEStar) */
/*            iCntl(2) is rank of interval length belonging to */
/*                     [iCntl(1)-iCntl(2), iCntl(1)+iCntl(2)] */
/*            iCntl(3) is the pointer to the middle of L2E */
/*            iCntl(4) is the output channel in case of debugging */
/*                     (=0 for no debugging) */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int lstmak_(integer *nestar, integer *ne, 
	integer *l1e, integer *l2e, integer *nl2, integer *icntl, integer *
	lhol)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer n, ipos2, nstep;

/* ================================================================ */
/* ================================================================ */
    /* Parameter adjustments */
    --lhol;
    --icntl;
    --l2e;
    l1e -= 3;

    /* Function Body */
    nstep = icntl[1];
    ipos2 = icntl[3];
    l1e[(l2e[1] << 1) + 1] = 0;
    l1e[(l2e[*ne] << 1) + 2] = 0;
    i__1 = *ne - 1;
    for (n = 1; n <= i__1; ++n) {
	l1e[(l2e[n] << 1) + 2] = l2e[n + 1];
    }
    i__1 = *ne;
    for (n = 2; n <= i__1; ++n) {
	l1e[(l2e[n] << 1) + 1] = l2e[n - 1];
    }
    *nl2 = 0;
    i__1 = *ne;
    i__2 = nstep;
    for (n = 1; i__2 < 0 ? n >= i__1 : n <= i__1; n += i__2) {
	++(*nl2);
	l2e[*nl2] = l2e[n];
	l2e[ipos2 + *nl2] = nstep;
    }
    l2e[ipos2 + *nl2] = *ne - nstep * (*nl2 - 1);
    lhol[1] = 0;
    return 0;
} /* lstmak_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int lstupd_(integer *ne, integer *l1e, integer *
	nl2, integer *l2e, integer *icntl, doublereal *qe, integer *ie, 
	doublereal *qie)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k, ia, ib, nthr, ipos2, ilast, nstep;
    static logical flagfp;
    static integer iprevl1, iprevl2, inextl1;

/* ================================================================ */
/*     Input:  iE,qiE */
/*     Output: Updated qE and L1E,L2E */
/* ================================================================ */
/* ================================================================ */
    /* Parameter adjustments */
    --qe;
    --icntl;
    --l2e;
    l1e -= 3;

    /* Function Body */
    nstep = icntl[1];
    nthr = icntl[2];
    ipos2 = icntl[3];
    if (l1e[(*ie << 1) + 1] == 0 && l1e[(*ie << 1) + 2] == 0) {
	errmes_(&c__5101, "lstUpd", "no element in the list", (ftnlen)6, (
		ftnlen)22);
    }
    fpcheck_(qie, &flagfp);
    if (flagfp) {
	errmes_(&c__5111, "lstUpd", "bad input value(NAN/INF)", (ftnlen)6, (
		ftnlen)24);
    }
    iprevl2 = prevl2ie_(nl2, &l2e[1], &qe[1], ie, &l1e[3]);
    if (l1e[(*ie << 1) + 1] != 0) {
	l1e[(l1e[(*ie << 1) + 1] << 1) + 2] = l1e[(*ie << 1) + 2];
    } else {
	l1e[(l1e[(*ie << 1) + 2] << 1) + 1] = 0;
    }
    if (l1e[(*ie << 1) + 2] != 0) {
	l1e[(l1e[(*ie << 1) + 2] << 1) + 1] = l1e[(*ie << 1) + 1];
    } else {
	l1e[(l1e[(*ie << 1) + 1] << 1) + 2] = 0;
    }
    --(*ne);
    if (iprevl2 == 0) {
	ia = 1;
	ib = 1;
    } else {
	ib = iprevl2;
	if (l2e[iprevl2] == *ie) {
	    ia = 0;
	} else {
	    ia = 1;
	}
	if (ib < *nl2) {
	    if (*ie == l2e[ib + 1]) {
		++ib;
	    }
	}
    }
    if (l2e[ipos2 + ib] - 1 >= nstep - nthr) {
	--l2e[ipos2 + ib];
	if (l2e[ib] == *ie) {
	    l2e[ib] = l1e[(l2e[ib] << 1) + 2];
	}
    } else {
	i__1 = *nl2;
	for (i__ = iprevl2 + ia; i__ <= i__1; ++i__) {
	    l2e[i__] = l1e[(l2e[i__] << 1) + 2];
	}
/*  ...    the last one may be just positive !!! */
	if (l2e[ipos2 + *nl2] - 1 >= 1) {
	    --l2e[ipos2 + *nl2];
	} else {
	    if (l2e[*nl2] != 0) {
		errmes_(&c__5102, "lstUpd", "L2E(nL2) must be 0 here", (
			ftnlen)6, (ftnlen)23);
	    }
	    --(*nl2);
	}
    }
    l1e[(*ie << 1) + 1] = 0;
    l1e[(*ie << 1) + 2] = 0;
/* ... add new element */
    iprevl2 = prevl2_(nl2, &l2e[1], &qe[1], qie);
    if (iprevl2 == 0) {
	iprevl1 = 0;
    } else {
	iprevl1 = prevl1_(&l1e[3], &l2e[1], &iprevl2, &qe[1], qie);
    }
    ++(*ne);
    qe[*ie] = *qie;
    if (iprevl1 != 0) {
	inextl1 = l1e[(iprevl1 << 1) + 2];
    } else {
	inextl1 = l2e[1];
    }
    if (inextl1 != 0) {
	l1e[(inextl1 << 1) + 1] = *ie;
    }
    if (iprevl1 != 0) {
	l1e[(iprevl1 << 1) + 2] = *ie;
    }
    l1e[(*ie << 1) + 1] = iprevl1;
    l1e[(*ie << 1) + 2] = inextl1;
    ib = iprevl2;
    if (iprevl2 == 0) {
	ib = 1;
    }
    if (l2e[ipos2 + ib] + 1 <= nstep + nthr) {
	++l2e[ipos2 + ib];
	if (iprevl2 == 0) {
	    l2e[ib] = l1e[(l2e[ib] << 1) + 1];
	}
    } else {
	i__1 = *nl2;
	for (i__ = iprevl2 + 1; i__ <= i__1; ++i__) {
	    l2e[i__] = l1e[(l2e[i__] << 1) + 1];
	}
	if (l2e[ipos2 + *nl2] + 1 <= nstep + nthr) {
	    ++l2e[ipos2 + *nl2];
	} else {
	    ilast = l2e[*nl2];
	    while(TRUE_) {
		k = l1e[(ilast << 1) + 2];
		if (k == 0) {
		    ++(*nl2);
		    l2e[*nl2] = ilast;
		    l2e[ipos2 + *nl2] = 1;
		    goto L1;
		}
		ilast = k;
	    }
	}
    }
L1:
    return 0;
} /* lstupd_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int lstdel_(integer *ne, integer *l1e, integer *
	nl2, integer *l2e, integer *icntl, integer *lhol, doublereal *qe, 
	integer *ie)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, ia, ib, nthr, ipos2, nstep, iprevl2;

/* ================================================================ */
/*     Input:  iE */
/*     Output: Updated L1E,L2E,nE,nL2,LHol */
/* ================================================================ */
/* ================================================================ */
    /* Parameter adjustments */
    --qe;
    --lhol;
    --icntl;
    --l2e;
    l1e -= 3;

    /* Function Body */
    nstep = icntl[1];
    nthr = icntl[2];
    ipos2 = icntl[3];
    if (l1e[(*ie << 1) + 1] == 0 && l1e[(*ie << 1) + 2] == 0) {
	errmes_(&c__5103, "lstDel", "no element in the list", (ftnlen)6, (
		ftnlen)22);
    }
    iprevl2 = prevl2ie_(nl2, &l2e[1], &qe[1], ie, &l1e[3]);
    if (l1e[(*ie << 1) + 1] != 0) {
	l1e[(l1e[(*ie << 1) + 1] << 1) + 2] = l1e[(*ie << 1) + 2];
    } else {
	l1e[(l1e[(*ie << 1) + 2] << 1) + 1] = 0;
    }
    if (l1e[(*ie << 1) + 2] != 0) {
	l1e[(l1e[(*ie << 1) + 2] << 1) + 1] = l1e[(*ie << 1) + 1];
    } else {
	l1e[(l1e[(*ie << 1) + 1] << 1) + 2] = 0;
    }
    --(*ne);
    if (iprevl2 == 0) {
	ia = 1;
	ib = 1;
    } else {
	ib = iprevl2;
	if (l2e[iprevl2] == *ie) {
	    ia = 0;
	} else {
	    ia = 1;
	}
	if (ib < *nl2) {
	    if (*ie == l2e[ib + 1]) {
		++ib;
	    }
	}
    }
    if (l2e[ipos2 + ib] - 1 >= nstep - nthr) {
	--l2e[ipos2 + ib];
	if (l2e[ib] == *ie) {
	    l2e[ib] = l1e[(l2e[ib] << 1) + 2];
	}
    } else {
	i__1 = *nl2;
	for (i__ = iprevl2 + ia; i__ <= i__1; ++i__) {
	    l2e[i__] = l1e[(l2e[i__] << 1) + 2];
	}
/*  ...   the last one may be just positive !!! */
	if (l2e[ipos2 + *nl2] - 1 >= 1) {
	    --l2e[ipos2 + *nl2];
	} else {
	    if (l2e[*nl2] != 0) {
		errmes_(&c__5104, "lstDel", "L2E(nL2) must be 0 here", (
			ftnlen)6, (ftnlen)23);
	    }
	    --(*nl2);
	}
    }
    l1e[(*ie << 1) + 1] = 0;
    l1e[(*ie << 1) + 2] = 0;
    ++lhol[1];
    lhol[lhol[1] + 1] = *ie;
    return 0;
} /* lstdel_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int lstadd_(integer *ne, integer *l1e, integer *
	nl2, integer *l2e, integer *icntl, integer *lhol, doublereal *qe, 
	doublereal *qie, integer *ie)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k, ib, nthr, ipos2, ilast, nstep;
    static logical flagfp;
    static integer iprevl1, iprevl2, inextl1;

/* ================================================================ */
/*     Input:  qiE,LHol */
/*     Output: iE, updated L1E,L2E,nE,nL2,LHol,qe */
/* ================================================================ */
/* ================================================================ */
    /* Parameter adjustments */
    --qe;
    --lhol;
    --icntl;
    --l2e;
    l1e -= 3;

    /* Function Body */
    nstep = icntl[1];
    nthr = icntl[2];
    ipos2 = icntl[3];
    fpcheck_(qie, &flagfp);
    if (flagfp) {
	errmes_(&c__5113, "lstAdd", "bad input value(NAN/INF)", (ftnlen)6, (
		ftnlen)24);
    }
    iprevl2 = prevl2_(nl2, &l2e[1], &qe[1], qie);
    if (iprevl2 == 0) {
	iprevl1 = 0;
    } else {
	iprevl1 = prevl1_(&l1e[3], &l2e[1], &iprevl2, &qe[1], qie);
    }
    ++(*ne);
    if (lhol[1] == 0) {
	*ie = *ne;
    } else {
	*ie = lhol[lhol[1] + 1];
	--lhol[1];
    }
    qe[*ie] = *qie;
    if (iprevl1 != 0) {
	inextl1 = l1e[(iprevl1 << 1) + 2];
    } else {
	inextl1 = l2e[1];
    }
    if (inextl1 != 0) {
	l1e[(inextl1 << 1) + 1] = *ie;
    }
    if (iprevl1 != 0) {
	l1e[(iprevl1 << 1) + 2] = *ie;
    }
    l1e[(*ie << 1) + 1] = iprevl1;
    l1e[(*ie << 1) + 2] = inextl1;
    ib = iprevl2;
    if (iprevl2 == 0) {
	ib = 1;
    }
    if (l2e[ipos2 + ib] + 1 <= nstep + nthr) {
	++l2e[ipos2 + ib];
	if (iprevl2 == 0) {
	    l2e[ib] = l1e[(l2e[ib] << 1) + 1];
	}
    } else {
	i__1 = *nl2;
	for (i__ = iprevl2 + 1; i__ <= i__1; ++i__) {
	    l2e[i__] = l1e[(l2e[i__] << 1) + 1];
	}
	if (l2e[ipos2 + *nl2] + 1 <= nstep + nthr) {
	    ++l2e[ipos2 + *nl2];
	} else {
	    ilast = l2e[*nl2];
	    while(TRUE_) {
		k = l1e[(ilast << 1) + 2];
		if (k == 0) {
		    ++(*nl2);
		    l2e[*nl2] = ilast;
		    l2e[ipos2 + *nl2] = 1;
		    goto L1;
		}
		ilast = k;
	    }
	}
    }
L1:
    return 0;
} /* lstadd_ */

/* ================================================================ */
/* @f2h@ */ integer prevl2_(integer *nl2, integer *l2e, doublereal *qe, 
	doublereal *qie)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer i1, i2, i3;

/* ================================================================ */
/* ================================================================ */
    /* Parameter adjustments */
    --qe;
    --l2e;

    /* Function Body */
    if (*qie < qe[l2e[1]]) {
	ret_val = 0;
	return ret_val;
    }
    if (*qie >= qe[l2e[*nl2]]) {
	ret_val = *nl2;
	return ret_val;
    }
    i1 = 1;
    i2 = *nl2 - 1;
    while(TRUE_) {
	i3 = (i1 + i2) / 2;
	if (i1 == i2) {
	    ret_val = i1;
	    return ret_val;
	}
	if (i1 == i2 - 1) {
	    if (*qie < qe[l2e[i2]]) {
		ret_val = i1;
	    } else {
		ret_val = i2;
	    }
	    return ret_val;
	}
	if (*qie < qe[l2e[i3]]) {
	    i2 = i3;
	} else {
	    i1 = i3;
	}
    }
    return ret_val;
} /* prevl2_ */

/* ================================================================ */
/* @f2h@ */ integer prevl2ie_(integer *nl2, integer *l2e, doublereal *qe, 
	integer *ie, integer *l1e)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer i__, j, k, i1, i2, i3, ii;
    static doublereal qie;

/* ================================================================ */
/* ================================================================ */
    /* Parameter adjustments */
    l1e -= 3;
    --qe;
    --l2e;

    /* Function Body */
    if (l2e[1] == *ie) {
	ret_val = 0;
	return ret_val;
    }
    qie = qe[*ie];
    if (qie < qe[l2e[1]]) {
	goto L55;
    }
    if (qie >= qe[l2e[*nl2]]) {
	goto L55;
    }
    i1 = 1;
    i2 = *nl2 - 1;
    while(TRUE_) {
	i3 = (i1 + i2) / 2;
	if (i1 == i2) {
	    i__ = i1 + 1;
	    goto L22;
	}
	if (i1 == i2 - 1) {
	    if (qie < qe[l2e[i2]]) {
		i__ = i1 + 1;
	    } else {
		i__ = i2 + 1;
	    }
	    goto L22;
	}
	if (qie < qe[l2e[i3]]) {
	    i2 = i3;
	} else {
	    i1 = i3;
	}
    }
L22:
    for (ii = i__; ii >= 2; --ii) {
	j = l2e[ii];
	while(TRUE_) {
	    if (j == *ie) {
		goto L33;
	    }
	    k = l1e[(j << 1) + 1];
	    if (k == l2e[ii - 1]) {
		goto L44;
	    }
	    j = k;
	}
L33:
	ret_val = ii - 1;
	return ret_val;
L44:
	;
    }
L55:
    ret_val = *nl2;
    if (qie == qe[l2e[*nl2]]) {
	for (ii = *nl2; ii >= 2; --ii) {
	    j = l2e[ii];
	    while(TRUE_) {
		if (j == *ie) {
		    goto L888;
		}
		k = l1e[(j << 1) + 1];
		if (k == l2e[ii - 1]) {
		    goto L999;
		}
		j = k;
	    }
L888:
	    ret_val = ii - 1;
	    return ret_val;
L999:
	    ;
	}
    }
    return ret_val;
} /* prevl2ie_ */

/* ================================================================ */
/* @f2h@ */ integer prevl1_(integer *l1e, integer *l2e, integer *iprevl2, 
	doublereal *qe, doublereal *qie)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer j, k;

/* ================================================================ */
/* ================================================================ */
    /* Parameter adjustments */
    --qe;
    --l2e;
    l1e -= 3;

    /* Function Body */
    if (*iprevl2 <= 0) {
	errmes_(&c__5005, "PrevL1", "wrong iPrevL2", (ftnlen)6, (ftnlen)13);
    }
    j = l2e[*iprevl2];
    while(TRUE_) {
	if (*qie < qe[j]) {
	    ret_val = l1e[(j << 1) + 1];
	    return ret_val;
	}
	k = l1e[(j << 1) + 2];
	if (k == 0) {
	    ret_val = j;
	    return ret_val;
	}
	j = k;
    }
    return ret_val;
} /* prevl1_ */

