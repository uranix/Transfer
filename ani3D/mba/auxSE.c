/* auxSE.f -- translated by f2c (version 20090411).
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

/* Table of constant values */

static integer c__1007 = 1007;
static doublereal c_b12 = 1.;

/* ====================================================================== */
/* @f2h@ */ /* Subroutine */ int makse_(integer *ie, integer *icp, integer *
	iep, integer *ipf, integer *ipe, integer *ife, integer *iee, 
	doublereal *qe, integer *maxs, integer *lf, integer *le, integer *ifs,
	 integer *ies, integer *ipfs, integer *ipes, doublereal *qes)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, m, n, k1, if__, ke, ip, ks, iet, keadd, nmaxs;

/* ====================================================================== */
/* ====================================================================== */
/* A superelement around element iE is created. */

/* Remark: The element iE and his face-heighboors should be included */
/*         independently to prevent a situation where vertices of iE */
/*         give a multy-connected (point-based) superelement. */
/*         In addition, we repeat algorithm whenever the "face closure" */
/*         of the superlement touches element iE. */

/* Remark: The 5th column of the IPE is temporary overloaded here. */
/*         The negative values of IPE(5, *) are used for a quick */
/*         search in the list iEs. */

/* Remark: Four elements of ICP are temporary overloded here. */
/*         The negative values of ICP(*) are used for a quick */
/*         search for vertices of iE. */

/* ====================================================================== */
/* group (M) */
/* group (S) */
/* ====================================================================== */
    /* Parameter adjustments */
    --qes;
    ipes -= 6;
    ipfs -= 5;
    --ies;
    --ifs;
    --qe;
    iee -= 5;
    ife -= 5;
    ipe -= 6;
    ipf -= 5;
    --iep;
    --icp;

    /* Function Body */
    ke = 1;
    ies[1] = *ie;
    for (i__ = 1; i__ <= 4; ++i__) {
	iet = iee[i__ + (*ie << 2)];
	if (iet > 0) {
	    ++ke;
	    ies[ke] = iet;
	}
	icp[ipe[i__ + *ie * 5]] = -icp[ipe[i__ + *ie * 5]];
    }
    nmaxs = *maxs - 1;
    for (i__ = 1; i__ <= 4; ++i__) {
	ip = ipe[i__ + *ie * 5];
	maksp_(&ip, &iep[1], &ipe[6], &iee[5], &nmaxs, &keadd, &ies[ke + 1]);
	ke += keadd;
	nmaxs -= keadd;
    }
    *le = 0;
    i__1 = ke;
    for (n = 1; n <= i__1; ++n) {
	iet = ies[n];
	if (ipe[iet * 5 + 5] < 0) {
	    goto L5;
	}
	++(*le);
	ies[*le] = iet;
	for (i__ = 1; i__ <= 5; ++i__) {
	    ipes[i__ + *le * 5] = ipe[i__ + iet * 5];
	}
	ipe[iet * 5 + 5] = -ipe[iet * 5 + 5];
L5:
	;
    }
/* ... 3D extention: adding of elements adajcent by a face */
    nmaxs = *maxs;
    k1 = 1;
    ke = *le;
    ks = 1;
L10:
    m = 0;
    i__1 = ke;
    i__2 = ks;
    for (n = k1; i__2 < 0 ? n >= i__1 : n <= i__1; n += i__2) {
	if (*le >= nmaxs - 4) {
	    goto L1000;
	}
	for (i__ = 1; i__ <= 4; ++i__) {
	    iet = iee[i__ + (ies[n] << 2)];
	    if (iet != 0) {
		if (ipe[iet * 5 + 5] < 0) {
		    goto L15;
		}
		++(*le);
		ies[*le] = iet;
		for (j = 1; j <= 5; ++j) {
		    ipes[j + *le * 5] = ipe[j + iet * 5];
		}
		ipe[iet * 5 + 5] = -ipe[iet * 5 + 5];
/*  ...  checking for intersection with tetrahedron iE */
		for (j = 1; j <= 4; ++j) {
		    if (icp[ipe[j + iet * 5]] < 0) {
			++m;
			ies[*maxs - m + 1] = iet;
			--nmaxs;
			goto L15;
		    }
		}
	    }
L15:
	    ;
	}
    }
    if (m > 0) {
	k1 = *maxs;
	ke = *maxs - m + 1;
	ks = -1;
	goto L10;
    }
/* ... restoring the 5th colunm of IPE */
    i__2 = *le;
    for (n = 1; n <= i__2; ++n) {
	ipe[ies[n] * 5 + 5] = -ipe[ies[n] * 5 + 5];
    }
/* ... restoring the ICP */
    for (i__ = 1; i__ <= 4; ++i__) {
	icp[ipe[i__ + *ie * 5]] = -icp[ipe[i__ + *ie * 5]];
    }
    *lf = 0;
    i__2 = *le;
    for (n = 1; n <= i__2; ++n) {
	ke = ies[n];
	qes[n] = qe[ke];
	for (i__ = 1; i__ <= 4; ++i__) {
	    if__ = ife[i__ + (ke << 2)];
	    if (if__ != 0) {
		i__1 = *lf;
		for (m = 1; m <= i__1; ++m) {
		    if (if__ == ifs[m]) {
			goto L20;
		    }
		}
		++(*lf);
		ifs[*lf] = if__;
		for (j = 1; j <= 4; ++j) {
		    ipfs[j + (*lf << 2)] = ipf[j + (if__ << 2)];
		}
	    }
L20:
	    ;
	}
    }
    return 0;
L1000:
    errmes_(&c__1007, "makSE", "local variable MaxS is small", (ftnlen)5, (
	    ftnlen)28);
    return 0;
} /* makse_ */

/* ====================================================================== */
/* @f2h@ */ /* Subroutine */ int calso_(doublereal *xyp, integer *ipe, 
	integer *le, integer *ies, integer *ios)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal d__;
    static integer n, ie, ip1, ip2, ip3, ip4;

/* ====================================================================== */
/* Oriented volumes of tetrahedra are saved in array iOs. */
/* ====================================================================== */
    /* Parameter adjustments */
    --ios;
    --ies;
    ipe -= 6;
    xyp -= 4;

    /* Function Body */
    i__1 = *le;
    for (n = 1; n <= i__1; ++n) {
	ie = ies[n];
	if (ie <= 0) {
	    goto L10;
	}
	ip1 = ipe[ie * 5 + 1];
	ip2 = ipe[ie * 5 + 2];
	ip3 = ipe[ie * 5 + 3];
	ip4 = ipe[ie * 5 + 4];
	d__ = calvol_(&xyp[ip1 * 3 + 1], &xyp[ip2 * 3 + 1], &xyp[ip3 * 3 + 1],
		 &xyp[ip4 * 3 + 1]);
	ios[n] = (integer) d_sign(&c_b12, &d__);
L10:
	;
    }
    return 0;
} /* calso_ */

/* ====================================================================== */
/* @f2h@ */ /* Subroutine */ int chkso_(integer *ip, doublereal *xyps, 
	doublereal *xyp, integer *ipe, integer *le, integer *ies, integer *
	ios, logical *flag__)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal d__;
    static integer n, ie, ip1, ip2, ip3, ip4, iot;

/* ====================================================================== */
/* flag = TRUE if orientation of new tetrahedra coincides with */
/*             the orientation given in array iOs */
/* ====================================================================== */
/* ====================================================================== */
    /* Parameter adjustments */
    --ios;
    --ies;
    ipe -= 6;
    xyp -= 4;
    --xyps;

    /* Function Body */
    *flag__ = FALSE_;
    i__1 = *le;
    for (n = 1; n <= i__1; ++n) {
	ie = ies[n];
	if (ie <= 0) {
	    goto L10;
	}
	ip1 = ipe[ie * 5 + 1];
	ip2 = ipe[ie * 5 + 2];
	ip3 = ipe[ie * 5 + 3];
	ip4 = ipe[ie * 5 + 4];
	if (ip1 == *ip) {
	    d__ = calvol_(&xyps[1], &xyp[ip2 * 3 + 1], &xyp[ip3 * 3 + 1], &
		    xyp[ip4 * 3 + 1]);
	} else if (ip2 == *ip) {
	    d__ = calvol_(&xyp[ip1 * 3 + 1], &xyps[1], &xyp[ip3 * 3 + 1], &
		    xyp[ip4 * 3 + 1]);
	} else if (ip3 == *ip) {
	    d__ = calvol_(&xyp[ip1 * 3 + 1], &xyp[ip2 * 3 + 1], &xyps[1], &
		    xyp[ip4 * 3 + 1]);
	} else if (ip4 == *ip) {
	    d__ = calvol_(&xyp[ip1 * 3 + 1], &xyp[ip2 * 3 + 1], &xyp[ip3 * 3 
		    + 1], &xyps[1]);
	}
	if (d__ == 0.) {
	    goto L1000;
	}
	iot = (integer) d_sign(&c_b12, &d__);
	if (iot != ios[n]) {
	    goto L1000;
	}
L10:
	;
    }
    *flag__ = TRUE_;
L1000:
    return 0;
} /* chkso_ */

/* ====================================================================== */
/*     Subroutine findSP(lP, iPs, iP, nPs) */
/*     Subroutine findSF(lF, iFs, iF, nFs) */
/* @f2h@ */ /* Subroutine */ int findse_(integer *le, integer *ies, integer *
	ie, integer *nes)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer n;

/* ====================================================================== */
/* Search for index i such that iEs(i) = iE. */
/* Zero is returned when the search fails. */
/* ====================================================================== */
    /* Parameter adjustments */
    --ies;

    /* Function Body */
    *nes = 0;
    i__1 = *le;
    for (n = 1; n <= i__1; ++n) {
	if (ies[n] == *ie) {
	    *nes = n;
	    goto L1000;
	}
    }
L1000:
    return 0;
} /* findse_ */

/* ====================================================================== */
/* @f2h@ */ /* Subroutine */ int copyse_(integer *lfu, integer *leu, integer *
	ifu, integer *ieu, integer *ipfu, integer *ipeu, doublereal *qeu, 
	integer *lf, integer *le, integer *ifs, integer *ies, integer *ipfs, 
	integer *ipes, doublereal *qes)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, n;

/* ====================================================================== */
/* Copy superelement data (u) into superelement data (s). */
/* ====================================================================== */
    /* Parameter adjustments */
    --qes;
    ipes -= 6;
    ipfs -= 5;
    --ies;
    --ifs;
    --qeu;
    ipeu -= 6;
    ipfu -= 5;
    --ieu;
    --ifu;

    /* Function Body */
    *lf = *lfu;
    i__1 = *lf;
    for (n = 1; n <= i__1; ++n) {
	ifs[n] = ifu[n];
	for (i__ = 1; i__ <= 4; ++i__) {
	    ipfs[i__ + (n << 2)] = ipfu[i__ + (n << 2)];
	}
    }
    *le = *leu;
    i__1 = *le;
    for (n = 1; n <= i__1; ++n) {
	ies[n] = ieu[n];
	qes[n] = qeu[n];
	for (i__ = 1; i__ <= 5; ++i__) {
	    ipes[i__ + n * 5] = ipeu[i__ + n * 5];
	}
    }
    return 0;
} /* copyse_ */

/* ====================================================================== */
/* @f2h@ */ /* Subroutine */ int copysemove_(integer *lfu, integer *leu, 
	integer *ifu, integer *ieu, doublereal *qeu, integer *lf, integer *le,
	 integer *ifs, integer *ies, doublereal *qes)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer n;

/* ====================================================================== */
/* A version of copySE optimized for local operation MOVE. */
/* ====================================================================== */
    /* Parameter adjustments */
    --qes;
    --ies;
    --ifs;
    --qeu;
    --ieu;
    --ifu;

    /* Function Body */
    *lf = *lfu;
    i__1 = *lf;
    for (n = 1; n <= i__1; ++n) {
	ifs[n] = ifu[n];
    }
    *le = *leu;
    i__1 = *le;
    for (n = 1; n <= i__1; ++n) {
	ies[n] = ieu[n];
	qes[n] = qeu[n];
    }
    return 0;
} /* copysemove_ */

/* ====================================================================== */
/* @f2h@ */ /* Subroutine */ int copysq_(integer *le, doublereal *qeu, 
	doublereal *xypu, doublereal *hespu, doublereal *detgu, doublereal *
	qes, doublereal *xyps, doublereal *hesps, doublereal *detgs)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, n;

/* ====================================================================== */
/* Routine copies data from u-arrays to s-arrays. */
/* ================================================================ */
    /* Parameter adjustments */
    --hesps;
    --xyps;
    --qes;
    --hespu;
    --xypu;
    --qeu;

    /* Function Body */
    i__1 = *le;
    for (n = 1; n <= i__1; ++n) {
	qes[n] = qeu[n];
    }
    for (i__ = 1; i__ <= 3; ++i__) {
	xyps[i__] = xypu[i__];
    }
    for (i__ = 1; i__ <= 6; ++i__) {
	hesps[i__] = hespu[i__];
    }
    *detgs = *detgu;
    return 0;
} /* copysq_ */

/* ================================================================ */
/* @f2h@ */ logical check1j_(integer *i__, integer *j)
{
    /* System generated locals */
    logical ret_val;

/* ================================================================ */
/* check1j = TRUE if i belongs to the set j(4)}. */
/* ================================================================ */
    /* Parameter adjustments */
    --j;

    /* Function Body */
    ret_val = TRUE_;
    if (*i__ == j[1]) {
	goto L1000;
    }
    if (*i__ == j[2]) {
	goto L1000;
    }
    if (*i__ == j[3]) {
	goto L1000;
    }
    if (*i__ == j[4]) {
	goto L1000;
    }
    ret_val = FALSE_;
L1000:
    return ret_val;
} /* check1j_ */

/* ================================================================ */
/* @f2h@ */ logical check2j_(integer *i1, integer *i2, integer *j)
{
    /* System generated locals */
    logical ret_val;

/* ================================================================ */
/* check2j = TRUE if {i1, i2} is a subset of j(3). */
/* ================================================================ */
    /* Parameter adjustments */
    --j;

    /* Function Body */
    ret_val = FALSE_;
    if (*i1 != j[1] && *i1 != j[2] && *i1 != j[3]) {
	goto L1000;
    }
    if (*i2 != j[1] && *i2 != j[2] && *i2 != j[3]) {
	goto L1000;
    }
    ret_val = TRUE_;
L1000:
    return ret_val;
} /* check2j_ */

/* ================================================================ */
/* @f2h@ */ logical check22_(integer *i1, integer *i2, integer *j1, integer *
	j2)
{
    /* System generated locals */
    logical ret_val;

/* ================================================================ */
/* check22 = TRUE if pair {i1, i2} coinsides with {j1, j2}. */
/* ================================================================ */
    ret_val = FALSE_;
    if (*i1 != *j1 && *i1 != *j2) {
	goto L1000;
    }
    if (*i2 != *j1 && *i2 != *j2) {
	goto L1000;
    }
    ret_val = TRUE_;
L1000:
    return ret_val;
} /* check22_ */

/* ================================================================ */
/* @f2h@ */ logical check3j_(integer *i1, integer *i2, integer *i3, integer *
	j)
{
    /* System generated locals */
    logical ret_val;

/* ================================================================ */
/* check3j = TRUE if {i1, i2, i3} is a subset of j(4). */
/* ================================================================ */
    /* Parameter adjustments */
    --j;

    /* Function Body */
    if (*i1 == j[1] || *i1 == j[2] || *i1 == j[3] || *i1 == j[4]) {
	if (*i2 == j[1] || *i2 == j[2] || *i2 == j[3] || *i2 == j[4]) {
	    if (*i3 == j[1] || *i3 == j[2] || *i3 == j[3] || *i3 == j[4]) {
		ret_val = TRUE_;
		goto L1000;
	    }
	}
    }
    ret_val = FALSE_;
L1000:
    return ret_val;
} /* check3j_ */

/* ================================================================ */
/* @f2h@ */ logical check33_(integer *i1, integer *i2, integer *i3, integer *
	j1, integer *j2, integer *j3)
{
    /* System generated locals */
    logical ret_val;

/* ================================================================ */
/* check33 = TRUE if triple {i1, i2, i3} coinsides with {j1, j2, j3} */
/* ================================================================ */
    if (*i1 == *j1 || *i1 == *j2 || *i1 == *j3) {
	if (*i2 == *j1 || *i2 == *j2 || *i2 == *j3) {
	    if (*i3 == *j1 || *i3 == *j2 || *i3 == *j3) {
		ret_val = TRUE_;
		goto L1000;
	    }
	}
    }
    ret_val = FALSE_;
L1000:
    return ret_val;
} /* check33_ */

/* ================================================================ */
/* @f2h@ */ logical check13_(integer *i1, integer *j1, integer *j2, integer *
	j3)
{
    /* System generated locals */
    logical ret_val;

/* ================================================================ */
/* check13 = TRUE if i1 belongs to the set {j1, j2, j3}. */
/* ================================================================ */
    ret_val = TRUE_;
    if (*i1 == *j1) {
	goto L1000;
    }
    if (*i1 == *j2) {
	goto L1000;
    }
    if (*i1 == *j3) {
	goto L1000;
    }
    ret_val = FALSE_;
L1000:
    return ret_val;
} /* check13_ */

/* ================================================================ */
/* @f2h@ */ logical check14_(integer *i1, integer *j1, integer *j2, integer *
	j3, integer *j4)
{
    /* System generated locals */
    logical ret_val;

/* ================================================================ */
/* check14 = TRUE if i1 belongs to the set {j1, j2, j3, j4}. */
/* ================================================================ */
    ret_val = TRUE_;
    if (*i1 == *j1) {
	goto L1000;
    }
    if (*i1 == *j2) {
	goto L1000;
    }
    if (*i1 == *j3) {
	goto L1000;
    }
    if (*i1 == *j4) {
	goto L1000;
    }
    ret_val = FALSE_;
L1000:
    return ret_val;
} /* check14_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int swapdd_(doublereal *d1, doublereal *d2)
{
    static doublereal d__;

/* ================================================================ */
/* Routine swaps two real*8 numbers. */
/* ================================================================ */
    d__ = *d1;
    *d1 = *d2;
    *d2 = d__;
    return 0;
} /* swapdd_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int swapii_(integer *i1, integer *i2)
{
    static integer i__;

/* ================================================================ */
/* Routine swaps two integer numbers. */
/* ================================================================ */
    i__ = *i1;
    *i1 = *i2;
    *i2 = i__;
    return 0;
} /* swapii_ */

