/* dsort.f -- translated by f2c (version 20090411).
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

/* @f2h@ */ /* Subroutine */ int dsort_(doublereal *dx, integer *ix, integer *
	n, integer *kflag, integer *ierr)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k, l, m;
    static doublereal r__, t;
    static integer ij, il[21], kk, nn, iu[21];
    static doublereal tt;
    static integer itx, ittx;

/* ***BEGIN PROLOGUE  DSORT */
/* ***PURPOSE  Sort an array and optionally make the same interchanges in */
/*            an auxiliary array.  The array may be sorted in increasing */
/*            or decreasing order.  A slightly modified QUICKSORT */
/*            algorithm is used. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  N6A2B */
/* ***TYPE      DOUBLE PRECISION (SSORT-S, DSORT-D, ISORT-I) */
/* ***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/*           Wisniewski, J. A., (SNLA) */
/* ***DESCRIPTION */

/*   DSORT sorts array DX and optionally makes the same interchanges in */
/*   array IX.  The array DX may be sorted in increasing order or */
/*   decreasing order.  A slightly modified quicksort algorithm is used. */

/*   Description of Parameters */
/*      DX - array of values to be sorted   (usually abscissas) */
/*      IX - array to be (optionally) carried along */
/*      N  - number of values in array DX to be sorted */
/*      KFLAG - control parameter */
/*            =  2  means sort DX in increasing order and carry IX along. */
/*            =  1  means sort DX in increasing order (ignoring IX) */
/*            = -1  means sort DX in decreasing order (ignoring IX) */
/*            = -2  means sort DX in decreasing order and carry IX along. */

/*      IERR - the error's indicator */
/*           = 0 - Okay */
/*           != 0 - errors */

/* ***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm */
/*                 for sorting with minimal storage, Communications of */
/*                 the ACM, 12, 3 (1969), pp. 185-187. */
/* ***ROUTINES CALLED  XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   761101  DATE WRITTEN */
/*   761118  Modified to use the Singleton quicksort algorithm.  (JAW) */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891009  Removed unreferenced statement labels.  (WRB) */
/*   891024  Changed category.  (WRB) */
/*   891024  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   901012  Declared all variables; changed X,Y to DX,IX; changed */
/*           code to parallel SSORT. (M. McClain) */
/*   920501  Reformatted the REFERENCES section.  (DWL, WRB) */
/*   920519  Clarified error messages.  (DWL) */
/*   920801  Declarations section rebuilt and code restructured to use */
/*           IF-THEN-ELSE-ENDIF.  (RWC, WRB) */
/* ***END PROLOGUE  DSORT */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. Intrinsic Functions .. */
/* ***FIRST EXECUTABLE STATEMENT  DSORT */
    /* Parameter adjustments */
    --ix;
    --dx;

    /* Function Body */
    *ierr = 0;
    nn = *n;
    if (nn < 1) {
	*ierr = 1201;
	return 0;
    }

    kk = dabs(*kflag);
    if (kk != 1 && kk != 2) {
	*ierr = 1202;
	return 0;
    }

/*     Initializing IX */

/*     DO 5 I=1,NN */
/*        IX(I)=I */
/* 5    CONTINUE */

/*     Alter array DX to get decreasing order if needed */

    if (*kflag <= -1) {
	i__1 = nn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dx[i__] = -dx[i__];
/* L10: */
	}
    }

    if (kk == 2) {
	goto L100;
    }

/*     Sort DX only */

    m = 1;
    i__ = 1;
    j = nn;
    r__ = .375;

L20:
    if (i__ == j) {
	goto L60;
    }
    if (r__ <= .5898437) {
	r__ += .0390625;
    } else {
	r__ += -.21875;
    }

L30:
    k = i__;

/*     Select a central element of the array and save it in location T */

    ij = i__ + (integer) ((j - i__) * r__);
    t = dx[ij];

/*     If first element of array is greater than T, interchange with T */

    if (dx[i__] > t) {
	dx[ij] = dx[i__];
	dx[i__] = t;
	t = dx[ij];
    }
    l = j;

/*     If last element of array is less than than T, interchange with T */

    if (dx[j] < t) {
	dx[ij] = dx[j];
	dx[j] = t;
	t = dx[ij];

/*        If first element of array is greater than T, interchange with T */

	if (dx[i__] > t) {
	    dx[ij] = dx[i__];
	    dx[i__] = t;
	    t = dx[ij];
	}
    }

/*     Find an element in the second half of the array which is smaller */
/*     than T */

L40:
    --l;
    if (dx[l] > t) {
	goto L40;
    }

/*     Find an element in the first half of the array which is greater */
/*     than T */

L50:
    ++k;
    if (dx[k] < t) {
	goto L50;
    }

/*     Interchange these elements */

    if (k <= l) {
	tt = dx[l];
	dx[l] = dx[k];
	dx[k] = tt;
	goto L40;
    }

/*     Save upper and lower subscripts of the array yet to be sorted */

    if (l - i__ > j - k) {
	il[m - 1] = i__;
	iu[m - 1] = l;
	i__ = k;
	++m;
    } else {
	il[m - 1] = k;
	iu[m - 1] = j;
	j = l;
	++m;
    }
    goto L70;

/*     Begin again on another portion of the unsorted array */

L60:
    --m;
    if (m == 0) {
	goto L190;
    }
    i__ = il[m - 1];
    j = iu[m - 1];

L70:
    if (j - i__ >= 1) {
	goto L30;
    }
    if (i__ == 1) {
	goto L20;
    }
    --i__;

L80:
    ++i__;
    if (i__ == j) {
	goto L60;
    }
    t = dx[i__ + 1];
    if (dx[i__] <= t) {
	goto L80;
    }
    k = i__;

L90:
    dx[k + 1] = dx[k];
    --k;
    if (t < dx[k]) {
	goto L90;
    }
    dx[k + 1] = t;
    goto L80;

/*     Sort DX and carry IX along */

L100:
    m = 1;
    i__ = 1;
    j = nn;
    r__ = .375;

L110:
    if (i__ == j) {
	goto L150;
    }
    if (r__ <= .5898437) {
	r__ += .0390625;
    } else {
	r__ += -.21875;
    }

L120:
    k = i__;

/*     Select a central element of the array and save it in location T */

    ij = i__ + (integer) ((j - i__) * r__);
    t = dx[ij];
    itx = ix[ij];

/*     If first element of array is greater than T, interchange with T */

    if (dx[i__] > t) {
	dx[ij] = dx[i__];
	dx[i__] = t;
	t = dx[ij];
	ix[ij] = ix[i__];
	ix[i__] = itx;
	itx = ix[ij];
    }
    l = j;

/*     If last element of array is less than T, interchange with T */

    if (dx[j] < t) {
	dx[ij] = dx[j];
	dx[j] = t;
	t = dx[ij];
	ix[ij] = ix[j];
	ix[j] = itx;
	itx = ix[ij];

/*        If first element of array is greater than T, interchange with T */

	if (dx[i__] > t) {
	    dx[ij] = dx[i__];
	    dx[i__] = t;
	    t = dx[ij];
	    ix[ij] = ix[i__];
	    ix[i__] = itx;
	    itx = ix[ij];
	}
    }

/*     Find an element in the second half of the array which is smaller */
/*     than T */

L130:
    --l;
    if (dx[l] > t) {
	goto L130;
    }

/*     Find an element in the first half of the array which is greater */
/*     than T */

L140:
    ++k;
    if (dx[k] < t) {
	goto L140;
    }

/*     Interchange these elements */

    if (k <= l) {
	tt = dx[l];
	dx[l] = dx[k];
	dx[k] = tt;
	ittx = ix[l];
	ix[l] = ix[k];
	ix[k] = ittx;
	goto L130;
    }

/*     Save upper and lower subscripts of the array yet to be sorted */

    if (l - i__ > j - k) {
	il[m - 1] = i__;
	iu[m - 1] = l;
	i__ = k;
	++m;
    } else {
	il[m - 1] = k;
	iu[m - 1] = j;
	j = l;
	++m;
    }
    goto L160;

/*     Begin again on another portion of the unsorted array */

L150:
    --m;
    if (m == 0) {
	goto L190;
    }
    i__ = il[m - 1];
    j = iu[m - 1];

L160:
    if (j - i__ >= 1) {
	goto L120;
    }
    if (i__ == 1) {
	goto L110;
    }
    --i__;

L170:
    ++i__;
    if (i__ == j) {
	goto L150;
    }
    t = dx[i__ + 1];
    itx = ix[i__ + 1];
    if (dx[i__] <= t) {
	goto L170;
    }
    k = i__;

L180:
    dx[k + 1] = dx[k];
    ix[k + 1] = ix[k];
    --k;
    if (t < dx[k]) {
	goto L180;
    }
    dx[k + 1] = t;
    ix[k + 1] = itx;
    goto L170;

/*     Clean up */

L190:
    if (*kflag <= -1) {
	i__1 = nn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dx[i__] = -dx[i__];
/* L200: */
	}
    }
    return 0;
} /* dsort_ */

