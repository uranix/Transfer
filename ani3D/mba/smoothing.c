/* smoothing.f -- translated by f2c (version 20090411).
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
#include "makQ.h"
#include "update.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c__4 = 4;

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int smoothingsol_(integer *np, integer *ne, 
	doublereal *xyp, integer *ipe, doublereal *sol, integer *maxwr, 
	integer *maxwi, doublereal *rw, integer *iw)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, n;
    static doublereal v;
    static integer i1, i2, ie, ip1, ip2, ip3, ip4, iend, iiep, inep, ierr, 
	    isol, ivol, isup;
    static char message[80];

    /* Fortran I/O blocks */
    static icilist io___6 = { 0, message, 0, "(A,I10)", 80, 1 };
    static icilist io___10 = { 0, message, 0, "(A,I10)", 80, 1 };


/* ================================================================ */
/* ================================================================ */
/* Routine smoothes the nodal piecewise linear solution Sol by */
/* averaging over superelements. The superelement for a point P is */
/* defined as the union of tets having vertex P. */
/* ================================================================ */
/* group (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --iw;
    --rw;
    --sol;
    ipe -= 5;
    xyp -= 4;

    /* Function Body */
    inep = 0;
    iiep = inep + *np;
    iend = iiep + (*ne << 2);
    if (iend > *maxwi) {
	ierr = 1001;
	s_wsfi(&io___6);
	do_fio(&c__1, "The approximate size of iW is ", (ftnlen)30);
	do_fio(&c__1, (char *)&iend, (ftnlen)sizeof(integer));
	e_wsfi();
	errmes_(&ierr, "smoothingSol", message, (ftnlen)12, (ftnlen)80);
    }
    ivol = 0;
    isup = ivol + *ne;
    isol = isup + *np;
    iend = isol + *np;
    if (iend > *maxwr) {
	ierr = 1002;
	s_wsfi(&io___10);
	do_fio(&c__1, "The approximate size of rW is ", (ftnlen)30);
	do_fio(&c__1, (char *)&iend, (ftnlen)sizeof(integer));
	e_wsfi();
	errmes_(&ierr, "smoothingSol", message, (ftnlen)12, (ftnlen)80);
    }
/* ... creating an auxiliary structure */
    backreferences_(np, ne, &c__4, &c__4, &ipe[5], &iw[inep + 1], &iw[iiep + 
	    1]);
/* ... computing volumes of elemnts and superelements */
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	ip1 = ipe[(n << 2) + 1];
	ip2 = ipe[(n << 2) + 2];
	ip3 = ipe[(n << 2) + 3];
	ip4 = ipe[(n << 2) + 4];
	v = (d__1 = calvol_(&xyp[ip1 * 3 + 1], &xyp[ip2 * 3 + 1], &xyp[ip3 * 
		3 + 1], &xyp[ip4 * 3 + 1]), dabs(d__1));
	rw[ivol + n] = v;
    }
    i2 = 0;
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	i1 = i2 + 1;
	i2 = iw[inep + n];
	v = 0.;
	i__2 = i2;
	for (i__ = i1; i__ <= i__2; ++i__) {
	    ie = iw[iiep + i__];
	    v += rw[ivol + ie];
	}
	rw[isup + n] = v;
    }
/* ... itegrating the piecewise linear function SOL */
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	rw[isol + n] = 0.;
    }
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	v = 0.;
	for (i__ = 1; i__ <= 4; ++i__) {
	    ip1 = ipe[i__ + (n << 2)];
	    v += sol[ip1];
	}
	v = v * rw[ivol + n] / 4;
	for (i__ = 1; i__ <= 4; ++i__) {
	    ip1 = ipe[i__ + (n << 2)];
	    rw[isol + ip1] += v / rw[isup + ip1];
	}
    }
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	sol[n] = rw[isol + n];
    }
    return 0;
} /* smoothingsol_ */

