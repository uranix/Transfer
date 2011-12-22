/* time_new.f -- translated by f2c (version 20090411).
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

/* ========================================================== */
/* @f2h@ */ doublereal anitime_(real *tmdata)
{
    /* System generated locals */
    real ret_val;

    /* Local variables */

/* ========================================================== */
/* The wrapper for etime(tmdata). ETIME returns a real number */
/* which is the total CPU time used by the program. ETIME */
/* uses TMDATA to return the user time (tmdata(1)) and the */
/* system time (tmdata(2)) separately. */

/* *** Remarks: */
/*        1. the user may set tmdata(1) = tmdata(2) = 0.0 and */
/*           ANItime = 0.0 without breaking the code. */
/* ========================================================== */
    /* Parameter adjustments */
    --tmdata;

    /* Function Body */
    ret_val = (real)etime_(&tmdata[1]);
    return ret_val;
} /* anitime_ */