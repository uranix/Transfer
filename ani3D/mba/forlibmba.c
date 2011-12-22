/* forlibmba.f -- translated by f2c (version 20090411).
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

/* ====================================================================== */
/* @f2h@ */ integer ani_metric_eucl__(doublereal *x, doublereal *y, 
	doublereal *z__, doublereal *metric)
{
    /* System generated locals */
    integer ret_val;

/* ====================================================================== */
/* ====================================================================== */
/*  This routine creates a metric at point (x,y,z). The */
/*  metric is a 3x3 positive definite symmetric tensor: */

/*                M11   M12   M13 */
/*      Metric =  M12   M22   M23 */
/*                M13   M23   M33 */

/*  Only the upper triangular part of Metric must be defined. */
/* ====================================================================== */
/* type of the metric: isotropic/scalar or full tensor */
    /* Parameter adjustments */

    /* Function Body */
    metric[0] = 1.;
    metric[4] = 1.;
    metric[8] = 1.;
    metric[3] = 0.;
    metric[6] = 0.;
    metric[7] = 0.;
    ret_val = 1;
    return ret_val;
} /* ani_metric_eucl__ */

