/* minim.f -- translated by f2c (version 20090411).
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
#include "update.h"

#include "makQ.h"
#include "nlnfnc.h"

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int minim_(integer *nu, doublereal *u, 
	doublereal *zz, doublereal *lamda1, doublereal *lamda2, doublereal *
	fmin, doublereal *u1, integer *icnt, doublereal *rmove, logical *
	flagini, logical *flagresult, doublereal *xyp, integer *ipe, 
	doublereal *detg, doublereal *hesp, doublereal *hstar, integer *ips, 
	integer *le, integer *ies, doublereal *xyps, integer *icps, integer *
	ipes, doublereal *detgs, doublereal *hesps, doublereal *qes, integer *
	npw, integer *new__, doublereal *xypw, doublereal *hespw, integer *
	ipew, I_fp metricfunction, logical *flaganalytic, integer *milintrp, 
	integer *mrlintrp, integer *ise, doublereal *rse, integer *icontrol, 
	logical *flagtm, integer *nbad)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__, n;
    static doublereal x, f0, f1;
    static logical flagorient;
    static doublereal qet[1000];
    static integer ios[1000], mbad;
    static doublereal step, xypt[3], xypu[3], lamda, detgt, astep, hespt[6], 
	    dlamda;

/* ================================================================ */
/* group(F) */
/* group (ANI) */
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
/* MINIMIZATION OF THE NONLINEAR FUNCTIONAL FUNC(U) */
/* IN GIVEN DIRECTION  ZZ */
/* ================================================================ */
/* group (F) */
/* ========================================================== */
/* The mesh generation can be controled in a few ways. */
/* The main way is through the metric field. However, some */
/* mesh features may be missing. In order to ensure additional */
/* mesh features, we introduce integer variable status. */
/* The non-zero bits of this variable turn on the */
/* following mesh features: */

/*  bit 01 - (ANIForbidBoundaryElements=1): the final mesh will */
/*           not contain tetrahedra with all vertices on the */
/*           domain boundary. */

/*  bit 02 - (ANIUse2ArmRule=2): each boundary point of the mesh */
/*           may be connected to an inner mesh point with at */
/*           most 2 mesh edges. */

/*           The feature allows more stable implementation of */
/*           Hessian recovering algorithms. */

/*  bit 03 - (ANIFixBoundaryFaces=4): the surface faces will be */
/*           added to the list of fixed faces. */

/*           The feature may be useful for generating multi-block */
/*           meshes. */

/*  bit 04 - (ANIDeleteTemporaryFaces=8): the new faces created by */
/*           the algorithm will be removed from the final mesh. */
/*           The new egdes may include material interfaces. */

/*  bit 05 - (ANIFixSurfacePoints=16): the surface points will be */
/*           added to the list of fixed points. */

/*           The feature may be useful for preserving complex */
/*           geometries in anisotropic metric fields. */

/*  bit 06 - (reserved=32) */

/* ========================================================== */
/* The following bits reflect mesh features discovered by */
/* the algorithm. We keep it for the user disposal. */

/*  bit 07 - (ANIMultiConnectedGeometry=64): the geometry is */
/*           multi-connected. */

/*  bit 08 - (reserved=128) */


/* ========================================================== */
/* The following bits are under development. They will be */
/* automatically removed from the initial value of status. */

/*  bit 09 - (reserved=256) */

/*  bit 10 - (ANITangledMesh=512): the tangled initial mesh */
/*           will be fixed. */

/*  bit 11 - (reserved=1024): */

/* ========================================================== */
/* group (ANI) */
/* ================================================================ */
/* group (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --rse;
    --ise;
    ipew -= 5;
    hespw -= 7;
    xypw -= 4;
    --qes;
    --hesps;
    ipes -= 6;
    xyps -= 4;
    --ies;
    hesp -= 7;
    --detg;
    ipe -= 6;
    xyp -= 4;
    --u1;
    --zz;
    --u;

    /* Function Body */
    *flagresult = FALSE_;
    x = 0.;
    i__1 = *nu;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = zz[i__];
	x += d__1 * d__1;
    }
    x = sqrt(x);
    if (x == 0.) {
	goto L9000;
    }
    i__1 = *nu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	zz[i__] /= x;
    }
    x = *lamda1;
    *lamda1 = min(x,*lamda2);
    *lamda2 = max(x,*lamda2);
    *icnt = 0;
    step = .005;
    astep = 0.;
    dlamda = *lamda2 - *lamda1;
    if (*flagini) {
	f0 = nlnfnc_(nu, &u[1], &xyp[4], &hesp[7], &detg[1], hstar, ips, le, &
		ies[1], &xyps[4], icps, &ipes[6], &hesps[1], detgs, &qes[1], 
		npw, new__, &xypw[4], &hespw[7], &ipew[5], (I_fp)
		metricfunction, flaganalytic, milintrp, mrlintrp, &ise[1], &
		rse[1], icontrol);
/* group (ANI) */
    } else {
	f0 = *fmin;
    }
    copysq_(le, &qes[1], &xyps[4], &hesps[1], detgs, qet, xypt, hespt, &detgt)
	    ;
L10:
    lamda = step * dlamda;
    if (lamda > (1. - astep - step) * dlamda) {
	f1 = f0 + 1.;
    } else {
	i__1 = *nu;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    u1[i__] = u[i__] + lamda * zz[i__];
	}
	++(*icnt);
	f1 = nlnfnc_(nu, &u1[1], &xyp[4], &hesp[7], &detg[1], hstar, ips, le, 
		&ies[1], &xyps[4], icps, &ipes[6], &hesps[1], detgs, &qes[1], 
		npw, new__, &xypw[4], &hespw[7], &ipew[5], (I_fp)
		metricfunction, flaganalytic, milintrp, mrlintrp, &ise[1], &
		rse[1], icontrol);
/* group (ANI) */
/*  ...  checking for inverted elements */
	if (*flagtm) {
	    for (i__ = 1; i__ <= 3; ++i__) {
		xypu[i__ - 1] = xyp[i__ + *ips * 3];
		xyp[i__ + *ips * 3] = xyps[i__ + 3];
	    }
	    i__1 = *le;
	    for (n = 1; n <= i__1; ++n) {
		if (ies[n] >= 0) {
		    updqb_(&n, le, &ies[1], &xyp[4], &ipes[6], &qes[1]);
		}
	    }
	    for (i__ = 1; i__ <= 3; ++i__) {
		xyp[i__ + *ips * 3] = xypu[i__ - 1];
	    }
	    mbad = 0;
	    i__1 = *le;
	    for (n = 1; n <= i__1; ++n) {
		if (qes[n] <= 0.) {
		    ++mbad;
		}
	    }
	    if (mbad > *nbad) {
		f1 = f0 + 1.;
	    }
/*  ...  checking for the tetrahedrons orientation */
	} else {
	    calso_(&xyp[4], &ipe[6], le, &ies[1], ios);
	    chkso_(ips, &xyps[4], &xyp[4], &ipe[6], le, &ies[1], ios, &
		    flagorient);
	    if (! flagorient) {
		f1 = f0 + 1.;
	    }
	}
    }
    if (f1 < f0) {
	f0 = f1;
	i__1 = *nu;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    u[i__] = u1[i__];
	}
	*flagresult = TRUE_;
	copysq_(le, &qes[1], &xyps[4], &hesps[1], detgs, qet, xypt, hespt, &
		detgt);
	astep += step;
	step *= 4.;
    } else {
	step /= 4.;
	if (step < .001) {
	    *fmin = f0;
	    *rmove = astep * dlamda / *lamda2;
	    goto L1000;
	}
    }
    goto L10;
L1000:
    copysq_(le, qet, xypt, hespt, &detgt, &qes[1], &xyps[4], &hesps[1], detgs)
	    ;
L9000:
    return 0;
} /* minim_ */

