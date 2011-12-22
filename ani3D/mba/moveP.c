/* moveP.f -- translated by f2c (version 20090411).
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

#include "movep.h"
#include "nlnfnc.h"
#include "minim.h"
#include "list.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c__128 = 128;
static integer c__512 = 512;
static integer c__2 = 2;
static integer c__8 = 8;
static integer c__64 = 64;
static doublereal c_b11 = 1.;
static integer c__3 = 3;
static logical c_false = FALSE_;

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int movep_(integer *iwp, integer *iwe, integer *
	ne, doublereal *xyp, integer *ipf, integer *ipe, doublereal *hstar, 
	integer *icp, integer *ife, integer *l1e, integer *l2e, integer *nl2, 
	integer *nstep, integer *status, doublereal *hesp, doublereal *
	rquality, doublereal *detg, doublereal *qe, I_fp metricfunction, 
	logical *flaganalytic, integer *lfu, integer *leu, integer *ifu, 
	integer *ieu, integer *ipfs, integer *ipes, doublereal *qeu, integer *
	npw, integer *new__, doublereal *xypw, doublereal *hespw, integer *
	ipew, integer *milintrp, integer *mrlintrp, integer *ise, doublereal *
	rse, integer *icontrol, doublereal *rmove, logical *flag__)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j, m, n;
    static doublereal q[4], s, u[3], v, z__[9]	/* was [3][3] */;
    static integer i1, i2, i3, i4;
    static doublereal t1, u1[3], t2;
    static integer ie, le, lf, ip[5], nu, nz, lc1, ip1, ip2, ip3, ip4, ics[
	    1000], ifs[1000], ies[1000];
    static doublereal qes[1000];
    static integer ipt[3], ips[1000], irs[1000], nbad, iref[4];
    static doublereal heit, hmin, hmax;
    static integer icps, icnt;
    static doublereal qmin, xyps[9]	/* was [3][3] */, xypt[3];
    static integer icrab;
    static doublereal detgs, hesps[6], xstep[3];
    static logical flagtm;

/* ================================================================ */
/* group (M) */
/* group (Q) */
/* group (S) */
/* !!! &           lFu, lEu, iFu, iEu, IPFu, IPEu, qEu, */
/* group (W) */
/* ================================================================ */
/* Remark: the routine was optimized (reference !!!): */
/*         IPFu => IPFs  since IPFu is not changed */
/*         IPEu => IPEs  since IPEu is not changed */
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
/* ========================================================== */
/* The local operations with grid: */
/*   F2E   : face is replaced by an edge. Two elements are */
/*           replaced by three elements. */

/*   E2F   : reverse to F2E operation. It works if 3 elements */
/*           have the common edge. */

/*   SWAP  : a generalization of E2F. The edge is destroyed */
/*           and a hole is filled with new tetrahedra. */

/*   DELET : delete a point. */

/*   MOVE  : movement of a point. */

/*   SPLTE : the elements is splitted into 4 elements by an */
/*           internal point. The operation is used only to */
/*           make grid to satisfy some FE restrictions. */
/* ========================================================== */
/* ...  available operations and their order (from 1 to 7 ONLY!) */
/* ...  number of operations */
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
/* ================================================================ */
/* group (M) */
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
/* group (Q) */
/* group (S) */
/* !!!  Integer iFu(*), iEu(*), IPFu(4, *), IPEu(5, *) */
/* group (W) */
/* group (Flag) */
/* ================================================================ */
/* group (Local Functions) */
/* group (Local variables) */
/* ... for nonlinear minimization procedure */
/* ================================================================ */
    /* Parameter adjustments */
    --rse;
    --ise;
    ipew -= 5;
    hespw -= 7;
    xypw -= 4;
    --qeu;
    ipes -= 6;
    ipfs -= 5;
    --ieu;
    --ifu;
    --qe;
    --detg;
    hesp -= 7;
    --nstep;
    --l2e;
    l1e -= 3;
    ife -= 5;
    --icp;
    ipe -= 6;
    ipf -= 5;
    xyp -= 4;

    /* Function Body */
    *flag__ = FALSE_;
    copysemove_(lfu, leu, &ifu[1], &ieu[1], &qeu[1], &lf, &le, ifs, ies, qes);
    iref[0] = 1;
    iref[1] = 2;
    iref[2] = 3;
    iref[3] = 1;
    ip[0] = 1;
    ip[1] = 2;
    ip[2] = 3;
    ip[3] = 4;
    ip[4] = 1;
    ip1 = ipe[*iwp + *iwe * 5];
    icps = icp[ip1];
/* ... checking the case when moving is impossible */
    if (ifxnode_(&icps, &c__1)) {
	goto L1000;
    }
    if (ifxnode_(&icps, &c__128)) {
	goto L1000;
    }
/* ... checking the number of inverted elements */
    flagtm = ifxnode_(status, &c__512);
    if (flagtm) {
	nbad = 0;
	i__1 = le;
	for (n = 1; n <= i__1; ++n) {
	    if (qes[n - 1] <= 0.) {
		++nbad;
	    }
	}
    }
    flagtm = flagtm && nbad > 0;
    heit = 1e12;
    i__1 = le;
    for (n = 1; n <= i__1; ++n) {
	ie = ies[n - 1];
	for (i1 = 1; i1 <= 4; ++i1) {
	    if (ipe[i1 + ie * 5] == ip1) {
		i2 = ip[i1];
		i3 = ip[i2];
		i4 = ip[i3];
		ip2 = ipe[i2 + ie * 5];
		ip3 = ipe[i3 + ie * 5];
		ip4 = ipe[i4 + ie * 5];
		v = calvol_(&xyp[ip1 * 3 + 1], &xyp[ip2 * 3 + 1], &xyp[ip3 * 
			3 + 1], &xyp[ip4 * 3 + 1]);
		s = calsqr_(&xyp[ip2 * 3 + 1], &xyp[ip3 * 3 + 1], &xyp[ip4 * 
			3 + 1]);
/* Computing MIN */
		d__1 = heit, d__2 = dabs(v) * 6. / s;
		heit = min(d__1,d__2);
		goto L10;
	    }
	}
	if (! flagtm) {
	    ies[n - 1] = -ie;
	}
L10:
	;
    }
/* ... calculating a search direction */
    nz = 0;
    lc1 = 0;
    if (ifxnode_(&icps, &c__2)) {
	nz = 1;
	i__1 = lf;
	for (n = 1; n <= i__1; ++n) {
	    ipt[0] = ipfs[(n << 2) + 1];
	    ipt[1] = ipfs[(n << 2) + 2];
	    ipt[2] = ipfs[(n << 2) + 3];
	    for (i1 = 1; i1 <= 3; ++i1) {
		if (ip1 == ipt[i1 - 1]) {
		    ++lc1;
		    ics[lc1 - 1] = n;
		    ips[lc1 - 1] = -1;
		    irs[lc1 - 1] = ipf[(ifs[n - 1] << 2) + 4];
		    i2 = i1;
		    for (j = 1; j <= 2; ++j) {
			i2 = iref[i2];
			clrsr_(&ip1, &ipt[i2 - 1], &icp[1], &ipf[5], &ife[5], 
				&lf, ifs, &le, ies, &icrab);
			if (icrab == icps) {
			    ips[lc1 - 1] = ipt[i2 - 1];
			    for (i__ = 1; i__ <= 3; ++i__) {
				z__[i__ - 1] = xyp[i__ + ipt[i2 - 1] * 3] - 
					xyp[i__ + ip1 * 3];
			    }
			    if (ips[lc1 - 1] == ips[0]) {
				goto L20;
			    }
			}
		    }
		}
	    }
L20:
	    ;
	}
    } else if (ifxnode_(&icps, &c__8)) {
	nz = 2;
	i__1 = lf;
	for (n = 1; n <= i__1; ++n) {
	    ipt[0] = ipfs[(n << 2) + 1];
	    ipt[1] = ipfs[(n << 2) + 2];
	    ipt[2] = ipfs[(n << 2) + 3];
	    for (i1 = 1; i1 <= 3; ++i1) {
		if (ip1 == ipt[i1 - 1]) {
		    i2 = iref[i1];
		    i3 = iref[i2];
		    for (i__ = 1; i__ <= 3; ++i__) {
			z__[i__ - 1] = xyp[i__ + ipt[i2 - 1] * 3] - xyp[i__ + 
				ip1 * 3];
			z__[i__ + 2] = xyp[i__ + ipt[i3 - 1] * 3] - xyp[i__ + 
				ip1 * 3];
		    }
		    vecmul_(z__, &z__[3], &z__[6]);
		    vecmul_(z__, &z__[6], &z__[3]);
		}
	    }
	}
    } else if (ifxnode_(&icps, &c__64)) {
	nz = 3;
	for (i__ = 1; i__ <= 3; ++i__) {
	    for (j = 1; j <= 3; ++j) {
		z__[i__ + j * 3 - 4] = 0.;
	    }
	    z__[i__ + i__ * 3 - 4] = 1.;
	}
    }
/* ... calculating the terminal points in the search directions */
    i__1 = nz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (nz == 1) {
	    xstep[i__ - 1] = heit;
	    for (j = 1; j <= 3; ++j) {
		xyps[j + i__ * 3 - 4] = xyp[j + ip1 * 3] + xstep[i__ - 1] * 
			z__[j + i__ * 3 - 4];
	    }
	} else if (nz == 2) {
	    t1 = distsr_(&ip1, &z__[i__ * 3 - 3], &z__[i__ * 3 - 2], &z__[i__ 
		    * 3 - 1], &lf, ifs, &ipfs[5], &xyp[4]);
	    d__1 = -z__[i__ * 3 - 3];
	    d__2 = -z__[i__ * 3 - 2];
	    d__3 = -z__[i__ * 3 - 1];
	    t2 = distsr_(&ip1, &d__1, &d__2, &d__3, &lf, ifs, &ipfs[5], &xyp[
		    4]);
	    if (t2 > t1) {
		t1 = -t2;
	    }
/* Computing MIN */
	    d__1 = dabs(t1);
	    xstep[i__ - 1] = min(d__1,heit) * .02 * d_sign(&c_b11, &t1);
	    for (j = 1; j <= 3; ++j) {
		xyps[j + i__ * 3 - 4] = xyp[j + ip1 * 3] + xstep[i__ - 1] * 
			z__[j + i__ * 3 - 4];
	    }
	} else if (nz == 3) {
	    t1 = distsf_(&ip1, &z__[i__ * 3 - 3], &z__[i__ * 3 - 2], &z__[i__ 
		    * 3 - 1], &le, ies, &ipes[6], &xyp[4]);
	    d__1 = -z__[i__ * 3 - 3];
	    d__2 = -z__[i__ * 3 - 2];
	    d__3 = -z__[i__ * 3 - 1];
	    t2 = distsf_(&ip1, &d__1, &d__2, &d__3, &le, ies, &ipes[6], &xyp[
		    4]);
/*  ...  the case when shut is along a face */
	    if (t1 < 0. && t2 < 0.) {
		goto L1000;
	    }
	    if (t2 > t1) {
		t1 = -t2;
	    }
/* Computing MIN */
	    d__1 = dabs(t1);
	    xstep[i__ - 1] = min(d__1,heit) * .02 * d_sign(&c_b11, &t1);
	    for (j = 1; j <= 3; ++j) {
		xyps[j + i__ * 3 - 4] = xyp[j + ip1 * 3];
	    }
	    xyps[i__ + i__ * 3 - 4] = xyp[i__ + ip1 * 3] + xstep[i__ - 1];
	}
    }
    q[0] = 1. - *rquality;
    i__1 = nz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	q[i__] = nlnfnc_(&c__3, &xyps[i__ * 3 - 3], &xyp[4], &hesp[7], &detg[
		1], hstar, &ip1, &le, ies, xypt, &icps, &ipes[6], hesps, &
		detgs, qes, npw, new__, &xypw[4], &hespw[7], &ipew[5], (I_fp)
		metricfunction, flaganalytic, milintrp, mrlintrp, &ise[1], &
		rse[1], icontrol);
/* group (F) */
/* group (ANI) */
/*  ...   updating the spoiled values of qEs */
	i__2 = le;
	for (n = 1; n <= i__2; ++n) {
	    qes[n - 1] = qeu[n];
	}
    }
    nu = 3;
    for (i__ = 1; i__ <= 3; ++i__) {
	u[i__ - 1] = xyp[i__ + ip1 * 3];
    }
    hmin = 0.;
    if (ifxnode_(&icps, &c__2)) {
	if (q[1] > q[0]) {
	    for (i__ = 1; i__ <= 3; ++i__) {
		z__[i__ - 1] = -z__[i__ - 1];
	    }
	}
	hmax = heit;
    } else if (ifxnode_(&icps, &c__8)) {
	t1 = (q[0] - q[1]) / xstep[0];
	t2 = (q[0] - q[2]) / xstep[1];
	for (i__ = 1; i__ <= 3; ++i__) {
	    z__[i__ - 1] = t1 * z__[i__ - 1] + t2 * z__[i__ + 2];
	}
	hmax = distsr_(&ip1, z__, &z__[1], &z__[2], &lf, ifs, &ipfs[5], &xyp[
		4]);
	if (hmax < 0.) {
	    goto L1000;
	}
	hmax = heit;
    } else if (ifxnode_(&icps, &c__64)) {
	for (i__ = 1; i__ <= 3; ++i__) {
	    z__[i__ - 1] = (q[0] - q[i__]) / xstep[i__ - 1];
	}
	hmax = distsf_(&ip1, z__, &z__[1], &z__[2], &le, ies, &ipes[6], &xyp[
		4]);
	if (hmax < 0.) {
	    goto L1000;
	}
    }
    hmax -= (hmax - hmin) * .2f;
    if (flagtm) {
	hmax *= 2;
    }
    qmin = q[0];
    minim_(&nu, u, z__, &hmin, &hmax, &qmin, u1, &icnt, rmove, &c_false, 
	    flag__, &xyp[4], &ipe[6], &detg[1], &hesp[7], hstar, &ip1, &le, 
	    ies, xyps, &icps, &ipes[6], &detgs, hesps, qes, npw, new__, &xypw[
	    4], &hespw[7], &ipew[5], (I_fp)metricfunction, flaganalytic, 
	    milintrp, mrlintrp, &ise[1], &rse[1], icontrol, &flagtm, &nbad);
/* group(F) */
/* group (ANI) */
    if (*rmove == 0.) {
	goto L1000;
    }
/* ... analysing information from the previous routine */
    m = 1;
    if (! (*flag__) && icps == 64) {
	if (q[1] < q[0]) {
	    m = 1;
	} else if (q[2] < q[0]) {
	    m = 2;
	} else if (q[3] < q[0]) {
	    m = 3;
	} else {
	    goto L1000;
	}
	qmin = nlnfnc_(&c__3, &xyps[m * 3 - 3], &xyp[4], &hesp[7], &detg[1], 
		hstar, &ip1, &le, ies, xypt, &icps, &ipes[6], hesps, &detgs, 
		qes, npw, new__, &xypw[4], &hespw[7], &ipew[5], (I_fp)
		metricfunction, flaganalytic, milintrp, mrlintrp, &ise[1], &
		rse[1], icontrol);
/* group (F) */
/* group (ANI) */
    } else if (! (*flag__)) {
	goto L1000;
    }
/* ... analyzing curvilinear and plane faces */
/*     Call copySQ(0, qEs, XYP(1, iP1), HesPs, detGs, */
/*    &               qEt, XYPt,        HesPt, detGt) */
    pntupd_(&ip1, &icp[1], &xyp[4], &hesp[7], &detg[1], &icps, &xyps[m * 3 - 
	    3], hesps, &detgs);
/* ... checking for inverted elements */
    if (flagtm) {
	i__1 = le;
	for (n = 1; n <= i__1; ++n) {
	    if (ies[n - 1] >= 0) {
		updqb_(&n, &le, ies, &xyp[4], &ipes[6], qes);
	    }
	}
    }
/* ... updating the grid */
    *flag__ = TRUE_;
    i__1 = le;
    for (n = 1; n <= i__1; ++n) {
	if (ies[n - 1] <= 0) {
	    goto L200;
	}
	lstupd_(ne, &l1e[3], nl2, &l2e[1], &nstep[1], &qe[1], &ies[n - 1], &
		qes[n - 1]);
L200:
	;
    }
L1000:
    return 0;
} /* movep_ */

/* ================================================================ */
/* @f2h@ */ doublereal distsr_(integer *ipo, doublereal *nx, doublereal *ny, 
	doublereal *nz, integer *lf, integer *ifs, integer *ipfs, doublereal *
	xyp)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer i__, n, i1, i2, i3;
    static doublereal v1, vv;
    static integer ipa, ipb;
    static doublereal xyc[3], xyd[3], veca[3], vecb[3];
    static integer iref[4];

/* ================================================================ */
/* ================================================================ */
/* group (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    xyp -= 4;
    ipfs -= 5;
    --ifs;

    /* Function Body */
    iref[0] = 1;
    iref[1] = 2;
    iref[2] = 3;
    iref[3] = 1;
    ret_val = 1e12;
    i__1 = *lf;
    for (n = 1; n <= i__1; ++n) {
	if (ifs[n] <= 0) {
	    goto L20;
	}
	for (i1 = 1; i1 <= 3; ++i1) {
	    if (ipfs[i1 + (n << 2)] == *ipo) {
		i2 = iref[i1];
		i3 = iref[i2];
		ipa = ipfs[i2 + (n << 2)];
		ipb = ipfs[i3 + (n << 2)];
		xyc[0] = xyp[*ipo * 3 + 1] + *nx;
		xyc[1] = xyp[*ipo * 3 + 2] + *ny;
		xyc[2] = xyp[*ipo * 3 + 3] + *nz;
		v1 = calvol_(&xyp[*ipo * 3 + 1], &xyp[ipa * 3 + 1], &xyp[ipb *
			 3 + 1], xyc);
		if (dabs(v1) > 1e-10) {
		    goto L20;
		}
		for (i__ = 1; i__ <= 3; ++i__) {
		    veca[i__ - 1] = xyp[i__ + ipa * 3] - xyp[i__ + *ipo * 3];
		    vecb[i__ - 1] = xyp[i__ + ipb * 3] - xyp[i__ + *ipo * 3];
		}
		vecmul_(veca, vecb, xyd);
		for (i__ = 1; i__ <= 3; ++i__) {
		    xyd[i__ - 1] = xyp[i__ + *ipo * 3] + xyd[i__ - 1];
		}
/*              v1 = calVol(XYP(1, iPo), XYD, XYC, XYP(1, iPa)) */
/*              v2 = calVol(XYP(1, iPo), XYD, XYC, XYP(1, iPb)) */
/*              If(v1 * v2.GT.0D0) Goto 20 */
		vv = mutualorientation_(&xyp[*ipo * 3 + 1], xyd, xyc, &xyp[
			ipa * 3 + 1], &xyp[ipb * 3 + 1]);
		if (vv > 0.) {
		    goto L20;
		}
		ret_val = 1.;
		goto L1000;
	    }
	}
L20:
	;
    }
    ret_val = -ret_val;
L1000:
    return ret_val;
} /* distsr_ */

/* ================================================================ */
/* @f2h@ */ doublereal distsf_(integer *ipo, doublereal *nx, doublereal *ny, 
	doublereal *nz, integer *le, integer *ies, integer *ipes, doublereal *
	xyp)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2, d__3;

    /* Local variables */
    static logical flagshut;
    static doublereal d__, h__;
    static integer n;
    static doublereal v;
    static integer i1, i2, i3, i4;
    static doublereal v1, v2, v3;
    static integer ipa, ipb, ipc;
    static doublereal xyd[3];
    static integer iref[5];

/* ================================================================ */
/* ================================================================ */
/* group (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    xyp -= 4;
    ipes -= 6;
    --ies;

    /* Function Body */
    iref[0] = 1;
    iref[1] = 2;
    iref[2] = 3;
    iref[3] = 4;
    iref[4] = 1;
/* Computing 2nd power */
    d__1 = *nx;
/* Computing 2nd power */
    d__2 = *ny;
/* Computing 2nd power */
    d__3 = *nz;
    d__ = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
    ret_val = 1e12;
    i__1 = *le;
    for (n = 1; n <= i__1; ++n) {
	if (ies[n] <= 0) {
	    goto L20;
	}
	for (i1 = 1; i1 <= 4; ++i1) {
	    if (ipes[i1 + n * 5] == *ipo) {
		i2 = iref[i1];
		i3 = iref[i2];
		i4 = iref[i3];
		ipa = ipes[i2 + n * 5];
		ipb = ipes[i3 + n * 5];
		ipc = ipes[i4 + n * 5];
		xyd[0] = xyp[*ipo * 3 + 1] + *nx;
		xyd[1] = xyp[*ipo * 3 + 2] + *ny;
		xyd[2] = xyp[*ipo * 3 + 3] + *nz;
		flagshut = TRUE_;
		shutf_(&xyp[*ipo * 3 + 1], xyd, &xyp[ipa * 3 + 1], &xyp[ipb * 
			3 + 1], &xyp[ipc * 3 + 1], &flagshut);
		if (! flagshut) {
		    goto L20;
		}
		v = calvol_(&xyp[*ipo * 3 + 1], &xyp[ipa * 3 + 1], &xyp[ipb * 
			3 + 1], &xyp[ipc * 3 + 1]);
		v1 = calvol_(&xyp[*ipo * 3 + 1], xyd, &xyp[ipa * 3 + 1], &xyp[
			ipb * 3 + 1]);
		v2 = calvol_(&xyp[*ipo * 3 + 1], xyd, &xyp[ipa * 3 + 1], &xyp[
			ipc * 3 + 1]);
		v3 = calvol_(&xyp[*ipo * 3 + 1], xyd, &xyp[ipb * 3 + 1], &xyp[
			ipc * 3 + 1]);
		h__ = dabs(v) / (dabs(v1) + dabs(v2) + dabs(v3));
/*  ...  checking for the orientation of points O and D */
		xyd[0] = xyp[*ipo * 3 + 1] + h__ * 2 * *nx;
		xyd[1] = xyp[*ipo * 3 + 2] + h__ * 2 * *ny;
		xyd[2] = xyp[*ipo * 3 + 3] + h__ * 2 * *nz;
		v1 = calvol_(xyd, &xyp[ipa * 3 + 1], &xyp[ipb * 3 + 1], &xyp[
			ipc * 3 + 1]);
		if (v * v1 < 0.) {
		    ret_val = h__ * d__;
		    goto L1000;
		}
	    }
	}
L20:
	;
    }
    ret_val = -ret_val;
L1000:
    return ret_val;
} /* distsf_ */

