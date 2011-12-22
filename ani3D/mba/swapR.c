/* swapR.f -- translated by f2c (version 20090411).
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


/* Table of constant values */

static integer c__128 = 128;
static integer c__2 = 2;
static integer c__1 = 1;
static integer c__64 = 64;
static integer c__4 = 4;
static integer c__3 = 3;
static integer c__0 = 0;
static integer c__512 = 512;
static integer c__1007 = 1007;

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int swapr_(integer *iwr, integer *iwe, integer *
	nf, integer *maxf, integer *ne, integer *maxe, doublereal *xyp, 
	integer *ipf, integer *ipe, doublereal *hstar, integer *icp, integer *
	iep, integer *ife, integer *iee, integer *l1e, integer *l2e, integer *
	nl2, integer *nstep, integer *iholf, integer *ihole, integer *status, 
	doublereal *hesp, doublereal *rquality, doublereal *qe, integer *lfu, 
	integer *leu, integer *ifu, integer *ieu, integer *ipfu, integer *
	ipeu, doublereal *qeu, logical *flag__)
{
    /* Initialized data */

    static integer ipr[12]	/* was [2][6] */ = { 1,2,1,3,1,4,2,3,2,4,3,4 }
	    ;

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static logical flagface, flagedge, flagbnds, flagloop, flagshut;
    static integer nitrswap, i__, j, k, m, n;
    static doublereal c1;
    static integer k1;
    static logical l1, l2;
    static integer n1, m1;
    static doublereal oldvolume, newvolume;
    static integer ic, le, lf, ip[5], lp, lr;
    static doublereal vv;
    static integer ne1, ip1, ip2, ip3, ip4, ide, ipa, ipb, ipc, ies[1000], 
	    ifs[1000], ipd;
    static doublereal qes[1000];
    static integer net, ipt[3], irs[3000]	/* was [3][1000] */, ips[1000]
	    ;
    static doublereal ves[1000];
    static integer idx, ipu, ift, iet, nbad, idel[1000], iref[4], icnt, icps[
	    1000], ipes[5000]	/* was [5][1000] */, ipfs[4000]	/* was [4][
	    1000] */, idps[2000]	/* was [2][1000] */, ilev[1000], nlev,
	     inps[2000]	/* was [2][1000] */;
    static doublereal xypd[3], xyps[3];
    static integer ipbad, icrab, iclrf, lpend, leold, lfold;
    static logical flagtm;

/* ================================================================ */
/* group (M) */
/* group (Q) */
/* group (S) */
/* group (W) */
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
/* group (Flag) */
/* ================================================================ */
/* group (Functions) */
/* group (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --qeu;
    ipeu -= 6;
    ipfu -= 5;
    --ieu;
    --ifu;
    --qe;
    hesp -= 7;
    --ihole;
    --iholf;
    --nstep;
    --l2e;
    l1e -= 3;
    iee -= 5;
    ife -= 5;
    --iep;
    --icp;
    ipe -= 6;
    ipf -= 5;
    xyp -= 4;

    /* Function Body */
/* ================================================================ */
    *flag__ = FALSE_;
    copyse_(lfu, leu, &ifu[1], &ieu[1], &ipfu[5], &ipeu[6], &qeu[1], &lf, &le,
	     ifs, ies, ipfs, ipes, qes);
    ip[0] = 1;
    ip[1] = 2;
    ip[2] = 3;
    ip[3] = 4;
    ip[4] = 1;
    iref[0] = 1;
    iref[1] = 2;
    iref[2] = 3;
    iref[3] = 1;
    ipa = ipe[ipr[(*iwr << 1) - 2] + *iwe * 5];
    ipb = ipe[ipr[(*iwr << 1) - 1] + *iwe * 5];
/* ... checking the case when swapping is impossible */
    if (ifxnode_(&icp[ipa], &c__128) && ifxnode_(&icp[ipb], &c__128)) {
	goto L1000;
    }
    clrsr_(&ipa, &ipb, &icp[1], &ipf[5], &ife[5], &lf, ifs, &le, ies, &icrab);
    if (ifxnode_(&icrab, &c__2)) {
	goto L1000;
    }
/* ... checking for surrounding points */
    maksr_(&ipa, &ipb, &le, ies, ipes, &lr, irs, &flagedge);
    if (! flagedge) {
	goto L1000;
    }
    if (ifxnode_(status, &c__2) || ifxnode_(status, &c__1)) {
	flagbnds = TRUE_;
	i__1 = lr;
	for (n = 1; n <= i__1; ++n) {
	    if (ifxnode_(&icp[irs[n * 3 - 2]], &c__64)) {
		flagbnds = FALSE_;
	    }
	}
    } else {
	flagbnds = FALSE_;
    }
    if (flagbnds) {
	idps[0] = ipa;
	idps[1] = ipb;
	if (ifxnode_(&icp[ipa], &c__4) && ifxnode_(&icp[ipb], &c__64)) {
	    chkspf_(&c__1, &ipa, &c__3, &icp[1], &iep[1], &ipe[6], &iee[5], &
		    lp, ips);
	    chkspb_(&c__2, &c__1, idps, &c__0, inps, &c__3, &icp[1], &iep[1], 
		    &ipe[6], &iee[5], &lp, ips, &ipbad, &flagbnds);
	    if (flagbnds) {
		goto L1000;
	    }
	} else if (ifxnode_(&icp[ipb], &c__4) && ifxnode_(&icp[ipa], &c__64)) 
		{
	    chkspf_(&c__1, &ipb, &c__3, &icp[1], &iep[1], &ipe[6], &iee[5], &
		    lp, ips);
	    chkspb_(&c__2, &c__1, idps, &c__0, inps, &c__3, &icp[1], &iep[1], 
		    &ipe[6], &iee[5], &lp, ips, &ipbad, &flagbnds);
	    if (flagbnds) {
		goto L1000;
	    }
	} else if (ifxnode_(&icp[ipb], &c__4) && ifxnode_(&icp[ipa], &c__4)) {
	    lp = 2;
	    ips[0] = ipa;
	    ips[1] = ipb;
	    chkspb_(&c__2, &c__1, idps, &c__0, inps, &c__3, &icp[1], &iep[1], 
		    &ipe[6], &iee[5], &lp, ips, &ipbad, &flagbnds);
	    if (flagbnds) {
		goto L1000;
	    }
	}
    }
    flagloop = TRUE_;
    if (irs[lr * 3 - 1] != irs[1]) {
	flagloop = FALSE_;
    }
    lpend = lr;
    i__1 = lr;
    for (n = 1; n <= i__1; ++n) {
	ips[n - 1] = irs[n * 3 - 2];
    }
    if (! flagloop) {
	lpend = lr + 1;
	ips[lpend - 1] = irs[lr * 3 - 1];
    }
/* ... search for boundary points adjacent to points A and B */
    flagface = FALSE_;
    icnt = 0;
    i__1 = lf;
    for (n = 1; n <= i__1; ++n) {
	ipt[0] = ipfs[(n << 2) - 4];
	ipt[1] = ipfs[(n << 2) - 3];
	ipt[2] = ipfs[(n << 2) - 2];
	l1 = check13_(&ipa, ipt, &ipt[1], &ipt[2]);
	l2 = check13_(&ipb, ipt, &ipt[1], &ipt[2]);
	if (l1 && l2) {
	    for (i__ = 1; i__ <= 3; ++i__) {
		if (ipt[i__ - 1] != ipa && ipt[i__ - 1] != ipb) {
		    flagface = TRUE_;
		    ++icnt;
		    if (icnt == 1) {
			ipc = ipt[i__ - 1];
		    }
		    if (icnt == 2) {
			ipd = ipt[i__ - 1];
		    }
		    iclrf = ipfs[(n << 2) - 1];
		    ifs[n - 1] = -ifs[n - 1];
		}
	    }
	}
    }
/* ... only two faces are allowed to have the same edge */
    if (flagface && icnt != 2) {
	goto L1000;
    }
/* ... checking piece-wise plane material boundaries */
    if (flagface) {
	c1 = angle2faces_(&xyp[ipa * 3 + 1], &xyp[ipb * 3 + 1], &xyp[ipc * 3 
		+ 1], &xyp[ipd * 3 + 1]);
	if (c1 > -.99999999999999989) {
	    goto L1000;
	}
    }
/* ... checking the number of inverted elements */
    flagtm = ifxnode_(status, &c__512);
    if (flagtm) {
	nbad = 0;
	i__1 = lr;
	for (n = 1; n <= i__1; ++n) {
	    if (qes[irs[n * 3 - 3] - 1] <= 0.) {
		++nbad;
	    }
	}
	if (nbad > 0) {
	    goto L1000;
	}
    }
/* ... coloring the superelement by element colors */
/* ... array iLev is overloaded */
    i__1 = lr + 1;
    for (n = 1; n <= i__1; ++n) {
	icps[n - 1] = 0;
    }
    i__1 = lr;
    for (n = 1; n <= i__1; ++n) {
	net = irs[n * 3 - 3];
	ide = ipes[net * 5 - 1];
	if (icps[n - 1] != ide) {
	    icps[n - 1] += ide;
	}
	if (icps[n] != ide) {
	    icps[n] += ide;
	}
    }
    if (flagloop && icps[0] != icps[lr]) {
	icps[0] += icps[lr];
    }
/* ... checking that we have at most 2 different colors */
    i__1 = lpend;
    for (n = 1; n <= i__1; ++n) {
	ilev[n - 1] = icps[n - 1];
    }
    ic = countcolors_(&lpend, ilev);
    if (ic > 3) {
	goto L1000;
    }
/* ... making a virtual evaluation of the quality */
    leold = le;
    oldvolume = 0.;
    i__1 = lr;
    for (n = 1; n <= i__1; ++n) {
	net = irs[n * 3 - 3];
	ip1 = ipes[net * 5 - 5];
	ip2 = ipes[net * 5 - 4];
	ip3 = ipes[net * 5 - 3];
	ip4 = ipes[net * 5 - 2];
	oldvolume += (d__1 = calvol_(&xyp[ip1 * 3 + 1], &xyp[ip2 * 3 + 1], &
		xyp[ip3 * 3 + 1], &xyp[ip4 * 3 + 1]), dabs(d__1));
    }
/* ... starting the main loop */
    idx = 1;
    nlev = lpend - 2;
    ilev[0] = 0;
    nitrswap = 0;
L100:
    ++ilev[idx - 1];
    i__1 = lpend;
    for (n1 = ilev[idx - 1]; n1 <= i__1; ++n1) {
	++nitrswap;
	if (nitrswap > 300) {
	    goto L1000;
	}
	if (ips[n1 - 1] <= 0) {
	    goto L400;
	}
	ip1 = ips[n1 - 1];
	i__2 = lpend;
	for (m = n1 + 1; m <= i__2; ++m) {
	    if (ips[m - 1] > 0) {
		m1 = m;
		goto L200;
	    }
	}
	goto L400;
L200:
	ip2 = ips[m1 - 1];
	i__2 = lpend;
	for (k = m1 + 1; k <= i__2; ++k) {
	    if (ips[k - 1] > 0) {
		k1 = k;
		goto L300;
	    }
	}
	goto L400;
L300:
	ip3 = ips[k1 - 1];
	if (idx < nlev) {
	    flagshut = FALSE_;
	    shutf_(&xyp[ip1 * 3 + 1], &xyp[ip3 * 3 + 1], &xyp[ipa * 3 + 1], &
		    xyp[ipb * 3 + 1], &xyp[ip2 * 3 + 1], &flagshut);
	    if (! flagshut) {
		goto L400;
	    }
	} else {
	    if (m1 != n1 + 1) {
		i__ = ip2;
		ip2 = ip3;
		ip3 = i__;
	    }
	    for (i__ = 1; i__ <= 3; ++i__) {
		xyps[i__ - 1] = (xyp[i__ + ip3 * 3] + xyp[i__ + ip1 * 3]) / 2;
		xypd[i__ - 1] = (xyp[i__ + ipa * 3] + xyp[i__ + ipb * 3]) / 2;
	    }
/*           v1 = calVol(XYP(1, iPa), XYP(1, iP1), XYP(1, iP2), XYPs) */
/*           v2 = calVol(XYP(1, iPa), XYP(1, iP1), XYP(1, iP2), XYPd) */

/*           If(v1 * v2.LE.0D0) Goto 400 */

/*           v1 = calVol(XYP(1, iPb), XYP(1, iP1), XYP(1, iP2), XYPs) */
/*           v2 = calVol(XYP(1, iPb), XYP(1, iP1), XYP(1, iP2), XYPd) */

/*           If(v1 * v2.LE.0D0) Goto 400 */
	    vv = mutualorientation_(&xyp[ipa * 3 + 1], &xyp[ip1 * 3 + 1], &
		    xyp[ip2 * 3 + 1], xyps, xypd);
	    if (vv <= 0.) {
		goto L400;
	    }
	    vv = mutualorientation_(&xyp[ipb * 3 + 1], &xyp[ip1 * 3 + 1], &
		    xyp[ip2 * 3 + 1], xyps, xypd);
	    if (vv <= 0.) {
		goto L400;
	    }
	}
	ne1 = irs[m1 * 3 - 3];
	calqe_(&hesp[ipa * 6 + 1], &xyp[ipa * 3 + 1], &hesp[ip1 * 6 + 1], &
		xyp[ip1 * 3 + 1], &hesp[ip2 * 6 + 1], &xyp[ip2 * 3 + 1], &
		hesp[ip3 * 6 + 1], &xyp[ip3 * 3 + 1], hstar, &qes[ne1 - 1], &
		ves[ne1 - 1]);
	if (qes[ne1 - 1] <= *rquality) {
	    goto L400;
	}
	++le;
	if (le > 1000) {
	    goto L9000;
	}
	ies[le - 1] = 0;
	calqe_(&hesp[ipb * 6 + 1], &xyp[ipb * 3 + 1], &hesp[ip1 * 6 + 1], &
		xyp[ip1 * 3 + 1], &hesp[ip2 * 6 + 1], &xyp[ip2 * 3 + 1], &
		hesp[ip3 * 6 + 1], &xyp[ip3 * 3 + 1], hstar, &qes[le - 1], &
		ves[le - 1]);
	if (qes[le - 1] <= *rquality) {
	    --le;
	    goto L400;
	}
/*  ...  go to the next level */
/* Computing MIN */
	i__2 = icps[n1 - 1], i__3 = icps[m1 - 1];
	ide = min(i__2,i__3);
/* Computing MIN */
	i__2 = ide, i__3 = icps[k1 - 1];
	ide = min(i__2,i__3);
/*        ide = IPEs(5, nE1)  keep the mistaken choice */
	ipes[ne1 * 5 - 5] = ipa;
	ipes[ne1 * 5 - 4] = ip1;
	ipes[ne1 * 5 - 3] = ip2;
	ipes[ne1 * 5 - 2] = ip3;
	ipes[ne1 * 5 - 1] = ide;
	ipes[le * 5 - 5] = ipb;
	ipes[le * 5 - 4] = ip1;
	ipes[le * 5 - 3] = ip2;
	ipes[le * 5 - 2] = ip3;
	ipes[le * 5 - 1] = ide;
	inps[(idx << 1) - 2] = ips[n1 - 1];
	inps[(idx << 1) - 1] = ips[k1 - 1];
	idel[idx - 1] = m1;
	ips[m1 - 1] = -ips[m1 - 1];
	if (idx == nlev) {
	    newvolume = 0.;
	    i__2 = le;
	    for (n = leold + 1; n <= i__2; ++n) {
		newvolume += (d__1 = ves[n - 1], dabs(d__1));
	    }
	    i__2 = nlev;
	    for (n = 1; n <= i__2; ++n) {
		net = irs[idel[n - 1] * 3 - 3];
		newvolume += (d__1 = ves[net - 1], dabs(d__1));
	    }
/* ... checking only in the case of plane surfaces */
	    if ((d__1 = oldvolume - newvolume, dabs(d__1)) > oldvolume * 1e-10)
		     {
/*              Write(*, 5000) oldVolume, newVolume */
		goto L400;
	    }
	    goto L500;
	} else if (idx == nlev - 1 && flagface) {
	    i__2 = nlev - 1;
	    for (n = 1; n <= i__2; ++n) {
		for (j = 1; j <= 2; ++j) {
		    if (inps[j + (n << 1) - 3] == ipc && inps[3 - j + (n << 1)
			     - 3] == ipd) {
			goto L350;
		    }
		}
	    }
	    goto L400;
	}
L350:
	++idx;
	ilev[idx - 1] = 0;
	goto L100;
L400:
	;
    }
    --idx;
    if (idx == 0) {
	goto L1000;
    }
    m1 = idel[idx - 1];
    ips[m1 - 1] = -ips[m1 - 1];
    --le;
    goto L100;
/* ... analyzing curvilinear and plane faces */
L500:
    lfold = lf;
    if (flagface) {
	lf += 2;
	if (lf > 1000) {
	    goto L9000;
	}
	i__1 = lfold;
	for (m = 1; m <= i__1; ++m) {
	    if (check33_(&ipa, &ipc, &ipd, &ipfs[(m << 2) - 4], &ipfs[(m << 2)
		     - 3], &ipfs[(m << 2) - 2])) {
		goto L1000;
	    }
	}
	ipfs[(lf - 1 << 2) - 4] = ipa;
	ipfs[(lf - 1 << 2) - 3] = ipc;
	ipfs[(lf - 1 << 2) - 2] = ipd;
	ipfs[(lf - 1 << 2) - 1] = iclrf;
	i__1 = lfold;
	for (m = 1; m <= i__1; ++m) {
	    if (check33_(&ipb, &ipc, &ipd, &ipfs[(m << 2) - 4], &ipfs[(m << 2)
		     - 3], &ipfs[(m << 2) - 2])) {
		goto L1000;
	    }
	}
	ipfs[(lf << 2) - 4] = ipb;
	ipfs[(lf << 2) - 3] = ipc;
	ipfs[(lf << 2) - 2] = ipd;
	ipfs[(lf << 2) - 1] = iclrf;
	facadd_(&ifs[lf - 2], nf, maxf, &iholf[1]);
	facadd_(&ifs[lf - 1], nf, maxf, &iholf[1]);
    }
/* ... checking for boundary elements */
    if (ifxnode_(status, &c__1)) {
	i__1 = lpend - 1;
	for (n = 2; n <= i__1; ++n) {
	    net = irs[n * 3 - 3];
	    for (i__ = 1; i__ <= 4; ++i__) {
		ipu = ipes[i__ + net * 5 - 6];
		if (ifxnode_(&icp[ipu], &c__64)) {
		    goto L700;
		}
	    }
	    goto L2000;
L700:
	    ;
	}
	i__1 = le;
	for (n = leold + 1; n <= i__1; ++n) {
	    for (i__ = 1; i__ <= 4; ++i__) {
		ipu = ipes[i__ + n * 5 - 6];
		if (ifxnode_(&icp[ipu], &c__64)) {
		    goto L800;
		}
	    }
	    goto L2000;
L800:
	    ;
	}
    }
/* ... updating the grid */
    *flag__ = TRUE_;
    if (flagface) {
	i__1 = lfold;
	for (n = 1; n <= i__1; ++n) {
	    ift = ifs[n - 1];
	    if (ift <= 0) {
		ift = -ift;
		facdel_(&ift, nf, &ipf[5], &iholf[1]);
	    }
	}
	i__1 = lf - 1;
	facupd_(&i__1, &ipf[5], ifs, ipfs);
	facupd_(&lf, &ipf[5], ifs, ipfs);
    }
    net = irs[0];
    iet = ies[net - 1];
    ies[net - 1] = -iet;
    lstdel_(ne, &l1e[3], nl2, &l2e[1], &nstep[1], &ihole[1], &qe[1], &iet);
    eledel_(&iet, &ipe[6], &iee[5]);
    if (flagloop) {
	net = irs[lr * 3 - 3];
	iet = ies[net - 1];
	ies[net - 1] = -iet;
	lstdel_(ne, &l1e[3], nl2, &l2e[1], &nstep[1], &ihole[1], &qe[1], &iet)
		;
	eledel_(&iet, &ipe[6], &iee[5]);
    }
    i__1 = lpend - 1;
    for (n = 2; n <= i__1; ++n) {
	net = irs[n * 3 - 3];
	iet = ies[net - 1];
	lstupd_(ne, &l1e[3], nl2, &l2e[1], &nstep[1], &qe[1], &iet, &qes[net 
		- 1]);
	eledel_(&iet, &ipe[6], &iee[5]);
    }
    i__1 = le;
    for (n = leold + 1; n <= i__1; ++n) {
	eleadd_(ne, maxe, &ihole[1]);
	lstadd_(ne, &l1e[3], nl2, &l2e[1], &nstep[1], &ihole[1], &qe[1], &qes[
		n - 1], &ies[n - 1]);
	eledel_(&ies[n - 1], &ipe[6], &iee[5]);
    }
    i__1 = lpend - 1;
    for (n = 2; n <= i__1; ++n) {
	net = irs[n * 3 - 3];
	eleupd_(&net, &iep[1], &ipe[6], &ife[5], &iee[5], &lf, &le, ifs, ies, 
		ipfs, ipes);
    }
    i__1 = le;
    for (n = leold + 1; n <= i__1; ++n) {
	eleupd_(&n, &iep[1], &ipe[6], &ife[5], &iee[5], &lf, &le, ifs, ies, 
		ipfs, ipes);
    }
/* L5000: */
L1000:
    return 0;
/* ... delete wrong faces */
L2000:
    i__1 = lf;
    for (n = lfold + 1; n <= i__1; ++n) {
	facdel_(&ifs[n - 1], nf, &ipf[5], &iholf[1]);
    }
    return 0;
L9000:
    errmes_(&c__1007, "swapR", "local parameter MaxS is small", (ftnlen)5, (
	    ftnlen)29);
    return 0;
} /* swapr_ */

