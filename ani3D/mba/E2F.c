/* E2F.f -- translated by f2c (version 20090411).
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

/* Table of constant values */

static integer c__2 = 2;
static integer c__512 = 512;
static integer c__4 = 4;
static integer c__64 = 64;
static integer c__1 = 1;
static integer c__0 = 0;
static integer c__1007 = 1007;

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int e2f_(integer *iwr, integer *iwe, integer *nf,
	 integer *maxf, integer *ne, doublereal *xyp, integer *ipf, integer *
	ipe, doublereal *hstar, integer *icp, integer *iep, integer *ife, 
	integer *iee, integer *l1e, integer *l2e, integer *nl2, integer *
	nstep, integer *iholf, integer *ihole, integer *status, doublereal *
	hesp, doublereal *rquality, doublereal *detg, doublereal *qe, integer 
	*lfu, integer *leu, integer *ifu, integer *ieu, integer *ipfu, 
	integer *ipeu, doublereal *qeu, logical *flag__)
{
    /* Initialized data */

    static integer ipr[12]	/* was [2][6] */ = { 1,2,1,3,1,4,2,3,2,4,3,4 }
	    ;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static logical flagface, flagedge, flagbnds;
    static integer i__, m, n;
    static doublereal c1;
    static logical l1, l2;
    static doublereal v1, v2;
    static integer le, lf, lp, lr;
    static doublereal vv;
    static integer ie3, ne1, ne2, ne3, ip1, ip2, ip3, ip4, ide, ipa, ipb, ipc,
	     ies[1000], ifs[1000], ids[2]	/* was [2][1] */;
    static doublereal qes[1000];
    static integer ins[2]	/* was [2][1] */, ips[1000], irs[3000]	/* 
	    was [3][1000] */, ipt[3], ipd, ift, nbad, icnt, ipes[5000]	/* 
	    was [5][1000] */, ipfs[4000]	/* was [4][1000] */, ipbad, 
	    icrab, iclrf, lfold;
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
    --detg;
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
/* ... checking operations which do not modify the initial superelement */
    ipa = ipe[ipr[(*iwr << 1) - 2] + *iwe * 5];
    ipb = ipe[ipr[(*iwr << 1) - 1] + *iwe * 5];
    maksr_(&ipa, &ipb, leu, &ieu[1], &ipeu[6], &lr, irs, &flagedge);
    if (! flagedge) {
	goto L1000;
    }
/* ... checking the case when edge -> face operation is impossible */
    if (lr != 3) {
	goto L1000;
    }
    ip1 = irs[1];
    ip2 = irs[2];
    ip3 = irs[5];
    ip4 = irs[8];
    if (ip1 != ip4) {
	goto L1000;
    }
    clrsr_(&ipa, &ipb, &icp[1], &ipf[5], &ife[5], lfu, &ifu[1], leu, &ieu[1], 
	    &icrab);
    if (ifxnode_(&icrab, &c__2)) {
	goto L1000;
    }
/* ... checking the orientation */
/*     v1 = calVol(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3), XYP(1, iPa)) */
/*     v2 = calVol(XYP(1, iP1), XYP(1, iP2),  XYP(1, iP3), XYP(1, iPb)) */

/*     If(v1 * v2.GE.0D0) Goto 1000 */
    vv = mutualorientation_(&xyp[ip1 * 3 + 1], &xyp[ip2 * 3 + 1], &xyp[ip3 * 
	    3 + 1], &xyp[ipa * 3 + 1], &xyp[ipb * 3 + 1]);
    if (vv >= 0.) {
	goto L1000;
    }
/* ... saving the initial superelement structure (it is changed in the sequel) */
    copyse_(lfu, leu, &ifu[1], &ieu[1], &ipfu[5], &ipeu[6], &qeu[1], &lf, &le,
	     ifs, ies, ipfs, ipes, qes);
/* ... number of inverted elements */
    flagtm = ifxnode_(status, &c__512);
    if (flagtm) {
	nbad = 0;
	if (qes[irs[0] - 1] <= 0.) {
	    ++nbad;
	}
	if (qes[irs[3] - 1] <= 0.) {
	    ++nbad;
	}
	if (qes[irs[6] - 1] <= 0.) {
	    ++nbad;
	}
	if (nbad > 0) {
	    goto L1000;
	}
    }
/* ... making a virtual evaluation of the quality */
    ne1 = irs[0];
    ne2 = irs[3];
    ne3 = irs[6];
    ide = ipes[ne1 * 5 - 1];
    calqe_(&hesp[ipa * 6 + 1], &xyp[ipa * 3 + 1], &hesp[ip1 * 6 + 1], &xyp[
	    ip1 * 3 + 1], &hesp[ip2 * 6 + 1], &xyp[ip2 * 3 + 1], &hesp[ip3 * 
	    6 + 1], &xyp[ip3 * 3 + 1], hstar, &qes[ne1 - 1], &v1);
    if (qes[ne1 - 1] <= *rquality) {
	goto L1000;
    }
    calqe_(&hesp[ipb * 6 + 1], &xyp[ipb * 3 + 1], &hesp[ip1 * 6 + 1], &xyp[
	    ip1 * 3 + 1], &hesp[ip2 * 6 + 1], &xyp[ip2 * 3 + 1], &hesp[ip3 * 
	    6 + 1], &xyp[ip3 * 3 + 1], hstar, &qes[ne2 - 1], &v2);
    if (qes[ne2 - 1] <= *rquality) {
	goto L1000;
    }
/* ... checking for surrounding points */
    if (ifxnode_(&icp[ip1], &c__4) && ifxnode_(&icp[ip2], &c__4) && ifxnode_(&
	    icp[ip3], &c__4)) {
	if (ifxnode_(status, &c__2)) {
	    ids[0] = ipa;
	    ids[1] = ipb;
	    if (ifxnode_(&icp[ipa], &c__4) && ifxnode_(&icp[ipb], &c__64)) {
		chkspf_(&c__1, &ipa, &c__2, &icp[1], &iep[1], &ipe[6], &iee[5]
			, &lp, ips);
		chkspb_(&c__2, &c__1, ids, &c__0, ins, &c__2, &icp[1], &iep[1]
			, &ipe[6], &iee[5], &lp, ips, &ipbad, &flagbnds);
		if (flagbnds) {
		    goto L1000;
		}
	    } else if (ifxnode_(&icp[ipb], &c__4) && ifxnode_(&icp[ipa], &
		    c__64)) {
		chkspf_(&c__1, &ipb, &c__2, &icp[1], &iep[1], &ipe[6], &iee[5]
			, &lp, ips);
		chkspb_(&c__2, &c__1, ids, &c__0, ins, &c__2, &icp[1], &iep[1]
			, &ipe[6], &iee[5], &lp, ips, &ipbad, &flagbnds);
		if (flagbnds) {
		    goto L1000;
		}
	    } else if (ifxnode_(&icp[ipb], &c__4) && ifxnode_(&icp[ipa], &
		    c__4)) {
		lp = 2;
		ips[0] = ipa;
		ips[1] = ipb;
		chkspb_(&c__2, &c__1, ids, &c__0, ins, &c__2, &icp[1], &iep[1]
			, &ipe[6], &iee[5], &lp, ips, &ipbad, &flagbnds);
		if (flagbnds) {
		    goto L1000;
		}
	    }
	}
    }
/* ... search for boundary points adjacent to point A and B */
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
/* ... analyzing curvilinear and plane faces */
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
	if (ifxnode_(&icp[ip1], &c__64)) {
	    goto L100;
	}
	if (ifxnode_(&icp[ip2], &c__64)) {
	    goto L100;
	}
	if (ifxnode_(&icp[ip3], &c__64)) {
	    goto L100;
	}
	if (ifxnode_(&icp[ipa], &c__4)) {
	    goto L1000;
	}
	if (ifxnode_(&icp[ipb], &c__4)) {
	    goto L1000;
	}
    }
/* ... updating the grid */
L100:
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
/* ... creating new elements */
    ipes[ne1 * 5 - 5] = ipa;
    ipes[ne1 * 5 - 4] = ip1;
    ipes[ne1 * 5 - 3] = ip2;
    ipes[ne1 * 5 - 2] = ip3;
    ipes[ne1 * 5 - 1] = ide;
    ipes[ne2 * 5 - 5] = ipb;
    ipes[ne2 * 5 - 4] = ip1;
    ipes[ne2 * 5 - 3] = ip2;
    ipes[ne2 * 5 - 2] = ip3;
    ipes[ne2 * 5 - 1] = ide;
    ie3 = ies[ne3 - 1];
    ies[ne3 - 1] = -ie3;
    lstdel_(ne, &l1e[3], nl2, &l2e[1], &nstep[1], &ihole[1], &qe[1], &ie3);
    eledel_(&ie3, &ipe[6], &iee[5]);
    lstupd_(ne, &l1e[3], nl2, &l2e[1], &nstep[1], &qe[1], &ies[ne1 - 1], &qes[
	    ne1 - 1]);
    lstupd_(ne, &l1e[3], nl2, &l2e[1], &nstep[1], &qe[1], &ies[ne2 - 1], &qes[
	    ne2 - 1]);
    eledel_(&ies[ne1 - 1], &ipe[6], &iee[5]);
    eledel_(&ies[ne2 - 1], &ipe[6], &iee[5]);
    eleupd_(&ne1, &iep[1], &ipe[6], &ife[5], &iee[5], &lf, &le, ifs, ies, 
	    ipfs, ipes);
    eleupd_(&ne2, &iep[1], &ipe[6], &ife[5], &iee[5], &lf, &le, ifs, ies, 
	    ipfs, ipes);
L1000:
    return 0;
L9000:
    errmes_(&c__1007, "E2F", "local parameter MaxS is small", (ftnlen)3, (
	    ftnlen)29);
    return 0;
} /* e2f_ */

