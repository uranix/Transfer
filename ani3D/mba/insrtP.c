/* insrtP.f -- translated by f2c (version 20090411).
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

static integer c__512 = 512;
static integer c__128 = 128;
static integer c__64 = 64;
static integer c__1 = 1;
static integer c__4 = 4;
static integer c__2 = 2;
static integer c__5 = 5;
static integer c__0 = 0;
static integer c__1007 = 1007;

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int insrtp_(integer *iwr, integer *iwe, integer *
	np, integer *maxp, integer *nf, integer *maxf, integer *ne, integer *
	maxe, doublereal *xyp, integer *ipf, integer *ipe, doublereal *hstar, 
	integer *icp, integer *iep, integer *ife, integer *iee, integer *l1e, 
	integer *l2e, integer *nl2, integer *nstep, integer *iholp, integer *
	iholf, integer *ihole, integer *status, doublereal *hesp, doublereal *
	rquality, doublereal *detg, doublereal *qe, I_fp metricfunction, 
	logical *flaganalytic, integer *lfu, integer *leu, integer *ifu, 
	integer *ieu, integer *ipfu, integer *ipeu, doublereal *qeu, integer *
	npw, integer *new__, doublereal *xypw, doublereal *hespw, integer *
	ipew, integer *milintrp, integer *mrlintrp, integer *ise, doublereal *
	rse, integer *icontrol, logical *flag__)
{
    /* Initialized data */

    static integer ipr[12]	/* was [2][6] */ = { 1,2,1,3,1,4,2,3,2,4,3,4 }
	    ;

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static logical flagedge, flagbnds, flagloop;
    static integer i__, m, n;
    static doublereal v1, v2;
    static integer le, lf, ip[5], lp, lr, ne1, nf1, ne2, ip1, ip2, ide, ldh, 
	    ipa, ipb, ipc, ifs[1000], ies[1000];
    static doublereal qes[1000];
    static integer irs[3000]	/* was [3][1000] */, ips[1000], ipt, nxy, 
	    nbad, icfs[1000], icps, ipes[5000]	/* was [5][1000] */, ipfs[
	    4000]	/* was [4][1000] */, idps[2]	/* was [2][1] */, 
	    inps[2]	/* was [2][1] */;
    static doublereal xypo[3], xyps[3];
    static integer ipbad, icrab, leold, lpend, lfold, iclrf;
    static doublereal detgs, hesps[6];
    static logical flagtm, flagfbe;

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
/* ================================================================ */
/*   Inserting a point in the middle of mesh edge. The new is */
/*   projected to a curvilinear surface or edge if necessary. */
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
/* group (W) */
/* group (Flag) */
/* ================================================================ */
/* group (Functions) */
/* group (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --rse;
    --ise;
    ipew -= 5;
    hespw -= 7;
    xypw -= 4;
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
    --iholp;
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
    ipa = ipe[ipr[(*iwr << 1) - 2] + *iwe * 5];
    ipb = ipe[ipr[(*iwr << 1) - 1] + *iwe * 5];
/* ... checking the case when insertion is impossible */
    maksr_(&ipa, &ipb, &le, ies, ipes, &lr, irs, &flagedge);
    if (! flagedge) {
	goto L1000;
    }
    flagloop = TRUE_;
    if (irs[lr * 3 - 1] != irs[1]) {
	flagloop = FALSE_;
    }
    lpend = lr;
    if (! flagloop) {
	lpend = lr + 1;
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
/* ... making a virtual evaluation of the quality */
    clrsr_(&ipa, &ipb, &icp[1], &ipf[5], &ife[5], &lf, ifs, &le, ies, &icrab);
    icps = icrab;
    if (ifxnode_(&icps, &c__128)) {
	goto L1000;
    }
    if (ifxnode_(&icp[ipa], &c__128) && ifxnode_(&icp[ipb], &c__128)) {
	goto L1000;
    }
    for (i__ = 1; i__ <= 3; ++i__) {
	xyps[i__ - 1] = (xyp[i__ + ipa * 3] + xyp[i__ + ipb * 3]) / 2;
    }
    if (ifxnode_(&icrab, &c__64)) {
	for (i__ = 1; i__ <= 6; ++i__) {
	    hesps[i__ - 1] = (hesp[i__ + ipa * 6] + hesp[i__ + ipb * 6]) / 2;
	}
	ldh = 6;
	nxy = 1;
	if (! (*flaganalytic)) {
	    lintrp3d_(new__, &ipew[5], npw, &xypw[4], &ldh, &hespw[7], &nxy, 
		    xyps, hesps, &ise[1], milintrp, &rse[1], mrlintrp, 
		    icontrol);
	} else {
	    scaleback_(xyps, xypo);
	    iniq_analytic__(&c__1, xypo, (I_fp)metricfunction, hesps);
	}
    } else {
	lp = lpend + 2;
	i__1 = lpend;
	for (n = 1; n <= i__1; ++n) {
	    if (n < lpend) {
		ips[n - 1] = irs[n * 3 - 2];
	    } else {
		ips[n - 1] = irs[(n - 1) * 3 - 1];
	    }
	}
	ips[lp - 2] = ipa;
	ips[lp - 1] = ipb;
	hesbnd_(&lp, ips, &icp[1], &hesp[7], hesps);
    }
    caldet_(hesps, &detgs);
    flagfbe = ifxnode_(status, &c__1);
    flagfbe = flagfbe && ifxnode_(&icps, &c__4);
    flagfbe = flagfbe && (ifxnode_(&icp[ipa], &c__4) || ifxnode_(&icp[ipb], &
	    c__4));
    leold = le;
    i__1 = lr;
    for (n = 1; n <= i__1; ++n) {
	ip1 = irs[n * 3 - 2];
	ip2 = irs[n * 3 - 1];
/*  ...  checking for boundary elements */
	if (flagfbe) {
	    if (ifxnode_(&icp[ip1], &c__4) && ifxnode_(&icp[ip2], &c__4)) {
		goto L1000;
	    }
	}
	ne1 = irs[n * 3 - 3];
	calqe_(&hesp[ipa * 6 + 1], &xyp[ipa * 3 + 1], &hesp[ip1 * 6 + 1], &
		xyp[ip1 * 3 + 1], &hesp[ip2 * 6 + 1], &xyp[ip2 * 3 + 1], 
		hesps, xyps, hstar, &qes[ne1 - 1], &v1);
	if (qes[ne1 - 1] <= *rquality) {
	    goto L1000;
	}
	++le;
	if (le > 1000) {
	    goto L9000;
	}
	ies[le - 1] = 0;
	calqe_(&hesp[ipb * 6 + 1], &xyp[ipb * 3 + 1], &hesp[ip1 * 6 + 1], &
		xyp[ip1 * 3 + 1], &hesp[ip2 * 6 + 1], &xyp[ip2 * 3 + 1], 
		hesps, xyps, hstar, &qes[le - 1], &v2);
	if (qes[le - 1] <= *rquality) {
	    goto L1000;
	}
    }
/* ... checking for surrounding points */
    flagbnds = FALSE_;
    if (ifxnode_(status, &c__2)) {
	i__1 = lr;
	for (n = 1; n <= i__1; ++n) {
	    ipt = irs[n * 3 - 2];
	    if (ifxnode_(&icp[ipt], &c__4)) {
		flagbnds = TRUE_;
	    }
	}
    }
    if (flagbnds) {
	idps[0] = ipa;
	idps[1] = ipb;
	if (ifxnode_(&icp[ipa], &c__4) && ifxnode_(&icp[ipb], &c__4)) {
	    chkspf_(&c__1, &ipa, &c__5, &icp[1], &iep[1], &ipe[6], &iee[5], &
		    lp, ips);
	    chkspb_(&c__2, &c__1, idps, &c__0, inps, &c__5, &icp[1], &iep[1], 
		    &ipe[6], &iee[5], &lp, ips, &ipbad, &flagbnds);
	    if (flagbnds) {
		goto L1000;
	    }
	    chkspf_(&c__1, &ipb, &c__5, &icp[1], &iep[1], &ipe[6], &iee[5], &
		    lp, ips);
	    chkspb_(&c__2, &c__1, idps, &c__0, inps, &c__5, &icp[1], &iep[1], 
		    &ipe[6], &iee[5], &lp, ips, &ipbad, &flagbnds);
	    if (flagbnds) {
		goto L1000;
	    }
	}
    }
/*  ...   checking for the new boundary point */
    if (ifxnode_(&icps, &c__4)) {
	if (ifxnode_(status, &c__2)) {
	    flagbnds = TRUE_;
	    i__1 = lr;
	    for (n = 1; n <= i__1; ++n) {
		ips[n - 1] = irs[n * 3 - 2];
		if (ifxnode_(&icp[ipt], &c__64)) {
		    flagbnds = FALSE_;
		}
	    }
	} else {
	    flagbnds = FALSE_;
	}
	if (flagbnds) {
	    ips[lr] = ipa;
	    ips[lr + 1] = ipb;
	    lp = lr + 2;
	    chkspb_(&c__1, &c__1, idps, &c__0, inps, &c__5, &icp[1], &iep[1], 
		    &ipe[6], &iee[5], &lp, ips, &ipbad, &flagbnds);
	    if (flagbnds) {
		goto L1000;
	    }
	}
/*  ...  checking for boundary elements */
	if (flagfbe && ifxnode_(&icp[ipa], &c__4) && ifxnode_(&icp[ipb], &
		c__4)) {
	    i__1 = lr;
	    for (n = 1; n <= i__1; ++n) {
		ip1 = irs[n * 3 - 2];
		ip2 = irs[n * 3 - 1];
		if (ifxnode_(&icp[ip1], &c__64)) {
		    goto L100;
		}
		if (ifxnode_(&icp[ip2], &c__64)) {
		    goto L100;
		}
		goto L1000;
L100:
		;
	    }
	}
    }
/* ... analyzing curvilinear and plane faces */
    pntadd_(&ipc, np, maxp, &icp[1], &xyp[4], &hesp[7], &detg[1], &iholp[1], &
	    icps, xyps, hesps, &detgs);
    lfold = lf;
    i__1 = lpend;
    for (n = 1; n <= i__1; ++n) {
	if (n < lpend) {
	    ip1 = irs[n * 3 - 2];
	} else {
	    ip1 = irs[(n - 1) * 3 - 1];
	}
	i__2 = lfold;
	for (m = 1; m <= i__2; ++m) {
	    if (check33_(&ipa, &ipb, &ip1, &ipfs[(m << 2) - 4], &ipfs[(m << 2)
		     - 3], &ipfs[(m << 2) - 2])) {
		nf1 = m;
		iclrf = ipf[(ifs[nf1 - 1] << 2) + 4];
		ipfs[(nf1 << 2) - 4] = ipa;
		ipfs[(nf1 << 2) - 3] = ip1;
		ipfs[(nf1 << 2) - 2] = ipc;
		ipfs[(nf1 << 2) - 1] = iclrf;
		++lf;
		ipfs[(lf << 2) - 4] = ipb;
		ipfs[(lf << 2) - 3] = ip1;
		ipfs[(lf << 2) - 2] = ipc;
		ipfs[(lf << 2) - 1] = iclrf;
		icfs[lf - 1] = nf1;
		facadd_(&ifs[lf - 1], nf, maxf, &iholf[1]);
	    }
	}
    }
/* ... updating the grid */
    *flag__ = TRUE_;
    i__1 = lf;
    for (n = lfold + 1; n <= i__1; ++n) {
	nf1 = icfs[n - 1];
	facupd_(&nf1, &ipf[5], ifs, ipfs);
	facupd_(&n, &ipf[5], ifs, ipfs);
    }
    i__1 = lr;
    for (n = 1; n <= i__1; ++n) {
	ip1 = irs[n * 3 - 2];
	ip2 = irs[n * 3 - 1];
	ne1 = irs[n * 3 - 3];
	ide = ipes[ne1 * 5 - 1];
	ipes[ne1 * 5 - 5] = ipa;
	ipes[ne1 * 5 - 4] = ip1;
	ipes[ne1 * 5 - 3] = ip2;
	ipes[ne1 * 5 - 2] = ipc;
	ipes[ne1 * 5 - 1] = ide;
	ne2 = leold + n;
	ipes[ne2 * 5 - 5] = ipb;
	ipes[ne2 * 5 - 4] = ip1;
	ipes[ne2 * 5 - 3] = ip2;
	ipes[ne2 * 5 - 2] = ipc;
	ipes[ne2 * 5 - 1] = ide;
	lstupd_(ne, &l1e[3], nl2, &l2e[1], &nstep[1], &qe[1], &ies[ne1 - 1], &
		qes[ne1 - 1]);
	eledel_(&ies[ne1 - 1], &ipe[6], &iee[5]);
	eleadd_(ne, maxe, &ihole[1]);
	lstadd_(ne, &l1e[3], nl2, &l2e[1], &nstep[1], &ihole[1], &qe[1], &qes[
		ne2 - 1], &ies[ne2 - 1]);
	eledel_(&ies[ne2 - 1], &ipe[6], &iee[5]);
    }
    i__1 = lr;
    for (n = 1; n <= i__1; ++n) {
	ne1 = irs[n * 3 - 3];
	eleupd_(&ne1, &iep[1], &ipe[6], &ife[5], &iee[5], &lf, &le, ifs, ies, 
		ipfs, ipes);
	ne2 = leold + n;
	eleupd_(&ne2, &iep[1], &ipe[6], &ife[5], &iee[5], &lf, &le, ifs, ies, 
		ipfs, ipes);
    }
L1000:
    return 0;
L9000:
    errmes_(&c__1007, "insrtP", "local parameter MaxS is small", (ftnlen)6, (
	    ftnlen)29);
    return 0;
} /* insrtp_ */

