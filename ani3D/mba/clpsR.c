/* clpsR.f -- translated by f2c (version 20090411).
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

static integer c__1 = 1;
static integer c__128 = 128;
static integer c__8 = 8;
static integer c__512 = 512;
static integer c__64 = 64;
static integer c__2 = 2;
static integer c__4 = 4;
static integer c__0 = 0;

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int clpsr_(integer *iwr, integer *iwe, integer *
	np, integer *ne, doublereal *xyp, integer *ipe, doublereal *hstar, 
	integer *icp, integer *iep, integer *ife, integer *iee, integer *l1e, 
	integer *l2e, integer *nl2, integer *nstep, integer *iholp, integer *
	ihole, integer *status, doublereal *hesp, doublereal *rquality, 
	doublereal *detg, doublereal *qe, I_fp metricfunction, logical *
	flaganalytic, integer *lfu, integer *leu, integer *ifu, integer *ieu, 
	integer *ipfu, integer *ipeu, doublereal *qeu, integer *npw, integer *
	new__, doublereal *xypw, doublereal *hespw, integer *ipew, integer *
	milintrp, integer *mrlintrp, integer *ise, doublereal *rse, integer *
	icontrol, logical *flag__)
{
    /* Initialized data */

    static integer ipr[12]	/* was [2][6] */ = { 1,2,1,3,1,4,2,3,2,4,3,4 }
	    ;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static logical flagbnds;
    static integer i__, n;
    static doublereal v;
    static integer i1, i2, i3, i4;
    static doublereal w1, w2;
    static integer ie, le, lf, ip[5], lp;
    static logical flagorient;
    static integer ip1, ip2, ldh, ipa, ipb, ipc, ipf[1000], ifs[1000], ies[
	    1000];
    static doublereal qes[1000];
    static integer ids[2], ins[2], iet, ipt, nxy, icp1, icp2, icps, ipes[5000]
	    	/* was [5][1000] */, ipfs[4000]	/* was [4][1000] */;
    static doublereal xypo[3], xyps[3];
    static integer ipbad, ifncs;
    static doublereal detgs, hesps[6];
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
/* group (Local functions) */
/* group (Local variables) */
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
    --iholp;
    --nstep;
    --l2e;
    l1e -= 3;
    iee -= 5;
    ife -= 5;
    --iep;
    --icp;
    ipe -= 6;
    xyp -= 4;

    /* Function Body */
/* ================================================================ */
    *flag__ = FALSE_;
    ip[0] = 1;
    ip[1] = 2;
    ip[2] = 3;
    ip[3] = 4;
    ip[4] = 1;
    i1 = ipr[(*iwr << 1) - 2];
    i2 = ipr[(*iwr << 1) - 1];
    ip1 = ipe[i1 + *iwe * 5];
    ip2 = ipe[i2 + *iwe * 5];
/* ... check the case when the edge can not be collapsed */
    icp1 = icp[ip1];
    icp2 = icp[ip2];
    icps = minclr_(&icp1, &icp2);
    if (ifxnode_(&icps, &c__1)) {
	goto L1000;
    }
    if (ifxnode_(&icps, &c__128)) {
	goto L1000;
    }
/* ... add miscaleneous restrictions including missing algorithms */
    if (ifxnode_(&icp1, &c__8) || ifxnode_(&icp2, &c__8)) {
	goto L1000;
    }
/* ... number of inverted elements (the code is missing) */
    flagtm = ifxnode_(status, &c__512);
    if (flagtm) {
	goto L1000;
/*        nBad = 0 */
/*        If(qEs(iE1).LE.0D0) nBad = nBad + 1 */
/*        If(qEs(iE2).LE.0D0) nBad = nBad + 1 */
    }
/* ... finding a point in which we collapse the edge */
    w1 = .5;
    w2 = .5;
    ifncs = 0;
    if (icp1 == 64 && icp2 == 64) {
	icps = 64;
	for (i__ = 1; i__ <= 3; ++i__) {
	    xyps[i__ - 1] = xyp[i__ + ip1 * 3] * w1 + xyp[i__ + ip2 * 3] * w2;
	}
	for (i__ = 1; i__ <= 6; ++i__) {
	    hesps[i__ - 1] = hesp[i__ + ip1 * 6] * w1 + hesp[i__ + ip2 * 6] * 
		    w2;
	}
	if (! (*flaganalytic)) {
	    ldh = 6;
	    nxy = 1;
	    lintrp3d_(new__, &ipew[5], npw, &xypw[4], &ldh, &hespw[7], &nxy, 
		    xyps, hesps, &ise[1], milintrp, &rse[1], mrlintrp, 
		    icontrol);
	} else {
	    scaleback_(xyps, xypo);
	    iniq_analytic__(&c__1, xypo, (I_fp)metricfunction, hesps);
	}
	caldet_(hesps, &detgs);
    } else if (ifxnode_(&icp1, &c__8) && icp2 == 64 || ifxnode_(&icp1, &
	    c__128) && icp2 == 64 || ifxnode_(&icp1, &c__1)) {
	icps = icp1;
	for (i__ = 1; i__ <= 3; ++i__) {
	    xyps[i__ - 1] = xyp[i__ + ip1 * 3];
	}
	for (i__ = 1; i__ <= 6; ++i__) {
	    hesps[i__ - 1] = hesp[i__ + ip1 * 6];
	}
	detgs = detg[ip1];
    } else if (ifxnode_(&icp2, &c__8) && icp1 == 64 || ifxnode_(&icp2, &
	    c__128) && icp1 == 64 || ifxnode_(&icp2, &c__1)) {
	icps = icp2;
	for (i__ = 1; i__ <= 3; ++i__) {
	    xyps[i__ - 1] = xyp[i__ + ip2 * 3];
	}
	for (i__ = 1; i__ <= 6; ++i__) {
	    hesps[i__ - 1] = hesp[i__ + ip2 * 6];
	}
	detgs = detg[ip2];
/* ...    changing order of points to be consistent with the previous case */
	i__ = ip2;
	ip2 = ip1;
	ip1 = i__;
    } else {
	goto L1000;
    }
/* ... saving the initial superelement structure */
    copyse_(lfu, leu, &ifu[1], &ieu[1], &ipfu[5], &ipeu[6], &qeu[1], &lf, &le,
	     ifs, ies, ipfs, ipes, qes);
/* ... making a virtual evaluation of the quality */
    i__1 = le;
    for (n = 1; n <= i__1; ++n) {
	ie = ies[n - 1];
	if (check1j_(&ip1, &ipe[ie * 5 + 1]) && check1j_(&ip2, &ipe[ie * 5 + 
		1])) {
	    ies[n - 1] = -ies[n - 1];
	    goto L10;
	}
	for (i1 = 1; i1 <= 4; ++i1) {
	    if (ipes[i1 + n * 5 - 6] == ip1 || ipes[i1 + n * 5 - 6] == ip2) {
		i2 = ip[i1];
		i3 = ip[i2];
		i4 = ip[i3];
		ipes[i1 + n * 5 - 6] = ip1;
		ipa = ipes[i2 + n * 5 - 6];
		ipb = ipes[i3 + n * 5 - 6];
		ipc = ipes[i4 + n * 5 - 6];
		calqe_(hesps, xyps, &hesp[ipa * 6 + 1], &xyp[ipa * 3 + 1], &
			hesp[ipb * 6 + 1], &xyp[ipb * 3 + 1], &hesp[ipc * 6 + 
			1], &xyp[ipc * 3 + 1], hstar, &qes[n - 1], &v);
		if (qes[n - 1] <= *rquality) {
		    goto L1000;
		}
/*  ...  check orientation of neigboors */
		iet = iee[i1 + (ie << 2)];
		if (iet > 0) {
		    chksotet_(&ip1, &ip2, &ipa, &ipb, xyps, &xyp[4], &ipe[iet 
			    * 5 + 1], &v, &flagorient);
		    if (! flagorient) {
			goto L1000;
		    }
		}
		iet = iee[i3 + (ie << 2)];
		if (iet > 0) {
		    chksotet_(&ip1, &ip2, &ipb, &ipc, xyps, &xyp[4], &ipe[iet 
			    * 5 + 1], &v, &flagorient);
		    if (! flagorient) {
			goto L1000;
		    }
		}
		iet = iee[i4 + (ie << 2)];
		if (iet > 0) {
		    chksotet_(&ip1, &ip2, &ipc, &ipa, xyps, &xyp[4], &ipe[iet 
			    * 5 + 1], &v, &flagorient);
		    if (! flagorient) {
			goto L1000;
		    }
		}
		goto L10;
	    }
	}
	ies[n - 1] = 0;
L10:
	;
    }
/* ... checking for boundary elements */
    if (ifxnode_(status, &c__1)) {
	i__1 = le;
	for (n = 1; n <= i__1; ++n) {
	    ie = ies[n - 1];
	    if (ie <= 0) {
		goto L20;
	    }
	    for (i__ = 1; i__ <= 4; ++i__) {
		ipt = ipes[i__ + n * 5 - 6];
		if (ifxnode_(&icp[ipt], &c__64)) {
		    goto L20;
		}
	    }
	    goto L1000;
L20:
	    ;
	}
    }
/* ... checking for surrounding points (not ICP1 but ICP(iP1)) */
/* ... We need to study only pathes which were ending at iP2. */
/* ... Thus, each boundary point connected to iP2 should have */
/* ... another interior point neighboor. */
    if (ifxnode_(status, &c__2)) {
	if (ifxnode_(&icp[ip1], &c__4) && icp[ip2] == 64) {
	    ids[0] = ip2;
	    ids[1] = ip1;
	    chkspf_(&c__2, &ip2, &c__8, &icp[1], &iep[1], &ipe[6], &iee[5], &
		    lp, ipf);
	    chkspb_(&c__2, &c__1, ids, &c__0, ins, &c__8, &icp[1], &iep[1], &
		    ipe[6], &iee[5], &lp, ipf, &ipbad, &flagbnds);
	    if (flagbnds) {
		goto L1000;
	    }
	}
    }
/* ... checking for inverted elements */
    if (flagtm) {
/*        Do n = 1, lE */
/*           If(iEs(n).GE.0 .AND. qEs(n).GT.0D0) Then */
/*              Call updQb(n, lE, iEs, XYP, IPEs, qEs) */
/*           End if */
/*        End do */

/*        mBad = 0 */
/*        Do n = 1, lE */
/*           If(iEs(n).GE.0 .AND. qEs(n).LE.0D0) mBad = mBad + 1 */
/*        End do */

/*        If(mBad.GE.nBad) Goto 1000 */
    }
/* ... updating the grid */
    *flag__ = TRUE_;
    pntupd_(&ip1, &icp[1], &xyp[4], &hesp[7], &detg[1], &icps, xyps, hesps, &
	    detgs);
    pntdel_(&ip2, np, &icp[1], &iholp[1]);
    i__1 = le;
    for (n = 1; n <= i__1; ++n) {
	iet = ies[n - 1];
	if (iet < 0) {
	    iet = -iet;
	    lstdel_(ne, &l1e[3], nl2, &l2e[1], &nstep[1], &ihole[1], &qe[1], &
		    iet);
	    eledel_(&iet, &ipe[6], &iee[5]);
	} else if (iet > 0) {
	    lstupd_(ne, &l1e[3], nl2, &l2e[1], &nstep[1], &qe[1], &iet, &qes[
		    n - 1]);
	    eleupd_(&n, &iep[1], &ipe[6], &ife[5], &iee[5], &lf, &le, ifs, 
		    ies, ipfs, ipes);
	}
    }
L1000:
    return 0;
} /* clpsr_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int chksotet_(integer *ip1, integer *ip2, 
	integer *ipa, integer *ipb, doublereal *xyps, doublereal *xyp, 
	integer *ipe, doublereal *v1, logical *flag__)
{
    static integer i__;
    static doublereal v2;
    static integer ipt;

/* ================================================================ */
/* Routine checks orientation of {iP1, iPa, iPb, *} and IPE. */

/* flag = .TRUE. if the tets are oriented correctly. */
/* flag = .TRUE. if the second tet has both iP1 and iP2 */
/* ================================================================ */
/* ================================================================ */
    /* Parameter adjustments */
    --ipe;
    xyp -= 4;
    --xyps;

    /* Function Body */
    for (i__ = 1; i__ <= 4; ++i__) {
	ipt = ipe[i__];
	if (! check13_(&ipt, ip1, ipa, ipb) && ipt != *ip2) {
	    v2 = calvol_(&xyps[1], &xyp[*ipa * 3 + 1], &xyp[*ipb * 3 + 1], &
		    xyp[ipt * 3 + 1]);
	    *flag__ = *v1 * v2 < 0.;
	    goto L9000;
	}
    }
    *flag__ = TRUE_;
L9000:
    return 0;
} /* chksotet_ */

