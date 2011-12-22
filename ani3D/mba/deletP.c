/* deletP.f -- translated by f2c (version 20090411).
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

#include "list.h"
#include "makQ.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c__128 = 128;
static integer c__1000 = 1000;
static integer c__512 = 512;
static doublereal c_b8 = 1.;
static integer c__2 = 2;
static integer c__64 = 64;
static integer c__4 = 4;
static integer c__1007 = 1007;

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int deletp_(integer *iwp, integer *iwe, integer *
	np, integer *nf, integer *maxf, integer *ne, integer *maxe, 
	doublereal *xyp, integer *ipf, integer *ipe, doublereal *hstar, 
	integer *icp, integer *iep, integer *ife, integer *iee, integer *l1e, 
	integer *l2e, integer *nl2, integer *nstep, integer *iholp, integer *
	iholf, integer *ihole, integer *status, doublereal *hesp, doublereal *
	rquality, doublereal *detg, doublereal *qe, integer *lfu, integer *
	leu, integer *ifu, integer *ieu, integer *ipfu, integer *ipeu, 
	doublereal *qeu, logical *flag__)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static logical flagface, flagedge, flagbnds, flagface2;
    static integer i__, k, m, n;
    static doublereal v;
    static integer i1, i2, i3;
    static logical l1;
    static integer j1, j2, j3;
    static doublereal oldvolume, newvolume;
    static integer le, lf, ip[5], lp, ls;
    static logical flagorient;
    static integer ip1, lp1, lp2, ip2, ip3, ip4, ied[1], lde, ldf, ned, ipa, 
	    ipb, ipd, ies[1000], ifs[1000], ldp;
    static doublereal qes[1000];
    static integer nft, ios[1000], ipt[3], net, lnp, ipu, npt[3], ift, iet, 
	    ip1s[1000], ip2s[1000], mbad, nbad, ides[1000], iref[4], idfs[
	    1000], icnt, ipes[5000]	/* was [5][1000] */, ipfs[4000]	/* 
	    was [4][1000] */, idps[2000]	/* was [2][1000] */, iess[
	    1000], inps[2000]	/* was [2][1000] */, ipss[3000]	/* was [3][
	    1000] */, icf1s, icf2s, icp2s[1000], leadd, ipbad, leold, lfold, 
	    lenew, narms;
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
/* group (Local functions) */
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
    ipd = ipe[*iwp + *iwe * 5];
/* ... checking the case when deleting is impossible */
    if (ifxnode_(&icp[ipd], &c__1)) {
	goto L1000;
    }
    if (ifxnode_(&icp[ipd], &c__128)) {
	goto L1000;
    }
/* ... analyzing the superelement entries */
    maksf_(&ipd, &lf, ifs, &le, ies, ipes, &ipf[5], &ife[5], &c__1000, &lp1, 
	    ip1s, &lp2, ip2s, icp2s, &ls, ipss, iess, &ldf, idfs, &lde, ides, 
	    &flagface, &flagedge);
/* ... intersection of two interfaces at the point has more than 3 edges */
    if (flagedge && lp1 >= 3) {
	goto L1000;
    }
    if (lp2 + 1 >= *np) {
	goto L1000;
    }
    ldp = lp2;
    i__1 = ldp;
    for (n = 1; n <= i__1; ++n) {
	idps[(n << 1) - 2] = ipd;
	idps[(n << 1) - 1] = ip2s[n - 1];
    }
    lfold = lf;
    i__1 = ldf;
    for (n = 1; n <= i__1; ++n) {
	nft = idfs[n - 1];
	ifs[nft - 1] = -ifs[nft - 1];
    }
    i__1 = lde;
    for (n = 1; n <= i__1; ++n) {
	net = ides[n - 1];
	ies[net - 1] = -ies[net - 1];
    }
/* ... checking the number of inverted elements */
    flagtm = ifxnode_(status, &c__512);
    if (flagtm) {
	nbad = 0;
	i__1 = lde;
	for (n = 1; n <= i__1; ++n) {
	    net = ides[n - 1];
	    if (qes[net - 1] <= 0.) {
		++nbad;
	    }
	}
    }
    flagtm = flagtm && nbad > 0;
/* ... making a virtual evaluation of the quality */
    oldvolume = 0.;
    i__1 = lde;
    for (n = 1; n <= i__1; ++n) {
	ned = ides[n - 1];
	ip1 = ipes[ned * 5 - 5];
	ip2 = ipes[ned * 5 - 4];
	ip3 = ipes[ned * 5 - 3];
	ip4 = ipes[ned * 5 - 2];
	v = calvol_(&xyp[ip1 * 3 + 1], &xyp[ip2 * 3 + 1], &xyp[ip3 * 3 + 1], &
		xyp[ip4 * 3 + 1]);
	oldvolume += dabs(v);
	ios[n - 1] = (integer) d_sign(&c_b8, &v);
/* L20: */
    }
    i__1 = lp1;
    for (n = 1; n <= i__1; ++n) {
	ipa = ip1s[n - 1];
/*  ...  checking for the orientation */
	i__2 = ls;
	for (k = 1; k <= i__2; ++k) {
	    ip1 = ipss[k * 3 - 3];
	    ip2 = ipss[k * 3 - 2];
	    ip3 = ipss[k * 3 - 1];
	    ied[0] = (i__3 = ies[ides[k - 1] - 1], dabs(i__3));
	    if (check13_(&ipa, &ip1, &ip2, &ip3)) {
/*  ...  checking in details this face */
		i__3 = ls;
		for (m = 1; m <= i__3; ++m) {
		    if (m == k) {
			goto L30;
		    }
		    icnt = 0;
		    for (i__ = 1; i__ <= 3; ++i__) {
			if (ipss[i__ + k * 3 - 4] != ipa) {
			    if (check13_(&ipss[i__ + k * 3 - 4], &ipss[m * 3 
				    - 3], &ipss[m * 3 - 2], &ipss[m * 3 - 1]))
				     {
				++icnt;
			    }
			}
		    }
		    if (icnt == 2) {
			goto L40;
		    }
L30:
		    ;
		}
		goto L80;
	    }
	    chkso_(&ipd, &xyp[ipa * 3 + 1], &xyp[4], &ipe[6], &c__1, ied, &
		    ios[k - 1], &flagorient);
	    if (! flagorient) {
		goto L80;
	    }
L40:
	    ;
	}
/*  ...  trying for the quality */
	leadd = 0;
	newvolume = 0.;
	i__2 = ls;
	for (k = 1; k <= i__2; ++k) {
	    ip1 = ipss[k * 3 - 3];
	    ip2 = ipss[k * 3 - 2];
	    ip3 = ipss[k * 3 - 1];
	    if (check13_(&ipa, &ip1, &ip2, &ip3)) {
		goto L50;
	    }
	    ++leadd;
	    lenew = le + leadd;
	    if (lenew > 1000) {
		goto L9000;
	    }
	    ies[lenew - 1] = 0;
	    calqe_(&hesp[ipa * 6 + 1], &xyp[ipa * 3 + 1], &hesp[ip1 * 6 + 1], 
		    &xyp[ip1 * 3 + 1], &hesp[ip2 * 6 + 1], &xyp[ip2 * 3 + 1], 
		    &hesp[ip3 * 6 + 1], &xyp[ip3 * 3 + 1], hstar, &qes[lenew 
		    - 1], &v);
	    if (qes[lenew - 1] <= *rquality) {
		goto L80;
	    }
	    ipes[lenew * 5 - 5] = ipa;
	    ipes[lenew * 5 - 4] = ip1;
	    ipes[lenew * 5 - 3] = ip2;
	    ipes[lenew * 5 - 2] = ip3;
	    ipes[lenew * 5 - 1] = ipes[iess[k - 1] * 5 - 1];
	    newvolume += dabs(v);
L50:
	    ;
	}
	leold = le;
	le = lenew;
	goto L100;
L80:
	;
    }
    goto L1000;
L100:
/* ... checking for the global volume */
    if (! flagtm && (d__1 = oldvolume - newvolume, dabs(d__1)) > oldvolume * 
	    1e-10) {
/*        Write(*, 5000) oldVolume, newVolume */
	goto L1000;
    }
/* ... checking for surrounding points */
    lnp = 0;
    i__1 = lp2;
    for (n = 1; n <= i__1; ++n) {
	ipb = ip2s[n - 1];
	if (ipa != ipb) {
	    ++lnp;
	    if (lnp > 1000) {
		goto L9000;
	    }
	    inps[(lnp << 1) - 2] = ipa;
	    inps[(lnp << 1) - 1] = ipb;
	}
    }
    narms = 0;
    if (ifxnode_(status, &c__2)) {
	if (ifxnode_(&icp[ipd], &c__64)) {
	    narms = 2;
	} else if (ifxnode_(&icp[ipd], &c__4)) {
	    narms = 1;
	}
    }
    if (narms > 0) {
	chkspf_(&narms, &ipd, &c__4, &icp[1], &iep[1], &ipe[6], &iee[5], &lp, 
		ios);
	chkspb_(&c__2, &ldp, idps, &lnp, inps, &c__4, &icp[1], &iep[1], &ipe[
		6], &iee[5], &lp, ios, &ipbad, &flagbnds);
	if (flagbnds) {
	    goto L1000;
	}
    }
/* ... checking for boundary elements */
    if (ifxnode_(status, &c__1)) {
	i__1 = le;
	for (n = leold + 1; n <= i__1; ++n) {
	    for (i__ = 1; i__ <= 4; ++i__) {
		ipu = ipes[i__ + n * 5 - 6];
		if (ifxnode_(&icp[ipu], &c__64)) {
		    goto L800;
		}
	    }
	    goto L1000;
L800:
	    ;
	}
    }
/* ... checking for inverted elements */
    if (flagtm) {
	i__1 = le;
	for (n = 1; n <= i__1; ++n) {
	    if (ies[n - 1] >= 0) {
		updqb_(&n, &le, ies, &xyp[4], ipes, qes);
	    }
	}
	mbad = 0;
	i__1 = le;
	for (n = leold; n <= i__1; ++n) {
	    if (ies[n - 1] >= 0 && qes[n - 1] <= 0.) {
		++mbad;
	    }
	}
	if (mbad >= nbad) {
	    goto L1000;
	}
    }
/* ... analyzing curvilinear and plane faces */
    if (flagface) {
	i__1 = le;
	for (n = leold + 1; n <= i__1; ++n) {
	    for (i1 = 1; i1 <= 4; ++i1) {
		i2 = ip[i1];
		i3 = ip[i2];
		ipt[0] = ipes[i1 + n * 5 - 6];
		ipt[1] = ipes[i2 + n * 5 - 6];
		ipt[2] = ipes[i3 + n * 5 - 6];
		for (i__ = 1; i__ <= 3; ++i__) {
		    findse_(&lp2, ip2s, &ipt[i__ - 1], &npt[i__ - 1]);
		    if (npt[i__ - 1] <= 0) {
			goto L200;
		    }
		}
/* Computing MIN */
		i__2 = icp2s[npt[0] - 1], i__3 = icp2s[npt[1] - 1];
		icf1s = min(i__2,i__3);
/* Computing MIN */
		i__2 = icf1s, i__3 = icp2s[npt[2] - 1];
		icf1s = min(i__2,i__3);
/* Computing MAX */
		i__2 = icp2s[npt[0] - 1], i__3 = icp2s[npt[1] - 1];
		icf2s = max(i__2,i__3);
/* Computing MAX */
		i__2 = icf2s, i__3 = icp2s[npt[2] - 1];
		icf2s = max(i__2,i__3);
		flagface2 = FALSE_;
		if (icf1s > 0) {
		    if (icf1s == icf2s) {
			flagface2 = TRUE_;
		    }
		    if (lp1 == 2) {
			icnt = 0;
			for (j1 = 1; j1 <= 3; ++j1) {
			    if (ipt[j1 - 1] == ip1s[0] || ipt[j1 - 1] == ip1s[
				    1]) {
				++icnt;
				j2 = iref[j1];
				j3 = iref[j2];
				l1 = icp2s[npt[j2 - 1] - 1] == icp2s[npt[j3 - 
					1] - 1];
			    }
			}
			if (icnt == 2) {
			    flagface2 = TRUE_;
			}
			if (icnt == 1 && l1) {
			    flagface2 = TRUE_;
			}
		    }
		}
		if (flagface2) {
/*                 Do m = lFold + 1, lF */
		    i__2 = lf;
		    for (m = 1; m <= i__2; ++m) {
			if (check33_(ipt, &ipt[1], &ipt[2], &ipfs[(m << 2) - 
				4], &ipfs[(m << 2) - 3], &ipfs[(m << 2) - 2]))
				 {
			    goto L200;
			}
		    }
		    ++lf;
		    if (le > 1000) {
			goto L9000;
		    }
		    facadd_(&ifs[lf - 1], nf, maxf, &iholf[1]);
		    ipfs[(lf << 2) - 4] = ipt[0];
		    ipfs[(lf << 2) - 3] = ipt[1];
		    ipfs[(lf << 2) - 2] = ipt[2];
		    ipfs[(lf << 2) - 1] = icf1s;
		}
L200:
		;
	    }
	}
    }
/* ... updating the grid */
    *flag__ = TRUE_;
    pntdel_(&ipd, np, &icp[1], &iholp[1]);
    if (flagface) {
	i__1 = lfold;
	for (n = 1; n <= i__1; ++n) {
	    ift = ifs[n - 1];
	    if (ift <= 0) {
		ift = -ift;
		facdel_(&ift, nf, &ipf[5], &iholf[1]);
	    }
	}
	i__1 = lf;
	for (n = lfold + 1; n <= i__1; ++n) {
	    facupd_(&n, &ipf[5], ifs, ipfs);
	}
    }
    i__1 = le;
    for (n = 1; n <= i__1; ++n) {
	iet = ies[n - 1];
	if (iet == 0) {
	    eleadd_(ne, maxe, &ihole[1]);
	    lstadd_(ne, &l1e[3], nl2, &l2e[1], &nstep[1], &ihole[1], &qe[1], &
		    qes[n - 1], &ies[n - 1]);
	    eledel_(&ies[n - 1], &ipe[6], &iee[5]);
	} else if (iet < 0) {
	    iet = -iet;
	    lstdel_(ne, &l1e[3], nl2, &l2e[1], &nstep[1], &ihole[1], &qe[1], &
		    iet);
	    eledel_(&iet, &ipe[6], &iee[5]);
	}
    }
    i__1 = le;
    for (n = 1; n <= i__1; ++n) {
	if (ies[n - 1] > 0) {
	    eleupd_(&n, &iep[1], &ipe[6], &ife[5], &iee[5], &lf, &le, ifs, 
		    ies, ipfs, ipes);
	}
    }
/* L5000: */
L1000:
    return 0;
L9000:
    errmes_(&c__1007, "deletP", "local parameter MaxS is small", (ftnlen)6, (
	    ftnlen)29);
    return 0;
} /* deletp_ */

