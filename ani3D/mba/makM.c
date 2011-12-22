/* makM.f -- translated by f2c (version 20090411).
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

#include "makM.h"
#include "makQ.h"

/* Common Block Declarations */

struct {
    doublereal refxyp[3], scaxyp[3];
} anicrv_;

#define anicrv_1 anicrv_

/* Table of constant values */

static integer c__4 = 4;
static integer c__5 = 5;
static integer c__3 = 3;
static integer c__0 = 0;
static integer c__1011 = 1011;
static integer c__1 = 1;
static integer c__128 = 128;
static integer c__1010 = 1010;
static integer c__8 = 8;
static integer c__64 = 64;
static integer c__2 = 2;
static integer c__16 = 16;

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int makm_(integer *np, integer *nf, integer *
	maxf, integer *ne, integer *maxe, doublereal *xyp, integer *ipf, 
	integer *ipe, integer *icp, integer *ipp, integer *iep, integer *ife, 
	integer *iee, integer *iholp, integer *iholf, integer *ihole, integer 
	*iepw, integer *nepw, integer *status, integer *npv, integer *nfv, 
	integer *nev, integer *ipv, integer *ifv, integer *iev, integer *ierr)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, k, n, i1, i2, i3, i4, ie, if__, le, mf, ip[5], ie2, 
	    ip1, ip2, ip3, ip4, iet, ift, ipt, mat1, mat2, kmax, ipfs[4], 
	    mlist[200]	/* was [2][100] */, icface;
    static logical flagfbf, flagbnd;

/* ================================================================ */
/* group (M) */
/* group (Dev) */
/* group (ERR) */
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
/* ================================================================ */
/* Routine analyses the initial geometry and creates auxiliary */
/* cross-refrences. */

/* Pre-conditions:  1. connectivity structure {IPE(5, *), XYP(3, *)} */

/* Post-conditions: 1. additional mesh structures, IEP, IFE, IEE */
/*                  2. coloring mesh point according to given */
/*                     surface and volume colors, and lists of */
/*                     fixed points, triangles, and tetrahedra. */
/* ================================================================ */
/* group (M) */
/*     Integer MaxF, MaxE */
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
/* group (Dev) */
/* ================================================================ */
/* group (Functions) */
/* group (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --iev;
    --ifv;
    --ipv;
    --nepw;
    --iepw;
    --ihole;
    --iholf;
    --iholp;
    iee -= 5;
    ife -= 5;
    --iep;
    --ipp;
    --icp;
    ipe -= 6;
    ipf -= 5;
    xyp -= 4;

    /* Function Body */
    *ierr = 0;
    ip[0] = 1;
    ip[1] = 2;
    ip[2] = 3;
    ip[3] = 4;
    ip[4] = 1;
    iholp[1] = 0;
    iholf[1] = 0;
    ihole[1] = 0;
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	icp[n] = 0;
    }
/* ... create an auxiliary structure */
    backreferences_(np, ne, &c__4, &c__5, &ipe[6], &nepw[1], &iepw[1]);
/* ... create IEE & IEP */
    i__1 = *maxe;
    for (n = 1; n <= i__1; ++n) {
	for (i__ = 1; i__ <= 4; ++i__) {
	    iee[i__ + (n << 2)] = 0;
	    ife[i__ + (n << 2)] = 0;
	}
    }
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	for (i1 = 1; i1 <= 4; ++i1) {
	    i2 = ip[i1];
	    i3 = ip[i2];
	    ip1 = ipe[i1 + n * 5];
	    ip2 = ipe[i2 + n * 5];
	    ip3 = ipe[i3 + n * 5];
	    iep[ip1] = n;
	    if (cmpe_(&ip1, &ip2, &ip3, &iepw[1], &nepw[1], &n, &ie2)) {
		iee[i1 + (n << 2)] = ie2;
	    }
	}
    }
/* ... create an auxiliary structure */
    backreferences_(np, nf, &c__3, &c__4, &ipf[5], &nepw[1], &iepw[1]);
/* ... create IFE: basic, fictitious, material */
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	for (i1 = 1; i1 <= 4; ++i1) {
	    i2 = ip[i1];
	    i3 = ip[i2];
	    ip1 = ipe[i1 + n * 5];
	    ip2 = ipe[i2 + n * 5];
	    ip3 = ipe[i3 + n * 5];
	    if (cmpe_(&ip1, &ip2, &ip3, &iepw[1], &nepw[1], &c__0, &if__)) {
		ife[i1 + (n << 2)] = if__;
	    }
	}
    }
/* ... create missing boundaries */
    i__1 = *nf;
    for (n = 1; n <= i__1; ++n) {
	if (ipf[(n << 2) + 4] >= 900) {
	    errmes_(&c__1011, "makM", "reserved boundary identificator is us"
		    "ed", (ftnlen)4, (ftnlen)39);
	}
    }
    flagfbf = ifxnode_(status, &c__4);
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	for (i1 = 1; i1 <= 4; ++i1) {
	    if (iee[i1 + (n << 2)] == 0 && ife[i1 + (n << 2)] == 0) {
		i2 = ip[i1];
		i3 = ip[i2];
		ipfs[0] = ipe[i1 + n * 5];
		ipfs[1] = ipe[i2 + n * 5];
		ipfs[2] = ipe[i3 + n * 5];
		ipfs[3] = 900;
		facadd_(&if__, nf, maxf, &iholf[1]);
		facupd_(&c__1, &ipf[5], &if__, ipfs);
		ife[i1 + (n << 2)] = if__;
		if (flagfbf) {
		    for (i__ = 1; i__ <= 3; ++i__) {
			ipt = ipfs[i__ - 1];
			addxnode_(&icp[ipt], &c__128);
		    }
		}
	    }
	}
    }
/* ... create material boundaries */
    k = 0;
    kmax = 99;
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	for (i1 = 1; i1 <= 4; ++i1) {
	    ift = ife[i1 + (n << 2)];
	    iet = iee[i1 + (n << 2)];
	    if (ift == 0 && iet != 0) {
		mat1 = ipe[n * 5 + 5];
		mat2 = ipe[iet * 5 + 5];
		if (mat1 != mat2) {
/*  ...  search for this pair in the list of material interfaces */
		    i__2 = k;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			if (mlist[(i__ << 1) - 2] == mat1 && mlist[(i__ << 1) 
				- 1] == mat2 || mlist[(i__ << 1) - 2] == mat2 
				&& mlist[(i__ << 1) - 1] == mat1) {
			    icface = i__ + 900;
			    goto L1;
			}
		    }
/*  ...  make a new material interface */
		    ++k;
		    if (k > kmax) {
			errmes_(&c__1010, "makM", "not enough memory for mat"
				"erial faces", (ftnlen)4, (ftnlen)36);
		    }
		    mlist[(k << 1) - 2] = mat1;
		    mlist[(k << 1) - 1] = mat2;
		    icface = k + 900;
L1:
		    i2 = ip[i1];
		    i3 = ip[i2];
		    ipfs[0] = ipe[i1 + n * 5];
		    ipfs[1] = ipe[i2 + n * 5];
		    ipfs[2] = ipe[i3 + n * 5];
		    ipfs[3] = icface;
		    facadd_(&if__, nf, maxf, &iholf[1]);
		    facupd_(&c__1, &ipf[5], &if__, ipfs);
		    ife[i1 + (n << 2)] = if__;
		    for (i__ = 1; i__ <= 4; ++i__) {
			if (iee[i__ + (iet << 2)] == n) {
			    ife[i__ + (iet << 2)] = if__;
			}
		    }
		}
	    }
	}
    }
/* ... color the points (boundary points) */
/* ... new auxiliary structure accumulates all faces */
    backreferences_(np, nf, &c__3, &c__4, &ipf[5], &nepw[1], &iepw[1]);
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	for (i1 = 1; i1 <= 4; ++i1) {
	    if__ = ife[i1 + (n << 2)];
	    ie = iee[i1 + (n << 2)];
	    if (if__ != 0) {
		i2 = i1;
		for (k = 1; k <= 3; ++k) {
		    ip1 = ipe[i2 + n * 5];
		    addxnode_(&icp[ip1], &c__8);
		    if (ie == 0) {
			addxnode_(&icp[ip1], &c__4);
		    }
		    i2 = ip[i2];
		}
	    }
	}
    }
/* ... coloring the points (edge points) */
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	if (! ifxnode_(&icp[n], &c__4)) {
	    addxnode_(&icp[n], &c__64);
	}
	if (cmpp_(&n, &iepw[1], &nepw[1], &ipf[5])) {
	    addxnode_(&icp[n], &c__2);
	}
	if (cmpr_(&n, &iepw[1], &nepw[1], &xyp[4], &ipf[5])) {
	    addxnode_(&icp[n], &c__2);
	}
    }
/* ... coloring the points (cross points) */
    i2 = 0;
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	i1 = i2 + 1;
	i2 = nepw[n];
	if (ifxnode_(&icp[n], &c__8)) {
	    i__2 = i2 - i1 + 1;
	    if (crosspoint_(&xyp[4], &i__2, &iepw[i1], &n, &ipf[5])) {
		addxnode_(&icp[n], &c__1);
	    }
	}
    }
/* ... coloring the points (fix vertices) */
    i__1 = *npv;
    for (n = 1; n <= i__1; ++n) {
	ip1 = ipv[n];
	addxnode_(&icp[ip1], &c__1);
    }
/* ... color the points of fix faces and elements */
    i__1 = *nfv;
    for (n = 1; n <= i__1; ++n) {
	if__ = ifv[n];
	for (i__ = 1; i__ <= 3; ++i__) {
	    ip1 = ipf[i__ + (if__ << 2)];
	    addxnode_(&icp[ip1], &c__128);
	}
    }
    i__1 = *nev;
    for (n = 1; n <= i__1; ++n) {
	ie = iev[n];
	for (i__ = 1; i__ <= 4; ++i__) {
	    ip1 = ipe[i__ + ie * 5];
	    addxnode_(&icp[ip1], &c__128);
	}
    }
/* ... color the points (fixed boundary points) */
    if (ifxnode_(status, &c__16)) {
	i__1 = *nf;
	for (n = 1; n <= i__1; ++n) {
	    for (i__ = 1; i__ <= 3; ++i__) {
		ip1 = ipf[i__ + (n << 2)];
		addxnode_(&icp[ip1], &c__1);
	    }
	}
    }
/* ... marking points with a multi-connected superelement */
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	nepw[n] = 0;
    }
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	for (i__ = 1; i__ <= 4; ++i__) {
	    i1 = ipe[i__ + n * 5];
	    ++nepw[i1];
	}
    }
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	if (ifxnode_(&icp[n], &c__8)) {
	    maksp_(&n, &iep[1], &ipe[6], &iee[5], maxe, &le, &iepw[1]);
	    if (nepw[n] != le) {
		addxnode_(status, &c__64);
		addxnode_(&icp[n], &c__1);
	    }
	}
    }
/* ... coloring the T-points if any (case of inner boundaries) */
/* ... array IPF(4, *) is overloaded to mark only boundary edges */
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	for (i__ = 1; i__ <= 4; ++i__) {
	    ift = ife[i__ + (n << 2)];
	    if (ift > 0) {
		ipf[(ift << 2) + 4] = -ipf[(ift << 2) + 4];
	    }
	}
    }
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	iepw[n] = 0;
    }
    i__1 = *nf;
    for (n = 1; n <= i__1; ++n) {
	if (ipf[(n << 2) + 4] > 0) {
	    goto L100;
	}
	ipf[(n << 2) + 4] = -ipf[(n << 2) + 4];
	flagbnd = FALSE_;
	for (i__ = 1; i__ <= 3; ++i__) {
	    if (ipp[ipf[i__ + (n << 2)]] == 0) {
		flagbnd = TRUE_;
	    }
	}
	if (flagbnd) {
	    for (i__ = 1; i__ <= 3; ++i__) {
		iepw[ipf[i__ + (n << 2)]] = 1;
	    }
	}
L100:
	;
    }
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	if (ipp[n] > 0) {
	    addxnode_(&icp[n], &c__128);
	    if (iepw[n] == 0) {
		addxnode_(&icp[n], &c__64);
		delxnode_(&icp[n], &c__4);
	    }
	}
    }
/* ... delete isolated points */
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	iepw[n] = 0;
    }
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	for (i__ = 1; i__ <= 4; ++i__) {
	    iepw[ipe[i__ + n * 5]] = 1;
	}
    }
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	if (iepw[n] == 0) {
	    pntdel_(&n, np, &icp[1], &iholp[1]);
	}
    }
/* ... delete isolated faces */
    i__1 = *nf;
    for (n = 1; n <= i__1; ++n) {
	iepw[n] = 0;
    }
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	for (i__ = 1; i__ <= 4; ++i__) {
	    ift = ife[i__ + (n << 2)];
	    if (ift > 0) {
		iepw[ift] = 1;
	    }
	}
    }
    mf = *nf;
    i__1 = mf;
    for (n = 1; n <= i__1; ++n) {
	if (iepw[n] == 0) {
	    facdel_(&n, nf, &ipf[5], &iholf[1]);
	    i__2 = *nfv;
	    for (k = 1; k <= i__2; ++k) {
		if (ifv[k] == n) {
		    ifv[k] = ifv[*nfv];
		    --(*nfv);
		    goto L500;
		}
	    }
L500:
	    ;
	}
    }
/* ... order boundary faces clockwise looking from inside domain */
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	for (i1 = 1; i1 <= 4; ++i1) {
	    if__ = ife[i1 + (n << 2)];
	    if (if__ > 0) {
		i2 = ip[i1];
		i3 = ip[i2];
		i4 = ip[i3];
		ip4 = ipe[i4 + n * 5];
		ip1 = ipf[(if__ << 2) + 1];
		ip2 = ipf[(if__ << 2) + 2];
		ip3 = ipf[(if__ << 2) + 3];
		if (calvol_(&xyp[ip1 * 3 + 1], &xyp[ip2 * 3 + 1], &xyp[ip3 * 
			3 + 1], &xyp[ip4 * 3 + 1]) > 0.) {
		    swapii_(&ipf[(if__ << 2) + 1], &ipf[(if__ << 2) + 2]);
		}
	    }
	}
    }
    return 0;
} /* makm_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int updm_(integer *np, integer *nf, integer *ne, 
	doublereal *xyp, integer *ipf, integer *ipe, integer *icp, integer *
	ipp, integer *ife, integer *iee, integer *iholp, integer *iholf, 
	integer *ihole, integer *status, integer *npv, integer *nfv, integer *
	nev, integer *ipv, integer *ifv, integer *iev, doublereal *hesp, 
	doublereal *qe, integer *ipw)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, m, n, ic, ie, if__, ke, le, lf, kf, ip, lp, mp, 
	    iet, icfree;
    static logical flagdtf;

/* ================================================================ */
/* group (M) */
/* group (Dev) */
/* group (Q) */
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
/* ================================================================ */
/* Routine removes holes (deleted elements) from the mesh data */
/* structures. The optimal complexity algorithms are implemented. */
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
/* group (Dev) */
/* group (Q) */
/* group (W) */
/* ================================================================ */
/* group (Functions) */
/* ================================================================ */
    /* Parameter adjustments */
    --ipw;
    --qe;
    hesp -= 7;
    --iev;
    --ifv;
    --ipv;
    --ihole;
    --iholf;
    --iholp;
    iee -= 5;
    ife -= 5;
    --ipp;
    --icp;
    ipe -= 6;
    ipf -= 5;
    xyp -= 4;

    /* Function Body */
    flagdtf = ifxnode_(status, &c__8);
/* ... recovering V-, VB- and I-points from TV-, TVB, and TB-points */
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	delxnode_(&icp[n], &c__128);
    }
/* ... delete references to material or fictitious faces */
    le = *ne + ihole[1];
    i__1 = le;
    for (n = 1; n <= i__1; ++n) {
	for (i__ = 1; i__ <= 4; ++i__) {
	    if__ = ife[i__ + (n << 2)];
	    if (if__ > 0) {
		if (ipf[(if__ << 2) + 4] == 900 && flagdtf || ipf[(if__ << 2) 
			+ 4] >= 901) {
		    ife[i__ + (n << 2)] = -1;
		}
	    }
	}
    }
/* ... delete all material faces */
    lf = *nf + iholf[1];
    i__1 = lf;
    for (n = 1; n <= i__1; ++n) {
	if (ipf[(n << 2) + 4] >= 901) {
	    facdel_(&n, nf, &ipf[5], &iholf[1]);
	}
    }
    lf = *nf + iholf[1];
    if (flagdtf) {
/*  ...  delete all fictitious faces */
	i__1 = lf;
	for (n = 1; n <= i__1; ++n) {
	    if (ipf[(n << 2) + 4] == 900) {
		facdel_(&n, nf, &ipf[5], &iholf[1]);
	    }
	}
    } else {
/*  ...  change the color of fictitious faces if it's possible */
	icfree = 900;
	for (ic = 1; ic <= 899; ++ic) {
	    i__1 = lf;
	    for (n = 1; n <= i__1; ++n) {
		if (ipf[(n << 2) + 4] == ic) {
		    goto L1;
		}
	    }
	    icfree = ic;
	    goto L2;
L1:
	    ;
	}
L2:
	i__1 = lf;
	for (n = 1; n <= i__1; ++n) {
	    if (ipf[(n << 2) + 4] == 900) {
		ipf[(n << 2) + 4] = icfree;
	    }
	}
    }
    *np += iholp[1];
    *nf += iholf[1];
    *ne += ihole[1];
/* ... fill-in holes in the list of mesh points */
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	ipw[n] = 0;
    }
    lp = iholp[1];
    i__1 = lp;
    for (n = 1; n <= i__1; ++n) {
	ip = iholp[n + 1];
	ipw[ip] = -1;
    }
    mp = 0;
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	if (ipw[n] == 0) {
	    ++mp;
	    ipw[n] = mp;
	    for (i__ = 1; i__ <= 3; ++i__) {
		xyp[i__ + mp * 3] = xyp[i__ + n * 3];
	    }
	    icp[mp] = icp[n];
	    ipp[mp] = ipp[n];
	    for (i__ = 1; i__ <= 6; ++i__) {
		hesp[i__ + mp * 6] = hesp[i__ + n * 6];
	    }
	}
    }
    i__1 = *np;
    for (n = mp + 1; n <= i__1; ++n) {
	icp[n] = 0;
	ipp[n] = 0;
    }
    i__1 = *npv;
    for (n = 1; n <= i__1; ++n) {
	ipv[n] = ipw[ipv[n]];
    }
    i__1 = *nf;
    for (n = 1; n <= i__1; ++n) {
	if (ipf[(n << 2) + 1] > 0) {
	    for (i__ = 1; i__ <= 3; ++i__) {
		ipf[i__ + (n << 2)] = ipw[ipf[i__ + (n << 2)]];
	    }
	}
    }
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	if (ipe[n * 5 + 1] > 0) {
	    for (i__ = 1; i__ <= 4; ++i__) {
		ipe[i__ + n * 5] = ipw[ipe[i__ + n * 5]];
	    }
	}
    }
    *np = mp;
/* ... fill-in holes in the list of mesh faces */
    lf = iholf[1];
    i__1 = lf;
    for (n = 1; n <= i__1; ++n) {
	if__ = iholf[n + 1];
	i__2 = if__ + 1;
	for (m = *nf; m >= i__2; --m) {
	    if (ipf[(m << 2) + 1] > 0) {
		kf = m;
		goto L20;
	    }
	}
	goto L200;
L20:
	for (i__ = 1; i__ <= 4; ++i__) {
	    ipf[i__ + (if__ << 2)] = ipf[i__ + (kf << 2)];
	}
	i__2 = *nfv;
	for (k = 1; k <= i__2; ++k) {
	    if (ifv[k] == kf) {
		ifv[k] = if__;
	    }
	}
/*  ...  auxiliary structures */
	i__2 = *ne;
	for (k = 1; k <= i__2; ++k) {
	    for (i__ = 1; i__ <= 4; ++i__) {
		if (ife[i__ + (k << 2)] == kf) {
		    ife[i__ + (k << 2)] = if__;
		}
	    }
	}
	ipf[(kf << 2) + 1] = 0;
L200:
	--(*nf);
    }
/* ... fill-in holes in the list of mesh elements */
    le = ihole[1];
    i__1 = le;
    for (n = 1; n <= i__1; ++n) {
	ie = ihole[n + 1];
	i__2 = ie + 1;
	for (m = *ne; m >= i__2; --m) {
	    if (ipe[m * 5 + 1] > 0) {
		ke = m;
		goto L30;
	    }
	}
	goto L300;
L30:
	for (i__ = 1; i__ <= 5; ++i__) {
	    ipe[i__ + ie * 5] = ipe[i__ + ke * 5];
	}
	qe[ie] = qe[ke];
	i__2 = *nev;
	for (k = 1; k <= i__2; ++k) {
	    if (iev[k] == ke) {
		iev[k] = ie;
	    }
	}
/*  ...  auxiliary structures */
	for (i__ = 1; i__ <= 4; ++i__) {
	    ife[i__ + (ie << 2)] = ife[i__ + (ke << 2)];
	    iee[i__ + (ie << 2)] = iee[i__ + (ke << 2)];
	    iet = iee[i__ + (ie << 2)];
	    if (iet > 0) {
		for (j = 1; j <= 4; ++j) {
		    if (iee[j + (iet << 2)] == ke) {
			iee[j + (iet << 2)] = ie;
		    }
		}
	    }
	}
	ipe[ke * 5 + 1] = 0;
L300:
	--(*ne);
    }
    return 0;
} /* updm_ */

/* ================================================================ */
/* @f2h@ */ logical cmpe_(integer *i1, integer *i2, integer *i3, integer *iep,
	 integer *nep, integer *ie1, integer *ie2)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    logical ret_val;

    /* Local variables */
    static integer i__, j, k, ib[3], ie[3], ip[3];

/* ================================================================ */
/* cmpE = TRUE if iE2 != iE1 and iE2 = {i1, i2, i3, *} */
/* ================================================================ */
/* group (Local variables) */
    /* Parameter adjustments */
    --nep;
    --iep;

    /* Function Body */
    ip[0] = *i1;
    ip[1] = *i2;
    ip[2] = *i3;
    for (i__ = 1; i__ <= 3; ++i__) {
	if (ip[i__ - 1] == 1) {
	    ib[i__ - 1] = 1;
	} else {
	    ib[i__ - 1] = nep[ip[i__ - 1] - 1] + 1;
	}
	ie[i__ - 1] = nep[ip[i__ - 1]];
    }
    i__1 = ie[0];
    for (i__ = ib[0]; i__ <= i__1; ++i__) {
	*ie2 = iep[i__];
	if (*ie2 == *ie1) {
	    goto L10;
	}
	i__2 = ie[1];
	for (j = ib[1]; j <= i__2; ++j) {
	    if (*ie2 == iep[j]) {
		i__3 = ie[2];
		for (k = ib[2]; k <= i__3; ++k) {
		    if (*ie2 == iep[k]) {
			ret_val = TRUE_;
			goto L1000;
		    }
		}
		goto L10;
	    }
	}
L10:
	;
    }
    ret_val = FALSE_;
L1000:
    return ret_val;
} /* cmpe_ */

/* ================================================================ */
/* @f2h@ */ logical cmpf_(integer *i1, integer *i2, integer *ifp, integer *
	nfp, integer *if1, integer *if2)
{
    /* System generated locals */
    integer i__1, i__2;
    logical ret_val;

    /* Local variables */
    static integer i__, j, ib[2], ie[2], ip[2];

/* ================================================================ */
/* group (Local variables) */
    /* Parameter adjustments */
    --nfp;
    --ifp;

    /* Function Body */
    ip[0] = *i1;
    ip[1] = *i2;
    for (i__ = 1; i__ <= 2; ++i__) {
	if (ip[i__ - 1] == 1) {
	    ib[i__ - 1] = 1;
	} else {
	    ib[i__ - 1] = nfp[ip[i__ - 1] - 1] + 1;
	}
	ie[i__ - 1] = nfp[ip[i__ - 1]];
    }
    i__1 = ie[0];
    for (i__ = ib[0]; i__ <= i__1; ++i__) {
	*if2 = ifp[i__];
	if (*if2 == *if1) {
	    goto L10;
	}
	i__2 = ie[1];
	for (j = ib[1]; j <= i__2; ++j) {
	    if (*if2 == ifp[j]) {
		ret_val = TRUE_;
		goto L1000;
	    }
	}
L10:
	;
    }
    ret_val = FALSE_;
L1000:
    return ret_val;
} /* cmpf_ */

/* ================================================================ */
/* @f2h@ */ logical cmpp_(integer *ip, integer *ifp, integer *nfp, integer *
	ipf)
{
    /* System generated locals */
    integer i__1, i__2;
    logical ret_val;

    /* Local variables */
    static integer i__, j, ib, ie, icf1, icf2;

/* ================================================================ */
/* cmpP = TRUE if point iP belongs to a common edge of two faces */
/* with different color. Otherwise cmpP = FALSE. */

/* Remark: the routine doesn't say anything about cross points. */
/* ================================================================ */
    /* Parameter adjustments */
    ipf -= 5;
    --nfp;
    --ifp;

    /* Function Body */
    if (*ip == 1) {
	ib = 1;
    } else {
	ib = nfp[*ip - 1] + 1;
    }
    ie = nfp[*ip];
    i__1 = ie;
    for (i__ = ib; i__ <= i__1; ++i__) {
	icf1 = ipf[(ifp[i__] << 2) + 4];
	i__2 = ie;
	for (j = i__ + 1; j <= i__2; ++j) {
	    icf2 = ipf[(ifp[j] << 2) + 4];
	    if (icf1 != icf2) {
		ret_val = TRUE_;
		goto L1000;
	    }
	}
    }
    ret_val = FALSE_;
L1000:
    return ret_val;
} /* cmpp_ */

/* ================================================================ */
/* @f2h@ */ logical cmpr_(integer *ip, integer *ifp, integer *nfp, doublereal 
	*xyp, integer *ipf)
{
    /* Initialized data */

    static integer iref[4] = { 1,2,3,1 };

    /* System generated locals */
    integer i__1, i__2;
    logical ret_val;

    /* Local variables */
    static doublereal c__;
    static integer i__, j, i1, i2, j1, j2, i3, j3, ib, ie, if__, jf, ip1, ip2,
	     jp1, jp2, ip3, jp3;

/* ================================================================ */
/* cmpR = TRUE if point iP belongs to a common edge of 2 flat faces */
/* and the angle between these faces is smaller than 120 degrees. */
/* ================================================================ */
/* ================================================================ */
    /* Parameter adjustments */
    ipf -= 5;
    xyp -= 4;
    --nfp;
    --ifp;

    /* Function Body */
/* ================================================================ */
    if (*ip == 1) {
	ib = 1;
    } else {
	ib = nfp[*ip - 1] + 1;
    }
    ie = nfp[*ip];
    i__1 = ie;
    for (i__ = ib; i__ <= i__1; ++i__) {
	if__ = ifp[i__];
	i__2 = ie;
	for (j = i__ + 1; j <= i__2; ++j) {
	    jf = ifp[j];
	    for (i1 = 1; i1 <= 3; ++i1) {
		i2 = iref[i1];
		ip1 = ipf[i1 + (if__ << 2)];
		ip2 = ipf[i2 + (if__ << 2)];
		for (j1 = 1; j1 <= 3; ++j1) {
		    j2 = iref[j1];
		    jp1 = ipf[j1 + (jf << 2)];
		    jp2 = ipf[j2 + (jf << 2)];
		    if (check22_(&ip1, &ip2, &jp1, &jp2)) {
			i3 = iref[i2];
			j3 = iref[j2];
			ip3 = ipf[i3 + (if__ << 2)];
			jp3 = ipf[j3 + (jf << 2)];
			c__ = angle2faces_(&xyp[ip1 * 3 + 1], &xyp[ip2 * 3 + 
				1], &xyp[ip3 * 3 + 1], &xyp[jp3 * 3 + 1]);
/*  ...               cos(120) = -0.5 */
			if (c__ > -.5) {
			    ret_val = TRUE_;
			    goto L1000;
			}
		    }
		}
	    }
/* L10: */
	}
/* L20: */
    }
    ret_val = FALSE_;
L1000:
    return ret_val;
} /* cmpr_ */

/* ================================================================ */
/* @f2h@ */ logical crosspoint_(doublereal *xyp, integer *nf, integer *ifp, 
	integer *ip, integer *ipf)
{
    /* System generated locals */
    integer i__1, i__2;
    logical ret_val;

    /* Local variables */
    static integer i__, m, n, lr, if1, if2;
    static doublereal ang;
    static integer ipt, irs[3], npt;

/* ================================================================ */
/* The number of edges with end point iP is evaluated. */
/* The array IFP is destroyed in out algorithm. */
/* ================================================================ */
/* ================================================================ */
    /* Parameter adjustments */
    ipf -= 5;
    --ifp;
    xyp -= 4;

    /* Function Body */
    ret_val = FALSE_;
/* ... analyze pairs of different faces */
    lr = 0;
    i__1 = *nf;
    for (n = 1; n <= i__1; ++n) {
	if1 = ifp[n];
	i__2 = *nf;
	for (m = n + 1; m <= i__2; ++m) {
	    if2 = ifp[m];
	    if (ipf[(if1 << 2) + 4] == ipf[(if2 << 2) + 4]) {
		goto L20;
	    }
	    for (i__ = 1; i__ <= 3; ++i__) {
		ipt = ipf[i__ + (if1 << 2)];
		if (ipt == *ip) {
		    goto L10;
		}
		if (check13_(&ipt, &ipf[(if2 << 2) + 1], &ipf[(if2 << 2) + 2],
			 &ipf[(if2 << 2) + 3])) {
		    findse_(&lr, irs, &ipt, &npt);
		    if (npt == 0) {
			++lr;
			if (lr >= 3) {
			    ret_val = TRUE_;
			    goto L1000;
			}
			irs[lr - 1] = ipt;
		    }
		}
L10:
		;
	    }
L20:
	    ;
	}
    }
    if (lr <= 1) {
	ret_val = FALSE_;
    } else {
	ang = angle2edges_(&xyp[*ip * 3 + 1], &xyp[irs[0] * 3 + 1], &xyp[irs[
		1] * 3 + 1]);
	if (ang > 0.) {
	    ret_val = TRUE_;
	}
    }
L1000:
    return ret_val;
} /* crosspoint_ */

/* ================================================================ */
/* @f2h@ */ integer countcolors_(integer *n, integer *ice)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static integer i__, k, m, ic;

/* ================================================================ */
/* Compute the amount of different colors in array ICE(N). */
/* The array is destroyed in our algorithm. The different colors */
/* are gathered in the first entrices of the array. */

/* Remark : only positive colors are counted */
/* ================================================================ */
    /* Parameter adjustments */
    --ice;

    /* Function Body */
    ret_val = 0;
    if (*n == 0) {
	return ret_val;
    }
    m = *n;
L1:
    ++ret_val;
    k = ret_val;
    ic = ice[k];
    i__1 = m;
    for (i__ = ret_val + 1; i__ <= i__1; ++i__) {
	if (ice[i__] != ic) {
	    ++k;
	    ice[k] = ice[i__];
	}
    }
    m = k;
    if (ret_val < m) {
	goto L1;
    }
/* ... delete negative colors */
    ret_val = 0;
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (ice[i__] > 0) {
	    ++ret_val;
	    ice[ret_val] = ice[i__];
	}
    }
    return ret_val;
} /* countcolors_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int scale2cube_(integer *np, doublereal *xyp, 
	logical *flag__)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, n;
    static doublereal scale, minxyp[3], maxxyp[3];

/* ================================================================ */
/* Routine scales the model to the square [0.1, 0.9]^2. We allow */
/* 10% freedom for curved edges. */
/* ================================================================ */
/* ================================================================ */
/* ================================================================ */
    /* Parameter adjustments */
    xyp -= 4;

    /* Function Body */
    if (*flag__) {
	for (i__ = 1; i__ <= 3; ++i__) {
	    minxyp[i__ - 1] = xyp[i__ + 3];
	    maxxyp[i__ - 1] = xyp[i__ + 3];
	}
	i__1 = *np;
	for (n = 2; n <= i__1; ++n) {
	    for (i__ = 1; i__ <= 3; ++i__) {
/* Computing MIN */
		d__1 = minxyp[i__ - 1], d__2 = xyp[i__ + n * 3];
		minxyp[i__ - 1] = min(d__1,d__2);
/* Computing MAX */
		d__1 = maxxyp[i__ - 1], d__2 = xyp[i__ + n * 3];
		maxxyp[i__ - 1] = max(d__1,d__2);
	    }
	}
/*  ...  add 5% for the bounding box */
/*        Do i = 1, 3 */
/*           size = (maxXYP(i) - minXYP(i)) / 20 */
/*           minXYP(i) = minXYP(i) - size */
/*           maxXYP(i) = maxXYP(i) + size */
/*        End do */
	for (i__ = 1; i__ <= 3; ++i__) {
	    anicrv_1.refxyp[i__ - 1] = minxyp[i__ - 1];
	    anicrv_1.scaxyp[i__ - 1] = 1. / (maxxyp[i__ - 1] - minxyp[i__ - 1]
		    );
	}
	scale = min(anicrv_1.scaxyp[0],anicrv_1.scaxyp[1]);
	scale = min(scale,anicrv_1.scaxyp[2]);
	for (i__ = 1; i__ <= 3; ++i__) {
	    anicrv_1.scaxyp[i__ - 1] = scale;
	}
	i__1 = *np;
	for (n = 1; n <= i__1; ++n) {
	    for (i__ = 1; i__ <= 3; ++i__) {
		xyp[i__ + n * 3] = (xyp[i__ + n * 3] - anicrv_1.refxyp[i__ - 
			1]) * anicrv_1.scaxyp[i__ - 1];
	    }
	}
    } else {
	for (i__ = 1; i__ <= 3; ++i__) {
	    anicrv_1.scaxyp[i__ - 1] = 1. / anicrv_1.scaxyp[i__ - 1];
	}
	i__1 = *np;
	for (n = 1; n <= i__1; ++n) {
	    for (i__ = 1; i__ <= 3; ++i__) {
		xyp[i__ + n * 3] = anicrv_1.refxyp[i__ - 1] + xyp[i__ + n * 3]
			 * anicrv_1.scaxyp[i__ - 1];
	    }
	}
    }
    return 0;
} /* scale2cube_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int scaleback_(doublereal *xypi, doublereal *
	xypo)
{
    static integer i__;

/* ================================================================ */
/*  Routine computes physical coordinates of point XYPi */
/* ================================================================ */
/* ================================================================ */
    /* Parameter adjustments */
    --xypo;
    --xypi;

    /* Function Body */
    for (i__ = 1; i__ <= 3; ++i__) {
	xypo[i__] = anicrv_1.refxyp[i__ - 1] + xypi[i__] / anicrv_1.scaxyp[
		i__ - 1];
    }
    return 0;
} /* scaleback_ */

/* ============================================================== */
/* @f2h@ */ /* Subroutine */ int randr_(doublereal *xy1, doublereal *xy2, 
	doublereal *xy3, doublereal *xy4, doublereal *rout, doublereal *rin)
{
    static doublereal c__[3], f[3];
    static integer i__, j;
    static doublereal v[12]	/* was [3][4] */, vol, sqr;

/* ============================================================== */
/* Computes curcumscribed and inscribed radii for the tetrahedron */
/* given by forth vertices. */
/* ============================================================== */
/* ============================================================== */
    /* Parameter adjustments */
    --xy4;
    --xy3;
    --xy2;
    --xy1;

    /* Function Body */
    for (i__ = 1; i__ <= 3; ++i__) {
	v[i__ * 3 - 3] = xy1[i__] - xy4[i__];
	v[i__ * 3 - 2] = xy2[i__] - xy4[i__];
	v[i__ * 3 - 1] = xy3[i__] - xy4[i__];
	v[i__ + 8] = 0.;
    }
    for (i__ = 1; i__ <= 3; ++i__) {
	f[i__ - 1] = 0.;
	for (j = 1; j <= 3; ++j) {
	    f[i__ - 1] += v[i__ + j * 3 - 4] * (xy1[j] + xy4[j]);
	}
	f[i__ - 1] /= 2;
    }
    vol = calvol_(v, &v[3], &v[6], &v[9]);
    c__[0] = calvol_(f, &v[3], &v[6], &v[9]) / vol;
    c__[1] = calvol_(v, f, &v[6], &v[9]) / vol;
    c__[2] = calvol_(v, &v[3], f, &v[9]) / vol;
    *rout = caledge_(&xy1[1], c__);
    sqr = calsqr_(&xy1[1], &xy2[1], &xy3[1]) + calsqr_(&xy2[1], &xy3[1], &xy4[
	    1]) + calsqr_(&xy3[1], &xy4[1], &xy1[1]) + calsqr_(&xy4[1], &xy1[
	    1], &xy2[1]);
    *rin = dabs(vol) * 6 / sqr;
    return 0;
} /* randr_ */

/* ============================================================== */
/* @f2h@ */ /* Subroutine */ int copymeshdata_(integer *np, integer *ne, 
	doublereal *xyp, doublereal *hesp, integer *ipe, doublereal *xypw, 
	doublereal *hespw, integer *ipew)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, n;

/* ============================================================== */
/* ============================================================== */
    /* Parameter adjustments */
    ipew -= 5;
    hespw -= 7;
    xypw -= 4;
    ipe -= 5;
    hesp -= 7;
    xyp -= 4;

    /* Function Body */
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	for (i__ = 1; i__ <= 3; ++i__) {
	    xypw[i__ + n * 3] = xyp[i__ + n * 3];
	}
    }
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	for (i__ = 1; i__ <= 6; ++i__) {
	    hespw[i__ + n * 6] = hesp[i__ + n * 6];
	}
    }
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	for (i__ = 1; i__ <= 4; ++i__) {
	    ipew[i__ + (n << 2)] = ipe[i__ + (n << 2)];
	}
    }
    return 0;
} /* copymeshdata_ */

/* ============================================================== */
/* @f2h@ */ doublereal surfacearea_(integer *nf, doublereal *xyp, integer *
	ipf, integer *ic)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer n;
    static doublereal s;
    static integer ip1, ip2, ip3;

/* ============================================================== */
/* The routine computes area of surface maked as ic. If ic <= 0, */
/* area of the total surface is computed. */
/* ============================================================== */
/* (Local variables) */
    /* Parameter adjustments */
    ipf -= 5;
    xyp -= 4;

    /* Function Body */
    s = 0.;
    i__1 = *nf;
    for (n = 1; n <= i__1; ++n) {
	if (ipf[(n << 2) + 1] <= 0) {
	    goto L10;
	}
	if (ipf[(n << 2) + 4] == *ic || *ic <= 0) {
	    ip1 = ipf[(n << 2) + 1];
	    ip2 = ipf[(n << 2) + 2];
	    ip3 = ipf[(n << 2) + 3];
	    s += calsqr_(&xyp[ip1 * 3 + 1], &xyp[ip2 * 3 + 1], &xyp[ip3 * 3 + 
		    1]);
	}
L10:
	;
    }
    ret_val = s;
    return ret_val;
} /* surfacearea_ */

/* ============================================================== */
/* @f2h@ */ doublereal fixedarea_(integer *nfv, doublereal *xyp, integer *ipf,
	 integer *ifv)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer n;
    static doublereal s;
    static integer if__, ip1, ip2, ip3;

/* ============================================================== */
/* The routine computes area of surface maked as ic. */
/* ============================================================== */
/* (Local variables) */
    /* Parameter adjustments */
    --ifv;
    ipf -= 5;
    xyp -= 4;

    /* Function Body */
    s = 0.;
    i__1 = *nfv;
    for (n = 1; n <= i__1; ++n) {
	if__ = ifv[n];
	ip1 = ipf[(if__ << 2) + 1];
	ip2 = ipf[(if__ << 2) + 2];
	ip3 = ipf[(if__ << 2) + 3];
	s += calsqr_(&xyp[ip1 * 3 + 1], &xyp[ip2 * 3 + 1], &xyp[ip3 * 3 + 1]);
/* L10: */
    }
    ret_val = s;
    return ret_val;
} /* fixedarea_ */

/* ============================================================== */
/* @f2h@ */ doublereal domainvolume_(integer *ne, doublereal *xyp, integer *
	ipe, integer *ic)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Local variables */
    static integer n;
    static doublereal s;
    static integer ip1, ip2, ip3, ip4;

/* ============================================================== */
/* The routine computes volume of subdomain maked as ic. */
/* If ic = 0, volume the whole domain is computed. */
/* ============================================================== */
/* (Local variables) */
    /* Parameter adjustments */
    ipe -= 6;
    xyp -= 4;

    /* Function Body */
    s = 0.;
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	if (ipe[n * 5 + 1] <= 0) {
	    goto L10;
	}
	if (ipe[n * 5 + 5] == *ic || *ic <= 0) {
	    ip1 = ipe[n * 5 + 1];
	    ip2 = ipe[n * 5 + 2];
	    ip3 = ipe[n * 5 + 3];
	    ip4 = ipe[n * 5 + 4];
	    s += (d__1 = calvol_(&xyp[ip1 * 3 + 1], &xyp[ip2 * 3 + 1], &xyp[
		    ip3 * 3 + 1], &xyp[ip4 * 3 + 1]), dabs(d__1));
	}
L10:
	;
    }
    ret_val = s;
    return ret_val;
} /* domainvolume_ */

/* ============================================================== */
/* @f2h@ */ doublereal fixedvolume_(integer *nev, doublereal *xyp, integer *
	ipe, integer *iev)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Local variables */
    static integer n;
    static doublereal s;
    static integer ie, ip1, ip2, ip3, ip4;

/* ============================================================== */
/* The routine computes volume of a fixed domain. */
/* ============================================================== */
/* (Local variables) */
    /* Parameter adjustments */
    --iev;
    ipe -= 6;
    xyp -= 4;

    /* Function Body */
    s = 0.;
    i__1 = *nev;
    for (n = 1; n <= i__1; ++n) {
	ie = iev[n];
	ip1 = ipe[ie * 5 + 1];
	ip2 = ipe[ie * 5 + 2];
	ip3 = ipe[ie * 5 + 3];
	ip4 = ipe[ie * 5 + 4];
	s += (d__1 = calvol_(&xyp[ip1 * 3 + 1], &xyp[ip2 * 3 + 1], &xyp[ip3 * 
		3 + 1], &xyp[ip4 * 3 + 1]), dabs(d__1));
/* L10: */
    }
    ret_val = s;
    return ret_val;
} /* fixedvolume_ */

