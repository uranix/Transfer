/* F2E.f -- translated by f2c (version 20090411).
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
static integer c__1 = 1;
static integer c__4 = 4;
static integer c__1007 = 1007;

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int f2e_(integer *iwf, integer *iwe, integer *ne,
	 integer *maxe, doublereal *xyp, integer *ipe, doublereal *hstar, 
	integer *icp, integer *iep, integer *ife, integer *iee, integer *l1e, 
	integer *l2e, integer *nl2, integer *nstep, integer *ihole, integer *
	status, doublereal *hesp, doublereal *rquality, doublereal *qe, 
	integer *lfu, integer *leu, integer *ifu, integer *ieu, integer *ipfu,
	 integer *ipeu, doublereal *qeu, logical *flag__)
{
    static logical flagshut;
    static doublereal v;
    static integer i1, i2, i3, i4, j1, j2, j3, j4, if__, le, lf, ip[5], ie1, 
	    ie2, ne1, ne2, ip1, ip2, ip3, ide, ipa, ipb, ies[1000], ifs[1000];
    static doublereal qes[1000];
    static integer mbad, nbad, ipes[5000]	/* was [5][1000] */, ipfs[
	    4000]	/* was [4][1000] */;
    static logical flag1, flag2, flag3, flagtm;

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
    hesp -= 7;
    --ihole;
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
    *flag__ = FALSE_;
    ip[0] = 1;
    ip[1] = 2;
    ip[2] = 3;
    ip[3] = 4;
    ip[4] = 1;
    if__ = ife[*iwf + (*iwe << 2)];
    i1 = *iwf;
    i2 = ip[i1];
    i3 = ip[i2];
    ip1 = ipe[i1 + *iwe * 5];
    ip2 = ipe[i2 + *iwe * 5];
    ip3 = ipe[i3 + *iwe * 5];
    ie1 = *iwe;
    ie2 = iee[*iwf + (ie1 << 2)];
/* ... checking the case when face -> edge operation is impossible */
    if (ie2 == 0) {
	goto L1000;
    }
    if (if__ != 0) {
	goto L1000;
    }
    i4 = ip[i3];
    ipa = ipe[i4 + *iwe * 5];
    for (j1 = 1; j1 <= 4; ++j1) {
	j2 = ip[j1];
	j3 = ip[j2];
	if (check33_(&i1, &i2, &i3, &j1, &j2, &j3)) {
	    j4 = ip[j3];
	    ipb = ipe[j4 + ie2 * 5];
	}
    }
/* ... checking the case when face -> edge operation is impossible */
    flagshut = TRUE_;
    shutf_(&xyp[ipa * 3 + 1], &xyp[ipb * 3 + 1], &xyp[ip1 * 3 + 1], &xyp[ip2 *
	     3 + 1], &xyp[ip3 * 3 + 1], &flagshut);
    if (! flagshut) {
	goto L1000;
    }
/* ... number of inverted elements */
    flagtm = ifxnode_(status, &c__512);
    if (flagtm) {
	nbad = 0;
	if (qes[ie1 - 1] <= 0.) {
	    ++nbad;
	}
	if (qes[ie2 - 1] <= 0.) {
	    ++nbad;
	}
    }
/* ... saving the initial superelement structure */
    copyse_(lfu, leu, &ifu[1], &ieu[1], &ipfu[5], &ipeu[6], &qeu[1], &lf, &le,
	     ifs, ies, ipfs, ipes, qes);
/* ... making a virtual evaluation of the quality */
    findse_(&le, ies, &ie1, &ne1);
    calqe_(&hesp[ipa * 6 + 1], &xyp[ipa * 3 + 1], &hesp[ipb * 6 + 1], &xyp[
	    ipb * 3 + 1], &hesp[ip1 * 6 + 1], &xyp[ip1 * 3 + 1], &hesp[ip3 * 
	    6 + 1], &xyp[ip3 * 3 + 1], hstar, &qes[ne1 - 1], &v);
    if (qes[ne1 - 1] <= *rquality) {
	goto L1000;
    }
    findse_(&le, ies, &ie2, &ne2);
    calqe_(&hesp[ipa * 6 + 1], &xyp[ipa * 3 + 1], &hesp[ipb * 6 + 1], &xyp[
	    ipb * 3 + 1], &hesp[ip2 * 6 + 1], &xyp[ip2 * 3 + 1], &hesp[ip3 * 
	    6 + 1], &xyp[ip3 * 3 + 1], hstar, &qes[ne2 - 1], &v);
    if (qes[ne2 - 1] <= *rquality) {
	goto L1000;
    }
    ++le;
    if (le > 1000) {
	goto L9000;
    }
    calqe_(&hesp[ipa * 6 + 1], &xyp[ipa * 3 + 1], &hesp[ipb * 6 + 1], &xyp[
	    ipb * 3 + 1], &hesp[ip1 * 6 + 1], &xyp[ip1 * 3 + 1], &hesp[ip2 * 
	    6 + 1], &xyp[ip2 * 3 + 1], hstar, &qes[le - 1], &v);
    if (qes[le - 1] <= *rquality) {
	goto L1000;
    }
/* ... checking for boundary elements */
    if (ifxnode_(status, &c__1) && ifxnode_(&icp[ipa], &c__4) && ifxnode_(&
	    icp[ipb], &c__4)) {
	flag1 = ifxnode_(&icp[ip1], &c__4);
	flag2 = ifxnode_(&icp[ip2], &c__4);
	if (flag1 && flag2) {
	    goto L1000;
	}
	flag3 = ifxnode_(&icp[ip3], &c__4);
	if (flag1 && flag3) {
	    goto L1000;
	}
	if (flag2 && flag3) {
	    goto L1000;
	}
    }
/* ... creating new elements */
    ide = ipes[ne1 * 5 - 1];
    ipes[ne1 * 5 - 5] = ipa;
    ipes[ne1 * 5 - 4] = ipb;
    ipes[ne1 * 5 - 3] = ip1;
    ipes[ne1 * 5 - 2] = ip3;
    ipes[ne1 * 5 - 1] = ide;
    ipes[ne2 * 5 - 5] = ipa;
    ipes[ne2 * 5 - 4] = ipb;
    ipes[ne2 * 5 - 3] = ip2;
    ipes[ne2 * 5 - 2] = ip3;
    ipes[ne2 * 5 - 1] = ide;
    ipes[le * 5 - 5] = ipa;
    ipes[le * 5 - 4] = ipb;
    ipes[le * 5 - 3] = ip1;
    ipes[le * 5 - 2] = ip2;
    ipes[le * 5 - 1] = ide;
/* ... checking for inverted elements */
    if (flagtm) {
	updqb_(&ne1, &le, ies, &xyp[4], ipes, qes);
	updqb_(&ne2, &le, ies, &xyp[4], ipes, qes);
	updqb_(&le, &le, ies, &xyp[4], ipes, qes);
	mbad = 0;
	if (qes[ne1 - 1] <= 0.) {
	    ++mbad;
	}
	if (qes[ne2 - 1] <= 0.) {
	    ++mbad;
	}
	if (qes[le - 1] <= 0.) {
	    ++mbad;
	}
	if (mbad >= nbad) {
	    goto L1000;
	}
    }
/* ... updating the grid */
    *flag__ = TRUE_;
    lstupd_(ne, &l1e[3], nl2, &l2e[1], &nstep[1], &qe[1], &ies[ne1 - 1], &qes[
	    ne1 - 1]);
    lstupd_(ne, &l1e[3], nl2, &l2e[1], &nstep[1], &qe[1], &ies[ne2 - 1], &qes[
	    ne2 - 1]);
    eledel_(&ies[ne1 - 1], &ipe[6], &iee[5]);
    eledel_(&ies[ne2 - 1], &ipe[6], &iee[5]);
    eleadd_(ne, maxe, &ihole[1]);
    lstadd_(ne, &l1e[3], nl2, &l2e[1], &nstep[1], &ihole[1], &qe[1], &qes[le 
	    - 1], &ies[le - 1]);
    eledel_(&ies[le - 1], &ipe[6], &iee[5]);
    eleupd_(&ne1, &iep[1], &ipe[6], &ife[5], &iee[5], &lf, &le, ifs, ies, 
	    ipfs, ipes);
    eleupd_(&ne2, &iep[1], &ipe[6], &ife[5], &iee[5], &lf, &le, ifs, ies, 
	    ipfs, ipes);
    eleupd_(&le, &iep[1], &ipe[6], &ife[5], &iee[5], &lf, &le, ifs, ies, ipfs,
	     ipes);
L1000:
    return 0;
L9000:
    errmes_(&c__1007, "F2E", "local parameter MaxS is small", (ftnlen)3, (
	    ftnlen)29);
    return 0;
} /* f2e_ */

