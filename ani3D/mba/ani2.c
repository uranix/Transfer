/* ani2.f -- translated by f2c (version 20090411).
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
#include "list.h"
#include "E2F.h"
#include "F2E.h"
#include "dsort.h"
#include "lintrp3d.h"
#include "insrtp.h"
#include "deletp.h"
#include "swapr.h"
#include "movep.h"
#include "output.h"
#include "clpsr.h"
#include "chkm.h"
#include "splite.h"

#include "time_mba.h"

/* Common Block Declarations */

struct {
    doublereal anixy0[3], anixya[3], anixyb[3];
} aniplane_;

#define aniplane_1 aniplane_

/* Table of constant values */

static integer c__4001 = 4001;
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__64 = 64;
static integer c__1000 = 1000;
static integer c__512 = 512;
static integer c__0 = 0;
static integer c__5201 = 5201;

/* ====================================================================== */
/* @f2h@ */ /* Subroutine */ int ani2_(integer *np, integer *maxp, integer *
	nf, integer *maxf, integer *ne, integer *maxe, doublereal *xyp, 
	integer *ipf, integer *ipe, integer *nestar, doublereal *hstar, 
	integer *icp, integer *ipp, integer *iep, integer *ife, integer *iee, 
	integer *l1e, integer *l2e, integer *iholp, integer *iholf, integer *
	ihole, integer *iepw, integer *nepw, integer *milintrp, integer *
	mrlintrp, integer *ipew, integer *ise, logical *flagauto, integer *
	status, integer *npv, integer *nfv, integer *nev, integer *ipv, 
	integer *ifv, integer *iev, integer *maxskipe, integer *maxqitr, 
	integer *maxbaskets, doublereal *hesp, doublereal *quality, 
	doublereal *rquality, doublereal *detg, doublereal *qe, integer *npw, 
	integer *new__, doublereal *xypw, doublereal *hespw, doublereal *rse, 
	I_fp metricfunction, logical *flaganalytic, logical *flagfile, 
	integer *chanelout, integer *nlines, char *output, integer *iprint, 
	integer *nqitr, integer *ierr, ftnlen output_len)
{
    /* Format strings */
    static char fmt_5009[] = "(\002Warning:\002,i6,\002 new faces have been "
	    "added\002)";
    static char fmt_5002[] = "(\002Domain:\002,i3,\002/\002,i2,\002,  clr"
	    "=\002,i4,\002,   vol=\002,e12.6)";
    static char fmt_5007[] = "(/,\002Fixed domain:\002,i6,\002  tets,   vol"
	    "=\002,e12.6)";
    static char fmt_5008[] = "(\002Fixed  area: \002,i6,\002 faces,  area"
	    "=\002,e12.6)";
    static char fmt_5003[] = "(\002Surface:\002,i3,\002/\002,i2,\002,  clr"
	    "=\002,i4,\002,  area=\002,e12.6)";
    static char fmt_5004[] = "(\002Maximal R/r =\002,e10.3,\002  (R/r = 3 fo"
	    "r equilateral tetrahedron),  status.fd:\002,i4)";
    static char fmt_5000[] = "(\002ITRs:\002,i7,\002 Q=\002,e10.4,\002  #V#F"
	    "#T:\002,i7,i8,i9,\002  tm=\002,f6.1,\002s\002)";
    static char fmt_5010[] = "(\002Avg Quality is\002,e11.4,\002,  Maximal R"
	    "/r =\002,e11.4,\002,  status.fd:\002,i5)";
    static char fmt_6001[] = "(i4,\002  elements in the \002,i3,\002-th bask"
	    "et,\002,i7,\002  is the current bad element\002)";
    static char fmt_6002[] = "(\002Face-Edge  Edge-Face  Gen.Swap  Insert  D"
	    "elete   Move  Collapse\002)";
    static char fmt_6003[] = "(i10,i12,i10,i8,i8,i7,i9)";
    static char fmt_5001[] = "(/,i4,\002  elements in the \002,i3,\002-th ba"
	    "sket,\002,i7,\002  is the current bad element\002,/,\002Face-Edg"
	    "e Edge-Face  Gen.Swap  Insert  Delete    Move  Collapse\002,/,2("
	    "i9,i10,i10,i8,i8,i8,i10/))";
    static char fmt_5006[] = "(\002Domain volume and total boundary area:"
	    " \002,2e14.6)";

    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2;
    doublereal d__1, d__2;

    /* Local variables */
    static real averageq;
    static integer nqitradd;
    static logical flagtest;
    static integer nqitrbig;
    static doublereal aquality;
    static integer icontrol, i__, n;
    static doublereal s;
    static integer ic, ke, le, me;
    static doublereal sa;
    static integer mf, lf;
    static doublereal rr;
    static integer ip1, ip2, nl2, ip3, ip4;
    static doublereal tm1, tm2, fao, aer;
    static integer ldh, ieu[1000], iwf[4], ifu[1000];
    static doublereal sao, qeu[1000], ver, dvo;
    static integer iwr[6];
    static doublereal fvo;
    static integer nfo, ipt, iwe, iop, ite, nxy, nf2e, ne2f, mf2e, me2f, ipeu[
	    5000]	/* was [5][1000] */, ipfu[4000]	/* was [4][1000] */;
    static doublereal xyps[3];
    static integer nfold, mclps, nclps;
    static doublereal hesps[6];
    static integer nmove, mmove, mswap;
    static doublereal rmove;
    static integer nstep[4], nswap;
    static logical flagtm;
    static integer ndelet, mdelet, nskipe, minsrt, ninsrt;
    static logical flagfbe, flaguar;
    static char message[100];

    /* Fortran I/O blocks */
    static icilist io___21 = { 0, message, 0, fmt_5009, 100, 1 };
    static cilist io___22 = { 0, 6, 0, fmt_5009, 0 };
    static icilist io___28 = { 0, message, 0, fmt_5002, 100, 1 };
    static cilist io___29 = { 0, 6, 0, fmt_5002, 0 };
    static icilist io___31 = { 0, message, 0, fmt_5007, 100, 1 };
    static icilist io___32 = { 0, message, 0, fmt_5008, 100, 1 };
    static cilist io___33 = { 0, 6, 0, fmt_5007, 0 };
    static cilist io___34 = { 0, 6, 0, fmt_5008, 0 };
    static icilist io___36 = { 0, message, 0, fmt_5003, 100, 1 };
    static cilist io___37 = { 0, 6, 0, fmt_5003, 0 };
    static icilist io___41 = { 0, message, 0, fmt_5004, 100, 1 };
    static cilist io___42 = { 0, 6, 0, fmt_5004, 0 };
    static icilist io___46 = { 0, message, 0, fmt_5000, 100, 1 };
    static cilist io___47 = { 0, 6, 0, fmt_5000, 0 };
    static icilist io___66 = { 0, message, 0, fmt_5010, 100, 1 };
    static icilist io___67 = { 0, message, 0, fmt_5000, 100, 1 };
    static cilist io___68 = { 0, 6, 0, fmt_5010, 0 };
    static cilist io___69 = { 0, 6, 0, fmt_5000, 0 };
    static icilist io___74 = { 0, message, 0, fmt_5000, 100, 1 };
    static cilist io___75 = { 0, 6, 0, fmt_5000, 0 };
    static icilist io___76 = { 0, message, 0, fmt_6001, 100, 1 };
    static icilist io___77 = { 0, message, 0, fmt_6002, 100, 1 };
    static icilist io___78 = { 0, message, 0, fmt_6003, 100, 1 };
    static icilist io___79 = { 0, message, 0, fmt_6003, 100, 1 };
    static cilist io___80 = { 0, 6, 0, fmt_5001, 0 };
    static icilist io___90 = { 0, message, 0, fmt_5000, 100, 1 };
    static cilist io___91 = { 0, 6, 0, fmt_5000, 0 };
    static icilist io___92 = { 0, message, 0, fmt_6001, 100, 1 };
    static icilist io___93 = { 0, message, 0, fmt_6002, 100, 1 };
    static icilist io___94 = { 0, message, 0, fmt_6003, 100, 1 };
    static icilist io___95 = { 0, message, 0, fmt_6003, 100, 1 };
    static cilist io___96 = { 0, 6, 0, fmt_5001, 0 };
    static icilist io___99 = { 0, message, 0, fmt_5010, 100, 1 };
    static cilist io___100 = { 0, 6, 0, fmt_5010, 0 };
    static icilist io___101 = { 0, message, 0, fmt_5002, 100, 1 };
    static cilist io___102 = { 0, 6, 0, fmt_5002, 0 };
    static icilist io___103 = { 0, message, 0, fmt_5006, 100, 1 };
    static cilist io___104 = { 0, 6, 0, fmt_5006, 0 };
    static icilist io___107 = { 0, message, 0, fmt_5004, 100, 1 };
    static cilist io___108 = { 0, 6, 0, fmt_5004, 0 };


/* ====================================================================== */
/* group (M) */
/* group (Dev) */
/* group (Q) */
/* group (ERR) */
/* ====================================================================== */
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
/* ========================================================== */
/* MaxLines   - the maximal number of output lines per processor */
/* LineLenght - the maximal lenght of a message */
/* ========================================================== */
/* ====================================================================== */
/* The core of the package where the algorithm structure is realized */

/* We mention here only essential parameters that differ from ones */
/* in MeshMetric and MeshSolution. */

/*     IPE(5, MaxE) - columns 1,2,3,4 are the connectivity */
/*                            list of elements; */
/*                    column  5 is the element identificator */

/*     IPF(4, MaxF) - columns 1,2,3 are the connectivity */
/*                            list of faces; */
/*                    column  4 is the face identificator */

/*     IPEw(4, nE)  - container for keeping original */
/*                    connectivity list for elements */
/* ====================================================================== */
/* group (M) */
/*     Integer MaxP, MaxF, MaxE, nPv, nEStar */
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
/* group (Dev) */
/* group (Q) */
/*     Integer MaxSkipE, MaxQItr, MaxBaskets */
/* group (output) */
/* ====================================================================== */
/* group (Local functions) */
/* group (Local variables) */
/* ====================================================================== */
/* ====================================================================== */
    /* Parameter adjustments */
    output -= 100;
    --rse;
    hespw -= 7;
    xypw -= 4;
    --qe;
    --detg;
    hesp -= 7;
    --iev;
    --ifv;
    --ipv;
    --ise;
    ipew -= 5;
    --nepw;
    --iepw;
    --ihole;
    --iholf;
    --iholp;
    --l2e;
    l1e -= 3;
    iee -= 5;
    ife -= 5;
    --iep;
    --ipp;
    --icp;
    ipe -= 6;
    ipf -= 5;
    xyp -= 4;

    /* Function Body */
    seconds_(&tm1);
    nf2e = 0;
    ne2f = 0;
    nswap = 0;
    ninsrt = 0;
    ndelet = 0;
    nmove = 0;
    nclps = 0;
    mf2e = 0;
    me2f = 0;
    mswap = 0;
    minsrt = 0;
    mdelet = 0;
    mmove = 0;
    mclps = 0;
    flagtm = TRUE_;
/* ... Loop initialization */
    nqitradd = 0;
    nqitrbig = 0;
    *nqitr = 0;
    *ierr = 0;
    nfo = *nf;
    makm_(np, nf, maxf, ne, maxe, &xyp[4], &ipf[5], &ipe[6], &icp[1], &ipp[1],
	     &iep[1], &ife[5], &iee[5], &iholp[1], &iholf[1], &ihole[1], &
	    iepw[1], &nepw[1], status, npv, nfv, nev, &ipv[1], &ifv[1], &iev[
	    1], ierr);
/* group (M) */
/* group (Dev) */
/* group (iERR) */
    if (*ierr != 0) {
	goto L9000;
    }
    if (nfo != *nf) {
	if (! (*flagauto)) {
	    errmes_(&c__4001, "ani2", "inconsistent input data", (ftnlen)4, (
		    ftnlen)23);
	}
	if (*iprint >= 1) {
	    if (*flagfile) {
		s_wsfi(&io___21);
		i__1 = *nf - nfo;
		do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
		e_wsfi();
		addout_(nlines, message, output + 100, (ftnlen)100, (ftnlen)
			100);
	    } else {
		s_wsfe(&io___22);
		i__1 = *nf - nfo;
		do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
		e_wsfe();
	    }
	}
    }
/* ... gather statistics about elements. We assume that IHolE(1) = 0 */
    if (*iprint >= 3) {
	me = 0;
	i__1 = *ne;
	for (n = 1; n <= i__1; ++n) {
	    if (ipe[n * 5 + 1] > 0) {
		++me;
		iepw[me] = ipe[n * 5 + 5];
	    }
	}
	ic = countcolors_(&me, &iepw[1]);
	i__1 = ic;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s = domainvolume_(ne, &xyp[4], &ipe[6], &iepw[i__]);
	    if (*flagfile) {
		s_wsfi(&io___28);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&iepw[i__], (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&s, (ftnlen)sizeof(doublereal));
		e_wsfi();
		addout_(nlines, message, output + 100, (ftnlen)100, (ftnlen)
			100);
	    } else {
		s_wsfe(&io___29);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&iepw[i__], (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&s, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	}
	s = fixedvolume_(nev, &xyp[4], &ipe[6], &iev[1]);
	sa = fixedarea_(nfv, &xyp[4], &ipf[5], &ifv[1]);
	if (*flagfile) {
	    s_wsfi(&io___31);
	    do_fio(&c__1, (char *)&(*nev), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&s, (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    addout_(nlines, message, output + 100, (ftnlen)100, (ftnlen)100);
	    s_wsfi(&io___32);
	    do_fio(&c__1, (char *)&(*nev), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&sa, (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    addout_(nlines, message, output + 100, (ftnlen)100, (ftnlen)100);
	} else {
	    s_wsfe(&io___33);
	    do_fio(&c__1, (char *)&(*nev), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&s, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    s_wsfe(&io___34);
	    do_fio(&c__1, (char *)&(*nev), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&sa, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }
/* ... gather statistics about faces. We allow IHolF(1) >= 0 */
    if (*iprint >= 4) {
	mf = 0;
	i__1 = *nf + iholf[1];
	for (n = 1; n <= i__1; ++n) {
	    if (ipf[(n << 2) + 4] > 0) {
		++mf;
		iepw[mf] = ipf[(n << 2) + 4];
	    }
	}
	ic = countcolors_(&mf, &iepw[1]);
	i__1 = ic;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *nf + iholf[1];
	    s = surfacearea_(&i__2, &xyp[4], &ipf[5], &iepw[i__]);
	    if (*flagfile) {
		s_wsfi(&io___36);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&iepw[i__], (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&s, (ftnlen)sizeof(doublereal));
		e_wsfi();
		addout_(nlines, message, output + 100, (ftnlen)100, (ftnlen)
			100);
	    } else {
		s_wsfe(&io___37);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&iepw[i__], (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&s, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	}
    }
/* ... check the mesh (level 1) */
    chkm_(np, nf, ne, &xyp[4], &ipf[5], &ipe[6], &icp[1], &ife[5], &iee[5], &
	    rr, status, &c__1, &nepw[1], &iepw[1]);
/* group (M) */
/* group (W) */
/* ... check the mesh (level 2) */
    fvo = fixedvolume_(nev, &xyp[4], &ipe[6], &iev[1]);
    fao = fixedarea_(nfv, &xyp[4], &ipf[5], &ifv[1]);
/* ... output the statistics */
    if (*iprint >= 2) {
	if (*flagfile) {
	    s_wsfi(&io___41);
	    do_fio(&c__1, (char *)&rr, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*status), (ftnlen)sizeof(integer));
	    e_wsfi();
	    addout_(nlines, message, output + 100, (ftnlen)100, (ftnlen)100);
	} else {
	    s_wsfe(&io___42);
	    do_fio(&c__1, (char *)&rr, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*status), (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }
/* ... modify the grid to satisfy some restrictios (usually with FEs) */
    flaguar = ifxnode_(status, &c__2);
    flagfbe = ifxnode_(status, &c__1);
    if (flaguar || flagfbe) {
	seconds_(&tm2);
	if (*iprint >= 2) {
	    *rquality = 1.;
	    i__1 = *ne;
	    for (n = 1; n <= i__1; ++n) {
/* Computing MIN */
		d__1 = *rquality, d__2 = qe[n];
		*rquality = min(d__1,d__2);
	    }
	    if (*flagfile) {
		s_wsfi(&io___46);
		do_fio(&c__1, (char *)&(*nqitr), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*rquality), (ftnlen)sizeof(doublereal)
			);
		do_fio(&c__1, (char *)&(*np), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*nf), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*ne), (ftnlen)sizeof(integer));
		d__1 = tm2 - tm1;
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		e_wsfi();
		addout_(nlines, message, output + 100, (ftnlen)100, (ftnlen)
			100);
	    } else {
		s_wsfe(&io___47);
		do_fio(&c__1, (char *)&(*nqitr), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*rquality), (ftnlen)sizeof(doublereal)
			);
		do_fio(&c__1, (char *)&(*np), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*nf), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*ne), (ftnlen)sizeof(integer));
		d__1 = tm2 - tm1;
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	}
	ke = *ne;
	i__1 = ke;
	for (n = 1; n <= i__1; ++n) {
	    for (i__ = 1; i__ <= 4; ++i__) {
		ipt = ipe[i__ + n * 5];
		if (ifxnode_(&icp[ipt], &c__64)) {
		    goto L20;
		}
	    }
	    ++(*nqitr);
	    makse_(&n, &icp[1], &iep[1], &ipf[5], &ipe[6], &ife[5], &iee[5], &
		    qe[1], &c__1000, &lf, &le, ifu, ieu, ipfu, ipeu, qeu);
	    splite_(&n, np, maxp, ne, maxe, &xyp[4], &ipe[6], hstar, &icp[1], 
		    &iep[1], &ife[5], &iee[5], &iholp[1], &ihole[1], &hesp[7],
		     &detg[1], &qe[1], &lf, &le, ifu, ieu, ipfu, ipeu, qeu, &
		    flagtest);
/* group (M) */
/* group (Q) */
/* group (S) */
L20:
	    ;
	}
    }
/* ... initialize the list */
    if (ifxnode_(status, &c__512)) {
	i__1 = *ne;
	for (n = 1; n <= i__1; ++n) {
	    updqa_(&n, &xyp[4], &ipe[6], &iee[5], &qe[1]);
	}
    }
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	rse[n] = qe[n];
	l2e[n] = n;
    }
    dsort_(&rse[1], &l2e[1], ne, &c__2, ierr);
    if (*ierr != 0) {
	goto L9000;
    }
/* Computing MAX */
    r__1 = 100.f, r__2 = (real)log((real) (*nestar));
    nstep[0] = (int)dmax(r__1,r__2);
    nstep[1] = nstep[0] / 4;
    nstep[2] = *maxe;
    nstep[3] = 0;
    lstmak_(nestar, ne, &l1e[3], &l2e[1], &nl2, nstep, &ihole[1]);
/* ... initialize the 8-tree */
    ldh = 6;
    nxy = 0;
    icontrol = *chanelout + 1100;
    if (! (*flaganalytic)) {
	lintrp3d_(new__, &ipew[5], npw, &xypw[4], &ldh, &hespw[7], &nxy, xyps,
		 hesps, &ise[1], milintrp, &rse[1], mrlintrp, &icontrol);
    }
/* ... switch-off initialization */
    icontrol += 1000;
/* ... output of initial qualities */
    if (*iprint >= 2) {
	countbadelements_(ne, &l1e[3], &l2e[1], &qe[1], quality, nlines, 
		output + 100, flagfile, (ftnlen)100);
    }
/* ... empty the basket and repeat the main loop */
L100:
    ++nqitrbig;
    if (nqitrbig > *maxbaskets) {
	*ierr = 1000;
	goto L1000;
    }
    seconds_(&tm2);
    if (*iprint >= 3 || *iprint >= 1 && nqitrbig == 1) {
	aquality = avgq_(ne, &qe[1], &l1e[3], &l2e[1]);
	if (*flagfile) {
	    s_wsfi(&io___66);
	    do_fio(&c__1, (char *)&aquality, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&rr, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*status), (ftnlen)sizeof(integer));
	    e_wsfi();
	    addout_(nlines, message, output + 100, (ftnlen)100, (ftnlen)100);
	    s_wsfi(&io___67);
	    do_fio(&c__1, (char *)&(*nqitr), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&qe[l2e[1]], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*np), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*nf), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*ne), (ftnlen)sizeof(integer));
	    d__1 = tm2 - tm1;
	    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    addout_(nlines, message, output + 100, (ftnlen)100, (ftnlen)100);
	} else {
	    s_wsfe(&io___68);
	    do_fio(&c__1, (char *)&aquality, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&rr, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*status), (ftnlen)sizeof(integer));
	    e_wsfe();
	    s_wsfe(&io___69);
	    do_fio(&c__1, (char *)&(*nqitr), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&qe[l2e[1]], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*np), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*nf), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*ne), (ftnlen)sizeof(integer));
	    d__1 = tm2 - tm1;
	    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }
    if (*maxqitr <= 0) {
	goto L1100;
    }
    nskipe = 0;
L300:
    ++(*nqitr);
    if (*nqitr > *maxqitr + nqitradd) {
	if (nqitrbig == 1) {
	    nqitradd = nskipe;
	    goto L100;
	}
	*ierr = 1000;
	goto L1000;
    }
/* ... check the mesh (level 3) */
    if (qe[l2e[1]] > 0. && flagtm) {
	flagtm = FALSE_;
	delxnode_(status, &c__512);
	i__1 = *ne + ihole[1];
	dvo = domainvolume_(&i__1, &xyp[4], &ipe[6], &c__0);
	i__1 = *nf + iholf[1];
	sao = surfacearea_(&i__1, &xyp[4], &ipf[5], &c__0);
    }
/* ... check for termination signals (dummy function for simulators) */
    if (*nqitr / 200 * 200 == *nqitr && probeany_()) {
	if (nqitrbig == 1) {
	    nqitradd = nskipe;
	    goto L100;
	}
	*ierr = 1000;
	goto L1000;
    }
    iwe = l2e[1];
    i__1 = nskipe;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iwe = l1e[(iwe << 1) + 2];
    }
    if (*nqitr / 1000 * 1000 == *nqitr) {
	seconds_(&tm2);
	if (*iprint >= 2) {
	    if (*flagfile) {
		s_wsfi(&io___74);
		do_fio(&c__1, (char *)&(*nqitr), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&qe[l2e[1]], (ftnlen)sizeof(doublereal))
			;
		do_fio(&c__1, (char *)&(*np), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*nf), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*ne), (ftnlen)sizeof(integer));
		d__1 = tm2 - tm1;
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		e_wsfi();
		addout_(nlines, message, output + 100, (ftnlen)100, (ftnlen)
			100);
	    } else {
		s_wsfe(&io___75);
		do_fio(&c__1, (char *)&(*nqitr), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&qe[l2e[1]], (ftnlen)sizeof(doublereal))
			;
		do_fio(&c__1, (char *)&(*np), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*nf), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*ne), (ftnlen)sizeof(integer));
		d__1 = tm2 - tm1;
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	}
    }
    if (*iprint >= 3 && *nqitr / 10000 * 10000 == *nqitr) {
	if (*flagfile) {
	    s_wsfi(&io___76);
	    do_fio(&c__1, (char *)&nskipe, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nqitrbig, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&iwe, (ftnlen)sizeof(integer));
	    e_wsfi();
	    addout_(nlines, message, output + 100, (ftnlen)100, (ftnlen)100);
	    s_wsfi(&io___77);
	    e_wsfi();
	    addout_(nlines, message, output + 100, (ftnlen)100, (ftnlen)100);
	    s_wsfi(&io___78);
	    do_fio(&c__1, (char *)&nf2e, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ne2f, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nswap, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ninsrt, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ndelet, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nmove, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nclps, (ftnlen)sizeof(integer));
	    e_wsfi();
	    addout_(nlines, message, output + 100, (ftnlen)100, (ftnlen)100);
	    s_wsfi(&io___79);
	    do_fio(&c__1, (char *)&mf2e, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&me2f, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&mswap, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&minsrt, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&mdelet, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&mmove, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&mclps, (ftnlen)sizeof(integer));
	    e_wsfi();
	    addout_(nlines, message, output + 100, (ftnlen)100, (ftnlen)100);
	} else {
	    s_wsfe(&io___80);
	    do_fio(&c__1, (char *)&nskipe, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nqitrbig, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&iwe, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nf2e, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ne2f, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nswap, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ninsrt, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ndelet, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nmove, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nclps, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&mf2e, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&me2f, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&mswap, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&minsrt, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&mdelet, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&mmove, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&mclps, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }
    *rquality = qe[iwe];
    if (*rquality > *quality) {
	if (nskipe == 0) {
	    goto L1000;
	}
	goto L100;
    }
    *rquality *= 1.001;
    ip1 = ipe[iwe * 5 + 1];
    ip2 = ipe[iwe * 5 + 2];
    ip3 = ipe[iwe * 5 + 3];
    ip4 = ipe[iwe * 5 + 4];
    calqf_(&hesp[ip1 * 6 + 1], &xyp[ip1 * 3 + 1], &hesp[ip2 * 6 + 1], &xyp[
	    ip2 * 3 + 1], &hesp[ip3 * 6 + 1], &xyp[ip3 * 3 + 1], &hesp[ip4 * 
	    6 + 1], &xyp[ip4 * 3 + 1], hstar, iwf, iwr);
    makse_(&iwe, &icp[1], &iep[1], &ipf[5], &ipe[6], &ife[5], &iee[5], &qe[1],
	     &c__1000, &lf, &le, ifu, ieu, ipfu, ipeu, qeu);
    for (iop = 1; iop <= 7; ++iop) {
/* ... analyze edges of the element */
	if (1 == iop) {
	    for (i__ = 1; i__ <= 6; ++i__) {
		insrtp_(&iwr[i__ - 1], &iwe, np, maxp, nf, maxf, ne, maxe, &
			xyp[4], &ipf[5], &ipe[6], hstar, &icp[1], &iep[1], &
			ife[5], &iee[5], &l1e[3], &l2e[1], &nl2, nstep, &
			iholp[1], &iholf[1], &ihole[1], status, &hesp[7], 
			rquality, &detg[1], &qe[1], (I_fp)metricfunction, 
			flaganalytic, &lf, &le, ifu, ieu, ipfu, ipeu, qeu, 
			npw, new__, &xypw[4], &hespw[7], &ipew[5], milintrp, 
			mrlintrp, &ise[1], &rse[1], &icontrol, &flagtest);
/* group (M) */
/* group (Q) */
/* group (S) */
/* group (W) */
		++minsrt;
		if (flagtest) {
		    ++ninsrt;
		    goto L400;
		}
	    }
/* ... analyze faces of the element */
	} else if (2 == iop) {
	    for (i__ = 1; i__ <= 4; ++i__) {
		f2e_(&iwf[i__ - 1], &iwe, ne, maxe, &xyp[4], &ipe[6], hstar, &
			icp[1], &iep[1], &ife[5], &iee[5], &l1e[3], &l2e[1], &
			nl2, nstep, &ihole[1], status, &hesp[7], rquality, &
			qe[1], &lf, &le, ifu, ieu, ipfu, ipeu, qeu, &flagtest)
			;
/* group (M) */
/* group (Q) */
/* group (S) */
/* group (W) */
		++mf2e;
		if (flagtest) {
		    ++nf2e;
		    goto L400;
		}
	    }
/* ... analyze edges of the element */
	} else if (3 == iop) {
	    for (i__ = 1; i__ <= 6; ++i__) {
		e2f_(&iwr[i__ - 1], &iwe, nf, maxf, ne, &xyp[4], &ipf[5], &
			ipe[6], hstar, &icp[1], &iep[1], &ife[5], &iee[5], &
			l1e[3], &l2e[1], &nl2, nstep, &iholf[1], &ihole[1], 
			status, &hesp[7], rquality, &detg[1], &qe[1], &lf, &
			le, ifu, ieu, ipfu, ipeu, qeu, &flagtest);
/* group (M) */
/* group (Q) */
/* group (S) */
/* group (W) */
		++me2f;
		if (flagtest) {
		    ++ne2f;
		    goto L400;
		}
	    }
/* ... analyze edges of the element */
	} else if (4 == iop) {
	    for (i__ = 1; i__ <= 6; ++i__) {
		clpsr_(&i__, &iwe, np, ne, &xyp[4], &ipe[6], hstar, &icp[1], &
			iep[1], &ife[5], &iee[5], &l1e[3], &l2e[1], &nl2, 
			nstep, &iholp[1], &ihole[1], status, &hesp[7], 
			rquality, &detg[1], &qe[1], (I_fp)metricfunction, 
			flaganalytic, &lf, &le, ifu, ieu, ipfu, ipeu, qeu, 
			npw, new__, &xypw[4], &hespw[7], &ipew[5], milintrp, 
			mrlintrp, &ise[1], &rse[1], &icontrol, &flagtest);
/* group (M) */
/* group (Q) */
/* group (S) */
/* group (W) */
		++mclps;
		if (flagtest) {
		    ++nclps;
		    goto L400;
		}
	    }
/* ... swap edge (general method) */
	} else if (5 == iop) {
	    for (i__ = 1; i__ <= 6; ++i__) {
		swapr_(&iwr[i__ - 1], &iwe, nf, maxf, ne, maxe, &xyp[4], &ipf[
			5], &ipe[6], hstar, &icp[1], &iep[1], &ife[5], &iee[5]
			, &l1e[3], &l2e[1], &nl2, nstep, &iholf[1], &ihole[1],
			 status, &hesp[7], rquality, &qe[1], &lf, &le, ifu, 
			ieu, ipfu, ipeu, qeu, &flagtest);
/* group (M) */
/* group (Q) */
/* group (S) */
/* group (W) */
		++mswap;
		if (flagtest) {
		    ++nswap;
		    goto L400;
		}
	    }
/* ... analyze points of the element */
	} else if (7 == iop) {
	    for (i__ = 1; i__ <= 4; ++i__) {
		deletp_(&i__, &iwe, np, nf, maxf, ne, maxe, &xyp[4], &ipf[5], 
			&ipe[6], hstar, &icp[1], &iep[1], &ife[5], &iee[5], &
			l1e[3], &l2e[1], &nl2, nstep, &iholp[1], &iholf[1], &
			ihole[1], status, &hesp[7], rquality, &detg[1], &qe[1]
			, &lf, &le, ifu, ieu, ipfu, ipeu, qeu, &flagtest);
/* group (M) */
/* group (Q) */
/* group (S) */
/* group (W) */
		++mdelet;
		if (flagtest) {
		    ++ndelet;
		    goto L400;
		}
	    }
/* ... moving point */
	} else if (6 == iop) {
	    for (i__ = 1; i__ <= 4; ++i__) {
		movep_(&i__, &iwe, ne, &xyp[4], &ipf[5], &ipe[6], hstar, &icp[
			1], &ife[5], &l1e[3], &l2e[1], &nl2, nstep, status, &
			hesp[7], rquality, &detg[1], &qe[1], (I_fp)
			metricfunction, flaganalytic, &lf, &le, ifu, ieu, 
			ipfu, ipeu, qeu, npw, new__, &xypw[4], &hespw[7], &
			ipew[5], milintrp, mrlintrp, &ise[1], &rse[1], &
			icontrol, &rmove, &flagtest);
/* group (M) */
/* group (Q) */
/* group (S) */
/* group (W) */
		++mmove;
		if (flagtest) {
		    ++nmove;
		    goto L400;
		}
	    }
	}
    }
/* ... does an operation increase the quality sufficiently enough? */
L400:
    ite = l2e[1];
    i__1 = nskipe;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ite = l1e[(ite << 1) + 2];
    }
    if (ite == iwe) {
	flagtest = FALSE_;
    }
    if (flagtest) {
	goto L300;
    }
    ++nskipe;
    if (nskipe > *maxskipe || nskipe >= *ne) {
	goto L100;
    }
    goto L300;
L1000:
    *rquality = qe[l2e[1]];
    seconds_(&tm2);
    if (*iprint >= 1) {
	if (*flagfile) {
	    s_wsfi(&io___90);
	    do_fio(&c__1, (char *)&(*nqitr), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&qe[l2e[1]], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*np), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*nf), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*ne), (ftnlen)sizeof(integer));
	    d__1 = tm2 - tm1;
	    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	    e_wsfi();
	    addout_(nlines, message, output + 100, (ftnlen)100, (ftnlen)100);
	} else {
	    s_wsfe(&io___91);
	    do_fio(&c__1, (char *)&(*nqitr), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&qe[l2e[1]], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*np), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*nf), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*ne), (ftnlen)sizeof(integer));
	    d__1 = tm2 - tm1;
	    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }
/* ... calculate the number of bad tetrahedrons */
L1100:
    if (*iprint >= 3) {
	if (*flagfile) {
	    s_wsfi(&io___92);
	    do_fio(&c__1, (char *)&nskipe, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nqitrbig, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&iwe, (ftnlen)sizeof(integer));
	    e_wsfi();
	    addout_(nlines, message, output + 100, (ftnlen)100, (ftnlen)100);
	    s_wsfi(&io___93);
	    e_wsfi();
	    addout_(nlines, message, output + 100, (ftnlen)100, (ftnlen)100);
	    s_wsfi(&io___94);
	    do_fio(&c__1, (char *)&nf2e, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ne2f, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nswap, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ninsrt, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ndelet, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nmove, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nclps, (ftnlen)sizeof(integer));
	    e_wsfi();
	    addout_(nlines, message, output + 100, (ftnlen)100, (ftnlen)100);
	    s_wsfi(&io___95);
	    do_fio(&c__1, (char *)&mf2e, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&me2f, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&mswap, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&minsrt, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&mdelet, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&mmove, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&mclps, (ftnlen)sizeof(integer));
	    e_wsfi();
	    addout_(nlines, message, output + 100, (ftnlen)100, (ftnlen)100);
	} else {
	    s_wsfe(&io___96);
	    do_fio(&c__1, (char *)&nskipe, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nqitrbig, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&iwe, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nf2e, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ne2f, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nswap, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ninsrt, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ndelet, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nmove, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&nclps, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&mf2e, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&me2f, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&mswap, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&minsrt, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&mdelet, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&mmove, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&mclps, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }
    if (*iprint >= 2) {
	countbadelements_(ne, &l1e[3], &l2e[1], &qe[1], quality, nlines, 
		output + 100, flagfile, (ftnlen)100);
    }
/* ... update the mesh */
    nfold = *nf;
    updm_(np, nf, ne, &xyp[4], &ipf[5], &ipe[6], &icp[1], &ipp[1], &ife[5], &
	    iee[5], &iholp[1], &iholf[1], &ihole[1], status, npv, nfv, nev, &
	    ipv[1], &ifv[1], &iev[1], &hesp[7], &qe[1], &iepw[1]);
/* group (M) */
/* group (Dev) */
/* group (Q) */
/* group (W) */
/* ... check the mesh (level 1) */
    chkm_(np, nf, ne, &xyp[4], &ipf[5], &ipe[6], &icp[1], &ife[5], &iee[5], &
	    rr, status, &c__0, &nepw[1], &iepw[1]);
/* group (M) */
/* group (W) */
/* ... print out details */
    if (*iprint >= 1) {
	averageq = (real)avgq_(ne, &qe[1], &l1e[3], &l2e[1]);
	if (*flagfile) {
	    s_wsfi(&io___99);
	    do_fio(&c__1, (char *)&averageq, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&rr, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*status), (ftnlen)sizeof(integer));
	    e_wsfi();
	    addout_(nlines, message, output + 100, (ftnlen)100, (ftnlen)100);
	} else {
	    s_wsfe(&io___100);
	    do_fio(&c__1, (char *)&averageq, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&rr, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*status), (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }
    if (*iprint >= 3) {
	me = 0;
	i__1 = *ne;
	for (n = 1; n <= i__1; ++n) {
	    if (ipe[n * 5 + 1] > 0) {
		++me;
		iepw[me] = ipe[n * 5 + 5];
	    }
	}
	ic = countcolors_(&me, &iepw[1]);
	i__1 = ic;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s = domainvolume_(ne, &xyp[4], &ipe[6], &iepw[i__]);
	    if (*flagfile) {
		s_wsfi(&io___101);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&iepw[i__], (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&s, (ftnlen)sizeof(doublereal));
		e_wsfi();
		addout_(nlines, message, output + 100, (ftnlen)100, (ftnlen)
			100);
	    } else {
		s_wsfe(&io___102);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&iepw[i__], (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&s, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	}
    }
/* ... check the mesh (level 2) */
    if (*nf > 0) {
	if (*iprint >= 3) {
	    if (*flagfile) {
		s_wsfi(&io___103);
		do_fio(&c__1, (char *)&dvo, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&sao, (ftnlen)sizeof(doublereal));
		e_wsfi();
		addout_(nlines, message, output + 100, (ftnlen)100, (ftnlen)
			100);
	    } else {
		s_wsfe(&io___104);
		do_fio(&c__1, (char *)&dvo, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&sao, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	}
	ver = (d__1 = dvo - domainvolume_(ne, &xyp[4], &ipe[6], &c__0), dabs(
		d__1)) / dvo;
	aer = (d__1 = sao - surfacearea_(nf, &xyp[4], &ipf[5], &c__0), dabs(
		d__1)) / sao;
	if (ver > 1e-10 || aer > 1e-10 && nfold == *nf) {
	    wrnmes_(&c__5201, "ani2", "total volume or/and surface is not pr"
		    "eserved", (ftnlen)4, (ftnlen)44);
	}
	ver = (d__1 = fvo - fixedvolume_(nev, &xyp[4], &ipe[6], &iev[1]), dabs(
		d__1)) / dvo;
	aer = (d__1 = fao - fixedarea_(nfv, &xyp[4], &ipf[5], &ifv[1]), dabs(
		d__1)) / sao;
	if (ver > 1e-10 || aer > 1e-10) {
	    errmes_(&c__5201, "ani2", "fixed volume or/and surface is not pr"
		    "eserved", (ftnlen)4, (ftnlen)44);
	}
    }
    if (*iprint >= 2) {
	if (*flagfile) {
	    s_wsfi(&io___107);
	    do_fio(&c__1, (char *)&rr, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*status), (ftnlen)sizeof(integer));
	    e_wsfi();
	    addout_(nlines, message, output + 100, (ftnlen)100, (ftnlen)100);
	} else {
	    s_wsfe(&io___108);
	    do_fio(&c__1, (char *)&rr, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*status), (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }
/* ... display output */
/* ... parallel output */
L9000:
    return 0;
} /* ani2_ */

