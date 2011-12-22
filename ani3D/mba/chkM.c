/* chkM.f -- translated by f2c (version 20090411).
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

/* Table of constant values */

static integer c__1 = 1;
static integer c__5001 = 5001;
static integer c__5002 = 5002;
static integer c__5003 = 5003;
static integer c__4 = 4;
static integer c__8 = 8;
static integer c__5103 = 5103;
static integer c__64 = 64;
static integer c__4103 = 4103;
static integer c__5006 = 5006;
static integer c__5022 = 5022;
static integer c__5023 = 5023;
static integer c__512 = 512;
static integer c__5007 = 5007;
static integer c__5 = 5;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__5008 = 5008;

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int chkm_(integer *np, integer *nf, integer *ne, 
	doublereal *xyp, integer *ipf, integer *ipe, integer *icp, integer *
	ife, integer *iee, doublereal *rr, integer *status, integer *icheck, 
	integer *nepw, integer *iepw)
{
    /* Format strings */
    static char fmt_5002[] = "(\002Error in checking face =\002,i6,\002   iE"
	    "RR=\002,i5,/,\002Points =\002,4i6,/,\002colors =\002,3i6)";
    static char fmt_4000[] = "(/,\002=============== ERROR details ========="
	    "========\002)";
    static char fmt_5000[] = "(\002Error in checking tet =\002,i7,\002  vo"
	    "l =\002,e16.9,\002   iERR=\002,i4,/,\002Points =\002,4i7,/,\002F"
	    "aces  =\002,4i7,/,\002Tetras =\002,4i7,/,\002colors =\002,4i7,/)";
    static char fmt_5004[] = "(\002Face \002,i7,\002  (\002,i6,\002)  of the"
	    " bad tetrahedron\002,/,\002Points =\002,3i7,/,\002colors =\002,3"
	    "i7,/)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, k, n, i1, i2, i3, j1, j2, j3, i4, j4;
    static doublereal v1;
    static integer ic, ie, if__, ip[5];
    static doublereal vv;
    static integer ie1, if1, if2, if3, if4, ie2, ie3, ie4, ip1, ip2, ip3, ip4,
	     jp1, jp2, jp3, jp4;
    static doublereal rin;
    static integer ipt, nfac;
    static logical flag__;
    static integer ierr, ntet;
    static doublereal rout;
    static integer iclrf;
    static logical flagfbe;
    static real crvprec;

    /* Fortran I/O blocks */
    static cilist io___40 = { 0, 6, 0, 0, 0 };
    static cilist io___43 = { 0, 6, 0, fmt_5002, 0 };
    static cilist io___44 = { 0, 6, 0, fmt_4000, 0 };
    static cilist io___45 = { 0, 6, 0, fmt_5000, 0 };
    static cilist io___48 = { 0, 6, 0, fmt_5000, 0 };
    static cilist io___49 = { 0, 6, 0, fmt_5004, 0 };


/* ================================================================ */
/* group (M) */
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
/* Routine checks topology of the input and output meshes. */

/* iCheck : 0 - check everything */
/*          1 - don't check for isolated faces and points marked as */
/*              destroyed structures */

/* nEPw  :  working array of size nP */
/* IEPw  :  working array of size 4*nE */

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
/* group (W) */
/* ================================================================ */
/* group (Local function) */
/* group (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --iepw;
    --nepw;
    iee -= 5;
    ife -= 5;
    --icp;
    ipe -= 6;
    ipf -= 5;
    xyp -= 4;

    /* Function Body */
    ip[0] = 1;
    ip[1] = 2;
    ip[2] = 3;
    ip[3] = 4;
    ip[4] = 1;
    *rr = 0.;
    crvprec = 2e-8f;
    flagfbe = ifxnode_(status, &c__1);
/* ... check for face markers */
    i__1 = *nf;
    for (n = 1; n <= i__1; ++n) {
	iclrf = ipf[(n << 2) + 4];
	if (iclrf <= 0 && *icheck == 0) {
	    errmes_(&c__5001, "chkM", "wrong face ID", (ftnlen)4, (ftnlen)13);
	}
	if (iclrf > 1000) {
	    errmes_(&c__5002, "chkM", "face ID is out of limits", (ftnlen)4, (
		    ftnlen)24);
	}
    }
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	if (*icheck == 0 && icp[n] == 0) {
	    errmes_(&c__5003, "chkM", "wrong point color", (ftnlen)4, (ftnlen)
		    17);
	}
	if (ifxnode_(&icp[n], &c__4)) {
	    if (! ifxnode_(&icp[n], &c__8)) {
		errmes_(&c__5103, "chkM", "wrong point color", (ftnlen)4, (
			ftnlen)17);
	    }
	    if (ifxnode_(&icp[n], &c__64)) {
		errmes_(&c__5103, "chkM", "wrong point color", (ftnlen)4, (
			ftnlen)17);
	    }
	}
    }
    i__1 = *nf;
    for (n = 1; n <= i__1; ++n) {
	nfac = n;
	for (i__ = 1; i__ <= 4; ++i__) {
	    if (*icheck == 0 && ipf[i__ + (n << 2)] <= 0) {
		ierr = 5012;
		goto L400;
	    }
	}
	if (*icheck == 0 && ipf[(n << 2) + 4] >= 900) {
	    errmes_(&c__4103, "chkM", "reserved boundary identificator is us"
		    "ed", (ftnlen)4, (ftnlen)39);
	}
    }
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	nepw[n] = 0;
    }
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	ntet = n;
	flag__ = flagfbe;
	for (i__ = 1; i__ <= 4; ++i__) {
	    ipt = ipe[i__ + n * 5];
	    if (ipt <= 0) {
		errmes_(&c__5006, "chkM", "wrong connectivity table", (ftnlen)
			4, (ftnlen)24);
	    }
	    ++nepw[ipt];
	    if (! ifxnode_(&icp[ipt], &c__4)) {
		flag__ = FALSE_;
	    }
	}
	if (*icheck == 0 && flag__) {
	    errmes_(&c__5022, "chkM", "boundary element", (ftnlen)4, (ftnlen)
		    16);
	}
	ip1 = ipe[n * 5 + 1];
	ip2 = ipe[n * 5 + 2];
	ip3 = ipe[n * 5 + 3];
	ip4 = ipe[n * 5 + 4];
	if (ip1 == ip2 || ip1 == ip3 || ip1 == ip4 || ip2 == ip3 || ip2 == 
		ip4 || ip3 == ip4) {
	    ierr = 5013;
	    goto L500;
	}
	if (ipe[n * 5 + 5] <= 0) {
	    errmes_(&c__5023, "chkM", "non-positive element label", (ftnlen)4,
		     (ftnlen)26);
	}
	if1 = ife[(n << 2) + 1];
	if2 = ife[(n << 2) + 2];
	if3 = ife[(n << 2) + 3];
	if4 = ife[(n << 2) + 4];
	if (if1 == if2 && if1 > 0 || if1 == if3 && if1 > 0 || if1 == if4 && 
		if1 > 0 || if2 == if3 && if2 > 0 || if2 == if4 && if2 > 0 || 
		if3 == if4 && if3 > 0) {
	    ierr = 5014;
	    goto L500;
	}
	ie1 = iee[(n << 2) + 1];
	ie2 = iee[(n << 2) + 2];
	ie3 = iee[(n << 2) + 3];
	ie4 = iee[(n << 2) + 4];
	if (ie1 == ie2 && ie1 != 0 || ie1 == ie3 && ie1 != 0 || ie1 == ie4 && 
		ie1 != 0 || ie2 == ie3 && ie2 != 0 || ie2 == ie4 && ie2 != 0 
		|| ie3 == ie4 && ie3 != 0) {
	    ierr = 5015;
	    goto L500;
	}
	for (i1 = 1; i1 <= 4; ++i1) {
	    if__ = ife[i1 + (n << 2)];
	    ie = iee[i1 + (n << 2)];
	    if (if__ == 0 && ie == 0) {
		ierr = 5016;
		goto L500;
	    }
	    i2 = ip[i1];
	    i3 = ip[i2];
	    ip1 = ipe[i1 + n * 5];
	    ip2 = ipe[i2 + n * 5];
	    ip3 = ipe[i3 + n * 5];
	    if (if__ != 0 && if__ != -1) {
		jp1 = ipf[(if__ << 2) + 1];
		jp2 = ipf[(if__ << 2) + 2];
		jp3 = ipf[(if__ << 2) + 3];
		if (! check33_(&ip1, &ip2, &ip3, &jp1, &jp2, &jp3)) {
		    ierr = 5017;
		    goto L500;
		}
	    }
	    if (ie != 0) {
		if (if__ != 0) {
		    for (j1 = 1; j1 <= 4; ++j1) {
			if (ife[j1 + (ie << 2)] == if__) {
			    goto L10;
			}
		    }
		    ierr = 5018;
		    goto L500;
		}
L10:
		for (j1 = 1; j1 <= 4; ++j1) {
		    j2 = ip[j1];
		    j3 = ip[j2];
		    jp1 = ipe[j1 + ie * 5];
		    jp2 = ipe[j2 + ie * 5];
		    jp3 = ipe[j3 + ie * 5];
		    if (check33_(&ip1, &ip2, &ip3, &jp1, &jp2, &jp3)) {
			if (iee[j1 + (ie << 2)] != n) {
			    ierr = 5019;
			    goto L500;
			}
			i4 = ip[i3];
			ip4 = ipe[i4 + n * 5];
			j4 = ip[j3];
			jp4 = ipe[j4 + ie * 5];
/*                    v1 = calVol(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3), XYP(1, iP4)) */
/*                    v2 = calVol(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3), XYP(1, jP4)) */
			vv = mutualorientation_(&xyp[ip1 * 3 + 1], &xyp[ip2 * 
				3 + 1], &xyp[ip3 * 3 + 1], &xyp[ip4 * 3 + 1], 
				&xyp[jp4 * 3 + 1]);
			if (vv >= 0. && ! ifxnode_(status, &c__512)) {
			    ierr = 5020;
			    goto L500;
			}
			goto L20;
		    }
		}
		ierr = 5021;
		goto L500;
	    }
L20:
	    ;
	}
    }
    if (*icheck == 0) {
	i__1 = *np;
	for (n = 1; n <= i__1; ++n) {
	    if (nepw[n] == 0) {
		errmes_(&c__5007, "chkM", "isolated point", (ftnlen)4, (
			ftnlen)14);
	    }
	}
    }
/* ... check color of edge points */
    backreferences_(np, ne, &c__4, &c__5, &ipe[6], &nepw[1], &iepw[1]);
    i__1 = *ne << 2;
    for (n = 1; n <= i__1; ++n) {
	iepw[n] = ipe[iepw[n] * 5 + 5];
    }
    i2 = 0;
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	i1 = i2 + 1;
	i2 = nepw[n];
	i__2 = i2 - i1 + 1;
	ic = countcolors_(&i__2, &iepw[i1]);
	if (ic >= 3 && ! (ifxnode_(&icp[n], &c__2) || ifxnode_(&icp[n], &c__1)
		)) {
	    s_wsle(&io___40);
	    do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
	    do_lio(&c__3, &c__1, (char *)&icp[n], (ftnlen)sizeof(integer));
	    do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(integer));
	    e_wsle();
	    errmes_(&c__5008, "chkM", "wrong color of edge point", (ftnlen)4, 
		    (ftnlen)25);
	}
    }
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	ip1 = ipe[n * 5 + 1];
	ip2 = ipe[n * 5 + 2];
	ip3 = ipe[n * 5 + 3];
	ip4 = ipe[n * 5 + 4];
	randr_(&xyp[ip1 * 3 + 1], &xyp[ip2 * 3 + 1], &xyp[ip3 * 3 + 1], &xyp[
		ip4 * 3 + 1], &rout, &rin);
/* Computing MAX */
	d__1 = *rr, d__2 = rout / rin;
	*rr = max(d__1,d__2);
    }
    return 0;
L400:
    s_wsfe(&io___43);
    do_fio(&c__1, (char *)&nfac, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ierr, (ftnlen)sizeof(integer));
    for (i__ = 1; i__ <= 4; ++i__) {
	do_fio(&c__1, (char *)&ipf[i__ + (nfac << 2)], (ftnlen)sizeof(integer)
		);
    }
    for (i__ = 1; i__ <= 3; ++i__) {
	do_fio(&c__1, (char *)&icp[ipf[i__ + (nfac << 2)]], (ftnlen)sizeof(
		integer));
    }
    e_wsfe();
    errmes_(&ierr, "chkM", "faces are wrong", (ftnlen)4, (ftnlen)15);
L500:
    s_wsfe(&io___44);
    e_wsfe();
    s_wsfe(&io___45);
    do_fio(&c__1, (char *)&ntet, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&v1, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&ierr, (ftnlen)sizeof(integer));
    for (i__ = 1; i__ <= 4; ++i__) {
	do_fio(&c__1, (char *)&ipe[i__ + ntet * 5], (ftnlen)sizeof(integer));
    }
    for (i__ = 1; i__ <= 4; ++i__) {
	do_fio(&c__1, (char *)&ife[i__ + (ntet << 2)], (ftnlen)sizeof(integer)
		);
    }
    for (i__ = 1; i__ <= 4; ++i__) {
	do_fio(&c__1, (char *)&iee[i__ + (ntet << 2)], (ftnlen)sizeof(integer)
		);
    }
    for (i__ = 1; i__ <= 4; ++i__) {
	do_fio(&c__1, (char *)&icp[ipe[i__ + ntet * 5]], (ftnlen)sizeof(
		integer));
    }
    e_wsfe();
    for (k = 1; k <= 4; ++k) {
	ie = iee[k + (ntet << 2)];
	if (ie > 0) {
	    ip1 = ipe[ie * 5 + 1];
	    ip2 = ipe[ie * 5 + 2];
	    ip3 = ipe[ie * 5 + 3];
	    ip4 = ipe[ie * 5 + 4];
	    v1 = calvol_(&xyp[ip1 * 3 + 1], &xyp[ip2 * 3 + 1], &xyp[ip3 * 3 + 
		    1], &xyp[ip4 * 3 + 1]);
	    s_wsfe(&io___48);
	    do_fio(&c__1, (char *)&ie, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&v1, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&ierr, (ftnlen)sizeof(integer));
	    for (i__ = 1; i__ <= 4; ++i__) {
		do_fio(&c__1, (char *)&ipe[i__ + ie * 5], (ftnlen)sizeof(
			integer));
	    }
	    for (i__ = 1; i__ <= 4; ++i__) {
		do_fio(&c__1, (char *)&ife[i__ + (ie << 2)], (ftnlen)sizeof(
			integer));
	    }
	    for (i__ = 1; i__ <= 4; ++i__) {
		do_fio(&c__1, (char *)&iee[i__ + (ie << 2)], (ftnlen)sizeof(
			integer));
	    }
	    for (i__ = 1; i__ <= 4; ++i__) {
		do_fio(&c__1, (char *)&icp[ipe[i__ + ie * 5]], (ftnlen)sizeof(
			integer));
	    }
	    e_wsfe();
	}
    }
    for (k = 1; k <= 4; ++k) {
	if__ = ife[k + (ntet << 2)];
	if (if__ > 0) {
	    s_wsfe(&io___49);
	    do_fio(&c__1, (char *)&if__, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*nf), (ftnlen)sizeof(integer));
	    for (i__ = 1; i__ <= 3; ++i__) {
		do_fio(&c__1, (char *)&ipf[i__ + (if__ << 2)], (ftnlen)sizeof(
			integer));
	    }
	    for (i__ = 1; i__ <= 3; ++i__) {
		do_fio(&c__1, (char *)&icp[ipf[i__ + (max(1,if__) << 2)]], (
			ftnlen)sizeof(integer));
	    }
	    e_wsfe();
	}
    }
/*     Do n = 1, nE */
/*        IEPw(n) = 1 */
/*     End do */
/*     Call saveMgmv(nP, nF, nE, */
/*    &              XYP, IPF, IPE, IEPw, IEPw, */
/*    &              'error.gmv', IEPw(nE+1)) */
    errmes_(&ierr, "chkM", "tetrahedra are wrong", (ftnlen)4, (ftnlen)20);
/* 5001 Format('Error =', E12.4, /, */
/*    &       'Curvilinear face=', I4, ' attributs=', 2I4) */
    return 0;
} /* chkm_ */

