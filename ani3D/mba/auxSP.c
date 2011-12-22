/* auxSP.f -- translated by f2c (version 20090411).
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

#include "error.h"
#include "makQ.h"
#include "auxSP.h"
#include "auxSE.h"
#include "auxSF.h"
#include "auxSR.h"

/* Table of constant values */

static integer c__1007 = 1007;
static integer c__6001 = 6001;
static integer c__1000 = 1000;
static integer c__4 = 4;
static integer c__1 = 1;
static integer c__64 = 64;
static integer c__512 = 512;
static integer c__16 = 16;
static integer c__2 = 2;
static integer c__8 = 8;

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int maksp_(integer *ip, integer *iep, integer *
	ipe, integer *iee, integer *maxs, integer *le, integer *ies)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k, n, n1, n2, ie;
    static logical repeat;

/* ================================================================ */
/* Remark: the 3rd column of IPE is overloaded in makSE. If it is */
/*         needed, the absolute value should be used. */

/*         We realize the same idea here by overloading the 3rd */
/*         column of IPE. */
/* ================================================================ */
/* group (M) */
/* group (S) */
/* group (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --ies;
    iee -= 5;
    ipe -= 6;
    --iep;

    /* Function Body */
    *le = 1;
    ies[*le] = iep[*ip];
    ie = ies[1];
    ipe[ie * 5 + 3] = -ipe[ie * 5 + 3];
    n2 = 0;
L1:
    repeat = FALSE_;
    n1 = n2 + 1;
    n2 = *le;
    i__1 = n2;
    for (n = n1; n <= i__1; ++n) {
	if (*le >= *maxs - 4) {
	    goto L1000;
	}
	for (i__ = 1; i__ <= 4; ++i__) {
	    ie = iee[i__ + (ies[n] << 2)];
	    if (ie == 0) {
		goto L2;
	    }
	    if (ipe[ie * 5 + 3] < 0) {
		goto L2;
	    }
	    if (*ip == ipe[ie * 5 + 1] || *ip == ipe[ie * 5 + 2] || *ip == 
		    ipe[ie * 5 + 3] || *ip == ipe[ie * 5 + 4]) {
/*              Do k = lE, 1, -1 */
/*                 If(iE.EQ.iEs(k)) Goto 2 */
/*              End do */
		repeat = TRUE_;
		++(*le);
		ies[*le] = ie;
		ipe[ie * 5 + 3] = -ipe[ie * 5 + 3];
	    }
L2:
	    ;
	}
/* L4: */
    }
    if (repeat) {
	goto L1;
    }
/* ... restoring the overloaded values */
    i__1 = *le;
    for (k = 1; k <= i__1; ++k) {
	ie = ies[k];
	ipe[ie * 5 + 3] = -ipe[ie * 5 + 3];
    }
    return 0;
L1000:
    errmes_(&c__1007, "makSP", "local variable MaxS is small", (ftnlen)5, (
	    ftnlen)28);
    return 0;
} /* maksp_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int chkspf_(integer *nm, integer *ip1, integer *
	ioperat, integer *icp, integer *iep, integer *ipe, integer *iee, 
	integer *lpf, integer *ipf)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, m, n, m2, ie, le, ies[1000], ipt, npt, nstep;

/* ================================================================ */
/*  Routine collects boundary points within NM mesh steps of iP1. */
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
/* group (M) */
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
/* group (S) */
/* group (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --ipf;
    iee -= 5;
    ipe -= 6;
    --iep;
    --icp;

    /* Function Body */
    if (*nm > 2) {
	errmes_(&c__6001, "chkSPf", "system error", (ftnlen)6, (ftnlen)12);
    }
    *lpf = 0;
    i__1 = *nm;
    for (nstep = 1; nstep <= i__1; ++nstep) {
	if (nstep == 1) {
	    m2 = 1;
	    ipf[1] = *ip1;
	} else {
	    m2 = *lpf;
	}
	i__2 = m2;
	for (m = 1; m <= i__2; ++m) {
	    maksp_(&ipf[m], &iep[1], &ipe[6], &iee[5], &c__1000, &le, ies);
	    i__3 = le;
	    for (n = 1; n <= i__3; ++n) {
		ie = ies[n - 1];
		for (i__ = 1; i__ <= 4; ++i__) {
		    ipt = ipe[i__ + ie * 5];
		    if (*ioperat == 4 && ipt == *ip1) {
			goto L10;
		    }
		    if (ifxnode_(&icp[ipt], &c__4)) {
			findse_(lpf, &ipf[1], &ipt, &npt);
			if (npt != 0) {
			    goto L10;
			}
			++(*lpf);
			if (*lpf > 1000) {
			    goto L1000;
			}
			ipf[*lpf] = ipt;
		    }
L10:
		    ;
		}
	    }
	}
    }
    return 0;
L1000:
    errmes_(&c__1007, "chkSPf", "local variable MaxS is small", (ftnlen)6, (
	    ftnlen)28);
    return 0;
} /* chkspf_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int info_(void)
{
    /* Initialized data */

    static integer idat[14] = { 83,116,111,110,101,32,70,108,111,119,101,114,
	    32,33 };

    /* System generated locals */
    char ch__1[1];

    /* Local variables */
    static integer i__;

    /* Fortran I/O blocks */
    static cilist io___19 = { 0, 6, 0, "(14A1)", 0 };


/* ================================================================ */
/* Routines prints ANI's logo. */
/* ================================================================ */
    s_wsfe(&io___19);
    for (i__ = 1; i__ <= 14; ++i__) {
	*(unsigned char *)&ch__1[0] = (char) idat[(0 + (0 + (i__ - 1 << 2))) /
		 4];
	do_fio(&c__1, ch__1, (ftnlen)1);
    }
    e_wsfe();
    return 0;
} /* info_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int chkspb_(integer *nm, integer *ld, integer *
	ids, integer *ln, integer *ins, integer *ioperat, integer *icp, 
	integer *iep, integer *ipe, integer *iee, integer *lpf, integer *ipf, 
	integer *ipbad, logical *flag__)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, j, k, l, m, n, m2, ie, le, ipb[1000], ipc, ies[1000], 
	    lpb, ipt, npt, nstep;

/* ================================================================ */
/* ================================================================ */
/*  Routine checks that each boundary point in array lPf can be */
/*  connected with an interior point (having color jInode) using */
/*  at most NM mesh edges. flag = T, when there is a point which */
/*  is surrounded by boundary points. */

/*  The first bad point is returned in iPbad. */
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
/* group (M) */
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
/* group (S) */
/* group (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --ipf;
    iee -= 5;
    ipe -= 6;
    --iep;
    --icp;
    ins -= 3;
    ids -= 3;

    /* Function Body */
    if (*nm > 2) {
	errmes_(&c__6001, "chkSPf", "system error", (ftnlen)6, (ftnlen)12);
    }
    *flag__ = TRUE_;
    i__1 = *lpf;
    for (l = 1; l <= i__1; ++l) {
	lpb = 1;
	ipb[0] = ipf[l];
	i__2 = *nm;
	for (nstep = 1; nstep <= i__2; ++nstep) {
	    m2 = lpb;
	    i__3 = m2;
	    for (m = 1; m <= i__3; ++m) {
		ipc = ipb[m - 1];
		maksp_(&ipc, &iep[1], &ipe[6], &iee[5], &c__1000, &le, ies);
		i__4 = le;
		for (n = 1; n <= i__4; ++n) {
		    ie = ies[n - 1];
		    for (i__ = 1; i__ <= 4; ++i__) {
			ipt = ipe[i__ + ie * 5];
			i__5 = *ld;
			for (k = 1; k <= i__5; ++k) {
			    if (ipc == ids[(k << 1) + 1] && ipt == ids[(k << 
				    1) + 2]) {
				goto L5;
			    }
			    if (ipc == ids[(k << 1) + 2] && ipt == ids[(k << 
				    1) + 1]) {
				goto L5;
			    }
			}
			if (ifxnode_(&icp[ipt], &c__64)) {
			    goto L10;
			}
			findse_(&lpb, ipb, &ipt, &npt);
			if (npt != 0) {
			    goto L5;
			}
			++lpb;
			if (lpb > 1000) {
			    goto L1000;
			}
			ipb[lpb - 1] = ipt;
			i__5 = *ln;
			for (k = 1; k <= i__5; ++k) {
			    for (j = 1; j <= 2; ++j) {
				if (ipc == ins[j + (k << 1)]) {
				    ipt = ins[3 - j + (k << 1)];
				    if (ifxnode_(&icp[ipt], &c__64)) {
					goto L10;
				    }
				    if (nstep < *nm) {
					findse_(&lpb, ipb, &ipt, &npt);
					if (npt != 0) {
					    goto L5;
					}
					++lpb;
					if (lpb > 1000) {
					    goto L1000;
					}
					ipb[lpb - 1] = ipt;
				    }
/*  ...  the next line restricted the search for inner nodes and was removed */
/*                              Goto 5 */
				}
			    }
			}
L5:
			;
		    }
		}
	    }
	}
	*ipbad = ipb[0];
	return 0;
L10:
	;
    }
    *flag__ = FALSE_;
    return 0;
L1000:
    errmes_(&c__1007, "chkSP", "local variable MaxS is small", (ftnlen)5, (
	    ftnlen)28);
    return 0;
} /* chkspb_ */

/* ================================================================ */
/* @f2h@ */ doublereal calnorm_(doublereal *xyz)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3;

/* ================================================================ */
/* Routines computes 2-norm of vector xyz. */
/* ================================================================ */
    /* Parameter adjustments */
    --xyz;

    /* Function Body */
/* Computing 2nd power */
    d__1 = xyz[1];
/* Computing 2nd power */
    d__2 = xyz[2];
/* Computing 2nd power */
    d__3 = xyz[3];
    ret_val = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
    return ret_val;
} /* calnorm_ */

/* ================================================================ */
/* @f2h@ */ doublereal sqrnorm_(doublereal *xyz)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3;

/* ================================================================ */
/* Routines computes square of 2-norm of vector xyz. */
/* ================================================================ */
    /* Parameter adjustments */
    --xyz;

    /* Function Body */
/* Computing 2nd power */
    d__1 = xyz[1];
/* Computing 2nd power */
    d__2 = xyz[2];
/* Computing 2nd power */
    d__3 = xyz[3];
    ret_val = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
    return ret_val;
} /* sqrnorm_ */

/* ================================================================ */
/*  All operations with colors are binary operations. */
/*  The same operations go with status. */
/* ================================================================ */
/* ================================================================ */
/*     Logical Function ifXstatus(status, iXstatus) */
/* @f2h@ */ logical ifxnode_(integer *clr, integer *ixnode)
{
    /* System generated locals */
    logical ret_val;

/* ================================================================ */
    ret_val = (*clr & *ixnode) == *ixnode;
    return ret_val;
} /* ifxnode_ */

/* ================================================================ */
/*     Subroutine addXstatus(status, iXstatus) */
/* @f2h@ */ /* Subroutine */ int addxnode_(integer *clr, integer *ixnode)
{
/* ================================================================ */
/* Binary operation +. */
/* ================================================================ */
    *clr |= *ixnode;
    return 0;
} /* addxnode_ */

/* ================================================================ */
/*     Subroutine delXstatus(status, iXstatus) */
/* @f2h@ */ /* Subroutine */ int delxnode_(integer *clr, integer *ixnode)
{
/* ================================================================ */
    *clr -= *clr & *ixnode;
    return 0;
} /* delxnode_ */

/* ================================================================ */
/* @f2h@ */ integer minclr_(integer *clr1, integer *clr2)
{
    /* System generated locals */
    integer ret_val;

/* ================================================================ */
/*  The function returns common color for both clr1 and clr2. */
/* ================================================================ */
    ret_val = *clr1 & *clr2;
    return ret_val;
} /* minclr_ */

/* ================================================================ */
/* @f2h@ */ integer maxclr_(integer *clr1, integer *clr2)
{
    /* System generated locals */
    integer ret_val;

/* ================================================================ */
/*  The function returns minimal color containing both clr1 and clr2. */
/* ================================================================ */
    ret_val = *clr1 | *clr2;
    return ret_val;
} /* maxclr_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int setstatus_(logical *flagauto, integer *
	status, integer *iprint)
{
    /* Format strings */
    static char fmt_5001[] = "(\002status.fd: -512  [ANITangledMesh]        "
	    "    [research]\002)";
    static char fmt_5003[] = "(\002status.fd: +1    [ANIForbidBoundaryElemen"
	    "ts] [user]\002)";
    static char fmt_5004[] = "(\002status.fd: +4    [ANIFixBoundaryFaces]   "
	    "    [user]\002)";
    static char fmt_5007[] = "(\002status.fd: +16   [ANIFixSurfacePoints]   "
	    "    [user]\002)";
    static char fmt_5002[] = "(\002status.fd: +64   [ANIMultiConnectedGeomet"
	    "ry] \002,a)";
    static char fmt_5005[] = "(\002status.fd: +2    [ANIUse2ArmRule]        "
	    "    \002,a)";
    static char fmt_5006[] = "(\002status.fd: +8    [ANIDeleteTemporaryFaces"
	    "]   \002,a)";

    /* Fortran I/O blocks */
    static cilist io___37 = { 0, 6, 0, fmt_5001, 0 };
    static cilist io___38 = { 0, 6, 0, fmt_5003, 0 };
    static cilist io___39 = { 0, 6, 0, fmt_5004, 0 };
    static cilist io___40 = { 0, 6, 0, fmt_5007, 0 };
    static cilist io___41 = { 0, 6, 0, fmt_5002, 0 };
    static cilist io___42 = { 0, 6, 0, fmt_5005, 0 };
    static cilist io___43 = { 0, 6, 0, fmt_5006, 0 };
    static cilist io___44 = { 0, 6, 0, fmt_5005, 0 };
    static cilist io___45 = { 0, 6, 0, fmt_5006, 0 };
    static cilist io___46 = { 0, 6, 0, 0, 0 };


/* ================================================================ */
/* ================================================================ */
/* Routine adds additional properties to the variable status */
/* ================================================================ */
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
/* ================================================================ */
    *status = max(0,*status);
/* ... remove obsolete and not-implemented input features */
    if (*status > 0) {
	if (ifxnode_(status, &c__512)) {
	    delxnode_(status, &c__512);
	    if (*iprint >= 2) {
		s_wsfe(&io___37);
		e_wsfe();
	    }
	}
    }
/* ... inform the user about requested features */
    if (*iprint >= 2) {
	if (ifxnode_(status, &c__1)) {
	    s_wsfe(&io___38);
	    e_wsfe();
	}
	if (ifxnode_(status, &c__4)) {
	    s_wsfe(&io___39);
	    e_wsfe();
	}
	if (ifxnode_(status, &c__16)) {
	    s_wsfe(&io___40);
	    e_wsfe();
	}
	if (ifxnode_(status, &c__64)) {
	    s_wsfe(&io___41);
	    do_fio(&c__1, "[user]", (ftnlen)6);
	    e_wsfe();
	}
	if (ifxnode_(status, &c__2)) {
	    s_wsfe(&io___42);
	    do_fio(&c__1, "[user]", (ftnlen)6);
	    e_wsfe();
	}
	if (ifxnode_(status, &c__8)) {
	    s_wsfe(&io___43);
	    do_fio(&c__1, "[user]", (ftnlen)6);
	    e_wsfe();
	}
    }
/* ... set up default features */
    if (*flagauto) {
	if (! ifxnode_(status, &c__2)) {
	    if (*iprint >= 2) {
		s_wsfe(&io___44);
		do_fio(&c__1, "[system]", (ftnlen)8);
		e_wsfe();
	    }
	    addxnode_(status, &c__2);
	}
	if (! ifxnode_(status, &c__8)) {
	    if (*iprint >= 2) {
		s_wsfe(&io___45);
		do_fio(&c__1, "[system]", (ftnlen)8);
		e_wsfe();
	    }
	    addxnode_(status, &c__8);
	}
    }
    if (*iprint >= 2) {
	s_wsle(&io___46);
	e_wsle();
    }
    return 0;
    return 0;
} /* setstatus_ */

