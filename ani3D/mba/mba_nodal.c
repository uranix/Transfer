/* mba_nodal.f -- translated by f2c (version 20090411).
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

#include "forlibmba.h"
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

#include "time_mba.h"

/* Common Block Declarations */

struct {
    doublereal refxyp[3], scaxyp[3];
} anicrv_;

#define anicrv_1 anicrv_

/* Table of constant values */

static integer c__4201 = 4201;
static integer c__1 = 1;
static logical c_true = TRUE_;
static integer c__60 = 60;
static logical c_false = FALSE_;

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int mbanodal_(integer *np, integer *maxp, 
	integer *nf, integer *maxf, integer *ne, integer *maxe, doublereal *
	xyp, integer *ipf, integer *ipe, integer *lbf, integer *lbe, integer *
	nestar, integer *npv, integer *nfv, integer *nev, integer *ipv, 
	integer *ifv, integer *iev, logical *flagauto, integer *status, 
	integer *maxskipe, integer *maxqitr, doublereal *metric, doublereal *
	quality, doublereal *rquality, integer *maxwr, integer *maxwi, 
	doublereal *rw, integer *iw, integer *iprint, integer *ierr)
{
    /* Format strings */
    static char fmt_5001[] = "(\002MBA: STONE FLOWER! (1997-2010), version 2"
	    ".3\002,/,5x,\002Target: Quality=\002,f4.2,\002 (nEStar:\002,i7"
	    ",\002, SkipE:\002,i5,\002, maxITR:\002,i8,\002)\002)";
    static char fmt_5003[] = "(\002Total:\002,i6,\002 Q=\002,e10.4,\002  #V#"
	    "F#E:\002,i7,i8,i9,\002  tm=\002,f6.1,\002s\002,/)";

    /* System generated locals */
    integer i__1, i__2;
    olist o__1;
    cllist cl__1;

    /* Local variables */
    static logical flagfile;
    static integer milintrp, mrlintrp, i__, m, n, chanelout;
    static doublereal tm0, tm1;
    static integer idg, iqe, statuswork, new__, npw, il1e, il2e, iiee, iife, 
	    iicp, iipe, iipf, iiep, iise, ilbp, iipp, irse;
    static logical flaganalytic;
    static integer ihesp, iiepw, iipew;
    static doublereal hstar;
    static integer inepw, nloop, nqitr, ixypw, iihole, iiholf, iiholp, nlines;
    static char output[100*100], message[80];
    static integer iprintl;

    /* Fortran I/O blocks */
    static cilist io___6 = { 0, 0, 0, fmt_5001, 0 };
    static cilist io___7 = { 0, 6, 0, fmt_5001, 0 };
    static icilist io___29 = { 0, message, 0, "(A,I10)", 80, 1 };
    static icilist io___36 = { 0, message, 0, "(A,I10)", 80, 1 };
    static cilist io___47 = { 0, 0, 0, "(A)", 0 };
    static cilist io___48 = { 0, 0, 0, fmt_5003, 0 };
    static cilist io___49 = { 0, 6, 0, fmt_5003, 0 };


/* ================================================================ */
/* group (M) */
/* group (Dev) */
/* group (Q) */
/* group (W) */
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
/*  VARIABLES & PARAMETER: */
/*     All input arguments are divided into two groups. The first */
/*     group (BASIC) containes the main parameters. The second */
/*     group (AUX) contains auxiliary parameters which can be */
/*     omitted. */

/*  ***BASIC 1: nP to IPF */
/*     nP, MaxP  - actual and maximal number of points (P) */
/*     nF, MaxF  - actual and maximal number of boundary and inner faces (F) */
/*     nE, MaxE  - actual and maximal number of elements (E) */

/*     XYP(3, MaxP) - Cartesian coordinates of mesh points */
/*     IPE(4, MaxE) - connectivity list of element */
/*     IPF(3, MaxF) - connectivity list of faces */

/*  ***AUX 1: lbF to lbE. */
/*     lbF(MaxF) - boundary identificator */
/*                   (a positive number strictly less than 1400=MaxS-100, */
/*                    colors from MaxS-99 to MaxS are reserved) */
/*                   (example: unit cube has 6 boundaries which */
/*                    may have or not different identificators) */

/*     lbE(MaxE) - element indentificator (a positive number) */

/*  ***BASIC 2: nEStar */
/*     nEstar    - target number of elements in the final mesh */

/*  ***AUX 2: nPv to MaxSkipE */
/*     nPv       - number of fixed points     (can be zero) */
/*     nFv       - number of fixed triangles  (can be zero) */
/*     nEv       - number of fixed tetrahedra (can be zero) */

/*     IPV(nPv)  - list of fixed points */
/*     IFV(nFv)  - list of fixed triangles */
/*     IEV(nEv)  - list of fixed tetrahedra */

/*     flagAuto  - flag controling the mesh generation: */
/*                 TRUE  - automatic recovering of missing mesh elements */
/*                 FALSE - rigorous checking of input data */

/*     status    - advanced control of mesh generation: */
/*                 0 or negative - no additional control */
/*                 positive      - enforce some of the mesh properties */
/*                    The detailed description of available properties */
/*                    is in file status.fd. Variable status is equal */
/*                    to the sum of positive numbers corresponding */
/*                    to the desired properties. Here is a short list */
/*                    of user-controled properties: */

/*                     1 - final mesh will not have boundary elements */
/*                     4 - subroutine will not change boundary faces */
/*                     8 - missing material interface and boundary */
/*                         faces created by the code will be removed */
/*                         from the final mesh */
/*                    16 - subroutine will not change surface points */
/*                   512 - tangled mesh will be fixed (testing) */

/*                    If only first two properties are required, set */
/*                    status = 5. Do not use the constant since status */
/*                    is the input/output parameter. */


/*     MaxSkipE  - the maximal number of skipped elements. */
/*                   The bad elements are put is a basket if their */
/*                   cannot be fixed. These bad elements will be */
/*                   skipped (temporary) from analysis. The basket */
/*                   size is MaxSkipE. When it is full, all elements */
/*                   are released back into the mesh. */

/*  ***BASIC 3: MaxQItr to iERR */
/*     MaxQItr   - the maximal number of local grid modifications */

/*     Metric(6,nP) - metric defined at mesh nodes. The metric is */
/*                    a 3x3 symmetric matrix M_ij. Each column of */
/*                    this array stores the upper diagonal entries */
/*                    in the following order: M_11, M_22, M_33, */
/*                    M_12, M_23, M_13. */

/*     Quality   - target quality for the final grid */
/*                    (a positive number between 0 and 1) */
/*     rQuality  - actual quality of the final grid */

/*     MaxWr     - size of the working real*8 array */
/*     MaxWi     - size of the wirking integer array */

/*     rW(MaxWr) - real*8  working array. On output: */
/*                   rW(1) - mesh generation time */
/*                   rW(2) - rQuality */
/*                   rW(3) - average size of elements with respect */
/*                           to the given metrici (hStar in theory) */
/*                   rW(4:nE+3) - quality of elements */


/*     iW(MaxWi) - integer working array. On output: */
/*                 iW(1 : nP) - colors of mesh points as described */
/*                               in file color.fd */
/*                     T - tought point:   cannot be changed anyway */
/*                                         T-edges cannot be changed */
/*                     V - vertex:         cannot be changed anyway */
/*                    VB - edge point:     has to live on the edge */
/*                     B - surface point:  has to live on the surface */
/*                     I - internal point: free point inside the domain */

/*                   VBI - internal edge point:     similar to VB */
/*                    BI - internal boundary point: similar to  B */

/*                 iW(nP + 1) - the number of performed iterations */

/*     iPrint - nonnegative number 10 * A + B: */
/*                B - level of output information: */
/*                      0 - nothing */
/*                      9 - maximum */

/*                A - output chanel number (between 0 and 99): */
/*                      A > 0 - the output goes to file aniMPI.log */
/*                      A = 0 - screen output is used */


/*     iERR   - error code: */
/*              0  - correct program termination */
/*           1000  - target quality was not reached (warning) */
/*           other - error as desribed in error.f */

/* ================================================================ */

/*  Note: */
/*       Input parameters:  MaxP, MaxF, MaxE, */
/*                          nEStar, nPv, nFv, nEv, IPV, IFV, IEV, */
/*                          flagAuto, MaxSkipE, MaxQItr, */
/*                          Metric, Quality, MaxWr, MaxWi, iPrint */

/*       Input / output:    nP, nF, nE, XYP, IPF, IPE, */
/*                          lbF, lbE, rW, iW */

/*       Output parameters: rQuality, iERR */

/* ================================================================ */

/*  A possible choice for input parameters: */
/*       MaxP > nP */
/*       MaxF > nF */
/*       MaxE > nE */
/*       MaxSkipE = 200 */
/*       MaxQItr = 25000 */
/*       Quality = 0.3 */
/*       MaxWr > 14 * MaxP + 2 * nP + MaxE  (approximate estimate) */
/*       MaxWi >  7 * MaxP + nP + 7 * MaxF + 18 * MaxE + 13 * nE */
/*                                                (approximately) */

/* ================================================================ */
/*  A short version (mbaNodalShort) can be used when there are */
/*  no curvilinear boundaries. The following parameters are set */
/*  by default and have to be changed if necessary: */

/*       MaxSkipE =  200 */
/*       MaxQItr  =  min(50000, nE) */
/*       Quality  =  0.3 */

/*  Remark: */
/*    The reached mesh quality is returned in rW(2). */

/* ================================================================ */

/* *** Authors: K. Lipnikov     (lipnikov@hotmail.com) */
/*              Yu. Vassilevski (vasilevs@dodo.inm.ras.ru) */
/* *** Date:    1997 - 2004 */
/* *** Updates: see ChangeLog */
/* *** External routines: DSORT, DSYEV */

/* ========================================================== */
/* group (M) */
/*     Integer MaxP, MaxF, MaxE */
/* ========================================================== */
/* MaxLines   - the maximal number of output lines per processor */
/* LineLenght - the maximal lenght of a message */
/* ========================================================== */
/* group (Dev) */
/* group (Q) */
/*     Integer MaxSkipE, MaxQItr */
/* group (W) */
/* group (Local variables) */
/* ================================================================ */
/* ================================================================ */
    /* Parameter adjustments */
    --iw;
    --rw;
    metric -= 7;
    --iev;
    --ifv;
    --ipv;
    --lbe;
    --lbf;
    ipe -= 5;
    ipf -= 4;
    xyp -= 4;

    /* Function Body */
    *ierr = 0;
    iprintl = *iprint % 10;
    chanelout = *iprint / 10;
    flagfile = *iprint >= 10;
    if (chanelout >= 100) {
	errmes_(&c__4201, "ani_solution", "output chanel number is wrong", (
		ftnlen)12, (ftnlen)29);
    }
    if (flagfile) {
	o__1.oerr = 0;
	o__1.ounit = chanelout;
	o__1.ofnmlen = 10;
	o__1.ofnm = "aniMPI.log";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
    }
    nlines = 0;
    nloop = 1;
/* ... print header */
/* hStar calculations */
    if (iprintl >= 1) {
	if (flagfile) {
	    io___6.ciunit = chanelout;
	    s_wsfe(&io___6);
	    do_fio(&c__1, (char *)&(*quality), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*nestar), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*maxskipe), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*maxqitr), (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    s_wsfe(&io___7);
	    do_fio(&c__1, (char *)&(*quality), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*nestar), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*maxskipe), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*maxqitr), (ftnlen)sizeof(integer));
	    e_wsfe();
	}
    }
    setstatus_(flagauto, status, iprint);
    statuswork = *status;
/* ... starting clocks */
    seconds_(&tm0);
/* ... memory for sequantial running (is cleaned) */
/*     iW(1) is overloaded to return colors & nQItr */
    iipf = 1;
    iipe = iipf + (*maxf << 2);
    iipp = iipe + *maxe * 5;
    iicp = iipp + *maxp;
    iiholp = iicp + *maxp;
    iiholf = iiholp + *maxp;
    iihole = iiholf + *maxf;
    iiep = iihole + *maxe;
    iife = iiep + *maxp;
    iiee = iife + (*maxe << 2);
    il1e = iiee + (*maxe << 2);
    il2e = il1e + (*maxe << 1);
    iiepw = il2e + *maxe;
/* Computing MAX */
    i__1 = *maxe << 2, i__2 = *maxf * 3;
    inepw = iiepw + max(i__1,i__2);
    iipew = inepw + *maxp;
    iise = iipew + (*ne << 2);
    ilbp = 1;
    milintrp = *maxwi - iise;
    if (milintrp <= *maxp + *maxf) {
	*ierr = 1001;
	s_wsfi(&io___29);
	do_fio(&c__1, "The approximate size of iW is ", (ftnlen)30);
	i__1 = iise + *maxp + *maxf;
	do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	e_wsfi();
	errmes_(ierr, "mbaNodal", message, (ftnlen)8, (ftnlen)80);
    }
/* ... memory for sequantial running (is cleaned) */
/*     rW(1) is overloaded to return total time */
    ihesp = 1;
    idg = ihesp + *maxp * 6;
    iqe = idg + *maxp;
    ixypw = iqe + *maxe;
    irse = ixypw + *np * 3;
    mrlintrp = *maxwr - irse;
    if (mrlintrp < *maxe) {
	*ierr = 1002;
	s_wsfi(&io___36);
	do_fio(&c__1, "The approximate size of rW is ", (ftnlen)30);
	i__1 = irse + *maxe;
	do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	e_wsfi();
	errmes_(ierr, "mbaNodal", message, (ftnlen)8, (ftnlen)80);
    }
    i__1 = *maxwr;
    for (n = 1; n <= i__1; ++n) {
	rw[n] = 0.;
    }
    i__1 = *maxwi;
    for (n = 1; n <= i__1; ++n) {
	iw[n] = 0;
    }
/* ... auxiliary structure of the data */
    m = iipf - 1;
    i__1 = *nf;
    for (n = 1; n <= i__1; ++n) {
	for (i__ = 1; i__ <= 3; ++i__) {
	    iw[m + i__] = ipf[i__ + n * 3];
	}
	iw[m + 4] = lbf[n];
	m += 4;
    }
    m = iipe - 1;
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	for (i__ = 1; i__ <= 4; ++i__) {
	    iw[m + i__] = ipe[i__ + (n << 2)];
	}
	iw[m + 5] = lbe[n];
	m += 5;
    }
/* ... scale geometry to unit cube */
    scale2cube_(np, &xyp[4], &c_true);
/* ... setting the fixed metric for the future loops */
    npw = *np;
    new__ = *ne;
    copymeshdata_(np, ne, &xyp[4], &metric[7], &ipe[5], &rw[ixypw], &rw[ihesp]
	    , &iw[iipew]);
/* ... quality of tetrahedra is computed */
    makq_(&nloop, np, ne, &xyp[4], &ipe[5], nev, &iev[1], nestar, &hstar, &rw[
	    ihesp], &rw[idg], &rw[iqe]);
/* group (M) */
/* group (Q) */
/* ... runing the basic algorithm for the global grid */
    flaganalytic = FALSE_;
    ani2_(np, maxp, nf, maxf, ne, maxe, &xyp[4], &iw[iipf], &iw[iipe], nestar,
	     &hstar, &iw[iicp], &iw[iipp], &iw[iiep], &iw[iife], &iw[iiee], &
	    iw[il1e], &iw[il2e], &iw[iiholp], &iw[iiholf], &iw[iihole], &iw[
	    iiepw], &iw[inepw], &milintrp, &mrlintrp, &iw[iipew], &iw[iise], 
	    flagauto, &statuswork, npv, nfv, nev, &ipv[1], &ifv[1], &iev[1], 
	    maxskipe, maxqitr, &c__60, &rw[ihesp], quality, rquality, &rw[idg]
	    , &rw[iqe], &npw, &new__, &rw[ixypw], &metric[7], &rw[irse], (
	    I_fp)ani_metric_eucl__, &flaganalytic, &flagfile, &chanelout, &
	    nlines, output, &iprintl, &nqitr, ierr, (ftnlen)100);
/* group (M) */
/* group (Dev) */
/* group (Q) */
/* group (ERR) */
    if (*ierr != 0 && *ierr != 1000) {
	errmes_(ierr, "mbaNodal", "memory problems", (ftnlen)8, (ftnlen)15);
    }
    seconds_(&tm1);
    tm1 -= tm0;
    if (iprintl >= 1) {
	if (flagfile) {
	    i__1 = nlines;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		io___47.ciunit = chanelout;
		s_wsfe(&io___47);
		do_fio(&c__1, output + (i__ - 1) * 100, (ftnlen)100);
		e_wsfe();
	    }
	    io___48.ciunit = chanelout;
	    s_wsfe(&io___48);
	    do_fio(&c__1, (char *)&nqitr, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*rquality), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*np), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*nf), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*ne), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&tm1, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else {
	    s_wsfe(&io___49);
	    do_fio(&c__1, (char *)&nqitr, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*rquality), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*np), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*nf), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*ne), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&tm1, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }
/* ... close output file */
    if (flagfile) {
	cl__1.cerr = 0;
	cl__1.cunit = chanelout;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }
/* ... rescale geometry back */
    scale2cube_(np, &xyp[4], &c_false);
/* ... original structure of the data */
    m = iipf - 1;
    i__1 = *nf;
    for (n = 1; n <= i__1; ++n) {
	for (i__ = 1; i__ <= 3; ++i__) {
	    ipf[i__ + n * 3] = iw[m + i__];
	}
	lbf[n] = iw[m + 4];
	m += 4;
    }
    m = iipe - 1;
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	for (i__ = 1; i__ <= 4; ++i__) {
	    ipe[i__ + (n << 2)] = iw[m + i__];
	}
	lbe[n] = iw[m + 5];
	m += 5;
    }
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	iw[n] = iw[iicp + n - 1];
    }
    iw[*np + 1] = nqitr;
    rw[1] = tm1;
    rw[2] = *rquality;
    rw[3] = hstar;
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	rw[n + 3] = rw[iqe + n - 1];
    }
    return 0;
} /* mbanodal_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int mbanodalshort_(integer *np, integer *maxp, 
	integer *nf, integer *maxf, integer *ne, integer *maxe, doublereal *
	xyp, integer *ipf, integer *ipe, integer *lbf, integer *nestar, 
	integer *status, doublereal *metric, integer *maxwr, integer *maxwi, 
	doublereal *rw, integer *iw, integer *iprint, integer *ierr)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static logical flagauto;
    static integer maxskipe;
    static doublereal rquality;
    static integer n, nev, nfv, npv, ilbe, iiev, iifv, iipv, iiend;
    static char message[80];
    static integer maxwiav, maxwrav, maxqitr;
    static doublereal quality;

    /* Fortran I/O blocks */
    static icilist io___53 = { 0, message, 0, "(A,I10)", 80, 1 };
    static icilist io___54 = { 0, message, 0, "(A,I10)", 80, 1 };


/* ================================================================ */
/* group (M) */
/* group (W) */
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
/* MaxLines   - the maximal number of output lines per processor */
/* LineLenght - the maximal lenght of a message */
/* ========================================================== */
/* group (Q) */
/* group (W) */
/* group (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --iw;
    --rw;
    metric -= 7;
    --lbf;
    ipe -= 5;
    ipf -= 4;
    xyp -= 4;

    /* Function Body */
    *ierr = 0;
    maxwiav = *maxp * 7 + *np + *maxf * 9 + *maxe * 19 + *ne * 13;
    maxwrav = *maxp * 14 + *np * 14 + *maxe;
    if (*maxwi <= maxwiav) {
	*ierr = 1001;
	s_wsfi(&io___53);
	do_fio(&c__1, "The approximate size of iW is ", (ftnlen)30);
	do_fio(&c__1, (char *)&maxwiav, (ftnlen)sizeof(integer));
	e_wsfi();
	errmes_(ierr, "mbaNodalShort", message, (ftnlen)13, (ftnlen)80);
    }
    if (*maxwr <= maxwrav) {
	*ierr = 1002;
	s_wsfi(&io___54);
	do_fio(&c__1, "The approximate size of rW is ", (ftnlen)30);
	do_fio(&c__1, (char *)&maxwrav, (ftnlen)sizeof(integer));
	e_wsfi();
	errmes_(ierr, "mbaNodalShort", message, (ftnlen)13, (ftnlen)80);
    }
/* ... setup of missing parameter */
    npv = 0;
    nfv = 0;
    nev = 0;
    flagauto = TRUE_;
    maxskipe = 200;
    maxqitr = min(20000,*ne);
    maxqitr = max(1000,maxqitr);
    quality = .3;
    iipv = 1;
    ilbe = iipv + *maxf;
    iifv = ilbe + *maxe;
    iiev = iifv + nfv + 1;
    iiend = iiev + nev + 1;
    i__1 = *ne - 1;
    for (n = 0; n <= i__1; ++n) {
	iw[ilbe + n] = 1;
    }
    i__1 = *maxwi - iiend;
    mbanodal_(np, maxp, nf, maxf, ne, maxe, &xyp[4], &ipf[4], &ipe[5], &lbf[1]
	    , &iw[ilbe], nestar, &npv, &nfv, &nev, &iw[iipv], &iw[iifv], &iw[
	    iiev], &flagauto, status, &maxskipe, &maxqitr, &metric[7], &
	    quality, &rquality, maxwr, &i__1, &rw[1], &iw[iiend], iprint, 
	    ierr);
/* group (M) */
/* group (Dev) */
/* group (Q) */
/* group (W) */
    return 0;
} /* mbanodalshort_ */

