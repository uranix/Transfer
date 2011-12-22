/* saveM.f -- translated by f2c (version 20090411).
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

#include "makQ.h"
#include "makM.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c__3 = 3;
static integer c__9 = 9;
static integer c__5 = 5;
static integer c__6004 = 6004;
static integer c__4 = 4;
static integer c__2 = 2;
static integer c__4001 = 4001;

/* ================================================================ */
/* Saves mesh in ANI format */
/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int savemani_(integer *np, integer *nf, integer *
	ne, doublereal *xyp, integer *ipf, integer *ipe, integer *lbf, 
	integer *lbe, integer *npv, integer *nfv, integer *nev, integer *ipv, 
	integer *ifv, integer *iev, logical *flagi, integer *ipp, integer *
	iff, char *fname, ftnlen fname_len)
{
    /* System generated locals */
    integer i__1;
    olist o__1;
    cllist cl__1;

    /* Local variables */
    static integer i__, j, n, kfe, kff, kcs, kfp, kes, kfs, kis, kps;

    /* Fortran I/O blocks */
    static cilist io___2 = { 0, 6, 0, "(A,A)", 0 };
    static cilist io___11 = { 0, 10, 0, "(3(A,I8),A)", 0 };
    static cilist io___12 = { 0, 10, 0, "(3(A,I8),A)", 0 };
    static cilist io___13 = { 0, 10, 0, "(3(A,I8),A)", 0 };
    static cilist io___14 = { 0, 10, 0, "(A)", 0 };
    static cilist io___15 = { 0, 10, 0, "(3(A,I8),A)", 0 };
    static cilist io___16 = { 0, 10, 0, "(A,I8)", 0 };
    static cilist io___17 = { 0, 10, 0, "(3(A,I8),A)", 0 };
    static cilist io___18 = { 0, 10, 0, "(A,I8)", 0 };
    static cilist io___19 = { 0, 10, 0, "(3(A,I8),A)", 0 };
    static cilist io___20 = { 0, 10, 0, "(A,I8)", 0 };
    static cilist io___21 = { 0, 10, 0, "(L1,A)", 0 };
    static cilist io___22 = { 0, 10, 0, 0, 0 };
    static cilist io___23 = { 0, 10, 0, 0, 0 };
    static cilist io___25 = { 0, 10, 0, "(3E23.15)", 0 };
    static cilist io___27 = { 0, 10, 0, 0, 0 };
    static cilist io___28 = { 0, 10, 0, 0, 0 };
    static cilist io___29 = { 0, 10, 0, 0, 0 };
    static cilist io___30 = { 0, 10, 0, 0, 0 };
    static cilist io___31 = { 0, 10, 0, 0, 0 };
    static cilist io___32 = { 0, 10, 0, 0, 0 };
    static cilist io___33 = { 0, 10, 0, 0, 0 };
    static cilist io___34 = { 0, 10, 0, 0, 0 };
    static cilist io___35 = { 0, 10, 0, 0, 0 };
    static cilist io___36 = { 0, 10, 0, 0, 0 };
    static cilist io___37 = { 0, 10, 0, 0, 0 };
    static cilist io___38 = { 0, 10, 0, 0, 0 };
    static cilist io___39 = { 0, 10, 0, 0, 0 };
    static cilist io___40 = { 0, 10, 0, 0, 0 };
    static cilist io___41 = { 0, 10, 0, 0, 0 };
    static cilist io___42 = { 0, 10, 0, 0, 0 };
    static cilist io___43 = { 0, 10, 0, 0, 0 };
    static cilist io___44 = { 0, 10, 0, 0, 0 };
    static cilist io___45 = { 0, 10, 0, 0, 0 };
    static cilist io___46 = { 0, 10, 0, 0, 0 };
    static cilist io___47 = { 0, 10, 0, 0, 0 };


/* ================================================================ */
/* group (M) */
/* group (Dev) */
/* group (I) */
/* ================================================================ */
/* group (M) */
/*     Integer MaxP, MaxF, MaxE */
/* group (Dev) */
/* group (I) */
/* ================================================================ */
    /* Parameter adjustments */
    --iff;
    --ipp;
    --iev;
    --ifv;
    --ipv;
    --lbe;
    --lbf;
    ipe -= 5;
    ipf -= 4;
    xyp -= 4;

    /* Function Body */
    i__ = 1;
    while(s_cmp(fname + (i__ - 1), ".ani", (ftnlen)4, (ftnlen)4) != 0) {
	++i__;
    }
    s_wsfe(&io___2);
    do_fio(&c__1, "Saving mesh ", (ftnlen)12);
    do_fio(&c__1, fname, i__ + 3);
    e_wsfe();
    kps = 11;
    kfs = kps + *np + 2;
    kes = kfs + *nf + 2;
    kcs = kes + *ne + 2;
    kfp = kcs + 2;
    kff = kfp + *npv + 2;
    kfe = kff + *nfv + 2;
    kis = kfe + *nev + 2;
    o__1.oerr = 0;
    o__1.ounit = 10;
    o__1.ofnmlen = i__ + 3;
    o__1.ofnm = fname;
    o__1.orl = 0;
    o__1.osta = "UNKNOWN";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_wsfe(&io___11);
    do_fio(&c__1, "T points:        ", (ftnlen)17);
    do_fio(&c__1, (char *)&(*np), (ftnlen)sizeof(integer));
    do_fio(&c__1, " (lines ", (ftnlen)8);
    do_fio(&c__1, (char *)&kps, (ftnlen)sizeof(integer));
    do_fio(&c__1, " - ", (ftnlen)3);
    i__1 = kps + *np - 1;
    do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
    do_fio(&c__1, ")", (ftnlen)1);
    e_wsfe();
    s_wsfe(&io___12);
    do_fio(&c__1, "T faces:         ", (ftnlen)17);
    do_fio(&c__1, (char *)&(*nf), (ftnlen)sizeof(integer));
    do_fio(&c__1, " (lines ", (ftnlen)8);
    do_fio(&c__1, (char *)&kfs, (ftnlen)sizeof(integer));
    do_fio(&c__1, " - ", (ftnlen)3);
    i__1 = kfs + *nf - 1;
    do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
    do_fio(&c__1, ")", (ftnlen)1);
    e_wsfe();
    s_wsfe(&io___13);
    do_fio(&c__1, "T elements:      ", (ftnlen)17);
    do_fio(&c__1, (char *)&(*ne), (ftnlen)sizeof(integer));
    do_fio(&c__1, " (lines ", (ftnlen)8);
    do_fio(&c__1, (char *)&kes, (ftnlen)sizeof(integer));
    do_fio(&c__1, " - ", (ftnlen)3);
    i__1 = kes + *ne - 1;
    do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
    do_fio(&c__1, ")", (ftnlen)1);
    e_wsfe();
    s_wsfe(&io___14);
    do_fio(&c__1, "T curved faces:     0", (ftnlen)21);
    e_wsfe();
    if (*npv != 0) {
	s_wsfe(&io___15);
	do_fio(&c__1, "T fixed points:  ", (ftnlen)17);
	do_fio(&c__1, (char *)&(*npv), (ftnlen)sizeof(integer));
	do_fio(&c__1, " (lines ", (ftnlen)8);
	do_fio(&c__1, (char *)&kfp, (ftnlen)sizeof(integer));
	do_fio(&c__1, " - ", (ftnlen)3);
	i__1 = kfp + *npv - 1;
	do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	do_fio(&c__1, ")", (ftnlen)1);
	e_wsfe();
    } else {
	s_wsfe(&io___16);
	do_fio(&c__1, "T fixed points:  ", (ftnlen)17);
	do_fio(&c__1, (char *)&(*npv), (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (*nfv != 0) {
	s_wsfe(&io___17);
	do_fio(&c__1, "T fixed faces:   ", (ftnlen)17);
	do_fio(&c__1, (char *)&(*nfv), (ftnlen)sizeof(integer));
	do_fio(&c__1, " (lines ", (ftnlen)8);
	do_fio(&c__1, (char *)&kff, (ftnlen)sizeof(integer));
	do_fio(&c__1, " - ", (ftnlen)3);
	i__1 = kff + *nfv - 1;
	do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	do_fio(&c__1, ")", (ftnlen)1);
	e_wsfe();
    } else {
	s_wsfe(&io___18);
	do_fio(&c__1, "T fixed faces:   ", (ftnlen)17);
	do_fio(&c__1, (char *)&(*nfv), (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (*nev != 0) {
	s_wsfe(&io___19);
	do_fio(&c__1, "T fixed elements:", (ftnlen)17);
	do_fio(&c__1, (char *)&(*nev), (ftnlen)sizeof(integer));
	do_fio(&c__1, " (lines ", (ftnlen)8);
	do_fio(&c__1, (char *)&kfe, (ftnlen)sizeof(integer));
	do_fio(&c__1, " - ", (ftnlen)3);
	i__1 = kfe + *nev - 1;
	do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	do_fio(&c__1, ")", (ftnlen)1);
	e_wsfe();
    } else {
	s_wsfe(&io___20);
	do_fio(&c__1, "T fixed elements:", (ftnlen)17);
	do_fio(&c__1, (char *)&(*nev), (ftnlen)sizeof(integer));
	e_wsfe();
    }
    s_wsfe(&io___21);
    do_fio(&c__1, (char *)&(*flagi), (ftnlen)sizeof(logical));
    do_fio(&c__1, " interface data:", (ftnlen)16);
    e_wsfe();
    s_wsle(&io___22);
    e_wsle();
    s_wsle(&io___23);
    do_lio(&c__3, &c__1, (char *)&(*np), (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, " # of nodes", (ftnlen)11);
    e_wsle();
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	s_wsfe(&io___25);
	for (j = 1; j <= 3; ++j) {
	    do_fio(&c__1, (char *)&xyp[j + n * 3], (ftnlen)sizeof(doublereal))
		    ;
	}
	e_wsfe();
    }
    s_wsle(&io___27);
    e_wsle();
    s_wsle(&io___28);
    do_lio(&c__3, &c__1, (char *)&(*nf), (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, " # of faces", (ftnlen)11);
    e_wsle();
    i__1 = *nf;
    for (n = 1; n <= i__1; ++n) {
	s_wsle(&io___29);
	for (j = 1; j <= 3; ++j) {
	    do_lio(&c__3, &c__1, (char *)&ipf[j + n * 3], (ftnlen)sizeof(
		    integer));
	}
	do_lio(&c__3, &c__1, (char *)&lbf[n], (ftnlen)sizeof(integer));
	e_wsle();
    }
    s_wsle(&io___30);
    e_wsle();
    s_wsle(&io___31);
    do_lio(&c__3, &c__1, (char *)&(*ne), (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, " # of elements", (ftnlen)14);
    e_wsle();
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	s_wsle(&io___32);
	for (j = 1; j <= 4; ++j) {
	    do_lio(&c__3, &c__1, (char *)&ipe[j + (n << 2)], (ftnlen)sizeof(
		    integer));
	}
	do_lio(&c__3, &c__1, (char *)&lbe[n], (ftnlen)sizeof(integer));
	e_wsle();
    }
/* ... saving curvilinear faces */
    s_wsle(&io___33);
    e_wsle();
    s_wsle(&io___34);
    do_lio(&c__9, &c__1, "0 # of curvilinear faces", (ftnlen)24);
    e_wsle();
/* ... saving the fixed mesh points, faces and elements */
    s_wsle(&io___35);
    e_wsle();
    s_wsle(&io___36);
    do_lio(&c__3, &c__1, (char *)&(*npv), (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, " # number of fixed points", (ftnlen)25);
    e_wsle();
    i__1 = *npv;
    for (n = 1; n <= i__1; ++n) {
	s_wsle(&io___37);
	do_lio(&c__3, &c__1, (char *)&ipv[n], (ftnlen)sizeof(integer));
	e_wsle();
    }
    s_wsle(&io___38);
    e_wsle();
    s_wsle(&io___39);
    do_lio(&c__3, &c__1, (char *)&(*nfv), (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, " # number of fixed faces", (ftnlen)24);
    e_wsle();
    i__1 = *nfv;
    for (n = 1; n <= i__1; ++n) {
	s_wsle(&io___40);
	do_lio(&c__3, &c__1, (char *)&ifv[n], (ftnlen)sizeof(integer));
	e_wsle();
    }
    s_wsle(&io___41);
    e_wsle();
    s_wsle(&io___42);
    do_lio(&c__3, &c__1, (char *)&(*nev), (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, " # number of fixed elements", (ftnlen)27);
    e_wsle();
    i__1 = *nev;
    for (n = 1; n <= i__1; ++n) {
	s_wsle(&io___43);
	do_lio(&c__3, &c__1, (char *)&iev[n], (ftnlen)sizeof(integer));
	e_wsle();
    }
    if (*flagi) {
	s_wsle(&io___44);
	e_wsle();
	s_wsle(&io___45);
	do_lio(&c__3, &c__1, (char *)&(*np), (ftnlen)sizeof(integer));
	do_lio(&c__3, &c__1, (char *)&(*nf), (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, " interface data", (ftnlen)15);
	e_wsle();
	i__1 = *np;
	for (n = 1; n <= i__1; ++n) {
	    s_wsle(&io___46);
	    do_lio(&c__3, &c__1, (char *)&ipp[n], (ftnlen)sizeof(integer));
	    e_wsle();
	}
	i__1 = *nf;
	for (n = 1; n <= i__1; ++n) {
	    s_wsle(&io___47);
	    do_lio(&c__3, &c__1, (char *)&iff[n], (ftnlen)sizeof(integer));
	    e_wsle();
	}
    }
    cl__1.cerr = 0;
    cl__1.cunit = 10;
    cl__1.csta = 0;
    f_clos(&cl__1);
/* L5000: */
    return 0;
} /* savemani_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int saves_(integer *np, doublereal *sol, char *
	fname, ftnlen fname_len)
{
    /* System generated locals */
    integer i__1;
    olist o__1;
    cllist cl__1;

    /* Local variables */
    static integer n;

    /* Fortran I/O blocks */
    static cilist io___48 = { 0, 6, 0, "(A,A)", 0 };
    static cilist io___49 = { 0, 10, 0, 0, 0 };
    static cilist io___50 = { 0, 10, 0, 0, 0 };
    static cilist io___52 = { 0, 10, 0, 0, 0 };


/* ================================================================ */
/* ================================================================ */
    /* Parameter adjustments */
    --sol;

    /* Function Body */
    s_wsfe(&io___48);
    do_fio(&c__1, "Saving solution ", (ftnlen)16);
    do_fio(&c__1, fname, fname_len);
    e_wsfe();
    o__1.oerr = 0;
    o__1.ounit = 10;
    o__1.ofnmlen = fname_len;
    o__1.ofnm = fname;
    o__1.orl = 0;
    o__1.osta = "UNKNOWN";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_wsle(&io___49);
    do_lio(&c__3, &c__1, (char *)&(*np), (ftnlen)sizeof(integer));
    e_wsle();
    s_wsle(&io___50);
    e_wsle();
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	s_wsle(&io___52);
	do_lio(&c__5, &c__1, (char *)&sol[n], (ftnlen)sizeof(doublereal));
	e_wsle();
    }
    cl__1.cerr = 0;
    cl__1.cunit = 10;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
} /* saves_ */

/* ===================================================================== */
/* @f2h@ */ /* Subroutine */ int savemgmv_(integer *np, integer *nf, integer *
	ne, doublereal *xyp, integer *ipf, integer *ipe, integer *lbf, 
	integer *lbe, char *fname, integer *iw, ftnlen fname_len)
{
    /* System generated locals */
    integer i__1;
    olist o__1;
    cllist cl__1;

    /* Local variables */
    static integer i__, k, n;
    static doublereal v;
    static integer ic, ip1, ip2, ip3, ip4;

    /* Fortran I/O blocks */
    static cilist io___54 = { 0, 6, 0, "(A,A)", 0 };
    static cilist io___55 = { 0, 10, 0, "(A)", 0 };
    static cilist io___56 = { 0, 10, 0, 0, 0 };
    static cilist io___57 = { 0, 10, 0, 0, 0 };
    static cilist io___59 = { 0, 10, 0, 0, 0 };
    static cilist io___60 = { 0, 10, 0, 0, 0 };
    static cilist io___61 = { 0, 10, 0, 0, 0 };
    static cilist io___67 = { 0, 10, 0, 0, 0 };
    static cilist io___68 = { 0, 10, 0, 0, 0 };
    static cilist io___69 = { 0, 10, 0, 0, 0 };
    static cilist io___70 = { 0, 10, 0, 0, 0 };
    static cilist io___72 = { 0, 10, 0, 0, 0 };
    static cilist io___73 = { 0, 10, 0, "(A,I1)", 0 };
    static cilist io___74 = { 0, 10, 0, "(A,I2)", 0 };
    static cilist io___75 = { 0, 10, 0, "(A,I3)", 0 };
    static cilist io___76 = { 0, 10, 0, "(A,I4)", 0 };
    static cilist io___77 = { 0, 10, 0, 0, 0 };
    static cilist io___78 = { 0, 10, 0, "(A)", 0 };
    static cilist io___79 = { 0, 10, 0, 0, 0 };
    static cilist io___81 = { 0, 10, 0, 0, 0 };
    static cilist io___82 = { 0, 10, 0, "(A)", 0 };
    static cilist io___83 = { 0, 10, 0, 0, 0 };
    static cilist io___84 = { 0, 10, 0, "(A)", 0 };


/* ===================================================================== */
/* ===================================================================== */
/* Routine saves mesh in the GMV file. (New version) */

/* *** Remarks: */
/*        1. The size of the working memory is nE */
/* ===================================================================== */
/* group (M) */
/* group (Local variables) */
/* ===================================================================== */
    /* Parameter adjustments */
    --iw;
    --lbe;
    --lbf;
    ipe -= 5;
    ipf -= 4;
    xyp -= 4;

    /* Function Body */
    i__ = 1;
    while(s_cmp(fname + (i__ - 1), ".gmv", (ftnlen)4, (ftnlen)4) != 0) {
	++i__;
    }
    s_wsfe(&io___54);
    do_fio(&c__1, "Saving GMV image ", (ftnlen)17);
    do_fio(&c__1, fname, i__ + 3);
    e_wsfe();
    o__1.oerr = 0;
    o__1.ounit = 10;
    o__1.ofnmlen = i__ + 3;
    o__1.ofnm = fname;
    o__1.orl = 0;
    o__1.osta = "UNKNOWN";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_wsfe(&io___55);
    do_fio(&c__1, "gmvinput ascii", (ftnlen)14);
    e_wsfe();
    s_wsle(&io___56);
    e_wsle();
    s_wsle(&io___57);
    do_lio(&c__9, &c__1, "nodev ", (ftnlen)6);
    do_lio(&c__3, &c__1, (char *)&(*np), (ftnlen)sizeof(integer));
    e_wsle();
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	s_wsle(&io___59);
	for (i__ = 1; i__ <= 3; ++i__) {
	    do_lio(&c__5, &c__1, (char *)&xyp[i__ + n * 3], (ftnlen)sizeof(
		    doublereal));
	}
	e_wsle();
    }
/* ... save cells */
    s_wsle(&io___60);
    e_wsle();
    s_wsle(&io___61);
    do_lio(&c__9, &c__1, "cells ", (ftnlen)6);
    do_lio(&c__3, &c__1, (char *)&(*ne), (ftnlen)sizeof(integer));
    e_wsle();
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	ip1 = ipe[(n << 2) + 1];
	ip2 = ipe[(n << 2) + 2];
	ip3 = ipe[(n << 2) + 3];
	ip4 = ipe[(n << 2) + 4];
	v = calvol_(&xyp[ip1 * 3 + 1], &xyp[ip2 * 3 + 1], &xyp[ip3 * 3 + 1], &
		xyp[ip4 * 3 + 1]);
	if (v > 0.) {
	    s_wsle(&io___67);
	    do_lio(&c__9, &c__1, " tet 4 ", (ftnlen)7);
	    do_lio(&c__3, &c__1, (char *)&ip1, (ftnlen)sizeof(integer));
	    do_lio(&c__3, &c__1, (char *)&ip2, (ftnlen)sizeof(integer));
	    do_lio(&c__3, &c__1, (char *)&ip3, (ftnlen)sizeof(integer));
	    do_lio(&c__3, &c__1, (char *)&ip4, (ftnlen)sizeof(integer));
	    e_wsle();
	} else {
	    s_wsle(&io___68);
	    do_lio(&c__9, &c__1, " tet 4 ", (ftnlen)7);
	    do_lio(&c__3, &c__1, (char *)&ip1, (ftnlen)sizeof(integer));
	    do_lio(&c__3, &c__1, (char *)&ip3, (ftnlen)sizeof(integer));
	    do_lio(&c__3, &c__1, (char *)&ip2, (ftnlen)sizeof(integer));
	    do_lio(&c__3, &c__1, (char *)&ip4, (ftnlen)sizeof(integer));
	    e_wsle();
	}
    }
    s_wsle(&io___69);
    e_wsle();
    s_wsle(&io___70);
    e_wsle();
/* ... save materials */
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	iw[n] = lbe[n];
    }
    ic = countcolors_(ne, &iw[1]);
    s_wsle(&io___72);
    do_lio(&c__9, &c__1, "material ", (ftnlen)9);
    do_lio(&c__3, &c__1, (char *)&ic, (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, " 0", (ftnlen)2);
    e_wsle();
    i__1 = ic;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (iw[i__] < 10) {
	    s_wsfe(&io___73);
	    do_fio(&c__1, "mat", (ftnlen)3);
	    do_fio(&c__1, (char *)&iw[i__], (ftnlen)sizeof(integer));
	    e_wsfe();
	} else if (iw[i__] < 100) {
	    s_wsfe(&io___74);
	    do_fio(&c__1, "mat", (ftnlen)3);
	    do_fio(&c__1, (char *)&iw[i__], (ftnlen)sizeof(integer));
	    e_wsfe();
	} else if (iw[i__] < 1000) {
	    s_wsfe(&io___75);
	    do_fio(&c__1, "mat", (ftnlen)3);
	    do_fio(&c__1, (char *)&iw[i__], (ftnlen)sizeof(integer));
	    e_wsfe();
	} else if (iw[i__] < 10000) {
	    s_wsfe(&io___76);
	    do_fio(&c__1, "mat", (ftnlen)3);
	    do_fio(&c__1, (char *)&iw[i__], (ftnlen)sizeof(integer));
	    e_wsfe();
	} else {
	    errmes_(&c__6004, "saveMgmv", "missing code", (ftnlen)8, (ftnlen)
		    12);
	}
    }
    s_wsle(&io___77);
    i__1 = *ne;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_lio(&c__3, &c__1, (char *)&lbe[i__], (ftnlen)sizeof(integer));
    }
    e_wsle();
/* ... save faces */
    s_wsfe(&io___78);
    do_fio(&c__1, "polygons", (ftnlen)8);
    e_wsfe();
    i__1 = *nf;
    for (n = 1; n <= i__1; ++n) {
	s_wsle(&io___79);
	do_lio(&c__3, &c__1, (char *)&lbf[n], (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, " 3 ", (ftnlen)3);
	e_wsle();
	for (k = 1; k <= 3; ++k) {
	    s_wsle(&io___81);
	    for (i__ = 1; i__ <= 3; ++i__) {
		do_lio(&c__5, &c__1, (char *)&xyp[k + ipf[i__ + n * 3] * 3], (
			ftnlen)sizeof(doublereal));
	    }
	    e_wsle();
	}
    }
    s_wsfe(&io___82);
    do_fio(&c__1, "endpoly", (ftnlen)7);
    e_wsfe();
    s_wsle(&io___83);
    e_wsle();
    s_wsfe(&io___84);
    do_fio(&c__1, "endgmv", (ftnlen)6);
    e_wsfe();
    cl__1.cerr = 0;
    cl__1.cunit = 10;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
} /* savemgmv_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int savemgeo_(integer *np, integer *maxp, 
	integer *nf, integer *maxf, integer *ne, integer *maxe, doublereal *
	xyp, integer *ipf, integer *ipe, integer *lbp, integer *lbf, integer *
	lbe, integer *npv, integer *nfv, integer *nev, integer *ipv, integer *
	ifv, integer *iev, char *fname, ftnlen fname_len)
{
    /* System generated locals */
    address a__1[2];
    integer i__1, i__2[2];
    olist o__1;
    cllist cl__1;

    /* Local variables */
    static char fnameext[70];
    static integer i__, j, k, n, nb;
    static logical flag__;

    /* Fortran I/O blocks */
    static cilist io___90 = { 0, 6, 0, "(A,A)", 0 };
    static cilist io___91 = { 0, 10, 0, 0, 0 };
    static cilist io___92 = { 0, 10, 0, 0, 0 };
    static cilist io___93 = { 0, 10, 0, 0, 0 };
    static cilist io___95 = { 0, 10, 0, 0, 0 };
    static cilist io___96 = { 0, 10, 0, 0, 0 };
    static cilist io___97 = { 0, 10, 0, 0, 0 };
    static cilist io___99 = { 0, 10, 0, 0, 0 };


/* ================================================================ */
/* group (M) */
/* group (Dev) */
/* ================================================================ */
/* ================================================================ */
/* group (M) */
/*     Integer MaxP, MaxF, MaxE */
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
/* group (Dev) */
/* group (Local variables) */
/* ================================================================ */
/* ... counting the boundary faces (array lbF is overloaded) */
    /* Parameter adjustments */
    --iev;
    --ifv;
    --ipv;
    --lbe;
    --lbf;
    --lbp;
    ipe -= 5;
    ipf -= 4;
    xyp -= 4;

    /* Function Body */
    nb = 0;
    i__1 = *nf;
    for (n = 1; n <= i__1; ++n) {
	flag__ = TRUE_;
	for (i__ = 1; i__ <= 3; ++i__) {
	    if (! ifxnode_(&lbp[ipf[i__ + n * 3]], &c__4)) {
		flag__ = FALSE_;
	    }
	}
	if (flag__) {
	    ++nb;
	    lbf[n] = -lbf[n];
	}
    }
/* Writing concatenation */
    i__2[0] = fname_len, a__1[0] = fname;
    i__2[1] = 4, a__1[1] = ".geo";
    s_cat(fnameext, a__1, i__2, &c__2, (ftnlen)70);
    s_wsfe(&io___90);
    do_fio(&c__1, "Saving mesh ", (ftnlen)12);
    do_fio(&c__1, fnameext, (ftnlen)70);
    e_wsfe();
    o__1.oerr = 0;
    o__1.ounit = 10;
    o__1.ofnmlen = 70;
    o__1.ofnm = fnameext;
    o__1.orl = 0;
    o__1.osta = "UNKNOWN";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_wsle(&io___91);
    do_lio(&c__3, &c__1, (char *)&(*np), (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&(*ne), (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&nb, (ftnlen)sizeof(integer));
    e_wsle();
    s_wsle(&io___92);
    e_wsle();
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	s_wsle(&io___93);
	do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(integer));
	for (j = 1; j <= 3; ++j) {
	    do_lio(&c__5, &c__1, (char *)&xyp[j + n * 3], (ftnlen)sizeof(
		    doublereal));
	}
	e_wsle();
    }
    s_wsle(&io___95);
    e_wsle();
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	s_wsle(&io___96);
	do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(integer));
	for (j = 1; j <= 4; ++j) {
	    do_lio(&c__3, &c__1, (char *)&ipe[j + (n << 2)], (ftnlen)sizeof(
		    integer));
	}
	do_lio(&c__3, &c__1, (char *)&lbe[n], (ftnlen)sizeof(integer));
	e_wsle();
    }
    s_wsle(&io___97);
    e_wsle();
    k = 0;
    i__1 = *nf;
    for (n = 1; n <= i__1; ++n) {
	if (lbf[n] < 0) {
	    lbf[n] = -lbf[n];
	    ++k;
	    s_wsle(&io___99);
	    do_lio(&c__3, &c__1, (char *)&k, (ftnlen)sizeof(integer));
	    for (j = 1; j <= 3; ++j) {
		do_lio(&c__3, &c__1, (char *)&ipf[j + n * 3], (ftnlen)sizeof(
			integer));
	    }
	    do_lio(&c__3, &c__1, (char *)&lbf[n], (ftnlen)sizeof(integer));
	    e_wsle();
	}
    }
    cl__1.cerr = 0;
    cl__1.cunit = 10;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
} /* savemgeo_ */

/* ================================================================ */
/* Saves mesh in AFT format */
/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int savemaft_(integer *np, integer *nf, integer *
	ne, doublereal *xyp, integer *ipf, integer *ipe, integer *lbf, 
	integer *lbe, char *fname, ftnlen fname_len)
{
    /* System generated locals */
    integer i__1;
    olist o__1;
    cllist cl__1;

    /* Local variables */
    static integer i__, n;

    /* Fortran I/O blocks */
    static cilist io___101 = { 0, 6, 0, "(A,A)", 0 };
    static cilist io___102 = { 0, 10, 0, 0, 0 };
    static cilist io___104 = { 0, 10, 0, 0, 0 };
    static cilist io___105 = { 0, 10, 0, 0, 0 };
    static cilist io___106 = { 0, 10, 0, 0, 0 };
    static cilist io___107 = { 0, 10, 0, 0, 0 };
    static cilist io___108 = { 0, 10, 0, 0, 0 };


/* ================================================================ */
/* ================================================================ */
/* ================================================================ */
    /* Parameter adjustments */
    --lbe;
    --lbf;
    ipe -= 5;
    ipf -= 4;
    xyp -= 4;

    /* Function Body */
    i__ = 1;
    while(s_cmp(fname + (i__ - 1), ".out", (ftnlen)4, (ftnlen)4) != 0) {
	++i__;
    }
    s_wsfe(&io___101);
    do_fio(&c__1, "Saving aft-mesh ", (ftnlen)16);
    do_fio(&c__1, fname, i__ + 3);
    e_wsfe();
    o__1.oerr = 1;
    o__1.ounit = 10;
    o__1.ofnmlen = i__ + 3;
    o__1.ofnm = fname;
    o__1.orl = 0;
    o__1.osta = "UNKNOWN";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__1 = f_open(&o__1);
    if (i__1 != 0) {
	goto L1000;
    }
/* ... write points */
    s_wsle(&io___102);
    do_lio(&c__3, &c__1, (char *)&(*np), (ftnlen)sizeof(integer));
    e_wsle();
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	s_wsle(&io___104);
	for (i__ = 1; i__ <= 3; ++i__) {
	    do_lio(&c__5, &c__1, (char *)&xyp[i__ + n * 3], (ftnlen)sizeof(
		    doublereal));
	}
	e_wsle();
    }
/* ... write elements */
    s_wsle(&io___105);
    do_lio(&c__3, &c__1, (char *)&(*ne), (ftnlen)sizeof(integer));
    e_wsle();
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	s_wsle(&io___106);
	for (i__ = 1; i__ <= 4; ++i__) {
	    do_lio(&c__3, &c__1, (char *)&ipe[i__ + (n << 2)], (ftnlen)sizeof(
		    integer));
	}
	do_lio(&c__3, &c__1, (char *)&lbe[n], (ftnlen)sizeof(integer));
	e_wsle();
    }
/* ... write faces */
    s_wsle(&io___107);
    do_lio(&c__3, &c__1, (char *)&(*nf), (ftnlen)sizeof(integer));
    e_wsle();
    i__1 = *nf;
    for (n = 1; n <= i__1; ++n) {
	s_wsle(&io___108);
	for (i__ = 1; i__ <= 3; ++i__) {
	    do_lio(&c__3, &c__1, (char *)&ipf[i__ + n * 3], (ftnlen)sizeof(
		    integer));
	}
	do_lio(&c__3, &c__1, (char *)&lbf[n], (ftnlen)sizeof(integer));
	e_wsle();
    }
    cl__1.cerr = 0;
    cl__1.cunit = 10;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
L1000:
    errmesio_(&c__4001, "saveMaft", "Input file name is wrong", (ftnlen)8, (
	    ftnlen)24);
    return 0;
} /* savemaft_ */

