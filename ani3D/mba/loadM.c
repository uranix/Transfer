/* loadM.f -- translated by f2c (version 20090411).
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

/* Table of constant values */

static integer c__1 = 1;
static integer c__8 = 8;
static integer c__3 = 3;
static integer c__1003 = 1003;
static integer c__5 = 5;
static integer c__1004 = 1004;
static integer c__1006 = 1006;
static integer c__4001 = 4001;
static integer c__2 = 2;
static integer c__9 = 9;
static integer c__1005 = 1005;
static integer c__1007 = 1007;

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int loadmani_(integer *maxp, integer *maxf, 
	integer *maxe, integer *np, integer *nf, integer *ne, doublereal *xyp,
	 integer *ipf, integer *ipe, integer *lbf, integer *lbe, integer *npv,
	 integer *nfv, integer *nev, integer *ipv, integer *ifv, integer *iev,
	 integer *ipp, integer *iff, char *fname, ftnlen fname_len)
{
    /* System generated locals */
    integer i__1;
    olist o__1;
    cllist cl__1;

    /* Local variables */
    static integer i__, j, n, nc, nfo, npo, icrv[2];
    static doublereal rcrv[6];
    static logical flagc, flage, flagf, flagi, flagp, flagev, flagfv, flagpv;

    /* Fortran I/O blocks */
    static cilist io___2 = { 0, 6, 0, "(A,A)", 0 };
    static cilist io___3 = { 0, 10, 0, 0, 0 };
    static cilist io___5 = { 0, 10, 0, 0, 0 };
    static cilist io___7 = { 0, 10, 0, 0, 0 };
    static cilist io___9 = { 0, 10, 0, 0, 0 };
    static cilist io___11 = { 0, 10, 0, 0, 0 };
    static cilist io___13 = { 0, 10, 0, 0, 0 };
    static cilist io___15 = { 0, 10, 0, 0, 0 };
    static cilist io___17 = { 0, 10, 0, 0, 0 };
    static cilist io___19 = { 0, 10, 0, 0, 0 };
    static cilist io___20 = { 0, 10, 0, 0, 0 };
    static cilist io___22 = { 0, 10, 0, 0, 0 };
    static cilist io___24 = { 0, 10, 0, 0, 0 };
    static cilist io___25 = { 0, 10, 0, 0, 0 };
    static cilist io___26 = { 0, 10, 0, 0, 0 };
    static cilist io___27 = { 0, 10, 0, 0, 0 };
    static cilist io___28 = { 0, 10, 0, 0, 0 };
    static cilist io___29 = { 0, 10, 0, 0, 0 };
    static cilist io___30 = { 0, 10, 0, 0, 0 };
    static cilist io___31 = { 0, 10, 0, 0, 0 };
    static cilist io___33 = { 0, 10, 0, 0, 0 };
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
    static cilist io___49 = { 0, 10, 0, 0, 0 };
    static cilist io___50 = { 0, 10, 0, 0, 0 };


/* ================================================================ */
/* group (M) */
/* group (Dev) */
/* group (I) */
/* ================================================================ */
/* group (M) */
/*     Integer MaxP, MaxF, MaxE */
/* group (Dev) */
/* group (I) */
/* group (Local variables) */
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
    do_fio(&c__1, "Loading mesh ", (ftnlen)13);
    do_fio(&c__1, fname, i__ + 3);
    e_wsfe();
    o__1.oerr = 1;
    o__1.ounit = 10;
    o__1.ofnmlen = i__ + 3;
    o__1.ofnm = fname;
    o__1.orl = 0;
    o__1.osta = "OLD";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__1 = f_open(&o__1);
    if (i__1 != 0) {
	goto L1000;
    }
    s_rsle(&io___3);
    do_lio(&c__8, &c__1, (char *)&flagp, (ftnlen)sizeof(logical));
    e_rsle();
    s_rsle(&io___5);
    do_lio(&c__8, &c__1, (char *)&flagf, (ftnlen)sizeof(logical));
    e_rsle();
    s_rsle(&io___7);
    do_lio(&c__8, &c__1, (char *)&flage, (ftnlen)sizeof(logical));
    e_rsle();
    s_rsle(&io___9);
    do_lio(&c__8, &c__1, (char *)&flagc, (ftnlen)sizeof(logical));
    e_rsle();
    s_rsle(&io___11);
    do_lio(&c__8, &c__1, (char *)&flagpv, (ftnlen)sizeof(logical));
    e_rsle();
    s_rsle(&io___13);
    do_lio(&c__8, &c__1, (char *)&flagfv, (ftnlen)sizeof(logical));
    e_rsle();
    s_rsle(&io___15);
    do_lio(&c__8, &c__1, (char *)&flagev, (ftnlen)sizeof(logical));
    e_rsle();
    s_rsle(&io___17);
    do_lio(&c__8, &c__1, (char *)&flagi, (ftnlen)sizeof(logical));
    e_rsle();
    s_rsle(&io___19);
    e_rsle();
    s_rsle(&io___20);
    do_lio(&c__3, &c__1, (char *)&(*np), (ftnlen)sizeof(integer));
    e_rsle();
    if (*np > *maxp) {
	errmesio_(&c__1003, "loadMani", "local parameter MaxP is small", (
		ftnlen)8, (ftnlen)29);
    }
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	s_rsle(&io___22);
	for (j = 1; j <= 3; ++j) {
	    do_lio(&c__5, &c__1, (char *)&xyp[j + n * 3], (ftnlen)sizeof(
		    doublereal));
	}
	e_rsle();
    }
    s_rsle(&io___24);
    e_rsle();
    s_rsle(&io___25);
    do_lio(&c__3, &c__1, (char *)&(*nf), (ftnlen)sizeof(integer));
    e_rsle();
    if (*nf > *maxf) {
	errmesio_(&c__1004, "loadMani", "local parameter MaxF is small", (
		ftnlen)8, (ftnlen)29);
    }
    i__1 = *nf;
    for (n = 1; n <= i__1; ++n) {
	s_rsle(&io___26);
	for (j = 1; j <= 3; ++j) {
	    do_lio(&c__3, &c__1, (char *)&ipf[j + n * 3], (ftnlen)sizeof(
		    integer));
	}
	do_lio(&c__3, &c__1, (char *)&lbf[n], (ftnlen)sizeof(integer));
	e_rsle();
    }
    s_rsle(&io___27);
    e_rsle();
    s_rsle(&io___28);
    do_lio(&c__3, &c__1, (char *)&(*ne), (ftnlen)sizeof(integer));
    e_rsle();
    if (*ne > *maxe) {
	errmesio_(&c__1006, "loadMani", "local parameter MaxE is small", (
		ftnlen)8, (ftnlen)29);
    }
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	s_rsle(&io___29);
	for (j = 1; j <= 4; ++j) {
	    do_lio(&c__3, &c__1, (char *)&ipe[j + (n << 2)], (ftnlen)sizeof(
		    integer));
	}
	do_lio(&c__3, &c__1, (char *)&lbe[n], (ftnlen)sizeof(integer));
	e_rsle();
    }
/* ... reading curvilinear faces */
    s_rsle(&io___30);
    e_rsle();
    s_rsle(&io___31);
    do_lio(&c__3, &c__1, (char *)&nc, (ftnlen)sizeof(integer));
    e_rsle();
    i__1 = nc;
    for (n = 1; n <= i__1; ++n) {
	s_rsle(&io___33);
	for (i__ = 1; i__ <= 2; ++i__) {
	    do_lio(&c__3, &c__1, (char *)&icrv[i__ - 1], (ftnlen)sizeof(
		    integer));
	}
	for (i__ = 1; i__ <= 6; ++i__) {
	    do_lio(&c__5, &c__1, (char *)&rcrv[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	}
	e_rsle();
    }
/* ... reading the fixed mesh points, faces and elements */
    s_rsle(&io___36);
    e_rsle();
    s_rsle(&io___37);
    do_lio(&c__3, &c__1, (char *)&(*npv), (ftnlen)sizeof(integer));
    e_rsle();
    i__1 = *npv;
    for (n = 1; n <= i__1; ++n) {
	s_rsle(&io___38);
	do_lio(&c__3, &c__1, (char *)&ipv[n], (ftnlen)sizeof(integer));
	e_rsle();
    }
    s_rsle(&io___39);
    e_rsle();
    s_rsle(&io___40);
    do_lio(&c__3, &c__1, (char *)&(*nfv), (ftnlen)sizeof(integer));
    e_rsle();
    i__1 = *nfv;
    for (n = 1; n <= i__1; ++n) {
	s_rsle(&io___41);
	do_lio(&c__3, &c__1, (char *)&ifv[n], (ftnlen)sizeof(integer));
	e_rsle();
    }
    s_rsle(&io___42);
    e_rsle();
    s_rsle(&io___43);
    do_lio(&c__3, &c__1, (char *)&(*nev), (ftnlen)sizeof(integer));
    e_rsle();
    i__1 = *nev;
    for (n = 1; n <= i__1; ++n) {
	s_rsle(&io___44);
	do_lio(&c__3, &c__1, (char *)&iev[n], (ftnlen)sizeof(integer));
	e_rsle();
    }
/* ... reading interface data */
    if (flagi) {
	s_rsle(&io___45);
	e_rsle();
	s_rsle(&io___46);
	do_lio(&c__3, &c__1, (char *)&npo, (ftnlen)sizeof(integer));
	do_lio(&c__3, &c__1, (char *)&nfo, (ftnlen)sizeof(integer));
	e_rsle();
	i__1 = *np;
	for (n = 1; n <= i__1; ++n) {
	    ipp[n] = 0;
	}
	i__1 = *nf;
	for (n = 1; n <= i__1; ++n) {
	    iff[n] = 0;
	}
	i__1 = min(*np,npo);
	for (n = 1; n <= i__1; ++n) {
	    s_rsle(&io___49);
	    do_lio(&c__3, &c__1, (char *)&ipp[n], (ftnlen)sizeof(integer));
	    e_rsle();
	}
	i__1 = min(*nf,nfo);
	for (n = 1; n <= i__1; ++n) {
	    s_rsle(&io___50);
	    do_lio(&c__3, &c__1, (char *)&iff[n], (ftnlen)sizeof(integer));
	    e_rsle();
	}
    }
    cl__1.cerr = 0;
    cl__1.cunit = 10;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
L1000:
    errmesio_(&c__4001, "loadMani", "Input file name is wrong ", (ftnlen)8, (
	    ftnlen)25);
    return 0;
} /* loadmani_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int loadmani_header__(integer *np, integer *nf, 
	integer *ne, integer *npv, integer *nfv, integer *nev, char *fname, 
	ftnlen fname_len)
{
    /* System generated locals */
    address a__1[2];
    integer i__1[2], i__2;
    olist o__1;
    cllist cl__1;

    /* Local variables */
    static char fnameext[30];
    static integer i__, nc;
    static logical flag__;
    static char text[30], text2[30];

    /* Fortran I/O blocks */
    static cilist io___53 = { 0, 6, 0, "(A,A)", 0 };
    static cilist io___54 = { 0, 10, 0, 0, 0 };
    static cilist io___57 = { 0, 10, 0, 0, 0 };
    static cilist io___58 = { 0, 10, 0, 0, 0 };
    static cilist io___59 = { 0, 10, 0, 0, 0 };
    static cilist io___62 = { 0, 10, 0, 0, 0 };
    static cilist io___63 = { 0, 10, 0, 0, 0 };
    static cilist io___64 = { 0, 10, 0, 0, 0 };


/* ================================================================ */
/* Routines reads the header from the input file */
/* ================================================================ */
/* ================================================================ */
    i__ = 1;
    while(s_cmp(fname + (i__ - 1), ".ani", (ftnlen)4, (ftnlen)4) != 0) {
	++i__;
    }
/* Writing concatenation */
    i__1[0] = i__, a__1[0] = fname;
    i__1[1] = 3, a__1[1] = "ani";
    s_cat(fnameext, a__1, i__1, &c__2, (ftnlen)30);
    s_wsfe(&io___53);
    do_fio(&c__1, "Reading header of mesh ", (ftnlen)23);
    do_fio(&c__1, fnameext, (ftnlen)30);
    e_wsfe();
    o__1.oerr = 1;
    o__1.ounit = 10;
    o__1.ofnmlen = 30;
    o__1.ofnm = fnameext;
    o__1.orl = 0;
    o__1.osta = "OLD";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__2 = f_open(&o__1);
    if (i__2 != 0) {
	goto L1000;
    }
    s_rsle(&io___54);
    do_lio(&c__8, &c__1, (char *)&flag__, (ftnlen)sizeof(logical));
    do_lio(&c__9, &c__1, text, (ftnlen)30);
    do_lio(&c__3, &c__1, (char *)&(*np), (ftnlen)sizeof(integer));
    e_rsle();
    s_rsle(&io___57);
    do_lio(&c__8, &c__1, (char *)&flag__, (ftnlen)sizeof(logical));
    do_lio(&c__9, &c__1, text, (ftnlen)30);
    do_lio(&c__3, &c__1, (char *)&(*nf), (ftnlen)sizeof(integer));
    e_rsle();
    s_rsle(&io___58);
    do_lio(&c__8, &c__1, (char *)&flag__, (ftnlen)sizeof(logical));
    do_lio(&c__9, &c__1, text, (ftnlen)30);
    do_lio(&c__3, &c__1, (char *)&(*ne), (ftnlen)sizeof(integer));
    e_rsle();
    s_rsle(&io___59);
    do_lio(&c__8, &c__1, (char *)&flag__, (ftnlen)sizeof(logical));
    do_lio(&c__9, &c__1, text, (ftnlen)30);
    do_lio(&c__9, &c__1, text2, (ftnlen)30);
    do_lio(&c__3, &c__1, (char *)&nc, (ftnlen)sizeof(integer));
    e_rsle();
    s_rsle(&io___62);
    do_lio(&c__8, &c__1, (char *)&flag__, (ftnlen)sizeof(logical));
    do_lio(&c__9, &c__1, text, (ftnlen)30);
    do_lio(&c__9, &c__1, text2, (ftnlen)30);
    do_lio(&c__3, &c__1, (char *)&(*npv), (ftnlen)sizeof(integer));
    e_rsle();
    s_rsle(&io___63);
    do_lio(&c__8, &c__1, (char *)&flag__, (ftnlen)sizeof(logical));
    do_lio(&c__9, &c__1, text, (ftnlen)30);
    do_lio(&c__9, &c__1, text2, (ftnlen)30);
    do_lio(&c__3, &c__1, (char *)&(*nfv), (ftnlen)sizeof(integer));
    e_rsle();
    s_rsle(&io___64);
    do_lio(&c__8, &c__1, (char *)&flag__, (ftnlen)sizeof(logical));
    do_lio(&c__9, &c__1, text, (ftnlen)30);
    do_lio(&c__9, &c__1, text2, (ftnlen)30);
    do_lio(&c__3, &c__1, (char *)&(*nev), (ftnlen)sizeof(integer));
    e_rsle();
    cl__1.cerr = 0;
    cl__1.cunit = 10;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
L1000:
    errmes_(&c__4001, "loadMani_header", "File name is wrong", (ftnlen)15, (
	    ftnlen)18);
    return 0;
} /* loadmani_header__ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int loads_(integer *np, doublereal *sol, char *
	fname, ftnlen fname_len)
{
    /* System generated locals */
    integer i__1;
    olist o__1;
    cllist cl__1;

    /* Local variables */
    static char fnameext[70];
    static integer n, npw;

    /* Fortran I/O blocks */
    static cilist io___65 = { 0, 6, 0, "(A,A)", 0 };
    static cilist io___67 = { 0, 10, 0, 0, 0 };
    static cilist io___69 = { 0, 10, 0, 0, 0 };
    static cilist io___71 = { 0, 10, 0, 0, 0 };


/* ================================================================ */
/* group (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    --sol;

    /* Function Body */
    s_wsfe(&io___65);
    do_fio(&c__1, "Loading solution ", (ftnlen)17);
    do_fio(&c__1, fname, fname_len);
    e_wsfe();
/* ... read the solution associated to the mesh */
    s_copy(fnameext, fname, (ftnlen)70, fname_len);
    o__1.oerr = 1;
    o__1.ounit = 10;
    o__1.ofnmlen = 70;
    o__1.ofnm = fnameext;
    o__1.orl = 0;
    o__1.osta = "OLD";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__1 = f_open(&o__1);
    if (i__1 != 0) {
	goto L1000;
    }
    s_rsle(&io___67);
    do_lio(&c__3, &c__1, (char *)&npw, (ftnlen)sizeof(integer));
    e_rsle();
    if (npw != *np) {
	errmesio_(&c__4001, "loadS", "the lines number differs from the poin"
		"ts number", (ftnlen)5, (ftnlen)47);
    }
    s_rsle(&io___69);
    e_rsle();
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	s_rsle(&io___71);
	do_lio(&c__5, &c__1, (char *)&sol[n], (ftnlen)sizeof(doublereal));
	e_rsle();
    }
    cl__1.cerr = 0;
    cl__1.cunit = 10;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
L1000:
    errmesio_(&c__4001, "loadS", "Input file name is wrong ", (ftnlen)5, (
	    ftnlen)25);
    return 0;
} /* loads_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int loadmaft_(integer *maxp, integer *maxf, 
	integer *maxe, integer *np, integer *nf, integer *ne, doublereal *xyp,
	 integer *ipf, integer *ipe, integer *lbf, integer *lbe, integer *npv,
	 integer *nfv, integer *nev, integer *ipv, integer *ifv, integer *iev,
	 integer *ipp, integer *iff, char *fname, ftnlen fname_len)
{
    /* System generated locals */
    integer i__1;
    olist o__1;
    cllist cl__1;

    /* Local variables */
    static integer i__, n;

    /* Fortran I/O blocks */
    static cilist io___73 = { 0, 6, 0, "(A,A)", 0 };
    static cilist io___74 = { 0, 10, 0, 0, 0 };
    static cilist io___76 = { 0, 10, 0, 0, 0 };
    static cilist io___77 = { 0, 10, 0, 0, 0 };
    static cilist io___78 = { 0, 10, 0, 0, 0 };
    static cilist io___79 = { 0, 10, 0, 0, 0 };
    static cilist io___80 = { 0, 10, 0, 0, 0 };


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
    while(s_cmp(fname + (i__ - 1), ".out", (ftnlen)4, (ftnlen)4) != 0) {
	++i__;
    }
    s_wsfe(&io___73);
    do_fio(&c__1, "Loading mesh ", (ftnlen)13);
    do_fio(&c__1, fname, i__ + 3);
    e_wsfe();
    o__1.oerr = 1;
    o__1.ounit = 10;
    o__1.ofnmlen = i__ + 3;
    o__1.ofnm = fname;
    o__1.orl = 0;
    o__1.osta = "OLD";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__1 = f_open(&o__1);
    if (i__1 != 0) {
	goto L1000;
    }
/* ... read points */
    s_rsle(&io___74);
    do_lio(&c__3, &c__1, (char *)&(*np), (ftnlen)sizeof(integer));
    e_rsle();
    if (*np > *maxp) {
	errmesio_(&c__1005, "loadMaft", "local parameter MaxP is small", (
		ftnlen)8, (ftnlen)29);
    }
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	s_rsle(&io___76);
	for (i__ = 1; i__ <= 3; ++i__) {
	    do_lio(&c__5, &c__1, (char *)&xyp[i__ + n * 3], (ftnlen)sizeof(
		    doublereal));
	}
	e_rsle();
    }
/* ... read elements */
    s_rsle(&io___77);
    do_lio(&c__3, &c__1, (char *)&(*ne), (ftnlen)sizeof(integer));
    e_rsle();
    if (*ne > *maxe) {
	errmesio_(&c__1006, "loadMaft", "local parameter MaxE is small", (
		ftnlen)8, (ftnlen)29);
    }
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	s_rsle(&io___78);
	for (i__ = 1; i__ <= 4; ++i__) {
	    do_lio(&c__3, &c__1, (char *)&ipe[i__ + (n << 2)], (ftnlen)sizeof(
		    integer));
	}
	do_lio(&c__3, &c__1, (char *)&lbe[n], (ftnlen)sizeof(integer));
	e_rsle();
    }
/* ... read faces */
    s_rsle(&io___79);
    do_lio(&c__3, &c__1, (char *)&(*nf), (ftnlen)sizeof(integer));
    e_rsle();
    if (*nf > *maxf) {
	errmesio_(&c__1007, "loadMaft", "local parameter MaxF is small", (
		ftnlen)8, (ftnlen)29);
    }
    i__1 = *nf;
    for (n = 1; n <= i__1; ++n) {
	s_rsle(&io___80);
	for (i__ = 1; i__ <= 3; ++i__) {
	    do_lio(&c__3, &c__1, (char *)&ipf[i__ + n * 3], (ftnlen)sizeof(
		    integer));
	}
	do_lio(&c__3, &c__1, (char *)&lbf[n], (ftnlen)sizeof(integer));
	e_rsle();
    }
    cl__1.cerr = 0;
    cl__1.cunit = 10;
    cl__1.csta = 0;
    f_clos(&cl__1);
/* ... set up the other parameters to default values */
    *npv = 0;
    *nfv = 0;
    *nev = 0;
    return 0;
L1000:
    errmesio_(&c__4001, "loadMaft", "Input file name is wrong", (ftnlen)8, (
	    ftnlen)24);
    return 0;
} /* loadmaft_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int loadmgmv_(integer *maxp, integer *maxf, 
	integer *maxe, integer *np, integer *nf, integer *ne, doublereal *xyp,
	 integer *ipf, integer *ipe, integer *lbf, integer *lbe, integer *npv,
	 integer *nfv, integer *nev, integer *ipv, integer *ifv, integer *iev,
	 integer *ipp, integer *iff, char *fname, ftnlen fname_len)
{
    /* System generated locals */
    integer i__1;
    olist o__1;
    cllist cl__1;

    /* Local variables */
    static integer i__, k, n, nc, nmat;
    static char text[30];

    /* Fortran I/O blocks */
    static cilist io___81 = { 0, 6, 0, "(A,A)", 0 };
    static cilist io___83 = { 0, 10, 0, 0, 0 };
    static cilist io___84 = { 0, 10, 0, 0, 0 };
    static cilist io___87 = { 0, 10, 0, 0, 0 };
    static cilist io___88 = { 0, 10, 0, 0, 0 };
    static cilist io___89 = { 0, 10, 0, 0, 0 };
    static cilist io___91 = { 0, 10, 0, 0, 0 };
    static cilist io___93 = { 0, 10, 0, 0, 0 };
    static cilist io___94 = { 0, 10, 0, 0, 0 };


/* ================================================================ */
/* group (M) */
/* group (Dev) */
/* group (I) */
/* ================================================================ */
/* group (M) */
/*     Integer MaxP, MaxF, MaxE */
/* group (Dev) */
/* group (I) */
/* group (Local variables */
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
    s_wsfe(&io___81);
    do_fio(&c__1, "Loading mesh ", (ftnlen)13);
    do_fio(&c__1, fname, fname_len);
    e_wsfe();
    o__1.oerr = 1;
    o__1.ounit = 10;
    o__1.ofnmlen = fname_len;
    o__1.ofnm = fname;
    o__1.orl = 0;
    o__1.osta = "OLD";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__1 = f_open(&o__1);
    if (i__1 != 0) {
	goto L1000;
    }
/* ... read header */
    for (i__ = 1; i__ <= 2; ++i__) {
	s_rsle(&io___83);
	e_rsle();
    }
/* ... read points */
    s_rsle(&io___84);
    do_lio(&c__9, &c__1, text, (ftnlen)30);
    do_lio(&c__3, &c__1, (char *)&(*np), (ftnlen)sizeof(integer));
    e_rsle();
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	s_rsle(&io___87);
	for (i__ = 1; i__ <= 3; ++i__) {
	    do_lio(&c__5, &c__1, (char *)&xyp[i__ + n * 3], (ftnlen)sizeof(
		    doublereal));
	}
	e_rsle();
    }
/* ... read elements */
    s_rsle(&io___88);
    do_lio(&c__9, &c__1, text, (ftnlen)30);
    do_lio(&c__3, &c__1, (char *)&(*ne), (ftnlen)sizeof(integer));
    e_rsle();
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	s_rsle(&io___89);
	do_lio(&c__9, &c__1, text, (ftnlen)30);
	do_lio(&c__3, &c__1, (char *)&k, (ftnlen)sizeof(integer));
	for (i__ = 1; i__ <= 4; ++i__) {
	    do_lio(&c__3, &c__1, (char *)&ipe[i__ + (n << 2)], (ftnlen)sizeof(
		    integer));
	}
	e_rsle();
    }
    s_rsle(&io___91);
    do_lio(&c__9, &c__1, text, (ftnlen)30);
    do_lio(&c__3, &c__1, (char *)&nmat, (ftnlen)sizeof(integer));
    e_rsle();
    i__1 = nmat;
    for (n = 1; n <= i__1; ++n) {
	s_rsle(&io___93);
	do_lio(&c__9, &c__1, text, (ftnlen)30);
	e_rsle();
    }
    s_rsle(&io___94);
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	do_lio(&c__3, &c__1, (char *)&lbe[n], (ftnlen)sizeof(integer));
    }
    e_rsle();
    cl__1.cerr = 0;
    cl__1.cunit = 10;
    cl__1.csta = 0;
    f_clos(&cl__1);
/* ... set up the other parameters to default values */
    *npv = 0;
    *nf = 0;
    *nfv = 0;
    nc = 0;
    *nev = 0;
    return 0;
L1000:
    errmesio_(&c__4001, "loadMgmv", "Input file name is wrong", (ftnlen)8, (
	    ftnlen)24);
    return 0;
} /* loadmgmv_ */

