/* output.f -- translated by f2c (version 20090411).
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

#include "output.h"

/* Table of constant values */

static integer c__1 = 1;

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int countbadelements_(integer *ne, integer *l1e, 
	integer *l2e, doublereal *qe, doublereal *quality, integer *nlines, 
	char *output, logical *flagfile, ftnlen output_len)
{
    /* Format strings */
    static char fmt_6005[] = "(\002 *** Distribution of tetrahedrons by the "
	    "quality\002)";
    static char fmt_6006[] = "(10f6.1)";
    static char fmt_6008[] = "(10i6)";
    static char fmt_6007[] = "(\002 *** Number of bad tetrahedrons =\002,i8)";
    static char fmt_5005[] = "(\002 *** Distribution of tetrahedrons by the "
	    "quality\002,/,10f8.1,/,10i8,/,\002 *** Number of bad tetrahedron"
	    "s =\002,i8,/)";

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__, k, n, ie, nbad[10], icnt;
    static char message[100];

    /* Fortran I/O blocks */
    static icilist io___7 = { 0, message, 0, fmt_6005, 100, 1 };
    static icilist io___8 = { 0, message, 0, fmt_6006, 100, 1 };
    static icilist io___10 = { 0, message, 0, fmt_6008, 100, 1 };
    static icilist io___11 = { 0, message, 0, fmt_6007, 100, 1 };
    static cilist io___12 = { 0, 6, 0, fmt_5005, 0 };


/* ================================================================ */
/* ================================================================ */
/* ========================================================== */
/* MaxLines   - the maximal number of output lines per processor */
/* LineLenght - the maximal lenght of a message */
/* ========================================================== */
/* group (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    output -= 100;
    --qe;
    --l2e;
    l1e -= 3;

    /* Function Body */
    for (n = 1; n <= 10; ++n) {
	nbad[n - 1] = 0;
    }
    icnt = 0;
    ie = l2e[1];
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	if (qe[ie] < *quality) {
	    ++icnt;
	}
	for (k = 1; k <= 10; ++k) {
	    if ((k - 1) * .1 <= qe[ie] && qe[ie] < k * .1) {
		++nbad[k - 1];
		goto L500;
	    }
	}
L500:
	ie = l1e[(ie << 1) + 2];
    }
    if (*flagfile) {
	s_wsfi(&io___7);
	e_wsfi();
	addout_(nlines, message, output + 100, (ftnlen)100, (ftnlen)100);
	s_wsfi(&io___8);
	for (i__ = 1; i__ <= 10; ++i__) {
	    d__1 = i__ * .1;
	    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	}
	e_wsfi();
	addout_(nlines, message, output + 100, (ftnlen)100, (ftnlen)100);
	s_wsfi(&io___10);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&nbad[i__ - 1], (ftnlen)sizeof(integer));
	}
	e_wsfi();
	addout_(nlines, message, output + 100, (ftnlen)100, (ftnlen)100);
	s_wsfi(&io___11);
	do_fio(&c__1, (char *)&icnt, (ftnlen)sizeof(integer));
	e_wsfi();
	addout_(nlines, message, output + 100, (ftnlen)100, (ftnlen)100);
    } else {
	s_wsfe(&io___12);
	for (i__ = 1; i__ <= 10; ++i__) {
	    d__1 = i__ * .1;
	    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	}
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&nbad[i__ - 1], (ftnlen)sizeof(integer));
	}
	do_fio(&c__1, (char *)&icnt, (ftnlen)sizeof(integer));
	e_wsfe();
    }
/* ... screen output */
/* ... parallel output */
    return 0;
} /* countbadelements_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int addout_(integer *nlines, char *message, char 
	*output, ftnlen message_len, ftnlen output_len)
{
    /* System generated locals */
    integer i__1;

/* ================================================================ */
/* ================================================================ */
/* ========================================================== */
/* MaxLines   - the maximal number of output lines per processor */
/* LineLenght - the maximal lenght of a message */
/* ========================================================== */
    /* Parameter adjustments */
    output -= output_len;

    /* Function Body */
/* Computing MIN */
    i__1 = *nlines + 1;
    *nlines = min(i__1,100);
    s_copy(output + *nlines * output_len, message, output_len, message_len);
    return 0;
} /* addout_ */

