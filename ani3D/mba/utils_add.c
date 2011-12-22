/* utils_add.f -- translated by f2c (version 20090411).
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
#include "makM.h"
#include "error.h"
#include "auxSE.h"
#include "makQ.h"


/* Table of constant values */

static integer c__4 = 4;
static integer c__3 = 3;
static integer c__0 = 0;
static integer c__5023 = 5023;
static integer c__1 = 1;

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int check_mesh__(integer *np, integer *nf, 
	integer *ne, doublereal *xyp, integer *ipf, integer *ipe, integer *
	lbf, integer *lbe, integer *nep, integer *iep, integer *ife, integer *
	iee)
{
    /* Format strings */
    static char fmt_5002[] = "(\002Error in face =\002,i6,\002,  Points ="
	    "\002,3i6,\002,  label=\002,i4)";
    static char fmt_5000[] = "(\002Error in checking tetrahedron =\002,i7"
	    ",\002   iERR=\002,i5,/,\002Points =\002,4i7,/,\002Faces  =\002,4"
	    "i7,/,\002Tetras =\002,4i7,/)";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k, n, i1, i2, i3, j1, j2, j3, i4, j4;
    static doublereal v1, v2;
    static integer ie, if__, ip[5], ie1, ie2, if1, if2, if3, if4, ie3, ie4, 
	    ip1, ip2, ip3, ip4, jp1, jp2, jp3, jp4, nfac, ierr, ntet, ifface;

    /* Fortran I/O blocks */
    static cilist io___36 = { 0, 6, 0, fmt_5002, 0 };
    static cilist io___37 = { 0, 6, 0, fmt_5000, 0 };
    static cilist io___39 = { 0, 6, 0, fmt_5000, 0 };


/* ================================================================ */
/* ================================================================ */
/* Routine checks topology of the input and output meshes. */

/* nEP :  working array of size nP */
/* IEP :  working array of size 4*nE */
/* IFE :  working array of size 4*nE */
/* IEE :  working array of size 4*nE */

/* ================================================================ */
/* group (M) */
/* group (W) */
/* ================================================================ */
/* group (Local variables) */
/* ================================================================ */
    /* Parameter adjustments */
    iee -= 5;
    ife -= 5;
    --iep;
    --nep;
    --lbe;
    --lbf;
    ipe -= 5;
    ipf -= 4;
    xyp -= 4;

    /* Function Body */
    ip[0] = 1;
    ip[1] = 2;
    ip[2] = 3;
    ip[3] = 4;
    ip[4] = 1;
/* ... check faces */
    i__1 = *nf;
    for (n = 1; n <= i__1; ++n) {
	nfac = n;
	for (i__ = 1; i__ <= 3; ++i__) {
	    if (ipf[i__ + n * 3] <= 0) {
		goto L400;
	    }
	}
	if (lbf[n] <= 0) {
	    goto L400;
	}
    }
/* ... create an auxiliary structure */
    backreferences_(np, ne, &c__4, &c__4, &ipe[5], &nep[1], &iep[1]);
/* ... create IEE */
    i__1 = *ne;
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
	    ip1 = ipe[i1 + (n << 2)];
	    ip2 = ipe[i2 + (n << 2)];
	    ip3 = ipe[i3 + (n << 2)];
	    if (cmpe_(&ip1, &ip2, &ip3, &iep[1], &nep[1], &n, &ie2)) {
		iee[i1 + (n << 2)] = ie2;
	    }
	}
    }
/* ... create an auxiliary structure */
    backreferences_(np, nf, &c__3, &c__3, &ipf[4], &nep[1], &iep[1]);
/* ... create IFE: basic, fictitious, material */
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	for (i1 = 1; i1 <= 4; ++i1) {
	    i2 = ip[i1];
	    i3 = ip[i2];
	    ip1 = ipe[i1 + (n << 2)];
	    ip2 = ipe[i2 + (n << 2)];
	    ip3 = ipe[i3 + (n << 2)];
	    if (cmpe_(&ip1, &ip2, &ip3, &iep[1], &nep[1], &c__0, &if__)) {
		ife[i1 + (n << 2)] = if__;
	    }
	}
    }
/* ... check for faces */
    i__1 = *np;
    for (n = 1; n <= i__1; ++n) {
	nep[n] = 0;
    }
    i__1 = *ne;
    for (n = 1; n <= i__1; ++n) {
	ntet = n;
	ip1 = ipe[(n << 2) + 1];
	ip2 = ipe[(n << 2) + 2];
	ip3 = ipe[(n << 2) + 3];
	ip4 = ipe[(n << 2) + 4];
	if (ip1 == ip2 || ip1 == ip3 || ip1 == ip4 || ip2 == ip3 || ip2 == 
		ip4 || ip3 == ip4) {
	    ierr = 5013;
	    goto L500;
	}
	if (lbe[n] <= 0) {
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
	    ip1 = ipe[i1 + (n << 2)];
	    ip2 = ipe[i2 + (n << 2)];
	    ip3 = ipe[i3 + (n << 2)];
	    if (if__ != 0 && if__ != ifface) {
		jp1 = ipf[if__ * 3 + 1];
		jp2 = ipf[if__ * 3 + 2];
		jp3 = ipf[if__ * 3 + 3];
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
		    jp1 = ipe[j1 + (ie << 2)];
		    jp2 = ipe[j2 + (ie << 2)];
		    jp3 = ipe[j3 + (ie << 2)];
		    if (check33_(&ip1, &ip2, &ip3, &jp1, &jp2, &jp3)) {
			if (iee[j1 + (ie << 2)] != n) {
			    ierr = 5019;
			    goto L500;
			}
			i4 = ip[i3];
			ip4 = ipe[i4 + (n << 2)];
			j4 = ip[j3];
			jp4 = ipe[j4 + (ie << 2)];
			v1 = calvol_(&xyp[ip1 * 3 + 1], &xyp[ip2 * 3 + 1], &
				xyp[ip3 * 3 + 1], &xyp[ip4 * 3 + 1]);
			v2 = calvol_(&xyp[ip1 * 3 + 1], &xyp[ip2 * 3 + 1], &
				xyp[ip3 * 3 + 1], &xyp[jp4 * 3 + 1]);
			if (v1 * v2 >= 0.) {
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
    return 0;
L400:
    s_wsfe(&io___36);
    do_fio(&c__1, (char *)&nfac, (ftnlen)sizeof(integer));
    for (i__ = 1; i__ <= 3; ++i__) {
	do_fio(&c__1, (char *)&ipf[i__ + nfac * 3], (ftnlen)sizeof(integer));
    }
    do_fio(&c__1, (char *)&lbf[nfac], (ftnlen)sizeof(integer));
    e_wsfe();
    return 0;
L500:
    s_wsfe(&io___37);
    do_fio(&c__1, (char *)&ntet, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ierr, (ftnlen)sizeof(integer));
    for (i__ = 1; i__ <= 4; ++i__) {
	do_fio(&c__1, (char *)&ipe[i__ + (ntet << 2)], (ftnlen)sizeof(integer)
		);
    }
    for (i__ = 1; i__ <= 4; ++i__) {
	do_fio(&c__1, (char *)&ife[i__ + (ntet << 2)], (ftnlen)sizeof(integer)
		);
    }
    for (i__ = 1; i__ <= 4; ++i__) {
	do_fio(&c__1, (char *)&iee[i__ + (ntet << 2)], (ftnlen)sizeof(integer)
		);
    }
    e_wsfe();
    for (k = 1; k <= 4; ++k) {
	ie = iee[k + (ntet << 2)];
	if (ie > 0) {
	    s_wsfe(&io___39);
	    do_fio(&c__1, (char *)&ie, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&ierr, (ftnlen)sizeof(integer));
	    for (i__ = 1; i__ <= 4; ++i__) {
		do_fio(&c__1, (char *)&ipe[i__ + (ie << 2)], (ftnlen)sizeof(
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
	    e_wsfe();
	}
    }
    errmes_(&ierr, "chkM", "tetrahedra are wrong", (ftnlen)4, (ftnlen)20);
/* L5004: */
    return 0;
} /* check_mesh__ */

