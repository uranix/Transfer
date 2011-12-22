/* error.f -- translated by f2c (version 20090411).
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

/* Table of constant values */

static integer c__1 = 1;

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int wrnmes_(integer *ierr, char *routine, char *
	message, ftnlen routine_len, ftnlen message_len)
{
    /* Format strings */
    static char fmt_5000[] = "(/,\002Error:\002,i7,/,\002Routine: \002,a,/"
	    ",\002Message: \002,a,/)";

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, fmt_5000, 0 };


/* ================================================================ */
/* Warnings:        1000 - the quality has been not reached */

/* Normal errors:   1000 - the quality has been not reached */
/* (memory errors)   1001 - not enough memory for Integer arrays */
/*                  1002 - not enough memory for Real*8 arrays */
/*                  1003 - local parameter MaxP is small */
/*                  1004 - local parameter MaxF is small */
/*                  1006 - local parameter MaxE is small */
/*                  1007 - local parameter MaxS is small */
/*                  1009 - local parameter MaxH is small */
/*                  1010 - not enough memory for material faces */
/*                  1011 - reserved boundary identificator is used */
/*                  1012 - local parameter MaxJP is small */
/*                  1013 - local parameter MaxIA is small */
/*                  1014 - local parameter MaxA  is small */
/*                  1101 - one of the local mesh parameters */
/*                        (MaxP, MaxF or MaxE) is small */

/* User errors:     2001 : 2003 - errors due to probably incorrect */
/*                                user routines */
/*                  2011 : 2011 - incorrect or out of bounds input data */

/* Library errors:  3001 : 3002 - errors in routine DSORT */
/*                  3011 : 3012 - errors in routine DSYEV */
/*                  3021        - errors in one of the AMG routines */

/* I/O errors:      4001 : 4002 - errors in input files */
/*                  4101 : 4103 - errors in input data */
/*                  4201        - output chanel number is wrong */

/* Internal errors: 5001 : 5022 - errors in mesh checking (1st level) */
/*                  5101 : 5102 - errors in list updating */
/*                  5201 : 5201 - errors in mesh checking (2nd level) */
/*                  6001 : 6005 - system errors */
/*                  6101 : 6102 - errorneous input for LINTRP3D */
/*                  6103 : 6105 - errors in LINTRP3D algorithm */
/*                  7001 : 7002 - errors in applications working with */
/*                                curvilinear boundaries */

/* Debug errors:    8001 : 8030 - debug errors in ani2_test.f */

/* ================================================================ */
/* ================================================================ */
    s_wsfe(&io___1);
    do_fio(&c__1, (char *)&(*ierr), (ftnlen)sizeof(integer));
    do_fio(&c__1, routine, routine_len);
    do_fio(&c__1, message, message_len);
    e_wsfe();
    return 0;
} /* wrnmes_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int errmes_(integer *ierr, char *routine, char *
	message, ftnlen routine_len, ftnlen message_len)
{

/* ================================================================ */
/* ================================================================ */
    wrnmes_(ierr, routine, message, routine_len, message_len);
    s_stop("", (ftnlen)0);
    return 0;
} /* errmes_ */

/* ================================================================ */
/* @f2h@ */ /* Subroutine */ int errmesio_(integer *ierr, char *routine, char 
	*message, ftnlen routine_len, ftnlen message_len)
{

/* ================================================================ */
/* ================================================================ */
    wrnmes_(ierr, routine, message, routine_len, message_len);
    s_stop("", (ftnlen)0);
    return 0;
} /* errmesio_ */

/* ========================================================== */
/* @f2h@ */ logical probeany_(void)
{
    /* System generated locals */
    logical ret_val;

/* ========================================================== */
/* This is the dummy functions for the simulation modes. */
/* ========================================================== */
    ret_val = FALSE_;
    return ret_val;
} /* probeany_ */

