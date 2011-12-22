/* f2h postprocessed file */
#ifndef __ERROR_H__
#define __ERROR_H__

#include "f2c.h"

/* @f2h@ */ /* Subroutine */ int wrnmes_(integer *ierr, char *routine, char *
	message, ftnlen routine_len, ftnlen message_len)
;

/* @f2h@ */ /* Subroutine */ int errmes_(integer *ierr, char *routine, char *
	message, ftnlen routine_len, ftnlen message_len)
;

/* @f2h@ */ /* Subroutine */ int errmesio_(integer *ierr, char *routine, char 
	*message, ftnlen routine_len, ftnlen message_len)
;

/* @f2h@ */ logical probeany_(void)
;

#endif
