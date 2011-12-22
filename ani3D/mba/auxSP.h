/* f2h postprocessed file */
#ifndef __AUXSP_H__
#define __AUXSP_H__

#include "f2c.h"

/* @f2h@ */ /* Subroutine */ int maksp_(integer *ip, integer *iep, integer *
	ipe, integer *iee, integer *maxs, integer *le, integer *ies)
;

/* @f2h@ */ /* Subroutine */ int chkspf_(integer *nm, integer *ip1, integer *
	ioperat, integer *icp, integer *iep, integer *ipe, integer *iee, 
	integer *lpf, integer *ipf)
;

/* @f2h@ */ /* Subroutine */ int info_(void)
;

/* @f2h@ */ /* Subroutine */ int chkspb_(integer *nm, integer *ld, integer *
	ids, integer *ln, integer *ins, integer *ioperat, integer *icp, 
	integer *iep, integer *ipe, integer *iee, integer *lpf, integer *ipf, 
	integer *ipbad, logical *flag__)
;

/* @f2h@ */ doublereal calnorm_(doublereal *xyz)
;

/* @f2h@ */ doublereal sqrnorm_(doublereal *xyz)
;

/* @f2h@ */ logical ifxnode_(integer *clr, integer *ixnode)
;

/* @f2h@ */ /* Subroutine */ int addxnode_(integer *clr, integer *ixnode)
;

/* @f2h@ */ /* Subroutine */ int delxnode_(integer *clr, integer *ixnode)
;

/* @f2h@ */ integer minclr_(integer *clr1, integer *clr2)
;

/* @f2h@ */ integer maxclr_(integer *clr1, integer *clr2)
;

/* @f2h@ */ /* Subroutine */ int setstatus_(logical *flagauto, integer *
	status, integer *iprint)
;

#endif
