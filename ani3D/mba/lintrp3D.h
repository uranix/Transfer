/* f2h postprocessed file */
#ifndef __LINTRP3D_H__
#define __LINTRP3D_H__

#include "f2c.h"

/* @f2h@ */ /* Subroutine */ int lintrp3d_(integer *nt, integer *tet, integer 
	*nv, doublereal *vrt, integer *ldf, doublereal *f, integer *nxyz, 
	doublereal *xyz, doublereal *g, integer *imem, integer *nimem, 
	doublereal *dmem, integer *ndmem, integer *icontrol)
;

/* @f2h@ */ /* Subroutine */ int initqt_(integer *nqt, integer *qt, 
	doublereal *h__, doublereal *xyzc)
;

/* @f2h@ */ integer newqt_(integer *nqt, integer *qt)
;

/* @f2h@ */ /* Subroutine */ int setij_(doublereal *xyzc, doublereal *xyz, 
	integer *i__, integer *j, integer *k)
;

/* @f2h@ */ /* Subroutine */ int drops_(integer *nqt, integer *qt, integer *
	nqtav, doublereal *h__, doublereal *xyzc, doublereal *xyz, integer *
	idx)
;

/* @f2h@ */ doublereal sqrdst_(doublereal *a, doublereal *b)
;

/* @f2h@ */ /* Subroutine */ int order2_(doublereal *xyzc, doublereal *h__, 
	doublereal *xyz, integer *ord, doublereal *sqrd)
;

/* @f2h@ */ integer nearst_(integer *qt, doublereal *xyzc, doublereal *xyz, 
	doublereal *point, doublereal *h__)
;

/* @f2h@ */ integer basetet_(integer *qt, doublereal *xyzc, doublereal *vrt, 
	integer *nt, integer *tet, integer *ref, integer *itet, doublereal *
	xyz, doublereal *h__, integer *buf, integer *nbuf, doublereal *prec)
;

/* @f2h@ */ /* Subroutine */ int restore_(integer *nt, integer *tet, integer *
	nv, doublereal *vrt, integer *ldf, doublereal *f, integer *nxyz, 
	doublereal *xyz, doublereal *g, integer *qt, doublereal *xyzc, 
	integer *ref, integer *itet, doublereal *h__, logical *flags, integer 
	*buf, integer *nbuf, integer *nunit)
;

/* @f2h@ */ logical enclose_(doublereal *xyz, doublereal *vrt, integer *tet, 
	doublereal *prec)
;

/* @f2h@ */ doublereal sizeqt_(doublereal *point, integer *imem, doublereal *
	dmem)
;

/* @f2h@ */ doublereal sizehost_(integer *qt, doublereal *xyzc, doublereal *
	point, doublereal *h__)
;

#endif
