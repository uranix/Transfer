#ifndef __ANGULARDATA_H__
#define __ANGULARDATA_H__

#include "util.h"

struct AngularData {

	idx slm;

	REAL *omega;
	idx *omega_pos;

	REAL *Ox, *Oy, *Oz;
	AngularData(int maxk);
	~AngularData();
};

#endif
