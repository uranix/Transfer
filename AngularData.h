#ifndef __ANGULARDATA_H__
#define __ANGULARDATA_H__

#include "util.h"

struct CudaContext;

struct AngularData {

	idx slm;

	REAL *omega;
	idx *omega_pos;

	REAL *Ox, *Oy, *Oz;
	AngularData(int maxk);
	~AngularData();
};

struct DeviceAngularDataRaw {
	idx slm;
	idx aslm;

	REAL *omega;
	idx *omega_pos;

	REAL *Ox, *Oy, *Oz;
};

struct DeviceAngularData : public DeviceAngularDataRaw {
	const CudaContext *ctx;
	DeviceAngularData(const CudaContext *ctx, const AngularData &);
	~DeviceAngularData();
private:
	DeviceAngularData();
};

#endif
