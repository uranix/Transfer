#include "AngularData.h"
#include "CudaContext.h"
#include "common.cuh"

#include <stdio.h>

DeviceAngularData::DeviceAngularData(const CudaContext *_ctx, const AngularData &host) {
	ctx = _ctx;
	slm = host.slm;
	aslm = align_power(host.slm, COALESCED_NUM(REAL));
	omega = (REAL *)ctx->deviceAlloc(3*aslm*aslm*sizeof(REAL));
	omega_pos = (idx *)ctx->deviceAlloc(3*aslm*aslm*sizeof(idx));
	Ox = (REAL *)ctx->deviceAlloc(aslm*aslm*sizeof(REAL));
	Oy = (REAL *)ctx->deviceAlloc(aslm*aslm*sizeof(REAL));
	Oz = (REAL *)ctx->deviceAlloc(aslm*aslm*sizeof(REAL));

	printf("DeviceAngularData:\n");
	printf("\tomega     = %p\n", omega);
	printf("\tomega_pos = %p\n", omega_pos);
	printf("\tOx        = %p\n", Ox);
	printf("\tOy        = %p\n", Oy);
	printf("\tOz        = %p\n", Oz);
	printf("\tmem used = %2.6f MB\n", 1e-6 * (6*sizeof(REAL) + 3*sizeof(idx)) * aslm * aslm);
	printf("\tslm = %d, aslm = %d\n", slm, aslm);

	REAL *omega_aligned = new REAL[3*aslm*aslm];
	idx *omega_pos_aligned = new idx[3*aslm*aslm];
	REAL *Ox_aligned = new REAL[aslm*aslm];
	REAL *Oy_aligned = new REAL[aslm*aslm];
	REAL *Oz_aligned = new REAL[aslm*aslm];

	for (idx i=0; i < aslm; i++)
		for (idx j=0; j < aslm; j++) 
			if (i < slm && j < slm) {
				omega_aligned[0*aslm*aslm + i*aslm + j] = host.omega[0*slm*slm + i*slm + j];
				omega_aligned[1*aslm*aslm + i*aslm + j] = host.omega[1*slm*slm + i*slm + j];
				omega_aligned[2*aslm*aslm + i*aslm + j] = host.omega[2*slm*slm + i*slm + j];
				omega_pos_aligned[0*aslm*aslm + i*aslm + j] = host.omega_pos[0*slm*slm + i*slm + j];
				omega_pos_aligned[1*aslm*aslm + i*aslm + j] = host.omega_pos[1*slm*slm + i*slm + j];
				omega_pos_aligned[2*aslm*aslm + i*aslm + j] = host.omega_pos[2*slm*slm + i*slm + j];
				Ox_aligned[i*aslm + j] = host.Ox[i*slm + j];
				Oy_aligned[i*aslm + j] = host.Oy[i*slm + j];
				Oz_aligned[i*aslm + j] = host.Oz[i*slm + j];
			} else {
				omega_aligned[0*aslm*aslm + i*aslm + j] = 0;
				omega_aligned[1*aslm*aslm + i*aslm + j] = 0;
				omega_aligned[2*aslm*aslm + i*aslm + j] = 0;
				omega_pos_aligned[0*aslm*aslm + i*aslm + j] = 0;
				omega_pos_aligned[1*aslm*aslm + i*aslm + j] = 0;
				omega_pos_aligned[2*aslm*aslm + i*aslm + j] = 0;
				Ox_aligned[i*aslm + j] = 0;
				Oy_aligned[i*aslm + j] = 0;
				Oz_aligned[i*aslm + j] = 0;
			}
		
	ctx->copyToDev(omega, omega_aligned, 3*aslm*aslm*sizeof(REAL));
	ctx->copyToDev(omega_pos, omega_pos_aligned, 3*aslm*aslm*sizeof(idx));
	ctx->copyToDev(Ox, Ox_aligned, aslm*aslm*sizeof(REAL));
	ctx->copyToDev(Oy, Oy_aligned, aslm*aslm*sizeof(REAL));
	ctx->copyToDev(Oz, Oz_aligned, aslm*aslm*sizeof(REAL));

	delete[] omega_aligned;
	delete[] omega_pos_aligned;
	delete[] Ox_aligned;
	delete[] Oy_aligned;
	delete[] Oz_aligned;
}

DeviceAngularData::~DeviceAngularData() {
	ctx->deviceFree(omega);
	ctx->deviceFree(omega_pos);
	ctx->deviceFree(Ox);
	ctx->deviceFree(Oy);
	ctx->deviceFree(Oz);
}
