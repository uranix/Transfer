#include "AngularData.h"

#include "common.cuh"

#include <stdio.h>

#ifndef _
#define _(x) do { \
	if ((x) != cudaSuccess) { \
	fprintf(stderr, "File %s line %d, %s failed with error `%s'\n", __FILE__, __LINE__, #x, cudaGetErrorString(cudaGetLastError())); \
	fflush(stderr); } \
} while (0)
#endif

DeviceAngularData::DeviceAngularData(const AngularData &host) {
	slm = host.slm;
	aslm = align_power(host.slm, COALESCED_NUM(REAL));
	_(cudaMalloc((void **)&omega, 3*aslm*aslm*sizeof(REAL)));
	_(cudaMalloc((void **)&omega_pos, 3*aslm*aslm*sizeof(idx)));
	_(cudaMalloc((void **)&Ox, aslm*aslm*sizeof(REAL)));
	_(cudaMalloc((void **)&Oy, aslm*aslm*sizeof(REAL)));
	_(cudaMalloc((void **)&Oz, aslm*aslm*sizeof(REAL)));

	printf("DeviceAngularData:\n");
	printf("\tomega     = %p\n", omega);
	printf("\tomega_pos = %p\n", omega_pos);
	printf("\tOx        = %p\n", Ox);
	printf("\tOy        = %p\n", Oy);
	printf("\tOz        = %p\n", Oz);

	dim3 block;
	REAL *omega_aligned = new REAL[3*aslm*aslm];
	idx *omega_pos_aligned = new idx[3*aslm*aslm];
	idx *Ox_aligned = new idx[aslm*aslm];
	idx *Oy_aligned = new idx[aslm*aslm];
	idx *Oz_aligned = new idx[aslm*aslm];

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
		
	_(cudaMemcpy(omega, omega_aligned, 3*aslm*aslm*sizeof(REAL), cudaMemcpyHostToDevice));
	_(cudaMemcpy(omega_pos, omega_pos_aligned, 3*aslm*aslm*sizeof(idx), cudaMemcpyHostToDevice));
	_(cudaMemcpy(Ox, Ox_aligned, aslm*aslm*sizeof(REAL), cudaMemcpyHostToDevice));
	_(cudaMemcpy(Oy, Oy_aligned, aslm*aslm*sizeof(REAL), cudaMemcpyHostToDevice));
	_(cudaMemcpy(Oz, Oz_aligned, aslm*aslm*sizeof(REAL), cudaMemcpyHostToDevice));
}

DeviceAngularData::~DeviceAngularData() {
	_(cudaFree(omega));
	_(cudaFree(omega_pos));
	_(cudaFree(Ox));
	_(cudaFree(Oy));
	_(cudaFree(Oz));
}
