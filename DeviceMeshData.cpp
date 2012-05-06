#include "MeshData.h"
#include "CudaContext.h"
#include "common.cuh"

#include <stdio.h>
#include <math.h>

DeviceMeshData::DeviceMeshData(const CudaContext *_ctx, const MeshData &host) {
	ctx = _ctx;
	int nps = (int)(sqrt((double)host.nP)+0.5000001);
	nPlow = align_power(nps, 16);
	nPhigh = align_power(nps, 16);
	nP = host.nP;

	idx tetcnt = host.tetstart[nP]; 
	idx facecnt = host.facestart[nP]; 

	tetstart	= (idx *)ctx->deviceAlloc((nP+1)*sizeof(idx));
	tetidx		= (idx *)ctx->deviceAlloc(tetcnt*sizeof(idx));
	tetpos		= (idx *)ctx->deviceAlloc(tetcnt*sizeof(idx));
	mesh		= (tetrahedron *)ctx->deviceAlloc(host.nT*sizeof(tetrahedron));

	facestart	= (idx *)ctx->deviceAlloc((nP+1)*sizeof(idx));
	faceidx		= (idx *)ctx->deviceAlloc(facecnt*sizeof(idx));
	facepos		= (idx *)ctx->deviceAlloc(facecnt*sizeof(idx));
	bnd			= (face *)ctx->deviceAlloc(host.nF*sizeof(face));

	idx totalmem = (2*nP + 2 + 2 * tetcnt + 2 * facecnt) * sizeof(idx) + 
		host.nT * sizeof(tetrahedron) + host.nF * sizeof(face);

	printf("DeviceMeshData:\n");
	printf("\tmesh      = %p\n",mesh);
	printf("\ttetstart  = %p\n",tetstart);
	printf("\ttetidx    = %p\n",tetidx);
	printf("\ttetpos    = %p\n",tetpos);
	printf("\tbnd       = %p\n",bnd);
	printf("\tfacestart = %p\n",facestart);
	printf("\tfaceidx   = %p\n",faceidx);
	printf("\tfacepos   = %p\n",facepos);
	printf("\tmem used = %2.6f MB\n", 1e-6 * totalmem);

	ctx->copyToDev(tetstart, host.tetstart, (nP+1)*sizeof(idx));
	ctx->copyToDev(tetidx, host.tetidx, tetcnt*sizeof(idx));
	ctx->copyToDev(tetpos, host.tetpos, tetcnt*sizeof(idx));
	ctx->copyToDev(mesh, host.mesh, host.nT*sizeof(tetrahedron));

	ctx->copyToDev(facestart, host.facestart, (nP+1)*sizeof(idx));
	ctx->copyToDev(faceidx, host.faceidx, facecnt*sizeof(idx));
	ctx->copyToDev(facepos, host.facepos, facecnt*sizeof(idx));
	ctx->copyToDev(bnd, host.bnd, host.nF*sizeof(face));
}

DeviceMeshData::~DeviceMeshData() {
	ctx->deviceFree(tetstart);
	ctx->deviceFree(tetidx);
	ctx->deviceFree(tetpos);
	ctx->deviceFree(mesh);

	ctx->deviceFree(facestart);
	ctx->deviceFree(faceidx);
	ctx->deviceFree(facepos);
	ctx->deviceFree(bnd);
}
