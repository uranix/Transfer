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

	tetstart	= (idx *)ctx->deviceAlloc((nP+1)*sizeof(idx));
	tetidx		= (idx *)ctx->deviceAlloc(host.tetstart[host.nP]*sizeof(idx));
	tetpos		= (idx *)ctx->deviceAlloc(host.tetstart[host.nP]*sizeof(idx));
	mesh		= (tetrahedron *)ctx->deviceAlloc(host.nT*sizeof(tetrahedron));

	facestart	= (idx *)ctx->deviceAlloc((nP+1)*sizeof(idx));
	faceidx		= (idx *)ctx->deviceAlloc(host.facestart[host.nP]*sizeof(idx));
	facepos		= (idx *)ctx->deviceAlloc(host.facestart[host.nP]*sizeof(idx));
	bnd			= (face *)ctx->deviceAlloc(host.nF*sizeof(face));

	printf("DeviceMeshData:\n");
	printf("\tmesh      = %p\n",mesh);
	printf("\ttetstart  = %p\n",tetstart);
	printf("\ttetidx    = %p\n",tetidx);
	printf("\ttetpos    = %p\n",tetpos);
	printf("\tbnd       = %p\n",bnd);
	printf("\tfacestart = %p\n",facestart);
	printf("\tfaceidx   = %p\n",faceidx);
	printf("\tfacepos   = %p\n",facepos);

	ctx->copyToDev(tetstart, host.tetstart, (host.nP+1)*sizeof(idx));
	ctx->copyToDev(tetidx, host.tetidx, host.tetstart[host.nP]*sizeof(idx));
	ctx->copyToDev(tetpos, host.tetpos, host.tetstart[host.nP]*sizeof(idx));
	ctx->copyToDev(mesh, host.mesh, host.nT*sizeof(tetrahedron));

	ctx->copyToDev(facestart, host.facestart, (host.nP+1)*sizeof(idx));
	ctx->copyToDev(faceidx, host.faceidx, host.facestart[host.nP]*sizeof(idx));
	ctx->copyToDev(facepos, host.facepos, host.facestart[host.nP]*sizeof(idx));
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
