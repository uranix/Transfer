#include "MeshData.h"
#include <stdio.h>

#ifndef _
#define _(x) do { \
	if ((x) != cudaSuccess) { \
	fprintf(stderr, "File %s line %d, %s failed with error `%s'\n", __FILE__, __LINE__, #x, cudaGetErrorString(cudaGetLastError())); \
	fflush(stderr); } \
} while (0)
#endif

DeviceMeshData::DeviceMeshData(const MeshData &host) {
	int nps = (int)(sqrt(host.nP)+0.5000001);
	nPlow = align_power(nps, 16);
	nPhigh = align_power(nps, 16);
	nP = host.nP;

	_(cudaMalloc((void **)&tetstart, (nP+1)*sizeof(idx)));
	_(cudaMalloc((void **)&tetidx, host.tetstart[host.nP]*sizeof(idx)));
	_(cudaMalloc((void **)&tetpos, host.tetstart[host.nP]*sizeof(idx)));
	_(cudaMalloc((void **)&mesh, host.nT*sizeof(tetrahedron)));

	_(cudaMalloc((void **)&facestart, (nP+1)*sizeof(idx)));
	_(cudaMalloc((void **)&faceidx, host.facestart[host.nP]*sizeof(idx)));
	_(cudaMalloc((void **)&facepos, host.facestart[host.nP]*sizeof(idx)));
	_(cudaMalloc((void **)&bnd, host.nF*sizeof(face)));

	printf("DeviceMeshData:\n");
	printf("\tmesh		= %p\n",mesh);
	printf("\ttetstart	= %p\n",tetstart);
	printf("\ttetidx	= %p\n",tetidx);
	printf("\ttetpos	= %p\n",tetpos);
	printf("\tbnd		= %p\n",bnd);
	printf("\tfacestart = %p\n",facestart);
	printf("\tfaceidx 	= %p\n",faceidx);
	printf("\tfacepos	= %p\n",facepos);

	_(cudaMemcpy(tetstart, host.tetstart, (host.nP+1)*sizeof(idx), cudaMemcpyHostToDevice));

	_(cudaMemcpy(facestart, host.facestart, (host.nP+1)*sizeof(idx), cudaMemcpyHostToDevice));

	_(cudaMemcpy(tetidx, host.tetidx, host.tetstart[host.nP]*sizeof(idx), cudaMemcpyHostToDevice));
	_(cudaMemcpy(tetpos, host.tetpos, host.tetstart[host.nP]*sizeof(idx), cudaMemcpyHostToDevice));
	_(cudaMemcpy(mesh, host.mesh, host.nT*sizeof(tetrahedron), cudaMemcpyHostToDevice));

	_(cudaMemcpy(faceidx, host.faceidx, host.facestart[host.nP]*sizeof(idx), cudaMemcpyHostToDevice));
	_(cudaMemcpy(facepos, host.facepos, host.facestart[host.nP]*sizeof(idx), cudaMemcpyHostToDevice));
	_(cudaMemcpy(bnd, host.bnd, host.nF*sizeof(face), cudaMemcpyHostToDevice));
}

DeviceMeshData::~DeviceMeshData() {
	printf("deleting\n");
	_(cudaFree(tetstart));
	_(cudaFree(tetidx));
	_(cudaFree(tetpos));
	_(cudaFree(mesh));

	_(cudaFree(facestart));
	_(cudaFree(faceidx));
	_(cudaFree(facepos));
	_(cudaFree(bnd));
}
