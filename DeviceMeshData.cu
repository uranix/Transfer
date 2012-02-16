#include <MeshData.h>

__global__ void copy_start(idx *dst, idx *src, idx nP) {
	idx i = blockIdx.x + blockIdx.y * gridDim.x;
	
	dst[i] = src[i>nP ? nP : i];
}

DeviceMeshData::DeviceMeshData(const MeshData host) {
	int nps = (int)(sqrt(nP)+0.5000001);
	nPlow = align_power(nps, 16);
	nPhigh = align_power(nps, 16);
	int nP = nPlow * nPhigh;

	idx *tmp;
	cudaMalloc((void **)&tmp, (host.nP+1)*sizeof(idx));

	cudaMalloc((void **)&tetstart, (nP+1)*sizeof(idx));
	cudaMalloc((void **)&tetidx, host.tetstart[host.nP]*sizeof(idx));
	cudaMalloc((void **)&tetpos, host.tetstart[host.nP]*sizeof(idx));
	cudaMalloc((void **)&mesh, host.nT*sizeof(tetrahedron));

	cudaMalloc((void **)&facestart, (nP+1)*sizeof(idx));
	cudaMalloc((void **)&faceidx, host.facestart[host.nP]*sizeof(idx));
	cudaMalloc((void **)&facepos, host.facestart[host.nP]*sizeof(idx));
	cudaMalloc((void **)&bnd, host.nF*sizeof(face));

	dim3 grid(nPlow, nPhigh);
	
	cudaMemcpy(tmp, host.tetstart, (host.nP+1)*sizeof(idx), cudaMemcpyHostToDevice);
	copy_start<<<grid,block>>>(tetstart, tmp, host.nP);

	cudaMemcpy(tmp, host.facestart, (host.nP+1)*sizeof(idx), cudaMemcpyHostToDevice);
	copy_start<<<grid,block>>>(facestart, tmp, host.nP);

	cudaMemcpy(tetidx, host.tetidx, host.tetstart[host.nP]*sizeof(idx));
	cudaMemcpy(tetpos, host.tetpos, host.tetstart[host.nP]*sizeof(idx));
	cudaMemcpy(mesh, host.mesh, host.nT*sizeof(tetrahedron));

	cudaMemcpy(faceidx, host.faceidx, host.facestart[host.nP]*sizeof(idx));
	cudaMemcpy(facepos, host.facepos, host.facestart[host.nP]*sizeof(idx));
	cudaMemcpy(bnd, host.bnd, host.nF*sizeof(face));

	cudaFree(tmp);
}

DeviceMeshData::~DeviceMeshData() {
	cudaFree(tetstart);
	cudaFree(tetidx);
	cudaFree(tetpos);
	cudaFree(mesh);

	cudaFree(facestart);
	cudaFree(faceidx);
	cudaFree(facepos);
	cudaFree(bnd);
}
