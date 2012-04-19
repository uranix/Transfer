#include <MeshData.h>

DeviceMeshData::DeviceMeshData(const MeshData &host) {
	int nps = (int)(sqrt(host.nP)+0.5000001);
	nPlow = align_power(nps, 16);
	nPhigh = align_power(nps, 16);
	nP = host.nP;

	cudaMalloc((void **)&tetstart, (nP+1)*sizeof(idx));
	cudaMalloc((void **)&tetidx, host.tetstart[host.nP]*sizeof(idx));
	cudaMalloc((void **)&tetpos, host.tetstart[host.nP]*sizeof(idx));
	cudaMalloc((void **)&mesh, host.nT*sizeof(tetrahedron));

	cudaMalloc((void **)&facestart, (nP+1)*sizeof(idx));
	cudaMalloc((void **)&faceidx, host.facestart[host.nP]*sizeof(idx));
	cudaMalloc((void **)&facepos, host.facestart[host.nP]*sizeof(idx));
	cudaMalloc((void **)&bnd, host.nF*sizeof(face));

	dim3 grid(nPlow, nPhigh);
	dim3 block(1);
	
	cudaMemcpy(tetstart, host.tetstart, (host.nP+1)*sizeof(idx), cudaMemcpyHostToDevice);

	cudaMemcpy(facestart, host.facestart, (host.nP+1)*sizeof(idx), cudaMemcpyHostToDevice);

	cudaMemcpy(tetidx, host.tetidx, host.tetstart[host.nP]*sizeof(idx), cudaMemcpyHostToDevice);
	cudaMemcpy(tetpos, host.tetpos, host.tetstart[host.nP]*sizeof(idx), cudaMemcpyHostToDevice);
	cudaMemcpy(mesh, host.mesh, host.nT*sizeof(tetrahedron), cudaMemcpyHostToDevice);

	cudaMemcpy(faceidx, host.faceidx, host.facestart[host.nP]*sizeof(idx), cudaMemcpyHostToDevice);
	cudaMemcpy(facepos, host.facepos, host.facestart[host.nP]*sizeof(idx), cudaMemcpyHostToDevice);
	cudaMemcpy(bnd, host.bnd, host.nF*sizeof(face), cudaMemcpyHostToDevice);
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
