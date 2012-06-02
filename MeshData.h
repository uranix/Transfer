#ifndef __MESHDATA_H__
#define __MESHDATA_H__

#include "util.h"
#include "tetrahedron.h"
#include "face.h"
#include "Config.h"
#include "AngularData.h"

struct CudaContext;
class Mesh;

struct MeshData {
	idx nP, nF, nT;

	idx *tetstart;
	idx *tetidx;
	idx *tetpos;
	tetrahedron *mesh;

	idx *facestart;
	idx *faceidx;
	idx *facepos;
	face *bnd;
	MeshData(const Config &);
	void ComputeMoments(const AngularData &ad, const REAL * const I[], REAL *U, REAL *T[3][3]);
	void ComputeFlux(	REAL *T[3][3], 
						REAL *W[3]);
	~MeshData();

	Mesh *_m;
};

struct DeviceMeshDataRaw {
	idx nPlow, nPhigh; // nP < nPlow * nPhigh
	idx nP;

	idx *tetstart;
	idx *tetidx;
	idx *tetpos;
	tetrahedron *mesh;

	idx *facestart;
	idx *faceidx;
	idx *facepos;
	face *bnd;
};

struct DeviceMeshData : public DeviceMeshDataRaw { 
	const CudaContext *ctx;
	DeviceMeshData(const CudaContext *ctx, const MeshData &);
	~DeviceMeshData();
private:
	DeviceMeshData();
};

#endif
