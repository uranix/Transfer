#ifndef __MESHDATA_H__
#define __MESHDATA_H__

#include "util.h"
#include "tetrahedron.h"
#include "face.h"

struct MeshData {
	idx nP, nT, nF;

	idx *tetstart;
	idx *tetidx;
	idx *tetpos;
	tetrahedron *mesh;

	idx *facestart;
	idx *faceidx;
	idx *facepos;
	face *bnd;
	MeshData(const MeshData &);
	~MeshData();
};

struct DeviceMeshData {
	idx nPlow, nPhigh; // nP < nPlow * nPhigh

	idx *tetstart;
	idx *tetidx;
	idx *tetpos;
	tetrahedron *mesh;

	idx *facestart;
	idx *faceidx;
	idx *facepos;
	face *bnd;
	DeviceMeshData(const MeshData &);
	~DeviceMeshData();
private:
	DeviceMeshData();
};

#endif
