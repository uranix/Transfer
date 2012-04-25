#ifndef __MESHDATA_H__
#define __MESHDATA_H__

#include "util.h"
#include "tetrahedron.h"
#include "face.h"
#include "Config.h"

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
	DeviceMeshData(const MeshData &);
	~DeviceMeshData();
private:
	DeviceMeshData();
};

#endif
