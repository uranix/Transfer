#ifndef __FRONT_H__
#define __FRONT_H__

#include "libscg.h"

#include "mesh.h"

#include <string.h>

class Front {
	int nVert, nFaces;
	int nVertMax, nFacesMax;
	double *vertex;
	int *face;
	Front *intersect_op(const Front *a, const Front *b, int op) const {
		Front *r = new Front(2*(a->nVert + b->nVert), 2*(a->nFaces + b->nFaces)); // i hope that would be enough
		scg_intersect_mesh(	a->nVert, a->nFaces, a->vertex, a->face, 
							b->nVert, b->nFaces, b->vertex, b->face, 
							&r->nVert, &r->nFaces, r->vertex, r->face, op);
		r->truncate();
		return r;
	}
	void truncate() {
		double *nv = new double[3*nVert];
		int *nf = new int[3*nFaces];
		memcpy(nv, vertex, 3*nVert*sizeof(double));
		memcpy(nf, face, 3*nFaces*sizeof(int));
		delete[] vertex;
		delete[] face;
		nFacesMax = nFaces;
		nVertMax = nVert;
		vertex = nv;
		face = nf;
	}
	friend class Mesh;
public:
	Front(int _nVertMax = 1000000, int _nFacesMax = 1000000 ) {
		nVertMax = _nVertMax;
		nFacesMax = _nFacesMax;
		vertex = new double[3*nVertMax];
		face = new int[3*nFacesMax];
		nVert = nFaces = 0;
	}
	~Front() {
		delete[] face;
		delete[] vertex;
	}
	void translate(double x, double y, double z) {
		scg_translate(nVert, vertex, x, y, z);
	}
	void scale(double mx, double my, double mz) {
		scg_scale(nVert, vertex, mx, my, mz);
	}
	void rotate(double yaw, double pitch, double roll) {
		scg_rotate(nVert, vertex, yaw, pitch, roll);
	}
	Front *unite(const Front *b) const {
		return intersect_op(this, b, 0);
	}
	Front *sub(const Front *b) const {
		return intersect_op(this, b, 2);
	}
	Front *isect(const Front *b) const {
		return intersect_op(this, b, 3);
	}
	void makeSphere(double radius, double meshsize) {
		scg_make_sphere(&nVert, &nFaces, vertex, face, radius, meshsize, nVertMax, nFacesMax);
		truncate();
	}
	void makeBox(double x, double y, double z, double meshsize) {
		scg_make_paral(&nVert, &nFaces, vertex, face, x, y, z, meshsize, nVertMax, nFacesMax);
		truncate();
	}
	void makeCylinder(double radius, double height, double meshsize) {
		scg_make_cylinder(&nVert, &nFaces, vertex, face, radius, height, meshsize, nVertMax, nFacesMax);
		truncate();
	}
};

#endif