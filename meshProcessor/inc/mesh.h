#ifndef __MESH_H__
#define __MESH_H__

#include <math.h>
#include <stdio.h>

#include "vector.h"
#include "list.h"
#include "front.h"

#include "libaft.h"
#include "libmba.h"

#if defined(_MSC_VER)
 #include <new.h>
#endif

enum ElementType /*: int*/ {
	EL_TETRAHEDRON = 1,
	EL_TYPE_END
};

enum FaceType /*: int*/ {
	FC_TRIANGLE = 1,
	FC_TYPE_END
};

struct Element;

struct Vertex {
	int index;
	Vector r;
	Node<Element *> *elems;
	Vertex(int _index, double _x, double _y, double _z): r(_x, _y, _z), index(_index) { }
	Vertex(int _index, const Vector &v): r(v), index(_index) { }
private:
	Vertex &operator=(const Vertex &p);
	Vertex(const Vertex &p);
};

struct Element {
	int index;
	Vector center;
	double volume;
	ElementType type;
	int region;
	double quality;

	virtual ~Element() {}
protected:
	Element(int _index) : index(_index) {}
private:
	Element();
	Element(const Element &);
};

struct Face {
	int index;
	Element *element;
	Face *flip;

	FaceType type;
	int borderType;
	Vector normal;
	Vector center;
	double surface;	
	virtual ~Face() {}
	virtual void setFlip(Face *_flip) {
		flip = _flip;
	}
protected:
	Face(int _index, Element *_element) : index(_index), element(_element) {}
private:
	Face();
	Face(const Face &);
};

class Mesh {
	int nVert, nFaces, nElems;
	Node<char *> *memory;
	Vertex **vertices;
	Element **elements;
	Face **faces;
	void fromAft(int nV, int nB, int nT, double *vert, int *bnd, int *tet, int *bndmat, int *tetmat);
	static double meshSize;
	static double usersize(double, double, double) {
		return meshSize;
	}
	static int euclid_metric(double *, double *, double *, double *metric) {
		metric[0] = 1.;
		metric[4] = 1.;
		metric[8] = 1.;
		metric[3] = 0.;
		metric[6] = 0.;
		metric[7] = 0.;
		return 1;
	}
public:
	/* Construct from ani format */
	/* @deprecated */
	Mesh(int nV, int nF, int nT, double *vert, int *bnd, int *tet, int *bndmat, int *tetmat);
	Mesh(char *fn);
	Mesh(const Front *f, double meshsize, bool fixShape, int nnV = 400000, int nnF = 1000000, int nnT = 400000);
	void saveVtk(char *fn);
	void saveBmf(char *fn);
	bool check();
	double quality();
	~Mesh();
};

#endif 