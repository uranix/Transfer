#ifndef __MESH3D__MESH_H__
#define __MESH3D__MESH_H__

#include <math.h>
#include <stdio.h>

#include "vector.h"
#include "list.h"

#include <new>

enum ElementType /*: int*/ {
	EL_TETRAHEDRON = 1,
	EL_TYPE_END
};

enum FaceType /*: int*/ {
	FC_TRIANGLE = 1,
	FC_TYPE_END
};

struct Element;
struct Face;

struct Vertex {
	int index;
	Vector r;
	Node<Element *> *elems;
	Node<int> *elidx;
	int elnum;
	Node<Face *> *bnds;
	Node<int> *bndidx;
	int bndnum;
	Vertex(int _index, double _x, double _y, double _z): index(_index), r(_x, _y, _z) { 
		elems = 0;
		elidx = 0;
		bnds = 0;
		bndidx = 0;
		elnum = 0;
		bndnum = 0;
	}
	Vertex(int _index, const Vector &v): index(_index), r(v) { }
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
	int bnd_index;
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
	Node<char *> *memory;
	Node<Element *> *vert2elem;
	Node<Face *> *vert2bnd;
	Node<int> *vert2elidx;
	Node<int> *vert2bndidx;
	void fromVol(int nV, int nB, int nT, double *vert, int *bnd, int *tet, int *bndmat, int *tetmat);
public:
	int nVert, nFaces, nElems, nBndFaces;
	Vertex **vertices;
	Element **elements;
	Face **faces;

	Mesh(const char *fn);
	void saveVtk(const char *fn, int realbytes, const char *cellfmt, const char *pointfmt, ...);
	bool check();
	double quality();
	~Mesh();
};

#endif 
