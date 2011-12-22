#include "mesh.h"

// OK. There is no way to make this non-static. Thanks to ani3D team
// Don't try to mesh something in parallel or so
double Mesh::meshSize = 1;

#include "tri_face.h"
#include "tetrahedron.h"
#include "list.h"

#include "stdint.h"
#include <stdlib.h>

void Mesh::fromAft(int nV, int nB, int nT, double *vert, int *bnd, int *tet, int *bndmat, int *tetmat) {
	nVert = nV;
	nElems = nT;
	nFaces = 4*nT + nB;
	
	memory = 0;

	Vertex *vertplace = (Vertex *)new char[sizeof(Vertex) * nV];
	Tetrahedron *tetraplace = (Tetrahedron *)new char[sizeof(Tetrahedron) * nT];
	TriFace *faceplace = (TriFace *)new char[sizeof(TriFace) * nFaces];

	addHead<char *>((char *)vertplace, &memory);
	addHead<char *>((char *)tetraplace, &memory);
	addHead<char *>((char *)faceplace, &memory);

	Node<int> *vert2face, *mem = list_allocate<int>(3*nFaces);
	vert2face = mem;

	vertices = new Vertex *[nVert];
	elements = new Element*[nElems];
	faces = new Face*[nFaces];

	Node<int> **lists = new Node<int> *[nV];
	
	for (int i = 0; i < nV; i++) {
		vertices[i] = new (vertplace + i) Vertex(i, vert[3*i], vert[3*i+1], vert[3*i+2]);
	}

	for (int i = 0; i < nV; i++)
		lists[i] = 0;
	for (int i = 0; i < nT; i++) {
		elements[i] = new (tetraplace + i) Tetrahedron(i, vertices[tet[4*i + 0]], vertices[tet[4*i + 1]], vertices[tet[4*i + 2]], vertices[tet[4*i + 3]], tetmat[i], faceplace+4*i);
		for (int j = 0; j < 4; j++) {
			faces[4*i+j] = faceplace + 4*i+j;
			for (int k = 0; k < 4; k++)
				 if (j != k)
					addHead(4*i + j, &lists[tet[4*i + k]], &mem);
		}
	}
	for (int i = 0; i < nB; i++) {
		faces[4*nT + i] = new(faceplace + 4*nT + i) TriFace(4*nT + i, vertices[bnd[3*i+0]], vertices[bnd[3*i+1]], vertices[bnd[3*i+2]], 0, bndmat[i]);
		for (int j = 0; j < 3; j++)
			addHead(4*nT + i, &lists[bnd[3*i+j]], &mem);
	}
	for (int i = 0; i < nFaces; i++) {
		TriFace *f = (TriFace *)faces[i];
		Vertex *p1 = f->p[0], *p2 = f->p[1], *p3 = f->p[2];
		int i1 = p1->index;
		int i2 = p2->index;
		int i3 = p3->index;
		/* TODO List -> Tree : O(m^2) -> O(m log m)*/
		int found = -1;
		for (Node<int> *n1 = lists[i1]; n1; n1 = n1->next)
			if (n1->data != i) 
				for (Node<int> *n2 = lists[i2]; n2; n2 = n2->next)
					if (n2->data == n1->data)
						for (Node<int> *n3 = lists[i3]; n3; n3 = n3->next)
							if (n3->data == n1->data) {
								found = n1->data;
								goto out;
							}
out:
		if (found == -1)
			throw "Flipped face not found";
		f->setFlip(faces[found]);
	}
	delete[] lists;
	list_deallocate(vert2face);
}

Mesh::Mesh(int nV, int nB, int nT, double *vert, int *bnd, int *tet, int *bndmat, int *tetmat) {
	fromAft(nV, nB, nT, vert, bnd, tet, bndmat, tetmat);
}

void Mesh::saveVtk(char *fn) {
	FILE *f;

	f=fopen(fn, "w");
	fprintf(f, "# vtk DataFile Version 2.0\n");
	fprintf(f, "Generated output\n");
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(f, "POINTS %9d double\n", nVert);
	for (int i=0; i<nVert; i++)
		fprintf(f, "% 2.10e % 2.10e % 2.10e\n", vertices[i]->r.x, vertices[i]->r.y, vertices[i]->r.z);
	fprintf(f, "\nCELLS %9d %9d\n", nElems, 5*nElems);
	for (int i=0; i<nElems; i++)
		fprintf(f, "4 %6d %6d %6d %6d\n",	((Tetrahedron *)elements[i])->p[0]->index, ((Tetrahedron *)elements[i])->p[1]->index, 
											((Tetrahedron *)elements[i])->p[2]->index, ((Tetrahedron *)elements[i])->p[3]->index);
	fprintf(f, "\nCELL_TYPES %9d\n", nElems);
	for (int i=0; i<nElems; i++)
		fprintf(f, "10 %s", (i&15)==15?"\n":"");

	fprintf(f, "\nCELL_DATA %9d\n", nElems);
	fprintf(f, "\nFIELD MerticData 5\n");
	fprintf(f, "\nVolume 1 %9d double\n", nElems);
	for (int i=0; i<nElems; i++)
		fprintf(f, "%2.10e %s", elements[i]->volume, (i&3)==3?"\n":"");
	fprintf(f, "\nQuality 1 %9d double\n", nElems);
	for (int i=0; i<nElems; i++)
		fprintf(f, "%2.10e %s", elements[i]->quality, (i&3)==3?"\n":"");
	fprintf(f, "\nCenterX 1 %9d double\n", nElems);
	for (int i=0; i<nElems; i++)
		fprintf(f, "%2.10e\n", elements[i]->center.x);
	fprintf(f, "\nCenterY 1 %9d double\n", nElems);
	for (int i=0; i<nElems; i++)
		fprintf(f, "%2.10e\n", elements[i]->center.y);
	fprintf(f, "\nCenterZ 1 %9d double\n", nElems);
	for (int i=0; i<nElems; i++)
		fprintf(f, "%2.10e\n", elements[i]->center.z);

	fprintf(f, "\nPOINT_DATA %9d\n", nVert);
	fprintf(f, "\nFIELD MerticData 1\n");
	fprintf(f, "\u 1 %9d double\n", nVert);
	for (int i=0; i<nVert; i++)
		fprintf(f, "%2.10e\n", vertices[i]->r.norm());
	fclose(f);
}

void Mesh::saveBmf(char *fn) {
	FILE *f = fopen(fn, "wb");
	bool bits32 = sizeof(int) == 4;
	if (bits32) {
		int magic = *(int *)"BMF";
		fwrite(&magic, sizeof(magic), 1, f);
	} else {
		int magic = *(int *)"BMF64\0\0";
		fwrite(&magic, sizeof(magic), 1, f);
	}
	fwrite(&nVert, sizeof(nVert), 1, f);
	fwrite(&nFaces, sizeof(nFaces), 1, f);
	fwrite(&nElems, sizeof(nElems), 1, f);
	for (int i=0; i<nVert; i++)
		vertices[i]->r.serialize(f);
	int facecount[FC_TYPE_END] = {0};
	for (int i=0; i<nFaces; i++)
		facecount[faces[i]->type]++;
	facecount[0] = FC_TYPE_END;
	fwrite(facecount, sizeof(int), FC_TYPE_END, f);
	for (int i=0; i<nFaces; i++) {
		fwrite(&faces[i]->type, sizeof(faces[i]->type), 1, f);
		switch (faces[i]->type) {
			case FC_TRIANGLE: {
				TriFace *face = (TriFace *)faces[i];
				int p1 = face->p[0]->index,
					p2 = face->p[1]->index,
					p3 = face->p[2]->index,
					flip = face->flip->index,
					elem = face->element ? face->element->index : -face->borderType;
				fwrite(&p1, sizeof(p1), 1, f);
				fwrite(&p2, sizeof(p2), 1, f);
				fwrite(&p3, sizeof(p3), 1, f);
				fwrite(&flip, sizeof(flip), 1, f);
				fwrite(&elem, sizeof(elem), 1, f);
			} break;
			default:
				throw "Face type not yet implemented";
		}
	}
	int elemcount[EL_TYPE_END] = {0};
	for (int i=0; i<nElems; i++)
		elemcount[elements[i]->type]++;
	elemcount[0] = EL_TYPE_END;
	fwrite(elemcount, sizeof(int), EL_TYPE_END, f);
	for (int i=0; i<nElems; i++) {
		fwrite(&elements[i]->type, sizeof(elements[i]->type), 1, f);
		switch (elements[i]->type) {
			case EL_TETRAHEDRON: {
				Tetrahedron *elem = (Tetrahedron *)elements[i];
				int p1 = elem->p[0]->index,
					p2 = elem->p[1]->index,
					p3 = elem->p[2]->index,
					p4 = elem->p[3]->index,
					f1 = elem->f[0]->index,
					f2 = elem->f[1]->index,
					f3 = elem->f[2]->index,
					f4 = elem->f[3]->index,
					region = elem->region;
				fwrite(&p1, sizeof(p1), 1, f);
				fwrite(&p2, sizeof(p2), 1, f);
				fwrite(&p3, sizeof(p3), 1, f);
				fwrite(&p4, sizeof(p4), 1, f);
				
				fwrite(&f1, sizeof(f1), 1, f);
				fwrite(&f2, sizeof(f2), 1, f);
				fwrite(&f3, sizeof(f3), 1, f);
				fwrite(&f4, sizeof(f4), 1, f);

				fwrite(&region, sizeof(region), 1, f);
			} break;
			default:
				throw "Element type not yet implemented";
		}
	}
	fclose(f);
}

Mesh::Mesh(char *fn) {
	FILE *f = fopen(fn, "rb");
	bool bits32 = sizeof(int) == 4;
	if (bits32) {
		char sig[4];
		fread(sig, 4, 1, f);
		if (*(int *)sig != *(int *)"BMF")
			throw "Mesh format error";
	} else {
		char sig[8];
		fread(sig, 8, 1, f);
		if (*(int *)sig != *(int *)"BMF64\0\0") {
			if (*(int32_t *)sig != *(int32_t *)"BMF")
				throw "32-bit mesh is unsupported on 64-bit platform";
			else
				throw "Mesh format error";
		}
	}
	fread(&nVert, sizeof(nVert), 1, f);
	fread(&nFaces, sizeof(nFaces), 1, f);
	fread(&nElems, sizeof(nElems), 1, f);

	vertices = new Vertex *[nVert];
	elements = new Element*[nElems];
	faces = new Face*[nFaces];

	memory = 0;

	Vertex *vertplace = (Vertex *)new char[sizeof(Vertex) * nVert];
	addHead<char *>((char *)vertplace, &memory);

	for (int i=0; i<nVert; i++) {
		Vector r(f);
		vertices[i] = new(vertplace + i) Vertex(i, r);
	}
	int facecount[FC_TYPE_END];
	fread(facecount, sizeof(int), 1, f);
	if (facecount[0] > FC_TYPE_END) 
		throw "Unknown face types";
	fread(facecount+1, sizeof(int), facecount[0]-1, f);

	char *faceplace[FC_TYPE_END] = {0};

	for (int i=0; i<nFaces; i++) {
		int type, flip;
		fread(&type, sizeof(type), 1, f);
		switch (type) {
			case FC_TRIANGLE: {
				if (!faceplace[type]) {
					faceplace[type] = new char[sizeof(TriFace) * facecount[type]];
					addHead<char *>((char *)faceplace[type], &memory);
				}
				TriFace *face = (TriFace *)faceplace[type];
				faceplace[type] += sizeof(TriFace);
				int p1, p2,	p3, elem;
				fread(&p1, sizeof(p1), 1, f);
				fread(&p2, sizeof(p2), 1, f);
				fread(&p3, sizeof(p3), 1, f);
				fread(&flip, sizeof(flip), 1, f);
				fread(&elem, sizeof(elem), 1, f);
				faces[i] = new(face) TriFace(i, vertices[p1], vertices[p2], vertices[p3], 0, elem<0?-elem:-1);
			} break;
			default:
				throw "Face type not yet implemented";
		}
		if (flip < i) {
			faces[i]->setFlip(faces[flip]);
			faces[flip]->setFlip(faces[i]);
		}
	}

	int elemcount[EL_TYPE_END];
	fread(elemcount, sizeof(int), 1, f);
	if (elemcount[0] > EL_TYPE_END) 
		throw "Unknown element types";
	fread(elemcount+1, sizeof(int), elemcount[0]-1, f);

	char *elemplace[EL_TYPE_END] = {0};

	for (int i=0; i<nElems; i++) {
		int type;
		fread(&type, sizeof(type), 1, f);

		switch (type) {
			case EL_TETRAHEDRON: {
				if (!elemplace[type]) {
					elemplace[type] = new char[sizeof(Tetrahedron) * elemcount[type]];
					addHead<char *>((char *)elemplace[type], &memory);
				}
				Tetrahedron *elem = (Tetrahedron *)elemplace[type];
				elemplace[type] += sizeof(Tetrahedron);
				int p1, p2, p3, p4, f1, f2, f3, f4, region;

				fread(&p1, sizeof(p1), 1, f);
				fread(&p2, sizeof(p2), 1, f);
				fread(&p3, sizeof(p3), 1, f);
				fread(&p4, sizeof(p4), 1, f);
				
				fread(&f1, sizeof(f1), 1, f);
				fread(&f2, sizeof(f2), 1, f);
				fread(&f3, sizeof(f3), 1, f);
				fread(&f4, sizeof(f4), 1, f);

				fread(&region, sizeof(region), 1, f);

				elements[i] = new(elem) Tetrahedron(i, vertices[p1], vertices[p2], vertices[p3], vertices[p4], region, 
													(TriFace *)faces[f1], (TriFace *)faces[f2], (TriFace *)faces[f3], (TriFace *)faces[f4]);
			} break;
			default:
				throw "Element type not yet implemented";
		}
	}
}

Mesh::~Mesh() {
	for (int i = 0; i < nVert; i++) {
		vertices[i]->~Vertex();
	}
	for (int i = 0; i < nElems; i++) {
		elements[i]->~Element();
	}
	for (int i = 0; i < nFaces; i++) {
		faces[i]->~Face();
	}
	delete[] vertices;
	delete[] faces;
	delete[] elements;
	Node<char *> *cnext;
	for (Node<char *> *chunk = memory; chunk; chunk = cnext) {
		delete[] chunk->data;
		cnext = chunk->next;
		delete[] chunk;
	}
}

double Mesh::quality() {
	double qual = 1, r;
	for (int i=0; i<nElems; i++) {
		if ((r = elements[i]->quality) < qual) 
			qual = r;
	}
	return qual;
}

bool Mesh::check() {
	bool ok = true, lastcheck;
	lastcheck = true;
	for (int i=0; i<nVert; i++) {
		lastcheck &= (vertices[i]->index == i);
	}
	ok &= lastcheck;
	if (!lastcheck)
		fprintf(stderr, "Vertex index is corrupted\n");

	lastcheck = true;
	for (int i=0; i<nFaces; i++) {
		lastcheck &= (faces[i]->index == i);
	}
	ok &= lastcheck;
	if (!lastcheck)
		fprintf(stderr, "Face index is corrupted\n");

	lastcheck = true;
	for (int i=0; i<nFaces; i++) {
		lastcheck &= (faces[i]->flip->flip->index == i);
	}
	ok &= lastcheck;
	if (!lastcheck)
		fprintf(stderr, "Face flip&flip != id\n");

	lastcheck = true;
	for (int i=0; i<nFaces; i++) {
		bool test = (fabs(faces[i]->flip->normal.dot(faces[i]->normal) + 1) < 1e-10);
		lastcheck &= test;
	}
	ok &= lastcheck;
	if (!lastcheck)
		fprintf(stderr, "Face flip not opposed to face (probably wrong oriented face, weird)\n");

	lastcheck = true;
	for (int i=0; i<nFaces; i++) {
		lastcheck &= ((faces[i]->element != 0) ^ (faces[i]->borderType != -1));
	}
	ok &= lastcheck;
	if (!lastcheck)
		fprintf(stderr, "Inner face with boundaryCondition or outer face without\n");

	lastcheck = true;
	for (int i=0; i<nFaces; i++) {
		lastcheck &= (faces[i]->surface > 0);
	}
	ok &= lastcheck;
	if (!lastcheck)
		fprintf(stderr, "Face with negative surface\n");

	lastcheck = true;
	for (int i=0; i<nElems; i++) {
		lastcheck &= (elements[i]->volume > 0);
	}
	ok &= lastcheck;
	if (!lastcheck)
		fprintf(stderr, "Element with negative volume\n");

	return ok;
}

Mesh::Mesh(const Front *front, double meshsize, bool fixShape, int nnV, int nnF, int nnT) {
	int    nF = front->nFaces, nV = front->nVert, nT = 0;
	int    *facefront = 0, *facedup = 0, *facematerial  = 0, *facematdup = 0;
	int    *tetra = 0, *tetramaterial = 0;
	double *vertexdup = 0;
	int    r;

	/* TODO malloc/free -> new/delete */

	vertexdup     = new double[3 * nnV];    // for running front & Mesh(...)
	facefront     = new int[3*nnF];    // for running front
	facematerial  = new int[  nnF];    // for running front
	facedup       = new int[3*front->nFaces]; // for Mesh(...)
	facematdup    = new int[  front->nFaces]; // for Mesh(...)
	tetra         = new int[4*nnT];    // for Mesh(...)
	tetramaterial = new int[  nnT];    // for Mesh(...)

	memcpy(vertexdup, front->vertex, 3*sizeof(double)*front->nVert);
	memcpy(facedup, front->face, 3*sizeof(int)*front->nFaces);

	for (int i=0; i < front->nFaces; i++) {
		facedup[3*i+0]--;
		facedup[3*i+1]--;
		facedup[3*i+2]--;
		facematdup[i] = 1;
	}
	memcpy(facefront, facedup, 3*sizeof(int)*front->nFaces);
	memcpy(facematerial, facematdup, sizeof(int)*front->nFaces);

	meshSize = meshsize;
	r = mesh_3d_aft_func(&nV, vertexdup, &nF, facefront, facematerial, &nT, tetra, tetramaterial, nnV, nnF, nnT, usersize);
	if (r!=0)  {
		// Meshing failed
		throw "mesh_3d_aft_func: failed";
	}

	if (fixShape) {
		/* Restore initial front */
		nF = front->nFaces;
		memcpy(facefront, facedup, 3*sizeof(int)*front->nFaces);
		memcpy(facematerial, facematdup, sizeof(int)*front->nFaces);

		int *ifv = new int[nF];
		for (int i=0; i<nF; i++) {
			ifv[i] = i+1;
			facefront[3*i+0]++;
			facefront[3*i+1]++;
			facefront[3*i+2]++;
		}
		for (int i=0; i<nT; i++) {
			tetra[4*i+0]++;
			tetra[4*i+1]++;
			tetra[4*i+2]++;
			tetra[4*i+3]++;
		}

		const int maxWr = 10000000, maxWi = 25000000;
		int *iW = new int[maxWi];
		double *rW = new double[maxWr];
		double rQual;

		r = mbaFixShape(&nV, nnV, &nF, nnF, &nT, nnT, 
			vertexdup, facefront, tetra, facematerial, tetramaterial,
			0, front->nFaces, 0, 0, ifv, 0,
			1, 0, 
			300, 500000, 
			euclid_metric, 1, &rQual, 
			maxWr, maxWi, rW, iW, 
			2);

		delete[] iW;
		delete[] rW;
		delete[] ifv;

		for (int i=0; i<nT; i++) {
			tetra[4*i+0]--;
			tetra[4*i+1]--;
			tetra[4*i+2]--;
			tetra[4*i+3]--;
		}

		/* just ignored r ... TODO checks */
	}

	fromAft(nV, front->nFaces, nT, vertexdup, facedup, tetra, facematdup, tetramaterial);

	delete[] vertexdup;
	delete[] facefront;
	delete[] facedup;
	delete[] facematerial;
	delete[] facematdup;
	delete[] tetra;
	delete[] tetramaterial;
}
