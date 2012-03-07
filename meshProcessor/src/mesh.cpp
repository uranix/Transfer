#include "mesh.h"

#include "tri_face.h"
#include "tetrahedron.h"
#include "list.h"

#include "stdint.h"
#include <stdlib.h>
#include <string.h>

void Mesh::fromVol(int nV, int nB, int nT, double *vert, int *bnd, int *tet, int *bndmat, int *tetmat) {
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

Mesh::Mesh(char *fn) {
	enum State {
		ST_NORM,
		ST_SURF_INIT,
		ST_SURF_DATA,
		ST_VOL_INIT,
		ST_VOL_DATA,
		ST_PTS_INIT,
		ST_PTS_DATA,
		ST_CURVE_INIT,
		ST_CURVE_DATA,
		ST_STOP,
	} state = ST_NORM;
	int i, cnt;

	FILE *f = fopen(fn, "r");
	if (!f) throw 0;
	char buf[1024];
	int nV, nB, nBx, nT;
	double *vert = 0;
	int *bnd = 0, *tet = 0, *bndmat = 0, *tetmat = 0;

	while (fgets(buf, 1024, f)) {
		/* printf("[%.2d] %s", state, buf); */
		if (buf[0]=='#')
			continue;
		if (state == ST_NORM) {
			if (!strncmp(buf, "endmesh", 7)) {
				state = ST_STOP;
				break;
			}
			if (!strncmp(buf, "surfaceelements", 15)) {
				state = ST_SURF_INIT;
				continue;
			}
			if (!strncmp(buf, "points", 6)) {
				state = ST_PTS_INIT;
				continue;
			}
			if (!strncmp(buf, "volumeelements", 14)) {
				state = ST_VOL_INIT;
				continue;
			}
			if (!strncmp(buf, "edgesegmentsgi2", 15)) {
				state = ST_CURVE_INIT;
				continue;
			}
		}
		if (state == ST_PTS_INIT) {
			nV = cnt = atoi(buf);
			vert = new double [3*nV];
			i = 0;
			state = ST_PTS_DATA;
			continue;
		}
		if (state == ST_VOL_INIT) {
			nT = cnt = atoi(buf);
			tet = new int [4*nT];
			tetmat = new int [nT];
			i = 0;
			state = ST_VOL_DATA;
			continue;
		}
		if (state == ST_SURF_INIT) {
			nB = cnt = atoi(buf);
			bnd = new int [3*nB];
			bndmat = new int [nB];
			nBx = 0;
			i = 0;
			state = ST_SURF_DATA;
			continue;
		}
		if (state == ST_CURVE_INIT) {
			cnt = atoi(buf);
			state = ST_CURVE_DATA;
			i = 0;
			/* Just ignore it */
			continue;
		};
		if (state == ST_PTS_DATA) {
			sscanf(buf, "%lf %lf %lf", vert + 3*i, vert + 3*i + 1, vert + 3*i + 2);
			if (++i == cnt)
				state = ST_NORM;
			continue;
		}
		if (state == ST_VOL_DATA) {
			int np;
			sscanf(buf, "%d %d %d %d %d %d", tetmat + i, &np, 
				tet + 4*i, tet + 4*i + 1, tet + 4*i + 2, tet + 4*i + 3);
			tet[4*i]--;
			tet[4*i+1]--;
			tet[4*i+2]--;
			tet[4*i+3]--;
			if (np != 4) {
				fprintf(stderr, "high-order tetrahedrons not implemented\n");
			}
			if (++i == cnt)
				state = ST_NORM;
			continue;
		}
		if (state == ST_SURF_DATA) {
			int sn, domin, domout, np;
			sscanf(buf, "%d %d %d %d %d %d %d %d", &sn, bndmat + nBx, &domin, &domout, &np, 
				bnd + 3*nBx, bnd + 3*nBx + 1, bnd + 3*nBx + 2);
			bnd[3*nBx]--;
			bnd[3*nBx+1]--;
			bnd[3*nBx+2]--;
			if (domout == 0) 
				nBx ++;
			if (np != 3) {
				fprintf(stderr, "high-order faces not implemented\n");
			}
			if (++i == cnt)
				state = ST_NORM;
			continue;
		}
		if (state == ST_CURVE_DATA) {
			if (++i == cnt)
				state = ST_NORM;
			continue;
		}
	}
	if (state == ST_STOP) {
		fromVol(nV, nBx, nT, vert, bnd, tet, bndmat, tetmat);
	}
	else {
		throw 0;
		/* Set Mesh object to valid state so ~Mesh() could safely destroy it */
	}
	if (vert) delete[] vert;
	if (bnd) delete[] bnd;
	if (bndmat) delete[] bndmat;
	if (tet) delete[] tet;
	if (tetmat) delete[] tetmat;
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
		fprintf(stderr, "Face flip not opposed to face (probably wrong oriented face)\n");

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
/*
	Element faces should be directed inside, so boundary faces could be directed outside
*/
/*
	lastcheck = true;
	for (int i=0; i<nElems; i++) {
		Tetrahedron *t = (Tetrahedron *)elements[i];
		for (int j=0; j<4; j++) {
			Vector r(t->f[j]->center);
			r.sub(t->center);
			lastcheck &= r.dot(t->f[j]->normal) < 0;
		}
	}
	ok &= lastcheck;
	if (!lastcheck)
		fprintf(stderr, "Element with face-normal directed outside\n");
*/

	return ok;
}