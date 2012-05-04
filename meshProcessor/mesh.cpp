#include "mesh.h"

#include "tri_face.h"
#include "tetrahedron.h"
#include "list.h"

#include <stdint.h>
#include <stdlib.h>
#include <stdarg.h>
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

	Node<int> *vert2face, *memv2f = list_allocate<int>(3*nFaces);
	Node<Element *> *memv2e = list_allocate<Element *>(4*nElems);
	Node<Face *> *memv2b = list_allocate<Face *>(3*nB);
	Node<int> *memv2ei = list_allocate<int>(4*nElems);
	Node<int> *memv2bi = list_allocate<int>(3*nB);
	vert2face = memv2f;
	vert2elem = memv2e;
	vert2bnd = memv2b;
	vert2elidx = memv2ei;
	vert2bndidx = memv2bi;

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
		elements[i] = new (tetraplace + i) Tetrahedron(i, 
				vertices[tet[4*i + 0]], 
				vertices[tet[4*i + 1]], 
				vertices[tet[4*i + 2]], 
				vertices[tet[4*i + 3]], 
				tetmat[i], faceplace+4*i);
		for (int j = 0; j < 4; j++) {
			addHead(elements[i], &vertices[tet[4*i+j]]->elems, &memv2e);
			addHead(j, &vertices[tet[4*i+j]]->elidx, &memv2ei);
			vertices[tet[4*i+j]]->elnum++;
			faces[4*i+j] = faceplace + 4*i+j;
			for (int k = 0; k < 4; k++)
				 if (j != k)
					addHead(4*i + j, &lists[tet[4*i + k]], &memv2f);
		}
	}
	for (int i = 0; i < nB; i++) {
		faces[4*nT + i] = new(faceplace + 4*nT + i) TriFace(4*nT + i, 
				vertices[bnd[3*i+0]], 
				vertices[bnd[3*i+1]], 
				vertices[bnd[3*i+2]], 
				0, bndmat[i]);
		for (int j = 0; j < 3; j++) {
			addHead(faces[4*nT + i], &vertices[bnd[3*i+j]]->bnds, &memv2b);
			addHead(j, &vertices[bnd[3*i+j]]->bndidx, &memv2bi);
			vertices[bnd[3*i+j]]->bndnum++;
			addHead(4*nT + i, &lists[bnd[3*i+j]], &memv2f);
		}
	}
	int bnd_idx = 0;
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
		if (!f->element) 
			f->bnd_index = bnd_idx++;
		else
			f->bnd_index = -1;
	}
	nBndFaces = bnd_idx;
	delete[] lists;
	list_deallocate(vert2face);
}

Mesh::Mesh(const char *fn) {
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
	int i = 0, cnt = 0;

	FILE *f = fopen(fn, "r");
	if (!f) {
		fprintf(stderr, "%s:%d error opening file `%s'\n", __FILE__, __LINE__, fn);
		fflush(stderr);
		throw 0;
	}
	char buf[1024];
	int nV = 0, nB = 0, nBx = 0, nT = 0;
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
			state = cnt?ST_PTS_DATA:ST_NORM;
			continue;
		}
		if (state == ST_VOL_INIT) {
			nT = cnt = atoi(buf);
			tet = new int [4*nT];
			tetmat = new int [nT];
			i = 0;
			state = cnt?ST_VOL_DATA:ST_NORM;
			continue;
		}
		if (state == ST_SURF_INIT) {
			nB = cnt = atoi(buf);
			bnd = new int [3*nB];
			bndmat = new int [nB];
			nBx = 0;
			i = 0;
			state = cnt?ST_SURF_DATA:ST_NORM;
			continue;
		}
		if (state == ST_CURVE_INIT) {
			cnt = atoi(buf);
			i = 0;
			state = cnt?ST_CURVE_DATA:ST_NORM;
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

void Mesh::saveVtk(const char *fn, int realbytes, int nExtraCellData, int nExtraPointData, ...) {
	typedef union {
		double *d;
		float *s;
	} realptr;

	FILE *f;
	va_list args;
	va_start(args, nExtraPointData);

	f=fopen(fn, "w");
	if (!f) {
		fprintf(stderr, "%s:%d error opening file `%s'\n", __FILE__, __LINE__, fn);
		fflush(stderr);
		throw 1;
	}
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
	fprintf(f, "\nFIELD MeshProcessorData 5\n");
	fprintf(f, "\n__mp_Volume 1 %9d double\n", nElems);
	for (int i=0; i<nElems; i++)
		fprintf(f, "%2.10e %s", elements[i]->volume, (i&3)==3?"\n":"");
	fprintf(f, "\n__mp_Quality 1 %9d double\n", nElems);
	for (int i=0; i<nElems; i++)
		fprintf(f, "%2.10e %s", elements[i]->quality, (i&3)==3?"\n":"");
	fprintf(f, "\n__mp_CenterX 1 %9d double\n", nElems);
	for (int i=0; i<nElems; i++)
		fprintf(f, "%2.10e\n", elements[i]->center.x);
	fprintf(f, "\n__mp_CenterY 1 %9d double\n", nElems);
	for (int i=0; i<nElems; i++)
		fprintf(f, "%2.10e\n", elements[i]->center.y);
	fprintf(f, "\n__mp_CenterZ 1 %9d double\n", nElems);
	for (int i=0; i<nElems; i++)
		fprintf(f, "%2.10e\n", elements[i]->center.z);

	fprintf(f, "\nFIELD ExtraCellData %d\n", nExtraCellData);
	for (int nextra = 0; nextra < nExtraCellData; nextra++) {
		realptr p = va_arg(args, realptr);
		fprintf(f, "\nextra_data%.4d 1 %9d double\n", nextra, nElems);
		for (int i=0; i<nElems; i++)
			fprintf(f, " %2.16e\n", (realbytes == 4) ? p.s[i] : p.d[i]);
	}

	fprintf(f, "\nPOINT_DATA %9d\n", nVert);
	fprintf(f, "\nFIELD ExtraPointData %d\n", nExtraPointData);
	for (int nextra = 0; nextra < nExtraPointData; nextra++) {
		realptr p = va_arg(args, realptr);
		fprintf(f, "\nextra_data%.4d 1 %9d double\n", nextra, nVert);
		for (int i=0; i<nVert; i++)
			fprintf(f, " %2.16e\n", (realbytes == 4) ? p.s[i] : p.d[i]);
	}
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
		delete chunk;
	}
	list_deallocate<Element *>(vert2elem);
	list_deallocate<Face *>(vert2bnd);
	list_deallocate<int>(vert2elidx);
	list_deallocate<int>(vert2bndidx);
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
	Check specific to tetrahedron meshes
*/

	lastcheck = true;
	for (int i=0; i<nElems; i++) {
		Tetrahedron *t = (Tetrahedron *)elements[i];
	/*	for (int j=0; j<4; j++) {
			Vector r(t->f[j]->center);
			r.sub(t->center);
			lastcheck &= r.dot(t->f[j]->normal) < 0;
		}*/
		for (int j=0; j<4; j++) {
			int k = (j+1) & 3;
			Vector r(t->p[j]->r);
			r.sub(t->p[k]->r);
			r.scale(t->f[j]->surface/3);
			lastcheck &= ( abs(r.dot(t->f[j]->normal) - t->volume) <= 1e-12 * abs(t->volume) );
		}
	}
	ok &= lastcheck;
	if (!lastcheck)
		fprintf(stderr, "Element with face-normal directed outside\n");

	int cnt1 = 0, cnt2 = 0;
	lastcheck = true;

	for (int i=0; i<nVert; i++) {
		Node<Element *> *q = vertices[i]->elems;
		Node<int> *qi = vertices[i]->elidx;
		Node<Face *> *p = vertices[i]->bnds;
		Node<int> *pi = vertices[i]->bndidx;
		for (;q;) {
			cnt1++;
			lastcheck &= ((Tetrahedron *)q->data)->p[qi->data]->index == i;
			q = q->next;
			qi = qi->next;
		}
		for (;p;) {
			cnt2++;
			lastcheck &= ((TriFace *)p->data)->p[pi->data]->index == i;
			p = p->next;
			pi = pi->next;
		}
	}

	ok &= lastcheck;
	if (!lastcheck)
		fprintf(stderr, "Incidental lists are corrupted\n");

	lastcheck = (cnt1 == 4*nElems);
	ok &= lastcheck;
	if (!lastcheck)
		fprintf(stderr, "Incidental elements list has wrong size\n");

	lastcheck = (cnt2 == 3*(nFaces - 4*nElems));
	ok &= lastcheck;
	if (!lastcheck)
		fprintf(stderr, "Incidental faces list has wrong size\n");

	return ok;
}
