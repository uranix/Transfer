#include "MeshData.h"

#include "meshProcessor/mesh.h"
#include "meshProcessor/tetrahedron.h"
#include "meshProcessor/tri_face.h"

MeshData::MeshData(const char *fn) {
	_m = new Mesh(fn);
	bool check = _m->check();
	tetstart = 0;
	printf("Mesh check: %s\n", check ? "OK": "Failed");
	if (!check)
		return;
	
	nP = _m->nVert;
	tetstart = new idx[nP+1];
	facestart = new idx[nP+1];
	tetstart[0] = facestart[0] = 0;
	idx st = 0, sb = 0, nt, nb;
	for (idx i=0; i<nP; i++) {
		Vertex *v = _m->vertices[i];
		nt = v->elnum;
		nb = v->bndnum;
		tetstart[i+1] = tetstart[i] + nt;
		facestart[i+1] = facestart[i] + nb;
	}
	nT = tetstart[nP];
	nF = facestart[nP];
	tetidx = new idx[nT];
	tetpos = new idx[nT];
	faceidx = new idx[nF];
	facepos = new idx[nF];
	sb = st = 0;
	for (idx i=0; i<nP; i++) {
		Vertex *v = _m->vertices[i];
		nt = v->elnum;
		nb = v->bndnum;
		Node<Element *> *e = v->elems;
		Node<int> *ei = v->elidx;
		for (idx j=0;j<nt;j++) {
			tetpos[st] = ei->data;
			tetidx[st] = e->data->index;
			st++;
			ei = ei->next;
			e = e->next;
		}
		Node<Face *> *b = v->bnds;
		Node<int> *bi = v->bndidx;
		for (idx j=0;j<nb;j++) {
			facepos[sb] = bi->data;
			faceidx[sb] = b->data->bnd_index;
			sb++;
			bi = bi->next;
			b = b->next;
		}
	}

	mesh = new tetrahedron[_m->nElems];
	bnd = new face[_m->nBndFaces];

	for (idx i = 0; i < _m->nElems; i++) {
		Tetrahedron *t = (Tetrahedron *)_m->elements[i];
		for (int j=0; j<4; j++)
			mesh[i].p[j] = t->p[j]->index;
		mesh[i].kappa_volume = t->volume * (t->region == 1? 10.: 1.);
		mesh[i].I_p = t->region == 1? 10.: 1.;
		for (int j=0; j<4; j++) {
			Vector s(t->f[j]->normal);
			s.scale( - t->f[j]->surface); /* note the minus */
			mesh[i].s[j][0] = s.x;
			mesh[i].s[j][1] = s.y;
			mesh[i].s[j][2] = s.z;
		}
	}

	for (idx i = 0, j = 0; i < _m->nFaces; i++) {
		TriFace *f = (TriFace *)_m->faces[i];
		if (f->bnd_index < 0)
			continue;
		for (int k=0; k<3; k++)
			bnd[j].p[k] = f->p[k]->index;
		Vector s(f->normal);
		s.scale( - f->surface); /* note the minus */
		bnd[j].s[0] = s.x;
		bnd[j].s[1] = s.y;
		bnd[j].s[2] = s.z;
		j++;
	}

	printf("Mesh has %d tetrahedra and %d boundary faces\n", _m->nElems, _m->nBndFaces);
}

MeshData::~MeshData() {
	delete _m;
	delete[] tetstart;
	delete[] tetidx;
	delete[] tetpos;

	delete[] facestart;
	delete[] faceidx;
	delete[] facepos;

	delete[] mesh;
	delete[] bnd;
}
