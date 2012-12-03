#include "MeshData.h"

#include "meshProcessor/mesh.h"
#include "meshProcessor/tetrahedron.h"
#include "meshProcessor/tri_face.h"

#include <iostream>

MeshData::MeshData(const Config &cfg) {
	_m = new Mesh(std::cerr, cfg.getMeshFilename());
	bool check = _m->check(std::cerr);
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
		nt = v->elems.size();
		nb = v->bnds.size();
		tetstart[i+1] = tetstart[i] + nt;
		facestart[i+1] = facestart[i] + nb;
		st += nt;
		sb += nb;
	}
	nT = _m->nElems;
	nF = _m->nBndFaces;
	tetidx = new idx[st];
	tetpos = new idx[st];
	faceidx = new idx[sb];
	facepos = new idx[sb];
	sb = st = 0;
	for (idx i=0; i<nP; i++) {
		Vertex *v = _m->vertices[i];
		nt = v->elems.size();
		nb = v->bnds.size();
		std::vector<std::pair<Element *, int> > &e = v->elems;
		for (std::vector<std::pair<Element *, int> >::iterator j = e.begin(); j != e.end(); j++) {
			tetpos[st] = j->second;
			tetidx[st] = j->first->index;
			st++;
		}
		std::vector<std::pair<Face *, int> > &b = v->bnds;
		for (std::vector<std::pair<Face *, int> >::iterator j = b.begin(); j != b.end(); j++) {
			facepos[sb] = j->second;
			faceidx[sb] = j->first->bnd_index;
			sb++;
		}
	}

	mesh = new tetrahedron[_m->nElems];
	bnd = new face[_m->nBndFaces];

	for (idx i = 0; i < (idx)_m->nElems; i++) {
		Tetrahedron *t = (Tetrahedron *)_m->elements[i];
		for (int j=0; j<4; j++)
			mesh[i].p[j] = t->p[j]->index;
		mesh[i].kappa_volume = t->volume * cfg.getKappa(t->region);
		mesh[i].I_p = cfg.getIp(t->region);
		Vector z(0,0,0);
		for (int j=0; j<4; j++) {
			Vector s(t->f[j]->normal);
			s.scale( - t->f[j]->surface); /* note the minus */
			z.add(s);
			mesh[i].s[j][0] = s.x;
			mesh[i].s[j][1] = s.y;
			mesh[i].s[j][2] = s.z;
		}
		if (z.norm() > 1e-14) {
			printf("Tet %d has sum of S_i = (%e,%e,%e)\n", i, z.x, z.y, z.z);
		}
	}

	for (idx i = 0, j = 0; i < (idx)_m->nFaces; i++) {
		TriFace *f = (TriFace *)_m->faces[i];
		if (f->bnd_index < 0)
			continue;
		for (int k=0; k<3; k++)
			bnd[j].p[k] = f->p[k]->index;
		Vector s(f->normal);
		s.scale(f->surface);
		bnd[j].s[0] = s.x;
		bnd[j].s[1] = s.y;
		bnd[j].s[2] = s.z;
		j++;
	}

	printf("Mesh has %d points, %d tetrahedra and %d boundary faces\n", _m->nVert, _m->nElems, _m->nBndFaces);
}

void MeshData::ComputeMoments(
		const AngularData &ad,
		const REAL * const I[],
		REAL *U,
		REAL *T[3][3]) 
{
	REAL fourpi = (REAL)12.5663706143591729538505735331;
	for (idx k = 0; k < nP; k++) {
		U[k] = fourpi * I[0][k];
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				T[i][j][k] = 0;
		for (idx lm = 0; lm < ad.slm; lm++) 
			for (int i = 0; i < 3; i++) {
				// lm' = 0
				int j = ad.omega_pos[i * ad.slm * ad.slm + lm * ad.slm];
				REAL v = ad.omega[i * ad.slm * ad.slm + lm * ad.slm];
				T[i][j][k] += fourpi * v * I[lm][k];
			}
	}
}

void MeshData::ComputeFlux(	
		REAL *T[3][3], 
		REAL *W[3])
{
	// W = -1/kappa div_j T_ij
	// assumed c = 1
	REAL w[3];

	for (idx i = 0; i < nT; i++) { 
		tetrahedron *t = mesh + i;	
		for (int j = 0; j < 3; j++) {
			w[j] = 0;
			for (int k = 0; k < 4; k++)
				for (int a = 0; a < 3; a++)
					w[j] += t->s[k][a] * T[a][j][t->p[k]];
			W[j][i] = ((REAL)(1.0 / 3.0)) * w[j] / t->kappa_volume;
		}
	}
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
