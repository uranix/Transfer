#include <stdio.h>
#include <errno.h>
#include <stdarg.h>
#include "tetra3.h"
#include "region3.h"
#include "tree3.h"
#include "user3.h"
#include "support.h"
#include "libfrtprm.h"


int set_surface_size_function(double (*f) (double, double, double)) {
	return setsizefunction(f);
}
int set_param_functions(int (*bounsurf) (int, double, double, double*, double*, double*), void (*v_u) (int, double, double*)) {
	return setbounsurffunction(bounsurf) + setvufunction(v_u);
}
int surface_boundary (
		int nVVert, double *VVertxyz,
		int nLine, int *LineD, int *LineP, double *LineT,
		int nSurface, int *SurfL, int *SurfI, double *SurfT,
		int *pnV, double *vertex,
		int *pnF, int *face, int *facecolor,
		int maxnV, int maxnF) {
	return aft3dboundary(nVVert, VVertxyz, nLine, LineD, LineP, LineT, nSurface, SurfL, SurfI, SurfT, pnV, vertex, pnF, face, facecolor, maxnV, maxnF, 0);
}
int surface_boundary_ (
		int *pnVVert, double *VVertxyz,
		int *pnLine, int *LineD, int *LineP, double *LineT,
		int *pnSurface, int *SurfL, int *SurfI, double *SurfT,
		int *pnV, double *vertex,
		int *pnF, int *face, int *facecolor,
		int *pmaxnV, int *pmaxnF) {
	return aft3dboundary(*pnVVert, VVertxyz, *pnLine, LineD, LineP, LineT, *pnSurface, SurfL, SurfI, SurfT, pnV, vertex, pnF, face, facecolor, *pmaxnV, *pmaxnF, 1);
}

int polyhedron_add_faceloop(int *pnE, int *edge, int nP, ...)
{
	va_list ap;
	int nE=*pnE;
	int i, v1, v2, v0;
	if (!nP) return 0;
	va_start(ap, nP);
	v0 = v2 = va_arg(ap, int);
	for (i=1; i<nP; i++)
	{
		v1 = v2, v2 = va_arg(ap, int);
		edge[2*nE+0] = v1;
		edge[2*nE+1] = v2;
		nE++;
	}
	edge[2*nE+0] = v2;
	edge[2*nE+1] = v0;
	nE++;

	*pnE = nE;
	return 0;
}

int polyhedron_make_front (int *pnV, double *vertex, int nS, int *nE, int **edge, int *color1, int *color2, 
		int *pnF, int *face, int *facematerial, int nnV, int nnF) 
{
	int nV = *pnV;
	int i, nEs=0;
	int *vls, *els;
	int v1, v2, inv, line, k, m;
	int    nVVert,    nLine,  nSurface;
	int    *LineD,    *LineP;
	double *LineT;
	int    *SurfL,    *SurfI;
	double *SurfT;
	double *VVert;

	nVVert = nV;
	*pnV = 0;

	// we compute the total number of edges -- nEs
	// and allocate memory for additional structures els and vls
	for (i=0; i<nS; i++) nEs += nE[i];
	els = (int*)malloc(sizeof(int) * (nEs+1) * 2);
	vls = (int*)malloc(sizeof(int) * (nV+1));
	for (i=0; i<nV+1; i++) vls[i] = -1;

	// allocate memory for boundary representation structure
	VVert = (double*)malloc(sizeof(double) *2* 3*nV);
	LineD = (int*)   malloc(sizeof(int)    *2* 3*nEs);
	LineP = (int*)   malloc(sizeof(int)    *2* 2*nEs);
	LineT = (double*)malloc(sizeof(double) *2* 2*nEs);
	nLine = 0;
	SurfL = (int*)   malloc(sizeof(int)    *2* 5*nS);
	SurfI = (int*)   malloc(sizeof(int)    *2* 2*nEs);
	SurfT = (double*)malloc(sizeof(double) *2* 4*nS);
	nSurface = 0;

	// copy coords from vertex to VVert
	for (i=0; i<nV; i++) 
	{
		VVert[3*i+0] = vertex[3*i+0];
		VVert[3*i+1] = vertex[3*i+1];
		VVert[3*i+2] = vertex[3*i+2];
	}

	m = 0;
	// for each surface
	for (nSurface=0; nSurface<nS; nSurface++) 
	{
		// fill boundary representation structure
		SurfL[5*nSurface+0] = nE[nSurface];
		SurfL[5*nSurface+1] = 0;
		SurfL[5*nSurface+2] = color1[nSurface];
		SurfL[5*nSurface+3] = (color2)?color2[nSurface]:0;
		SurfL[5*nSurface+4] = 0;
		SurfT[4*nSurface+0] = 0.0;
		SurfT[4*nSurface+1] = 0.0;
		SurfT[4*nSurface+2] = 0.0;
		SurfT[4*nSurface+3] = 0.0;
		// for each edge in surface
		for (k=0; k<nE[nSurface]; k++) 
		{
			// (v1,v2) is the current edge
			v1 = edge[nSurface][2*k+0];
			v2 = edge[nSurface][2*k+1];
			inv = -1; line = -1;
			// check if edge (v1,v2) is already added
			for (i=vls[v1]; i>=0; i=els[2*i+1]) 
			{
				if (els[2*i+0]==v2) 
				{
					// edge (v1,v2) found
					inv = 0;
					line = i + 1;
				}
			}
			if (inv<0) 
			{
				// check if edge (v2,v1) is already added
				for (i=vls[v2]; i>=0; i=els[2*i+1])
				{
					if (els[2*i+0]==v1) 
					{
						// edge (v2,v1) found
						inv = 1;
						line = i + 1;
					}
				}
			}
			if (inv<0) 
			{
				// add edge (v1,v2) into boundary representation structure
				els[2*nLine+0] = v2;
				els[2*nLine+1] = vls[v1];
				vls[v1] = nLine;
				line = 1 + nLine;
				inv = 0;
				LineD[3*nLine+0] = v1;
				LineD[3*nLine+1] = v2;
				LineD[3*nLine+2] = 1;
				LineP[2*nLine+0] = 0;
				LineP[2*nLine+1] = 0;
				LineT[2*nLine+0] = 0.0;
				LineT[2*nLine+1] = 0.0;
				nLine++;
			}
			// add edge into surface with respect of orientation
			SurfI[2*m+0] = line;
			SurfI[2*m+1] = inv;
			m++;
		}
		//		printf("\n");
	}
	// create surface mesh
	i = surface_boundary_(&nVVert, VVert, &nLine, LineD, LineP, LineT, &nSurface, SurfL, SurfI, SurfT,
			pnV, vertex, pnF, face, facematerial, &nnV, &nnF);
	free(VVert), free(LineD), free(LineP), free(LineT), free(SurfL), free(SurfI), free(SurfT);
	free(els), free(vls);
	return i;
}
