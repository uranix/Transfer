#include <stdlib.h>
#include <string.h>
#include "tetra3.h"
#include "region3.h"
#include "tree3.h"
#include "support.h"
#include "memory3.h"

extern StrucMesh3 mesh3;
extern StrucTree3 tree3;

static int niceit(int *pnVout, double *vertexout,
	int *pnFout, int *faceout, int *facecolor,
	int maxnV, int maxnF, int is) {
    int i, nf, nv;
    (void) maxnV, (void) maxnF;
    nv = *pnVout;
    nf = *pnFout;
    *pnVout += mesh3.nPoint;
    for (i=nv; i<*pnVout; i++) {
	vertexout[3*i+0] = mesh3.vert[i-nv].x;
	vertexout[3*i+1] = mesh3.vert[i-nv].y;
	vertexout[3*i+2] = mesh3.vert[i-nv].z;
    }
    *pnFout += tree3.nFace;
    for (i=nf; i<*pnFout; i++) {
	faceout[3*i+0] = tree3.face[i-nf]->v1 + nv + is;
	faceout[3*i+1] = tree3.face[i-nf]->v2 + nv + is;
	faceout[3*i+2] = tree3.face[i-nf]->v3 + nv + is;
	facecolor[i]   = tree3.face[i-nf]->color;
    }
    return 0;
}

int aft3dboundary (
	int nVVert, double *VVertxyz,
	int nLine, int *LineD, int *LineP, double *LineT,
	int nSurface, int *SurfL, int *SurfI, double *SurfT,
	int *pnVout, double *vertexout,
	int *pnFout, int *faceout, int *facecolor,
	int maxnV, int maxnF, int indexshift) {
    int r;
    initAFS_(&nVVert, VVertxyz, &nLine, LineD, LineP, LineT, &nSurface, SurfL, SurfI, SurfT);
    r = niceit(pnVout, vertexout, pnFout, faceout, facecolor, maxnV, maxnF, indexshift);
    freeMemory();
    return r;
}

