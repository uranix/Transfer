#include "aft.h"
#include "tetra3.h"
#include "tria32.h"
#include "region3.h"
#include "error3.h"
#include "memory3.h"
#include "tree3.h"
#include "tree32.h"
#include "section3.h"
#include "sectio32.h"
#include "refine3.h"
#include "user3.h"


/* extern  variables */
extern StrucMesh3  mesh3;
extern StrucSect3  sect3;
extern StrucTree3  tree32;
extern StrucTree3  tree3;
extern StrucRegion3  reg3;
extern double  minNear;
extern int     nBadVert,badVert[100],nSphereVert,sphereVert[100];


double CF2_5D = 0.0;

int tria_dump_front = 0;
int tria_debug_front = 0;
int tria_final_front = 0;

#define max(a,b)  (((a) > (b)) ? (a) : (b))
#define min(a,b)  (((a) < (b)) ? (a) : (b))

#define N_SPLIT 16
#define N_DEEP 24
#define N_ANGLE 20

double      xInv,yInv,zInv;
StrucVert3  surfNorm;

static int dump();
static int dump3();

void addSphereVert(int v);
int isSphereVert(int v);
int goodPlane32(int v1, int v2, int v);
void surfNormal(int iSurf, int iNorm, double u, double v, StrucVert3 *norm);
static int choseWorkVert(void);


static double ddet3(double a, double b, double c,
		    double d, double e, double f,
		    double g, double h, double i) {
    return a*e*i - c*e*g + b*f*g - a*f*h + c*d*h - b*d*i;
}
static int idet4(int a, int b, int c, int d) {
    double len = (distance(mesh3.vert[a].x, mesh3.vert[a].y, mesh3.vert[a].z,  mesh3.vert[d].x, mesh3.vert[d].y, mesh3.vert[d].z) +
		  distance(mesh3.vert[b].x, mesh3.vert[b].y, mesh3.vert[b].z,  mesh3.vert[d].x, mesh3.vert[d].y, mesh3.vert[d].z) +
		  distance(mesh3.vert[c].x, mesh3.vert[c].y, mesh3.vert[c].z,  mesh3.vert[d].x, mesh3.vert[d].y, mesh3.vert[d].z)) / 3.0;
    double ddet = ddet3(mesh3.vert[a].x-mesh3.vert[d].x, mesh3.vert[b].x-mesh3.vert[d].x, mesh3.vert[c].x-mesh3.vert[d].x,
	    		mesh3.vert[a].y-mesh3.vert[d].y, mesh3.vert[b].y-mesh3.vert[d].y, mesh3.vert[c].y-mesh3.vert[d].y,
			mesh3.vert[a].z-mesh3.vert[d].z, mesh3.vert[b].z-mesh3.vert[d].z, mesh3.vert[c].z-mesh3.vert[d].z);
    ddet /= len*len*len;
    if (ddet > 1e-5)  return +1;
    else if (ddet < -1e-5)  return -1;
    else return 0;
}

static int faceIntersectNew(int a, int b, int c, int u, int v, int d) {
    int uv, dup=0;

    if ((u==a)||(u==b)||(u==c)) dup++;
    if ((v==a)||(v==b)||(v==c)) dup++;
    if (dup==0) {
	uv = idet4(u, v, a, d) + idet4(u, v, b, d) + idet4(u, v, c, d);
	if ((uv==3) || (uv==-3)) return 0;
	if (idet4(b, c, u, d) + idet4(b, c, v, d) == -2) return 0;
	if (idet4(c, a, u, d) + idet4(c, a, v, d) == -2) return 0;
	if (idet4(a, b, u, d) + idet4(a, b, v, d) == -2) return 0;
    } else if (dup==1) {
	if (idet4(b, c, u, d) + idet4(b, c, v, d) == -1) return 0;
	if (idet4(c, a, u, d) + idet4(c, a, v, d) == -1) return 0;
	if (idet4(a, b, u, d) + idet4(a, b, v, d) == -1) return 0;
    } else return 0;
    return 1;
}

static int checkIntersect(PStrucFace3  face, int v) {
    int i, iFace, v1, v2, workVert, bln;
    PStrucFace3 f;

    v1 = face->v1;  v2 = face->v2;

    for (i=0; i<3; i++) {
	sect3.neigSin[i] = 1.;
	sect3.neigBool[i]= 0;
	sect3.neigInt[i] = 0;
    }
    iFace = 0;
    while (iFace < tree32.nVicinityFace) {
	f = tree32.vicinityFace[iFace].face;
	if ( (v1==f->v2) && (v==f->v1) ) {
	    sect3.neigBool[0] = 1;
	    sect3.neigFace[0] = f;
	}

	if ( (v2==f->v1) && (v==f->v2) ) {
	    sect3.neigBool[1] = 1;
	    sect3.neigFace[1] = f;
	}
	iFace++;
    }

    iFace = bln = 0;
    while (iFace < tree32.nVicinityFace) {
	if (1) {
	    f = tree32.vicinityFace[iFace].face;
	    {/* proection on perpendicular to surfNormal */
		StrucVert3  norm,e1,e2;
		double  x,y,z,p1,p2,p,len;

		x = (mesh3.vert[f->v1].x + mesh3.vert[f->v2].x)/2.0;
		y = (mesh3.vert[f->v1].y + mesh3.vert[f->v2].y)/2.0;
		z = (mesh3.vert[f->v1].z + mesh3.vert[f->v2].z)/2.0;
		p1 = 1.0/(1e-6 + distance(x,y,z, mesh3.vert[v1].x,mesh3.vert[v1].y,mesh3.vert[v1].z));
		surfNormal(reg3.iSurf, reg3.iNorm, mesh3.vert[v1].u, mesh3.vert[v1].v, &e1);
		p2 = 1.0/(1e-6 + distance(x,y,z, mesh3.vert[v2].x,mesh3.vert[v2].y,mesh3.vert[v2].z));
		surfNormal(reg3.iSurf, reg3.iNorm, mesh3.vert[v2].u, mesh3.vert[v2].v, &e2);
		p =  1.0/(1e-6 + distance(x,y,z, mesh3.vert[v].x,mesh3.vert[v].y,mesh3.vert[v].z));
		surfNormal(reg3.iSurf, reg3.iNorm, mesh3.vert[v].u, mesh3.vert[v].v, &norm);
		len = p + p1 + p2;
		sect3.v[0].x = (p1*e1.x+p2*e2.x+p*norm.x)/len;
		sect3.v[0].y = (p1*e1.y+p2*e2.y+p*norm.y)/len;
		sect3.v[0].z = (p1*e1.z+p2*e2.z+p*norm.z)/len;
		x += sect3.v[0].x*3.0/len,  y += sect3.v[0].y*3.0/len,  z += sect3.v[0].z*3.0/len;
		mesh3.nPoint++;
		addPoint(x,y,z);
		mesh3.nPoint--;
		mesh3.nPoint--;
	    }

/*	    if (faceIntersectNew(v1, v2, v, f->v1, f->v2, mesh3.nPoint+1) != faceIntersect32(f, v1, v2, v)) {
		printf("\nintersection test failed\nold: %d, new: %d\na=%d, b=%d, c=%d, u=%d, v=%d\n",
			faceIntersect32(f, v1, v2, v), faceIntersectNew(v1, v2, v, f->v1, f->v2, mesh3.nPoint+1),
			v1+1, v2+1, v+1, f->v1+1, f->v2+1);
		mesh3.nPoint+=2;
		dump();
		mesh3.nPoint-=2;
	    }
*/
	    if (1 && faceIntersectNew(v1, v2, v, f->v1, f->v2, mesh3.nPoint+1)) {
		bln++;
		sect3.sect = tree32.vicinityFace[iFace].face;
		workVert = choseWorkVert();
		if (workVert >= 0) {
		    return  bln;
		}
	    }
	    if (0 && faceIntersect32(f, v1, v2, v)) {
		bln++;
		sect3.sect = tree32.vicinityFace[iFace].face;
		workVert = choseWorkVert();
		if (workVert >= 0) {
		    return  bln;
		}
	    }
	}
	iFace++;
    }/* while  iFace */

    return  bln;
}/*checkIntersect*/


static int choseWorkVert(void) {
    int i, best=-1, v[2], b[2]={1,1};
    int v1, v2;
    double p, len=-1.0;

    v[0] = sect3.sect->v1;
    v[1] = sect3.sect->v2;

    v1 = sect3.work->v1;   v2 = sect3.work->v2;
    for (i=0; i<2; i++) {
	if ( (v[i] == v1) || (v[i] == v2) || isBadVert(v[i]) || !goodPlane32(v1, v2, v[i]) )  b[i] = 0;
    }
    for (i=0; i<2; i++) {
	if (b[i]) {
	    p = distanceS(mesh3.vert[v1].x,mesh3.vert[v1].y,mesh3.vert[v1].z,  mesh3.vert[v[i]].x,mesh3.vert[v[i]].y,mesh3.vert[v[i]].z) +
		distanceS(mesh3.vert[v2].x,mesh3.vert[v2].y,mesh3.vert[v2].z,  mesh3.vert[v[i]].x,mesh3.vert[v[i]].y,mesh3.vert[v[i]].z);
	    if ((best<0) || (p < len)) {
		len = p;
		best = v[i];
	    }
	}
    }
    return best;
} /*choseWorkVert*/


static void checkNearEdge(void) {
    int i, iFace, v1, v2, v, bestFace=-1, best;
    PStrucFace3 f;
    StrucVert3  e1,e2;
    double c, p, dist=-1.0, x, y, z, len;
    int vs[2], b[2]={1,1};

    v = mesh3.nPoint;

    for (iFace = 0; iFace < tree32.nVicinityFace; iFace++) {
	f = tree32.vicinityFace[iFace].face;
	v1 = f->v1,  v2 = f->v2;
	makeVector(&e1, v,  v1);
	makeVector(&e2, v2, v1);
	c = e1.x*e2.x + e1.y*e2.y + e1.z*e2.z;
	c /= e2.x*e2.x + e2.y*e2.y + e2.z*e2.z;
	if ((c>=0.0) && (c<=1.0)) {
	    x = (1.0-c)*mesh3.vert[v1].x + c*mesh3.vert[v2].x;
	    y = (1.0-c)*mesh3.vert[v1].y + c*mesh3.vert[v2].y;
	    z = (1.0-c)*mesh3.vert[v1].z + c*mesh3.vert[v2].z;
	    p = distance(x,y,z, mesh3.vert[v].x,mesh3.vert[v].y,mesh3.vert[v].z);
	} else if (c<0.5) {
	    x = mesh3.vert[v1].x;
	    y = mesh3.vert[v1].y;
	    z = mesh3.vert[v1].z;
	    p = distance(x,y,z, mesh3.vert[v].x,mesh3.vert[v].y,mesh3.vert[v].z);
	} else {
	    x = mesh3.vert[v2].x;
	    y = mesh3.vert[v2].y;
	    z = mesh3.vert[v2].z;
	    p = distance(x,y,z, mesh3.vert[v].x,mesh3.vert[v].y,mesh3.vert[v].z);
	}
	if ((bestFace < 0) || (p < dist)) {
	    dist = p;
	    bestFace = iFace;
	}
    }
    if (bestFace >= 0) {
	f = tree32.vicinityFace[bestFace].face;
	vs[0] = f->v1;
	vs[1] = f->v2;
	v1 = sect3.work->v1;   v2 = sect3.work->v2;
	for (i=0; i<2; i++) {
	    if ( (vs[i] == v1) || (vs[i] == v2) || isBadVert(vs[i]) || !goodPlane32(v1, v2, vs[i]) )  b[i] = 0;
	}
	best = -1,  len = -1.0;
	for (i=0; i<2; i++) {
	    if (b[i]) {
		p = distance(mesh3.vert[v].x,mesh3.vert[v].y,mesh3.vert[v].z,  mesh3.vert[vs[i]].x,mesh3.vert[vs[i]].y,mesh3.vert[vs[i]].z);
		if ((best<0) || (p < len)) {
		    len = p;
		    best = vs[i];
		}
	    }
	}
	if (best>=0)  sect3.workVert = best;
    }
} /*checkNearEdge*/


double determ(double a11, double a12, double a21, double a22) {
    return a11*a22 - a12*a21;
} /*determ*/


void surfNormal (int iSurf, int iNorm, double u, double v, StrucVert3 *norm) {
    double x0, y0, z0,  x, y, z;
    double dxdu, dydu, dzdu, dxdv, dydv, dzdv;
    double du, dv, er;

    bounSurf(iSurf, u, v, &x0, &y0, &z0);
    er = 0.1*sizeFace(x0, y0, z0);

    du = (u < (reg3.uMax + reg3.uMin)/2.0) ? reg3.uMax - reg3.uMin : reg3.uMin - reg3.uMax;
    do {
	du /= 2.0;
	bounSurf(iSurf, u+du, v, &x, &y, &z);
    } while (distance(x0,y0,z0, x,y,z) > er);
    dxdu = (x - x0)/du;
    dydu = (y - y0)/du;
    dzdu = (z - z0)/du;

    dv = (v < (reg3.vMax + reg3.vMin)/2.0) ? reg3.vMax - reg3.vMin : reg3.vMin - reg3.vMax;
    do {
	dv /= 2.0;
	bounSurf(iSurf, u, v+dv, &x, &y, &z);
    } while (distance(x0,y0,z0, x,y,z) > er);
    dxdv = (x - x0)/dv;
    dydv = (y - y0)/dv;
    dzdv = (z - z0)/dv;

    norm->x = -determ(dydu, dzdu, dydv, dzdv);
    norm->y = -determ(dzdu, dxdu, dzdv, dxdv);
    norm->z = -determ(dxdu, dydu, dxdv, dydv);

    du = sqrt(norm->x*norm->x + norm->y*norm->y + norm->z*norm->z);
    if (du < 1e-12) {
	if (dxdu*dxdu + dydu*dydu + dzdu*dzdu > dxdv*dxdv + dydv*dydv + dzdv*dzdv) {
	    du = (u < (reg3.uMax + reg3.uMin)/2.0) ? reg3.uMax - reg3.uMin : reg3.uMin - reg3.uMax;
	    dv = (v < (reg3.vMax + reg3.vMin)/2.0) ? reg3.vMax - reg3.vMin : reg3.vMin - reg3.vMax;
	    do {
		du /= 2.0;
		dv /= 2.0;
		bounSurf(iSurf, u+du, v+dv, &x, &y, &z);
	    } while (distance(x0,y0,z0, x,y,z) > er);
	    if (du*dv > 0.0)  iNorm = (iNorm > 0) ? 0 : 1;
	    dxdv = (x - x0)/sqrt(du*du + dv*dv);
	    dydv = (y - y0)/sqrt(du*du + dv*dv);
	    dzdv = (z - z0)/sqrt(du*du + dv*dv);
	    if (v-dv >= reg3.vMin && v-dv <= reg3.vMax)  bounSurf(iSurf, u+du, v-dv, &x, &y, &z),  dv = sqrt(du*du + dv*dv);
	    else /*if (u-du >= reg3.uMin && u-du <= reg3.uMax) */ bounSurf(iSurf, u-du, v+dv, &x, &y, &z),  dv = -sqrt(du*du + dv*dv);
	    //	else bounSurf(iSurf, u+du, v, &x, &y, &z),  dv = du;
	    dxdu = (x - x0)/dv;
	    dydu = (y - y0)/dv;
	    dzdu = (z - z0)/dv;
	} else {
	    du = (u < (reg3.uMax + reg3.uMin)/2.0) ? reg3.uMax - reg3.uMin : reg3.uMin - reg3.uMax;
	    dv = (v < (reg3.vMax + reg3.vMin)/2.0) ? reg3.vMax - reg3.vMin : reg3.vMin - reg3.vMax;
	    do {
		du /= 2.0;
		dv /= 2.0;
		bounSurf(iSurf, u+du, v+dv, &x, &y, &z);
	    } while (distance(x0,y0,z0, x,y,z) > er);
	    if (du*dv > 0.0)  iNorm = (iNorm > 0) ? 0 : 1;
	    dxdu = (x - x0)/sqrt(du*du + dv*dv);
	    dydu = (y - y0)/sqrt(du*du + dv*dv);
	    dzdu = (z - z0)/sqrt(du*du + dv*dv);
	    if (u-du >= reg3.uMin && u-du <= reg3.uMax)  bounSurf(iSurf, u-du, v+dv, &x, &y, &z),  du = -sqrt(du*du + dv*dv);
	    else /*if (u-du >= reg3.uMin && u-du <= reg3.uMax) */ bounSurf(iSurf, u+du, v-dv, &x, &y, &z),  du = sqrt(du*du + dv*dv);
	    //	else bounSurf(iSurf, u+du, v, &x, &y, &z),  dv = du;
	    dxdv = (x - x0)/du;
	    dydv = (y - y0)/du;
	    dzdv = (z - z0)/du;
	}
	norm->x = -determ(dydu, dzdu, dydv, dzdv);
	norm->y = -determ(dzdu, dxdu, dzdv, dxdv);
	norm->z = -determ(dxdu, dydu, dxdv, dydv);

	du = sqrt(norm->x*norm->x + norm->y*norm->y + norm->z*norm->z);
	if (du < 1e-12) errorExit3(3," du == 0.   in  surfNormal ");
    }
    norm->x /= du;
    norm->y /= du;
    norm->z /= du;

    if (iNorm > 0) {
	norm->x = -norm->x;
	norm->y = -norm->y;
	norm->z = -norm->z;
    }

    return;
} /*surfNormal*/


void rePlane(int vv) {
    double x, y, z, u, v;

    x = mesh3.vert[vv].x - reg3.x0;
    y = mesh3.vert[vv].y - reg3.y0;
    z = mesh3.vert[vv].z - reg3.z0;

    u = x*reg3.x1 + y*reg3.y1 + z*reg3.z1;
    v = x*reg3.x2 + y*reg3.y2 + z*reg3.z2;

    mesh3.vert[vv].u = u;
    mesh3.vert[vv].v = v;

    return;
} /*rePlane*/


int goodPlane32(int v1, int v2, int v) {
    double p;
    StrucVert3 norm, e1, e2;

    makeVector(&e1, v, v2);
    makeVector(&e2, v2, v1);
    makeNormal(&norm, &e1, &e2);
    surfNormal(reg3.iSurf, reg3.iNorm, mesh3.vert[v1].u, mesh3.vert[v1].v, &e1);
    surfNormal(reg3.iSurf, reg3.iNorm, mesh3.vert[v2].u, mesh3.vert[v2].v, &e2);
    e1.x += e2.x;
    e1.y += e2.y;
    e1.z += e2.z;
    p = norm.x*e1.x + norm.y*e1.y + norm.z*e1.z;
    if (p < 1.e-7) return 0;
    return 1;
} /*goodPlane32*/        

static int shiftperiodic(int k) {
    double pu, pv;
    double x, y, z;
    int i, q, j;
//    if (reg3.iSurf >= 0)  return 0;
    pu = periodic(0),  pv = periodic(1);
    if (pu > 0.0) {
	reg3.uMax = mesh3.vert[k].u + pu/2.0;
	reg3.uMin = mesh3.vert[k].u - pu/2.0;
	for (j=0; j<tree32.nFace; j++) {
	    i = tree32.face[j]->v1;
	    if (i==k)  continue;
	    q = 0;
	    while (mesh3.vert[i].u > reg3.uMax)  mesh3.vert[i].u -= pu,  q--;
	    while (mesh3.vert[i].u < reg3.uMin)  mesh3.vert[i].u += pu,  q++;
	    if (q) {
		bounSurf(reg3.iSurf, mesh3.vert[i].u, mesh3.vert[i].v, &x, &y, &z);
		if (distance(mesh3.vert[i].x, mesh3.vert[i].y, mesh3.vert[i].z,  x, y, z) > 1e-3)
		    printf("U period: %lf. Shift: %d. d = %le. old: %lf, %lf, %lf. new: %lf, %lf, %lf [%lf, %lf].\n",
			    pu, q, distance(mesh3.vert[i].x, mesh3.vert[i].y, mesh3.vert[i].z,  x, y, z),
			    mesh3.vert[i].x, mesh3.vert[i].y, mesh3.vert[i].z,  x, y, z, mesh3.vert[i].u, mesh3.vert[i].v);
//			mesh3.vert[i].u -= pu*q;
	    }
	}
    }
    if (pv > 0.0) {
	reg3.vMax = mesh3.vert[k].v + pv/2.0;
	reg3.vMin = mesh3.vert[k].v - pv/2.0;
	for (j=0; j<tree32.nFace; j++) {
	    i = tree32.face[j]->v1;
	    if (i==k)  continue;
	    q = 0;
	    while (mesh3.vert[i].v > reg3.vMax)  mesh3.vert[i].v -= pv,  q--;
	    while (mesh3.vert[i].v < reg3.vMin)  mesh3.vert[i].v += pv,  q++;
	    if (q) {
		bounSurf(reg3.iSurf, mesh3.vert[i].u, mesh3.vert[i].v, &x, &y, &z);
		if (distance(mesh3.vert[i].x, mesh3.vert[i].y, mesh3.vert[i].z,  x, y, z) > 1e-3)
		    printf("V period: %lf. Shift: %d. d = %le. old: %lf, %lf, %lf. new: %lf, %lf, %lf [%lf, %lf].\n",
			    pv, q, distance(mesh3.vert[i].x, mesh3.vert[i].y, mesh3.vert[i].z,  x, y, z),
			    mesh3.vert[i].x, mesh3.vert[i].y, mesh3.vert[i].z,  x, y, z, mesh3.vert[i].u, mesh3.vert[i].v);
//			mesh3.vert[i].v -= pv*q;
	    }
	}
    }
    return 0;
}

#define NQ 256
static int newPoint(PStrucFace3 face) {
    int i, j, v1, v2, nearVert;
    int nTest=0, sec=0;
    int bound=0;
    double weakness = 1.0, realsize;
    double x, y, z,  x1, y1, z1,  x2, y2, z2, size, dist, d1, d2;
    StrucVert3 norm, e3, e4;

    double fi1, p;
    double u, v, du, dv;

    double fi, fiMin, fiMax;
    double ra, raMin, raMax;
    double sin0, cos0, si, co;


    v1 = face->v1;  v2 = face->v2;
    shiftperiodic(v1);
    x1 = mesh3.vert[v1].x;  y1 = mesh3.vert[v1].y;  z1 = mesh3.vert[v1].z;
    x2 = mesh3.vert[v2].x;  y2 = mesh3.vert[v2].y;  z2 = mesh3.vert[v2].z;

    u = mesh3.vert[v1].u;
    v = mesh3.vert[v1].v;
    bounSurf(reg3.iSurf, u, v, &x, &y, &z);
    surfNormal(reg3.iSurf, reg3.iNorm, u, v, &surfNorm);
    reg3.uSave = mesh3.vert[v1].u;
    reg3.vSave = mesh3.vert[v1].v;
    dist = distance(x1, y1, z1,  x2, y2, z2);
    if (CF2_5D <= 0.0)  size = sizeFace(0.5*(x1+x2), 0.5*(y1+y2), 0.5*(z1+z2));
    else  size = CF2_5D * dist;
    //	size = dist*1.3;
    if (size < 0.5*dist)  {
	size = 0.6*dist;
	printf("\n...increasing local mesh size...\n");
    }
    realsize = size;

    du = mesh3.vert[v2].u - mesh3.vert[v1].u;
    dv = mesh3.vert[v2].v - mesh3.vert[v1].v;
    p = sqrt(du*du + dv*dv);
    if (p == 0.)  {
	printf("\nv1=%d, v2=%d, dist=%lf\n", v1, v2, dist);
	printf("%lf, %lf  :: %lf, %lf\n", mesh3.vert[v1].u, mesh3.vert[v1].v, mesh3.vert[v2].u, mesh3.vert[v2].v);
	printf("%lf, %lf, %lf\n", x1, y1, z1);
	printf("%lf, %lf, %lf\n", x2, y2, z2);
	printf("%lf, %lf, %lf\n", x, y, z);
	if (0)  dump();
	printf("p == 0. in newPoint \n");
	return -21;
//	errorExit3(2, "p == 0. in newPoint ");
    }
    cos0 = du/p;
    sin0 = dv/p;
    sec = -1;
    while (1) {
	sec++;
	if (sec>2*N_SPLIT) {
	    bound++;
//	    fprintf(stderr, "\nwarning! decreasing size!\n");
	    size = (dist/2.0 + 10.0*size)/11.0;
	    if (fabs(dist-2.0*size)/dist < 1e-5)  {
		weakness *= 2;
		size = realsize;
//		if (weakness > 5.0)  fprintf(stderr, "\nwarning! weakness = %lf!\n", weakness);
		sec = -1;
		size = realsize;
		if (weakness > 20.0) {
		    if (0)  dump();
		    return -20;
		    co = cos0;
		    si = sin0;
		    ra = dist*0.5;
		    break;
		}
		continue;
/*		if (!weak) {
		    fprintf(stderr, "\nwarning! weak conditions!!!\n");
		    weak = 1;
		    sec = -1;
		    continue;
		} else {
		    co = cos0;
		    si = sin0;
		    ra = dist*0.5;
		    break;
		}*/
	    } else {
		sec = -1;
		continue;
	    }
	}
	u = reg3.uSave;
	v = reg3.vSave;
	if (sec % 2 == 0) {
	    fiMin = M_PI/N_SPLIT * (sec/2);
	    fiMax = M_PI/N_SPLIT * (sec/2 + 1);
	} else {
	    fiMin = 2*M_PI - M_PI/N_SPLIT * (sec/2);
	    fiMax = 2*M_PI - M_PI/N_SPLIT * (sec/2 + 1);
	}
	for (i=0; i<N_ANGLE; i++) {
	    fi = 0.5*(fiMin + fiMax);
	    co = cos0*cos(fi) - sin0*sin(fi);
	    si = sin0*cos(fi) + cos0*sin(fi);
	    raMin = 0;
	    raMax = 2.0*sqrt((reg3.uMax-reg3.uMin)*(reg3.uMax-reg3.uMin) + (reg3.vMax-reg3.vMin)*(reg3.vMax-reg3.vMin));
	    for (j=0; j<N_DEEP; j++) {/*calculate ra for fi*/
		ra = 0.5*(raMin + raMax);
		if ((u+ra*co > reg3.uMax) || (u+ra*co < reg3.uMin) || (v+ra*si > reg3.vMax) || (v+ra*si < reg3.vMin)) {
		    raMax = ra;
		    continue;
		}
		bounSurf(reg3.iSurf,u+ra*co,v+ra*si,&x,&y,&z);
		d1 = distance(x1,y1,z1,x,y,z);
		if (fabs(d1-size)/size < 0.0001)
		    break;
		else if (d1 > size)
		    raMax = ra;
		else
		    raMin = ra;
	    }
	    if ((u+ra*co > reg3.uMax) || (u+ra*co < reg3.uMin) || (v+ra*si > reg3.vMax) || (v+ra*si < reg3.vMin)) {
		fiMax = fi;
	    }
	    d1 = distance(x1,y1,z1,x,y,z);
	    d2 = distance(x2,y2,z2,x,y,z);
	    if (fabs(d2-d1)/d1 < 0.0001)
		break;
	    else if (d1 < d2)
		fiMax = fi;
	    else
		fiMin = fi;
	}
	if ((u+ra*co > reg3.uMax) || (u+ra*co < reg3.uMin) || (v+ra*si > reg3.vMax) || (v+ra*si < reg3.vMin))  continue;
//	if ((!weak)&&(( fabs(size-distance(x1,y1,z1,x,y,z))/size > 1e-3 ) || (fabs(size-distance(x2,y2,z2,x,y,z))/size > 1e-3 )))  continue;
	if (( fabs(size-distance(x1,y1,z1,x,y,z))/size > 1e-3 * weakness ) || (fabs(size-distance(x2,y2,z2,x,y,z))/size > 1e-3 * weakness ))  continue;
	u += ra*co;
	v += ra*si;
	if (bounSurf(reg3.iSurf,u,v,&x,&y,&z) == -1)  continue;
	mesh3.vert[mesh3.nPoint].u = u;
	mesh3.vert[mesh3.nPoint].v = v;
	addPoint(x,y,z);
	mesh3.nPoint--;

	/*   printf("v= %5d   %5.5lf  %5.5lf  x=%5.5lf  y=%5.5lf  z=%5.5lf \n",mesh3.nPoint,u,v,x0,y0,z0); */
	makeVector(&e3, mesh3.nPoint, v2);
	makeVector(&e4, v2, v1);
	makeNormal(&norm, &e3, &e4);
	surfNormal(reg3.iSurf, reg3.iNorm, mesh3.vert[v1].u, mesh3.vert[v1].v, &e3);
	fi1 = norm.x*e3.x + norm.y*e3.y + norm.z*e3.z;
	if (fi1 < 1.e-7)  continue;
	break;
    }

    if (bounSurf(reg3.iSurf,mesh3.vert[mesh3.nPoint].u,mesh3.vert[mesh3.nPoint].v,&x,&y,&z) == -1)  printf("\nbounsurf failed!\n");
//    if ((fabs(size-distance(x1,y1,z1,x,y,z))/size > 1e-3) || (fabs(size-distance(x2,y2,z2,x,y,z))/size > 1e-3))
    if (( fabs(size-distance(x1,y1,z1,x,y,z))/size > 1e-2 ) || (fabs(size-distance(x2,y2,z2,x,y,z))/size > 1e-2 )) {
	printf("	size = %lf, s1 = %lf (%lf), s2 = %lf (%lf)\n", size, distance(x1,y1,z1,x,y,z), (size-distance(x1,y1,z1,x,y,z))/size,
		distance(x2,y2,z2,x,y,z), (size-distance(x2,y2,z2,x,y,z))/size);
	if (0)  dump();
	return -10;
    }
    sect3.workVert = mesh3.nPoint;
    sect3.work = face;

    /******************  TEST  ****************/
    x = (mesh3.vert[v1].x+mesh3.vert[v2].x+mesh3.vert[mesh3.nPoint].x)/3.;
    y = (mesh3.vert[v1].y+mesh3.vert[v2].y+mesh3.vert[mesh3.nPoint].y)/3.;
    z = (mesh3.vert[v1].z+mesh3.vert[v2].z+mesh3.vert[mesh3.nPoint].z)/3.;
    size = max(distance(x,y,z, mesh3.vert[v1].x,mesh3.vert[v1].y,mesh3.vert[v1].z),
	    distance(x,y,z, mesh3.vert[v2].x,mesh3.vert[v2].y,mesh3.vert[v2].z));
    size = max(size, distance(x,y,z,mesh3.vert[mesh3.nPoint].x,mesh3.vert[mesh3.nPoint].y,mesh3.vert[mesh3.nPoint].z));
    size *= 1.25;
    vicinityFaces32(x, y, z, size);
    size /= 1.25;

    if (1) {  /* sphereVert */
	int          i,iFace,v[2];
	double       min,p,rad;
	PStrucFace3  face;

/*	x = 0.5*(mesh3.vert[v1].x+mesh3.vert[v2].x);
	y = 0.5*(mesh3.vert[v1].y+mesh3.vert[v2].y);
	z = 0.5*(mesh3.vert[v1].z+mesh3.vert[v2].z);
	radS = distanceS(x,y,z,mesh3.vert[sect3.workVert].x,mesh3.vert[sect3.workVert].y,mesh3.vert[sect3.workVert].z);  
*/
	rad = 1.25*size;
	nSphereVert = 0;
	for(iFace=0;iFace<tree32.nVicinityFace;iFace++){
	    face = tree32.vicinityFace[iFace].face;
	    v[0] = face->v1;
	    v[1] = face->v2;
	    for (i=0; i<2; i++) {
		if ((v[i]==v1) || (v[i]==v2))  continue;
		if (isSphereVert(v[i]))  continue;
		if (!goodPlane32(v1,v2,v[i]))  continue;
		p = distance(x,y,z, mesh3.vert[v[i]].x,mesh3.vert[v[i]].y,mesh3.vert[v[i]].z);  
		if (p < rad)  addSphereVert(v[i]);
	    }
	}
	min = -1.0;
	v[0] = -1;
	for (i=0; i<nSphereVert; i++) {
	    p = distance(x,y,z, mesh3.vert[sphereVert[i]].x,mesh3.vert[sphereVert[i]].y,mesh3.vert[sphereVert[i]].z);
	    if ((v[0]<0) || (p < min)) {
		min = p;
		v[0] = sphereVert[i];
	    }  
	}
	if ( v[0] != -1 )
	    sect3.workVert = v[0];
    }  /* sphereVert */     

    x = mesh3.vert[mesh3.nPoint].x;
    y = mesh3.vert[mesh3.nPoint].y;
    z = mesh3.vert[mesh3.nPoint].z;
    if (sect3.workVert == mesh3.nPoint) {
	dist = nearest32(&nearVert,x,y,z);
	if ( (dist < minNear*size) && (nearVert != v1) && (nearVert != v2) ) {
	    if (goodPlane32(v1,v2,nearVert))  sect3.workVert = nearVert;
	}
    }
    if ( (sect3.workVert == mesh3.nPoint) && (bound>0) )  checkNearEdge();
    nBadVert = 0;
    nTest = 0;

    x = (mesh3.vert[v1].x+mesh3.vert[v2].x+mesh3.vert[mesh3.nPoint].x)/3.;
    y = (mesh3.vert[v1].y+mesh3.vert[v2].y+mesh3.vert[mesh3.nPoint].y)/3.;
    z = (mesh3.vert[v1].z+mesh3.vert[v2].z+mesh3.vert[mesh3.nPoint].z)/3.;

    while (1) {
	nTest++;
	if( nTest > NQ )
		errorExit3(2," nTest > NQ ");

	if (size < distance(x,y,z,mesh3.vert[sect3.workVert].x,mesh3.vert[sect3.workVert].y,mesh3.vert[sect3.workVert].z)) {
	    size = distance(x,y,z,mesh3.vert[sect3.workVert].x,mesh3.vert[sect3.workVert].y,mesh3.vert[sect3.workVert].z);
	    vicinityFaces32(x, y, z, size*1.25);
	}
	if (checkIntersect(face, sect3.workVert) == 0)  break;

	addBadVert(sect3.workVert);
	sect3.workVert = choseWorkVert();

	if (sect3.workVert < 0) {
	    double  p, len=-1.0;
	    int best=-1;

	    for (i=0; i<nSphereVert; i++) {
		if (isBadVert(sphereVert[i]))  continue;
		p = distanceS(mesh3.vert[v1].x,mesh3.vert[v1].y,mesh3.vert[v1].z,
			mesh3.vert[sphereVert[i]].x,mesh3.vert[sphereVert[i]].y,mesh3.vert[sphereVert[i]].z) +
		    distanceS(mesh3.vert[v2].x,mesh3.vert[v2].y,mesh3.vert[v2].z,
			    mesh3.vert[sphereVert[i]].x,mesh3.vert[sphereVert[i]].y,mesh3.vert[sphereVert[i]].z);
		if ((best<0) || (p<len)) {
		    len = p;
		    best = sphereVert[i];
		}  
	    }
	    if (best >= 0)  sect3.workVert = best;
	    else {
		printf("\nintersection in TRIA32\n");
		printf("Please, try with smaller mesh size\n");
		if (0)  dump3();
		return -1;
	    }
	}
    }
    return 1;
} /*newPoint*/

static void newTria(void) {
    static int fresh = 0;
    int i;
    int v[3];
    PStrucFace3 face;

    for (i=0; i<tree32.nFace; i++) {
	face = tree32.face[i];
	if (face->fail)  continue;
	if (newPoint(face) > 0)  break;
	face->fail = 1;
    }
    if (i>=tree32.nFace)  {
	if (!fresh) {
	printf("\nFront restart\n");
	for (i=0; i<tree32.nFace; i++)  face->fail = 0;
	for (i=0; i<tree32.nFace; i++) {
	    face = tree32.face[i];
	    if (face->fail)  continue;
	    if (newPoint(face) > 0)  break;
	    face->fail = 1;
	}
	fresh = 1;
	} else {
	    if (0)  dump();
	    errorExit3(3,"newTria failed");
	}
    }

    v[0] = face->v1;
    v[1] = sect3.workVert;
    v[2] = face->v2;
    /*printf("%5d: %5d %5d %5d\n",reg3.nTria,v[0],v[2],v[1]);*/
    remFace32(face);
    addTria(v[0], v[2], v[1]);

    if (sect3.workVert == mesh3.nPoint) {
	if (  (fabs(mesh3.vert[mesh3.nPoint].x-tree3.boxcx)>tree3.boxsize*1.0000001) ||
		(fabs(mesh3.vert[mesh3.nPoint].y-tree3.boxcy)>tree3.boxsize*1.0000001) ||
		(fabs(mesh3.vert[mesh3.nPoint].z-tree3.boxcz)>tree3.boxsize*1.0000001)  ) {
	    printf("\nbox: x: %le, y: %le, z: %le, size: %le\n", tree3.boxcx, tree3.boxcy, tree3.boxcz, tree3.boxsize);
	    printf("x: %le, y: %le, z: %le\n", mesh3.vert[mesh3.nPoint].x, mesh3.vert[mesh3.nPoint].y, mesh3.vert[mesh3.nPoint].z);
	    errorExit3(3,"x_y_z_Point");
	}
	mesh3.nPoint++;
    }
    for (i=0; i<2; i++) {
	if (sect3.neigBool[i]) {
	    remFace32(sect3.neigFace[i]);
	} else {
	    addFace32(v[i], v[i+1]);
	}
    }
    fresh = 0;
    return;
} /*newTria*/


static int dump() {
	static int num=0;
	char fname[1024];
	FILE *f;
	int i;
	sprintf(fname, "_dump_frt2_smv.%03d", num);
	f = fopen(fname, "w");

	fprintf(f, "%d %d %d %d 1\n", mesh3.nPoint, reg3.nTria, tree32.nFace, tree32.nFace+2);
	for (i=0; i<mesh3.nPoint; i++) fprintf(f, "%20.15lf %20.15lf %20.15lf\n", mesh3.vert[i].x, mesh3.vert[i].y, mesh3.vert[i].z);
	for (i=0; i<reg3.nTria; i++) fprintf(f, "%d %d %d\n", reg3.v1[i]+1, reg3.v2[i]+1, reg3.v3[i]+1);
	for (i=0; i<tree32.nFace; i++) fprintf(f, "%d %d\n", tree32.face[i]->v1+1, tree32.face[i]->v2+1);
	for (i=0; i<tree32.nFace; i++) fprintf(f, "%d\n", tree32.face[i]->v1+1);
	fprintf(f, "%d A\n", tree32.face[0]->v1+1);
	fprintf(f, "%d B\n", tree32.face[0]->v2+1);
	fclose(f);
	/*
	sprintf(fname, "frt2_gmv.%03d", num);
	f = fopen(fname, "w");
	fprintf(f, "gmvinput ascii\n\n nodes %5d\n", mesh3.nPoint);
	for (i=0; i<mesh3.nPoint; i++) fprintf(f, " %20.15lf", mesh3.vert[i].x);
	fprintf(f, "\n");
	for (i=0; i<mesh3.nPoint; i++) fprintf(f, " %20.15lf", mesh3.vert[i].y);
	fprintf(f, "\n");
	for (i=0; i<mesh3.nPoint; i++) fprintf(f, " %20.15lf", mesh3.vert[i].z);
	fprintf(f, "\n");
	fprintf(f, "\n cells %5d\n", tree32.nFace);
	for (i=0; i<tree32.nFace; i++) {
		fprintf(f, "  line 2\n  %4d %4d\n", tree32.face[i]->v1+1, tree32.face[i]->v2+1);
	}
	fprintf(f, "\nendgmv");
	fclose(f);
	*/
	return num++;
}

static int dump3() {
	static int num=0;
	char fname[1024];
	FILE *f;
	int i;
	sprintf(fname, "_dump_frt3_smv.%03d", num);
	f = fopen(fname, "w");

	fprintf(f, "%d %d %d %d 1\n", mesh3.nPoint+1, reg3.nTria, tree32.nFace, tree32.nFace+3);
	for (i=0; i<=mesh3.nPoint; i++) fprintf(f, "%20.15lf %20.15lf %20.15lf\n", mesh3.vert[i].x, mesh3.vert[i].y, mesh3.vert[i].z);
	for (i=0; i<reg3.nTria; i++) fprintf(f, "%d %d %d\n", reg3.v1[i]+1, reg3.v2[i]+1, reg3.v3[i]+1);
	for (i=0; i<tree32.nFace; i++) fprintf(f, "%d %d\n", tree32.face[i]->v1+1, tree32.face[i]->v2+1);
	for (i=0; i<tree32.nFace; i++) fprintf(f, "%d\n", tree32.face[i]->v1+1);
	fprintf(f, "%d A\n", tree32.face[0]->v1+1);
	fprintf(f, "%d B\n", tree32.face[0]->v2+1);
	fprintf(f, "%d C\n", (sect3.workVert >= 0) ? sect3.workVert+1 : mesh3.nPoint+1);
	fclose(f);
	return num++;
}


void makeTria(void) {
	if (tria_dump_front) dump();
	reg3.suMax = reg3.uMax,  reg3.suMin = reg3.uMin;
	reg3.svMax = reg3.vMax,  reg3.svMin = reg3.vMin;
	while (tree32.nFace > 0) {
		newTria();
#ifdef SHOWPROGRESS
		printf("\r nP = %5d  nT = %5d  nE = %5d", mesh3.nPoint, reg3.nTria, tree32.nFace);
		fflush(stdout);
#endif
		if (tria_debug_front) dump();
	}
#ifdef SHOWPROGRESS
	printf("\r nP = %5d  nT = %5d  done.        \n", mesh3.nPoint, reg3.nTria);
#endif
	if (tria_final_front) dump();
	reg3.uMax = reg3.suMax,  reg3.uMin = reg3.suMin;
	reg3.vMax = reg3.svMax,  reg3.vMin = reg3.svMin;
	return;
} /*makeTria*/
