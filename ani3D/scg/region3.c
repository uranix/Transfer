#include "error3.h"
#include "memory3.h"
#include "tree3.h"
#include "tree32.h"
#include "tria32.h"
#include "refine3.h"
#include "region3.h"
#include "user3.h"
#include "det.h"

int region_dump_face = 0;

StrucRegion3  reg3;
StrucCrvT     CrvT;
int           countCrvT;
static double precision = 1e-3;

double  surfacesizeratio = 0.05;

/* extern  variables */
extern StrucMesh3 mesh3;
extern StrucTree3 tree3;
extern StrucTree3 tree32;

extern double S0;

double sizeFace (double x, double y, double z) {
	return userSizeFace(x,y,z);
} /* sizeFace  */


int bounSurf0 (double u, double v, double *x, double *y, double *z) {
   x[0] = reg3.x0 + u*reg3.x1 + v*reg3.x2;
   y[0] = reg3.y0 + u*reg3.y1 + v*reg3.y2;
   z[0] = reg3.z0 + u*reg3.z1 + v*reg3.z2;
   return  1;
} /*bounSurf0*/


void V_U0 (double u, double *v) {
	(void)u;
   v[0] = 0.;
   return;
} /*V_U0*/


static double normvec(double x1, double y1, double z1, double x2, double y2, double z2, double *px, double *py, double *pz) {
	double w1[4], w2[4], w3[4], x, y, z;
	w1[0] = x1, w2[0] = x2, w3[0] = 0.0;
	w1[1] = y1, w2[1] = y2, w3[1] = 0.0;
	w1[2] = z1, w2[2] = z2, w3[2] = 0.0;
	w1[3] = x1, w2[3] = x2, w3[3] = 0.0;
	*px = x = orient2d(w3+1, w2+1, w1+1);
	*py = y = orient2d(w3+2, w2+2, w1+2);
	*pz = z = orient2d(w3+0, w2+0, w1+0);
	return x*x + y*y + z*z;
}


void initAFT_ (int *nF, int *faces, int *nV, double *vertices) {
	int v1, v2, v3, i;
#ifdef WARN
	double vol=0.0;
#endif

	init();
		
	for (i=0; i<*nV; i++) {
		addPoint(vertices[3*i+0], vertices[3*i+1], vertices[3*i+2]);
	}
	mesh3.nBoundPoint = *nV;
	for (i=0; i<*nF; i++) {
		v1 = faces[3*i+0] - 1;
		v2 = faces[3*i+1] - 1;
		v3 = faces[3*i+2] - 1;
		addFace(v1, v2, v3, 0, 1);
#ifdef WARN
		vol+=(mesh3.vert[v1].x * (mesh3.vert[v2].y*mesh3.vert[v3].z - mesh3.vert[v2].z*mesh3.vert[v3].y) +
				mesh3.vert[v1].y * (mesh3.vert[v2].z*mesh3.vert[v3].x - mesh3.vert[v2].x*mesh3.vert[v3].z) +
				mesh3.vert[v1].z * (mesh3.vert[v2].x*mesh3.vert[v3].y - mesh3.vert[v2].y*mesh3.vert[v3].x));
#endif
	}
#ifdef WARN
	if (vol < 0.0)
		printf("initAFT: wrong orientation of faces.\n");
#endif
	return;
}

double nextU (int iSurf, int iV_U, double t) {
	int i=0;
	double v, t0, t1, x0, y0, z0, x1, y1, z1, x, y, z, s1=0., s, er;

	t0 = t;
	t1 = reg3.uEnd;
	V_U(iV_U, t0, &v);
	bounSurf(iSurf, t0, v, &x0, &y0, &z0);
	V_U(iV_U, t1, &v);
	bounSurf(iSurf, t1, v, &x1, &y1, &z1);

	s = sizeFace(0.5*(x0+x1), 0.5*(y0+y1), 0.5*(z0+z1));
	er = 0.01*s;
	if (distance(x0, y0, z0, x1, y1, z1) <= s)
		return  -1.0e11;
	while (fabs(s1-s) >= er) {
		t = 0.5*(t0+t1);
		V_U(iV_U, t, &v);
		bounSurf(iSurf, t, v, &x, &y, &z);
		s1 = distance(x0, y0, z0, x, y, z);
		s = sizeFace(0.5*(x0+x), 0.5*(y0+y), 0.5*(z0+z));
		er = 0.01*s;
		if (s1>s) t1=t; else t0=t;
		i++;
		if (i > 100) {
			if (fabs(s1-s) >= 10.*er)
				errorExit3(3, "i>100 in  nextU");
			break;
		}
	}

	return t;
} /* nextU */

void smoothingLine (int *vert, double *pU, double *pV, int iSurf, int iV_U, int n) {
	int i, j=0, jj;
	double t, tp, tn, tt, sp, sn;
	double x, y, z, x0, y0, z0, x1, y1, z1, x2, y2, z2;
	double du, k=0.0, k0, k1, fi0_, fi1_, f1, p;
	int v, vp, vn;

	while (j < 5) {
		for (i=n-2; i>0; i--) {
			v  = vert[i];
			vp = vert[i-1];  tp = pU[i-1];
			vn = vert[i+1];  tn = pU[i+1];
			if (tp > tn) {
				t  = tp;
				tp = tn;
				tn = t;
			}
			x1 = mesh3.vert[vp].x;
			y1 = mesh3.vert[vp].y;
			z1 = mesh3.vert[vp].z;
			x2 = mesh3.vert[vn].x;
			y2 = mesh3.vert[vn].y;
			z2 = mesh3.vert[vn].z;

			sp = sizeFace(0.5*(mesh3.vert[v].x+x1), 0.5*(mesh3.vert[v].y+y1), 0.5*(mesh3.vert[v].z+z1));
			sn = sizeFace(0.5*(mesh3.vert[v].x+x2), 0.5*(mesh3.vert[v].y+y2), 0.5*(mesh3.vert[v].z+z2));

			/* find  point  on  surface */
			t = pU[i];
			V_U(iV_U, t, &tt);
			bounSurf(iSurf, t, tt, &x, &y, &z);
			f1 = (sp/sn - distance(x1, y1, z1, x, y, z)/distance(x2, y2, z2, x, y, z))*
				(sp/sn - distance(x1, y1, z1, x, y, z)/distance(x2, y2, z2, x, y, z));
			for (jj=0; jj<222; jj++) { /* NEW METH */
				if (f1 < 0.01) break;
				du = 0.2; /*???*/
				do {
					du /= 2.0;
					if (iSurf>=0) {
						V_U(iV_U, t+du, &tt);
						bounSurf(iSurf, t+du, tt, &x0, &y0, &z0);
					} else {
#ifdef CGM
						CGM_GetEdgeCoordsFromU(edge, t+du, xyz);
						x0 = xyz[0]; y0 = xyz[1]; z0 = xyz[2];
#endif
					}
					p = distanceS(x, y, z, x0, y0, z0);
				} while (p > 0.1*sp*sn);

				k0 = -0.9;
				k1 = 1.0;
				if (du == 0.0)
					errorExit3(3, "du == 0.0 ");
				if (t + k0*du < tp)
					k0 = (tp-t)/du;
				if (t + k1*du > tn)
					k1 = (tn-t)/du;
				while (k1 - k0 > 0.00001) {/*min of funk*/

/*					if( iSurf > 0 ){   double  k;
						printf("fun : \n");
						for(i=0;i<10;i++){
							k = k0 + i*(k1-k0)/10;
							V_U(iV_U,t+k*du,&tt);
							bounSurf(iSurf,t+k*du,tt,&x,&y,&z);
							fi0_ = (sp/sn-distance(x1,y1,z1,x,y,z)/distance(x2,y2,z2,x,y,z))*
								(sp/sn-distance(x1,y1,z1,x,y,z)/distance(x2,y2,z2,x,y,z));
							printf("k = %5.5f  fun = %9.9f  \n",k,fi0_);
						}
						printf("\n\n\n");
					}*/

					k = 0.5*(k0+k1);
					V_U(iV_U, t+k*du, &tt);
					bounSurf(iSurf, t+k*du, tt, &x, &y, &z);
					fi0_ = (sp/sn - distance(x1, y1, z1, x, y, z)/distance(x2, y2, z2, x, y, z))*
						(sp/sn - distance(x1, y1, z1, x, y, z)/distance(x2, y2, z2, x, y, z));
					V_U(iV_U, t+(k+fabs(k)*0.01)*du, &tt);
					bounSurf(iSurf, t+(k+fabs(k)*0.01)*du, tt, &x, &y, &z);
					fi1_ = (sp/sn - distance(x1, y1, z1, x, y, z)/distance(x2, y2, z2, x, y, z))*
						(sp/sn - distance(x1, y1, z1, x, y, z)/distance(x2, y2, z2, x, y, z));
					if (fi1_ - fi0_ > 0.)
						k1 = k;
					else
						k0 = k;
					}
				t += k*du;
				V_U(iV_U, t, &tt);
				bounSurf(iSurf, t, tt, &x, &y, &z);
				f1 = (sp/sn - distance(x1, y1, z1, x, y, z)/distance(x2, y2, z2, x, y, z))*
					(sp/sn - distance(x1, y1, z1, x, y, z)/distance(x2, y2, z2, x, y, z));
			} /*NEW METH*/
			pU[i] = t;
			V_U(iV_U, t, &tt);
			bounSurf(iSurf, t, tt, &x, &y, &z);
			pV[i] = tt;
			mesh3.vert[v].x = x;
			mesh3.vert[v].y = y;
			mesh3.vert[v].z = z;
		}
		j++;
	}
	return;
} /* smoothingLine */

void adjustingLine (int *vert, double *pU, double *pV, int iSurf, int iV_U, int n) {
	int i, j;
	double t, tp, tn, ta, tb, tt;
	double x, y, z, x1, y1, z1, x2, y2, z2, xa, ya, za, xb, yb, zb, xt, yt, zt;
	int v, vp, vn;

	vp = vert[0];   tp = pU[0];
	vn = vert[n-1]; tn = pU[n-1];
	x1 = mesh3.vert[vp].x;
	y1 = mesh3.vert[vp].y;
	z1 = mesh3.vert[vp].z;
	x2 = mesh3.vert[vn].x;
	y2 = mesh3.vert[vn].y;
	z2 = mesh3.vert[vn].z;

	if (distance(x2, y2, z2, x1, y1, z1) == 0)
		errorExit3(2, " bounds for V_U parametrization point to the same node ");

	t = pU[0];
	V_U(iV_U, t, &tt);
	pV[0] = tt;
	bounSurf(iSurf, t, tt, &x, &y, &z);
	if (distance(x1, y1, z1, x, y, z) > precision) {
		printf("x1 = %12.11lf y1 = %12.11lf z1 = %12.11lf \n", x1, y1, z1);
		printf("x2 = %12.11lf y2 = %12.11lf z2 = %12.11lf \n", x2, y2, z2);
		printf("x = %12.11lf y = %12.11lf z = %12.11lf \n", x, y, z);
		errorExit3(2, " bounds for V_U parametrization do not correspond, pU[0] "); 
	}
	t = pU[n-1];
	V_U(iV_U, t, &tt);
	pV[n-1] = tt;
	bounSurf(iSurf, t, tt, &x, &y, &z);
	if (distance(x2, y2, z2, x, y, z) > precision) {
		printf("x1 = %12.11lf y1 = %12.11lf z1 = %12.11lf \n", x1, y1, z1);
		printf("x2 = %12.11lf y2 = %12.11lf z2 = %12.11lf \n", x2, y2, z2);
		printf("x = %12.11lf y = %12.11lf z = %12.11lf \n", x, y, z);
		errorExit3(2, " bounds for V_U parametrization do not correspond, pU[last] ");
	} 

	for (i=n-2; i>0; i--) {
		v  = vert[i];
		x  = mesh3.vert[v].x;
		y  = mesh3.vert[v].y;
		z  = mesh3.vert[v].z;
		xa = x1;
		ya = y1;
		za = z1;
		xb = x2;
		yb = y2;
		zb = z2;
		ta = tp;
		tb = tn;

		for (j=0; j<40; j++) {
			if (distance(xa, ya, za, x, y, z) < precision) {
				t = ta; break;
			}
			if (distance(xb, yb, zb, x, y, z) < precision) {
				t = tb; break;
			}
			t = (ta + tb)/2.0;
			V_U(iV_U, t, &tt);
			bounSurf(iSurf, t, tt, &xt, &yt, &zt);
			if (distance(xa, ya, za, x, y, z) < distance(xb, yb, zb, x, y, z)) {
				if (  (distance(xa, ya, za, xt, yt, zt) > distance(xa, ya, za, x, y, z)) &&
						(distance(xa, ya, za, xt, yt, zt) > distance(xt, yt, zt, x, y, z))  ) {
					tb = t; xb = xt; yb = yt; zb = zt;
				} else {
					ta = t; xa = xt; ya = yt; za = zt;
				}
			} else {
				if (  (distance(xt, yt, zt, xb, yb, zb) > distance(xt, yt, zt, x, y, z)) &&
						(distance(xt, yt, zt, xb, yb, zb) > distance(xb, yb, zb, x, y, z))  ) {
					ta = t; xa = xt; ya = yt; za = zt;
				} else {
					tb = t; xb = xt; yb = yt; zb = zt;
				}
			}
		}
		pU[i] = t;
		V_U(iV_U, t, &tt);
		pV[i] = tt;
	}
	return;
} /* adjustingLine */


void makeAFLine (int *vert, int nVert, int bInverse) {
	int i;

	for (i=0; i<nVert-1; i++) {
		if (vert[i]==vert[i+1]) continue; // cheating ;-)
		if (bInverse!=0) {
			addFace32(vert[i+1], vert[i]);
//			printf("edge %d %d\n", vert[i+1], vert[i]);
		}
		if (bInverse!=1) {
			addFace32(vert[i], vert[i+1]);
//			printf("edge %d %d\n", vert[i], vert[i+1]);
		}
	}
	return;
} /*makeAFLine*/

static int dump(void) {
	static int num=0;
	char fname[1024];
	FILE *f;
	int i;
	sprintf(fname, "_dump_region_smv.%03d", num);
	f = fopen(fname, "w");

	fprintf(f, "%d %d 0 0\n", mesh3.nPoint, reg3.nTria);
	for (i=0; i<mesh3.nPoint; i++) fprintf(f, "%20.15lf %20.15lf %20.15lf\n", mesh3.vert[i].x, mesh3.vert[i].y, mesh3.vert[i].z);
	for (i=0; i<reg3.nTria; i++) fprintf(f, "%d %d %d\n", reg3.v1[i]+1, reg3.v2[i]+1, reg3.v3[i]+1);
	fclose(f);
	return num++;
}

int writeBound (int n, int label) {
	int  i, j;
	double xv, yv, zv, xw, yw, zw, t, tt;

	for (i=0, j=n; i<reg3.nTria; i++, j++) {
		if (j > mesh3.maxTetra)
			errorExit3(3, "j > mesh3.maxTetra ");
		mesh3.tetra[j].v1 = reg3.v1[i];
		mesh3.tetra[j].v2 = reg3.v2[i];
		mesh3.tetra[j].v3 = reg3.v3[i];
		mesh3.tetra[j].v4 = label;
	}

	if (region_dump_face) dump();

	if (reg3.iSurf != 0 && 0) { // FIXME check this logic
		for (i=0, j=n; i<reg3.nTria; i++, j++) {
			CrvT.u1[countCrvT] = mesh3.vert[mesh3.tetra[j].v1].u;
			CrvT.u2[countCrvT] = mesh3.vert[mesh3.tetra[j].v2].u;
			CrvT.u3[countCrvT] = mesh3.vert[mesh3.tetra[j].v3].u;
			CrvT.v1[countCrvT] = mesh3.vert[mesh3.tetra[j].v1].v;
			CrvT.v2[countCrvT] = mesh3.vert[mesh3.tetra[j].v2].v;
			CrvT.v3[countCrvT] = mesh3.vert[mesh3.tetra[j].v3].v;
			CrvT.bnd[countCrvT] = j;
			CrvT.iSurf[countCrvT] = reg3.iSurf;
			countCrvT = countCrvT + 1;
			xv = mesh3.vert[mesh3.tetra[j].v1].x;
			yv = mesh3.vert[mesh3.tetra[j].v1].y;
			zv = mesh3.vert[mesh3.tetra[j].v1].z;
			t  = mesh3.vert[mesh3.tetra[j].v1].u;
			tt = mesh3.vert[mesh3.tetra[j].v1].v;
			bounSurf(reg3.iSurf, t, tt, &xw, &yw, &zw);
			t = distance(xv, yv, zv, xw, yw, zw);
			if (t > precision) {
				printf("dist= %12.11lf \n", t);
//				errorExit3(2, "Wrong parametrization!");
			}
			xv = mesh3.vert[mesh3.tetra[j].v2].x;
			yv = mesh3.vert[mesh3.tetra[j].v2].y;
			zv = mesh3.vert[mesh3.tetra[j].v2].z;
			t  = mesh3.vert[mesh3.tetra[j].v2].u;
			tt = mesh3.vert[mesh3.tetra[j].v2].v;
			bounSurf(reg3.iSurf, t, tt, &xw, &yw, &zw);
			t = distance(xv, yv, zv, xw, yw, zw);
			if (t > precision) {
				printf("dist= %12.11lf \n", t);
//				errorExit3(2, "Wrong parametrization!");
			}
			xv = mesh3.vert[mesh3.tetra[j].v3].x;
			yv = mesh3.vert[mesh3.tetra[j].v3].y;
			zv = mesh3.vert[mesh3.tetra[j].v3].z;
			t  = mesh3.vert[mesh3.tetra[j].v3].u;
			tt = mesh3.vert[mesh3.tetra[j].v3].v;
			bounSurf(reg3.iSurf, t, tt, &xw, &yw, &zw);
			t = distance(xv, yv, zv, xw, yw, zw);
			if (t > precision) {
				printf("dist= %12.11lf \n", t);
//				errorExit3(2, "Wrong parametrization!");
			}
		}
	}
	return reg3.nTria;
} /*writeBound*/


void initAFS_ (
		int *pnVVert, double *VVertxyz,
		int *pnLine, int *LineD, int *LineP, double *LineT,
		int *pnSurface, int *SurfL, int *SurfI, double *SurfT
		) {
	int i, j, k, iLine, iSub, m;

	int nVVert = *pnVVert;
	int nLine = *pnLine;
	int nSurface = *pnSurface;

	int vBegin, vEnd, iSurf, iNorm, iV_U;
	double uBegin, uEnd;
	double uMin, uMax, vMin, vMax;
	double pmin[3]={0,0,0}, pmax[3]={0,0,0};

	int bInverse, bCut, iLabel, nBoundTria=0;

	int nSub, iSurfLine;
	int nVert;
	double x, y, z, x1, y1, z1, u, v;
	int *boundVert, *vVert;
	double *boundU, *boundV;
	StrucLine3 *line;
	StrucSurface *surface;

	init();
	detinit();

	tree3.boxcx = 0.5;
	tree3.boxcy = 0.5;
	tree3.boxcz = 0.5;
	tree3.boxsize = 0.5;

	boundVert = (int *)myAlloc(MAX1 * sizeof(int));
	boundU = (double *)myAlloc(MAX1 * sizeof(double));
	boundV = (double *)myAlloc(MAX1 * sizeof(double));

	vVert = (int *)myAlloc(nVVert * sizeof(int));

	for (i=0; i<nVVert; i++) {
		x = VVertxyz[3*i+0];
		y = VVertxyz[3*i+1];
		z = VVertxyz[3*i+2];
		if (i==0) {
			pmin[0] = pmax[0] = x;
			pmin[1] = pmax[1] = y;
			pmin[2] = pmax[2] = z;
		} else {
			if (pmin[0]>x) pmin[0]=x;
			if (pmax[0]<x) pmax[0]=x;
			if (pmin[1]>y) pmin[1]=y;
			if (pmax[1]<y) pmax[1]=y;
			if (pmin[2]>z) pmin[2]=z;
			if (pmax[2]<z) pmax[2]=z;
		}
		vVert[i] = mesh3.nPoint;
		addPoint(x, y, z);
	}
	
//	printf("BOX:%lf %lf %lf -- %lf %lf %lf\n", pmin[0], pmin[1], pmin[2], pmax[0], pmax[1], pmax[2]);
	S0 = (pmax[0]-pmin[0]);
	if (pmax[1]-pmin[1] > S0) S0 = pmax[1]-pmin[1];
	if (pmax[2]-pmin[2] > S0) S0 = pmax[2]-pmin[2];
	S0 *= surfacesizeratio;
//	printf("S0: %lf\n", S0);
	
	tree3.boxcx = 0.5*(pmin[0]+pmax[0]);
	tree3.boxcy = 0.5*(pmin[1]+pmax[1]);
	tree3.boxcz = 0.5*(pmin[2]+pmax[2]);
	tree3.boxsize = pmax[0]-pmin[0];
	if (pmax[1]-pmin[1]>tree3.boxsize) tree3.boxsize = pmax[1]-pmin[1];
	if (pmax[2]-pmin[2]>tree3.boxsize) tree3.boxsize = pmax[2]-pmin[2];


	line = (StrucLine3 *)myAlloc(2*nLine * sizeof(StrucLine3) );

	for (iSub=0, iLine=0; iLine<nLine; iLine++) {
		vBegin = LineD[3*iLine+0];
		vEnd   = LineD[3*iLine+1];
		nSub   = LineD[3*iLine+2];
		line[iLine].vBegin = vBegin - 1;
		line[iLine].vEnd   = vEnd - 1;
		line[iLine].nSub   = nSub;
		line[iLine].sub = (StrucSub *)myAlloc(nSub * sizeof(StrucSub));
		for (i=0; i<nSub; i++, iSub++) {
			iSurf = LineP[2*iSub+0];
			iV_U  = LineP[2*iSub+1];
			uBegin = LineT[2*iSub+0];
			uEnd   = LineT[2*iSub+1];
			line[iLine].sub[i].iSurf = iSurf;
			line[iLine].sub[i].iV_U = iV_U;
			line[iLine].sub[i].uBegin = uBegin;
			line[iLine].sub[i].uEnd = uEnd;
		}
	}

	/* for  mooving  vVert ... */
	for (iLine=0; iLine<nLine; iLine++) {
		iSurf = line[iLine].sub[0].iSurf;
		if (iSurf == 0)
			continue;
		vBegin = line[iLine].vBegin;
		vEnd = line[iLine].vEnd;
		uBegin = line[iLine].sub[0].uBegin;
		uEnd = line[iLine].sub[0].uEnd;
		iV_U = line[iLine].sub[0].iV_U;

		V_U(iV_U, uBegin, &v);
		bounSurf(iSurf, uBegin, v, &x, &y, &z);
		mesh3.vert[vBegin].x = x;
		mesh3.vert[vBegin].y = y;
		mesh3.vert[vBegin].z = z;
		mesh3.vert[vBegin].u = uBegin;
		mesh3.vert[vBegin].v = v;

		V_U(iV_U, uEnd, &v);
		bounSurf(iSurf, uEnd, v, &x, &y, &z);
		mesh3.vert[vEnd].x = x;
		mesh3.vert[vEnd].y = y;
		mesh3.vert[vEnd].z = z;
		mesh3.vert[vEnd].u = uEnd;
		mesh3.vert[vEnd].v = v;
	}

	surface = (StrucSurface *)myAlloc(nSurface * sizeof(StrucSurface));

	for (iSurfLine=0, iLine=0; iLine<nSurface; iLine++) {
		i      = SurfL[5*iLine+0];
		iSurf  = SurfL[5*iLine+1];
		iLabel = SurfL[5*iLine+2];
		bCut   = SurfL[5*iLine+3];
		iNorm  = SurfL[5*iLine+4];
		uMin = SurfT[4*iLine+0];
		uMax = SurfT[4*iLine+1];
		vMin = SurfT[4*iLine+2];
		vMax = SurfT[4*iLine+3];
		surface[iLine].nLine = i;
		surface[iLine].line = (int *)myAlloc(i * sizeof(int));
		surface[iLine].inverse = (int *)myAlloc(i * sizeof(int));

		for (j=0; j<i; j++, iSurfLine++) {
			k        = SurfI[2*iSurfLine+0];
			bInverse = SurfI[2*iSurfLine+1];
			surface[iLine].line[j] = k - 1;
			surface[iLine].inverse[j] = bInverse;
		}
		surface[iLine].iSurf = iSurf;
		surface[iLine].iLabel = iLabel;
		surface[iLine].bCut = bCut;
		surface[iLine].iNorm = iNorm;
		surface[iLine].uMin = uMin;
		surface[iLine].uMax = uMax;
		surface[iLine].vMin = vMin;
		surface[iLine].vMax = vMax;
	}


	/* lines into  sides */
	for (iLine=0; iLine<nLine; iLine++) {
		vBegin = line[iLine].vBegin;
		vEnd = line[iLine].vEnd;
		iSurf = line[iLine].sub[0].iSurf;
		if (iSurf > 0) {
			reg3.uBegin = uBegin = line[iLine].sub[0].uBegin;
			reg3.uEnd = uEnd = line[iLine].sub[0].uEnd;
			iV_U = line[iLine].sub[0].iV_U;
		} else {/*iSurf == 0*/
			reg3.uBegin = uBegin = 0.;
			reg3.uEnd = uEnd = 1.;
			iV_U = 0;
			reg3.x0 = mesh3.vert[vBegin].x;
			reg3.y0 = mesh3.vert[vBegin].y;
			reg3.z0 = mesh3.vert[vBegin].z;
			reg3.x1 = mesh3.vert[vEnd].x - reg3.x0;
			reg3.y1 = mesh3.vert[vEnd].y - reg3.y0;
			reg3.z1 = mesh3.vert[vEnd].z - reg3.z0;
			reg3.x2 = 0.;
			reg3.y2 = 0.;
			reg3.z2 = 0.;
		}

		nVert = 0;
		boundU[nVert] = uBegin;
		V_U(iV_U,uBegin,&v);
		boundV[nVert] = v;
		boundVert[nVert++] = vVert[vBegin];
		u = uBegin;
		for (j=0; ; j++) {
			if (nVert >= MAX1)
				errorExit3(4, "MAX1");
			u = nextU(iSurf, iV_U, u);
			if (u < -1.e10)
				break;

			V_U(iV_U, u, &v);
			mesh3.vert[mesh3.nPoint].u = u;
			mesh3.vert[mesh3.nPoint].v = v;
			bounSurf(iSurf, u, v, &x, &y, &z);
			boundU[nVert] = u;
			boundV[nVert] = v;
			boundVert[nVert++] = mesh3.nPoint;
			addPoint(x, y, z);
		}/* for(j) */
		boundU[nVert] = uEnd;
		V_U(iV_U, uEnd, &v);
		boundV[nVert] = v;
		boundVert[nVert++] = vVert[vEnd];
		smoothingLine(boundVert, boundU, boundV, iSurf, iV_U, nVert); 

		line[iLine].nVert = nVert;
		line[iLine].vert = (int *)myAlloc(nVert * sizeof(int));
		line[iLine].sub[0].u = (double *)myAlloc(nVert * sizeof(double));
		line[iLine].sub[0].v = (double *)myAlloc(nVert * sizeof(double));
		for (j=0; j<nVert; j++) {
			line[iLine].vert[j] = boundVert[j];
			line[iLine].sub[0].u[j] = boundU[j];
			line[iLine].sub[0].v[j] = boundV[j];
		}

		for (iSub=1; iSub<line[iLine].nSub; iSub++) {
			vBegin = line[iLine].vBegin;
			vEnd = line[iLine].vEnd;
			reg3.uBegin = uBegin = line[iLine].sub[iSub].uBegin;
			reg3.uEnd = uEnd = line[iLine].sub[iSub].uEnd;
			iV_U = line[iLine].sub[iSub].iV_U;
			iSurf = line[iLine].sub[iSub].iSurf;


			boundU[0] = uBegin;
			boundU[nVert-1] = uEnd;

			adjustingLine(boundVert, boundU, boundV, iSurf, iV_U, nVert);
			line[iLine].sub[iSub].u = (double *)myAlloc(nVert * sizeof(double));
			line[iLine].sub[iSub].v = (double *)myAlloc(nVert * sizeof(double));
			for (j=0; j<nVert; j++) {
				line[iLine].sub[iSub].u[j] = boundU[j];
				line[iLine].sub[iSub].v[j] = boundV[j];
			}
		}/*for iSub < nSub*/
	}/* for iLine */
	mesh3.nLinePoint = mesh3.nPoint;

	reg3.maxTria = tree3.maxFace;
	reg3.v1 = (int *)myAlloc(reg3.maxTria * sizeof(int));
	reg3.v2 = (int *)myAlloc(reg3.maxTria * sizeof(int));
	reg3.v3 = (int *)myAlloc(reg3.maxTria * sizeof(int));

	countCrvT = 0;
	CrvT.bnd = (int *)myAlloc(reg3.maxTria * sizeof(int));
	CrvT.iSurf = (int *)myAlloc(reg3.maxTria * sizeof(int));
	CrvT.u1 = (double *)myAlloc(reg3.maxTria * sizeof(double));
	CrvT.u2 = (double *)myAlloc(reg3.maxTria * sizeof(double));
	CrvT.u3 = (double *)myAlloc(reg3.maxTria * sizeof(double));
	CrvT.v1 = (double *)myAlloc(reg3.maxTria * sizeof(double));
	CrvT.v2 = (double *)myAlloc(reg3.maxTria * sizeof(double));
	CrvT.v3 = (double *)myAlloc(reg3.maxTria * sizeof(double));


	for (iSurf=0; iSurf<nSurface; iSurf++) {
		tree32.nFace = 0;
		tree32.root->flag = EMPTY;
		nLine = surface[iSurf].nLine;
		iLabel = surface[iSurf].iLabel;
		bCut = surface[iSurf].bCut;

		reg3.nTria = 0;
		reg3.uMin = surface[iSurf].uMin;
		reg3.uMax = surface[iSurf].uMax;
		reg3.vMin = surface[iSurf].vMin;
		reg3.vMax = surface[iSurf].vMax;
		reg3.iSurf = surface[iSurf].iSurf;
		if (reg3.iSurf > 0)
			reg3.iNorm = surface[iSurf].iNorm;
		else
			reg3.iNorm = 0;
		for (iLine=0; iLine<nLine; iLine++) {
			i = surface[iSurf].line[iLine];
			j = line[i].nSub;
			if ( (reg3.iSurf>0) && (j>=1) ) {/*reload u&v*/
				iSub = -1;
				for (k=0; k<j; k++) {
					if (line[i].sub[k].iSurf == reg3.iSurf) {
						iSub = k;
						break;
					}
				}
				if (iSub == -1)
					errorExit3(3, " iSub == 0");
				for (k=0; k<line[i].nVert; k++) {
					mesh3.vert[ line[i].vert[k] ].u = line[i].sub[iSub].u[k];
					mesh3.vert[ line[i].vert[k] ].v = line[i].sub[iSub].v[k];
				}
			}
			bInverse = surface[iSurf].inverse[iLine];
			makeAFLine(line[i].vert, line[i].nVert, bInverse);
		}/*for(iLine<nLine)*/

		if (reg3.iSurf == 0) {/*init  for  plane*/
			reg3.uMin = 0.;
			reg3.uMax = 1.;
			reg3.vMin = 0.;
			reg3.vMax = 1.;
			i = surface[iSurf].line[0];
			x = 0.0;
			y = 0.0;
			z = 0.0;
			for (k=0; k<tree32.nFace; k++) {
			    reg3.vv1 = tree32.face[k]->v1;
			    reg3.vv2 = tree32.face[k]->v2;
			    normvec(mesh3.vert[reg3.vv1].x, mesh3.vert[reg3.vv1].y, mesh3.vert[reg3.vv1].z,  mesh3.vert[reg3.vv2].x, mesh3.vert[reg3.vv2].y, mesh3.vert[reg3.vv2].z, &x1, &y1, &z1);
			    x+=x1, y+=y1, z+=z1;
			}

/*			for (iLine=0; iLine<nLine; iLine++) {
				i = surface[iSurf].line[iLine];
				x1 = 0.0,  y1 = 0.0,  z1 = 0.0;
				for (k=0; k<line[i].vert[line[i].nVert]; k++) {
				}
				reg3.vv1 = line[i].vert[0];
				reg3.vv2 = line[i].vert[ line[i].nVert - 1 ];
				normvec(mesh3.vert[reg3.vv1].x, mesh3.vert[reg3.vv1].y, mesh3.vert[reg3.vv1].z,  mesh3.vert[reg3.vv2].x, mesh3.vert[reg3.vv2].y, mesh3.vert[reg3.vv2].z, &x1, &y1, &z1);
				if (surface[iSurf].inverse[iLine]) {
					x-=x1, y-=y1, z-=z1;
				} else {
					x+=x1, y+=y1, z+=z1;
				}
			}*/
			u = sqrt(x*x + y*y + z*z);
			if (u == 0.0)
				errorExit3(3, " u == 0.0   in  surfNormal ");
			x /= u;
			y /= u;
			z /= u;
			reg3.vv0 = line[i].vert[0];
			reg3.vv1 = line[i].vert[ line[i].nVert - 1 ];
			if (surface[iSurf].inverse[0] > 0) {/*for  right  normal*/
				k = reg3.vv1;
				reg3.vv1 = reg3.vv0;
				reg3.vv0 = k;
			}
			x1 = mesh3.vert[reg3.vv1].x - mesh3.vert[reg3.vv0].x;
			y1 = mesh3.vert[reg3.vv1].y - mesh3.vert[reg3.vv0].y;
			z1 = mesh3.vert[reg3.vv1].z - mesh3.vert[reg3.vv0].z;
			u = sqrt(x1*x1 + y1*y1 + z1*z1);
			if (u == 0.0)
				errorExit3(3, " u == 0.0   in  surfNormal ");
			x1 /= u;
			y1 /= u;
			z1 /= u;

			reg3.x0 = mesh3.vert[reg3.vv0].x;
			reg3.y0 = mesh3.vert[reg3.vv0].y;
			reg3.z0 = mesh3.vert[reg3.vv0].z;
			reg3.x1 = x1;
			reg3.y1 = y1;
			reg3.z1 = z1;
			normvec(x, y, z,  x1, y1, z1,  &reg3.x2, &reg3.y2, &reg3.z2);

			for (iLine=0; iLine<nLine; iLine++) {/*reload u&v*/
				i = surface[iSurf].line[iLine];
				for (j=0; j<line[i].nVert; j++) {
					m = line[i].vert[j];
					rePlane(m);
					if (mesh3.vert[m].u > reg3.uMax) reg3.uMax = mesh3.vert[m].u;
					if (mesh3.vert[m].u < reg3.uMin) reg3.uMin = mesh3.vert[m].u;
					if (mesh3.vert[m].v > reg3.vMax) reg3.vMax = mesh3.vert[m].v;
					if (mesh3.vert[m].v < reg3.vMin) reg3.vMin = mesh3.vert[m].v;
				}
			}
		}/*if( reg3.iSurf == 0 )*/

		makeTria();
		for (i=0; i<5; i++) {
			smoothingSurf(); 
		}
		nBoundTria += writeBound(nBoundTria, iLabel);

		/* make AF for surface  */
		for (i=0; i<reg3.nTria; i++) {
			addFace(reg3.v1[i], reg3.v2[i], reg3.v3[i], 0, iLabel);
			if (bCut)
				addFace(reg3.v1[i], reg3.v3[i], reg3.v2[i], 1, bCut);
		}
	}/*for(iSurf<nSurface)*/

	mesh3.nBoundPoint = mesh3.nPoint;
//	outBound(nBoundTria);
	free(boundVert);
	free(boundU);
	free(boundV);
	free(vVert);
	for (iLine=0; iLine<nLine; iLine++) {
		for (i=0; i<line[iLine].nSub; i++) {
			free(line[iLine].sub[i].u);
			free(line[iLine].sub[i].v);
		}
		free(line[iLine].vert);
		free(line[iLine].sub);
	}
	free(line);
	for (iLine=0; iLine<nSurface; iLine++) {
		free(surface[iLine].line);
		free(surface[iLine].inverse);
	}
	free(surface);
	free(reg3.v1);
	free(reg3.v2);
	free(reg3.v3);
	free(CrvT.bnd);
	free(CrvT.iSurf);
	free(CrvT.u1);
	free(CrvT.u2);
	free(CrvT.u3);
	free(CrvT.v1);
	free(CrvT.v2);
	free(CrvT.v3);
	return;
} /* makeAFSurface */


void addTria (int v1, int v2, int v3) {
	if (reg3.nTria >= reg3.maxTria)
		errorExit3(2, "reg3.nTria >= reg3.maxTria");
	reg3.v1[reg3.nTria] = v1;
	reg3.v2[reg3.nTria] = v2;
	reg3.v3[reg3.nTria] = v3;
	reg3.nTria++;

	return;
} /*addTria*/

