#ifndef H_REGION3_MESH3D
#define H_REGION3_MESH3D

#ifndef CGMWRAP_H
#define CGMface void*
#endif

typedef struct{
	double  x0,y0,z0;
	double  x1,y1,z1;
	double  x2,y2,z2;
	double  uBegin,uEnd;
	double  uMin,uMax,vMin,vMax,uSave,vSave;
	double  suMin,suMax,svMin,svMax;
	int     *v1,*v2,*v3,nTria,maxTria;
	int     vv0,vv1,vv2;
	int     iSurf,iNorm;
	CGMface face;
}  StrucRegion3;


typedef struct{
	int     iSurf,iV_U;
	double  uBegin,uEnd,*u,*v;
}  StrucSub;


typedef struct{
	int      nVert,*vert;
	int      vBegin,vEnd;
	int      nSub;
	StrucSub *sub;
}  StrucLine3;


typedef struct{
	int     nLine,*line,*inverse;
	int     iSurf,iNorm,iLabel,bCut;
	double  uMin,uMax,vMin,vMax;
}  StrucSurface;


typedef struct{
	int     *bnd,*iSurf;
	double  *u1,*u2,*u3,*v1,*v2,*v3;
}  StrucCrvT;



#define  MAX1   50000


/* exported  functions */
int  bounSurf0( double u, double v, double *x, double *y, double *z );
void V_U0( double u, double *v );
double sizeFace( double x, double y, double z );
void addTria( int v1, int v2, int v3 );
void initAFT_(int *nF, int *faces, int *nV, double *vertices);
void initAFS_ (int *pnVVert, double *VVertxyz, int *pnLine, int *LineD, int *LineP, double *LineT, int *pnSurface, int *SurfL, int *SurfI, double *SurfT);

#endif

