#ifndef H_SECTION3_MESH3D
#define H_SECTION3_MESH3D


typedef struct{
	int     vert;    /*   */
	int     next;    /*   */
	int     nea;
	double  sin;     /*   */
	double  len;     /*   */
}  StrucConList3;


typedef struct{
	int          bool;            /*   */
	StrucVert3   v[6];            /*   */
	StrucVert3   vv[5];           /*   */
	int          workVert;        /*   */
	int          neigAdd[4];      /*   */
	int          neigBool[3];     /*   */
	int          neigInt[3];      /*   */
	int          neigVert[3];     /*   */
	double       neigSin[4];      /*   */
	PStrucFace3  neigFace[3];     /*   */
	PStrucFace3  work;            /*   */
	PStrucFace3  sect;            /*   */
	StrucConList3  list[55];     /*   */
	double       xSame,ySame;
	double       x,y,z,x0,y0,z0;
	double       xMin,xMax,yMin,yMax,zMin,zMax;
	int          nVert;
}  StrucSect3;


/* exported  function */
void  shift( PStrucVert3  v, double  dx, double  dy, double  dz );
void  rotate( PStrucVert3  v, unsigned  char  exix, double  cos, double  sin );
void  makeVector( PStrucVert3 v, int v1, int v2 );
int  makeNormal(PStrucVert3 v, PStrucVert3 v1, PStrucVert3 v2 );
int  intersectOne( void );


#endif

