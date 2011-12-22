#ifndef H_STRUCT3_MESH3D
#define H_STRUCT3_MESH3D


#include<ctype.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>


#define TEST_
#define WARN
#define DEBUG_
#define BADSIT


#ifdef TEST
#include<time.h>
#endif


#define  M_PI           3.14159265358979323846
#define  M_PI_2         1.57079632679489661923
#define  M_1_3          0.33333333333333333333
#define  M_ERROR        1.e-13
#define  S_MAX           1.1   
#define  S_MIN          -0.1   

//#define  MAX_NEIGBOR  55    //by Danilov for tests
#define  MAX_NEIGBOR  33


typedef struct{
	double    x,y,z;           /*   */
	double    u,v;             /*   */
}  StrucVert3;
typedef  StrucVert3  *PStrucVert3;


typedef struct{
	int       v1,v2,v3,v4;           /*   */
}  StrucTetra3;
typedef  StrucTetra3  *PStrucTetra3;


typedef struct{
	int       n;
	int       neig[MAX_NEIGBOR];    /*   */
}  StrucNeigbor;
typedef  StrucNeigbor  *PStrucNeigbor;


typedef struct{
	int      v1;
	int      v2;
	int      v3;
	int      t1;
	int      t2;
}  StrucEdge3;
typedef  StrucEdge3  *PStrucEdge3;


typedef struct{
	int      t1;
	int      t2;
	int      t3;
	int      t4;
}  StrucTetTet3;
typedef  StrucTetTet3  *PStrucTetTet3;


typedef struct{
	int     v1,v2,v3;   /**/
	int     f;          /* number  of  face       */
	double  x,y,z,s;    /**/
	int     color;
	int fail;
}  StrucFace3;
typedef  StrucFace3  *PStrucFace3;


typedef struct{
	int         flag;
	PStrucFace3 faceNode;
	double      maxs;
}  StrucNode3;
typedef  StrucNode3  *PStrucNode3;


typedef  PStrucNode3  StrucList3[8];
typedef  StrucList3   *PStrucList3;


typedef struct{
	char    surf;

	int     maxPoint,maxTetra;             /**/
	int     nPoint,nTetra;    
	PStrucVert3   vert;
	PStrucTetra3  tetra;

	int     nLinePoint;
	int     nBoundPoint;
	PStrucNeigbor  neigbor; /**/
	PStrucNeigbor  neigTetra; /**/

	char    *bPoint;
	char    *bTetra;
	int     nEdge;
	PStrucEdge3    edge;
	PStrucTetTet3  tetTet;

	double  size;
	double  sminsize;
}  StrucMesh3;


#define  S_StrucVert3     sizeof(StrucVert3)
#define  S_StrucTetra3    sizeof(StrucTetra3)
#define  S_StrucNeigbor   sizeof(StrucNeigbor)
#define  S_StrucEdge3     sizeof(StrucEdge3)
#define  S_StrucTetTet3   sizeof(StrucTetTet3)
#define  S_StrucList3     sizeof(StrucList3)
#define  S_StrucNode3     sizeof(StrucNode3)
#define  S_StrucFace3     sizeof(StrucFace3)


/* exported  function  function */
void init( void );
void reInit( void );
void addPoint( double x, double y, double z );
void addTetra( int v1, int v2, int v3, int v4 );
void addNeigbor( int v1, int v2 );
void outMesh( void );


#endif

