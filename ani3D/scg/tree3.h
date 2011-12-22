#ifndef H_TREE3_MESH3D
#define H_TREE3_MESH3D

#include "struct3.h"

typedef struct {
	PStrucFace3	face;	/**/
	double		dist;	/**/
} StrucFace4;

typedef struct {
	PStrucNode3	root;	/*  root  of  the  quadtree  */
	PStrucFace3	*face;	/*  array of  the  faces  in  advanced  front */
	int			nFace;
	long		maxFace;	/*  number  &  max  ...  of  faces  in  ... */

	StrucFace4	*vicinityFace;	/*  array of  faces  the  vicinity  */
	int			nVicinityFace,maxVicinityFace;	/*  number  &  max  ...  of  faces  in  ... */
	double		xVicinity,yVicinity,zVicinity;	/*  center  of  vicinity */
	double		sVicinity,ssVicinity;	/*  size  of  vicinity */

	int			fill,empty;	/*  global  for  recursive  remove  function  */
	double		xc,yc,zc,side;	/*  global  for  recursive  remove  &  insert  function  */
	double		x,y,z;	/*  global  for  recursive  remove  &  insert  function  */
	double      boxcx, boxcy, boxcz;
	double      boxsize;
} StrucTree3;

#define	EMPTY	0
#define	FACE	1
#define	TREE	2

/* exported functions */
int direction(double x,double y,double z,double xc,double yc,double zc);
void center(double *x,double *y,double *z,double side,int d);
double distance(double x,double y,double z,double xc,double yc,double zc);
double distanceS(double x,double y,double z,double xc,double yc,double zc);
double distance9(double x,double y,double z,double xc,double yc,double zc);
PStrucFace3 addFace(int v1, int v2, int v3, int twin, int color );
void remFace( PStrucFace3  face );
double nearest3( int *vert, double x, double y, double z );
int compare( const void *p1, const void *p2 );
void vicinityFaces( double x, double y, double z, double size );
int findMaxTetra( void );
PStrucFace3 findFace (int v1, int v2, int v3);

#endif


