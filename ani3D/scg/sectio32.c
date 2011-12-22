#include "struct3.h"
#include "error3.h"
#include "section3.h"
#include "sectio32.h"
#include "tree3.h"


/* extern  variables */
extern StrucMesh3 mesh3;
extern StrucSect3 sect3;
extern int sameVert;/* same  vert */

/* static  function */
static int isSamePoint(PStrucFace3 face, int v1, int v2, int v3);


static int isWorkFace(PStrucFace3 face) {
	if (face == sect3.work)
		return 1;
	if ( (sect3.work->v1==face->v2) && (sect3.work->v2==face->v1) )
		return 1;
	return 0;
} /*isWorkFace*/


static int isCommonEdge (PStrucFace3 face, int v1, int v2, int v3) {
	if ( (v1==face->v1) && (v3==face->v2) )
		return 1;
	if ( (v1==face->v2) && (v3==face->v1) ) {
		sect3.neigBool[0] = 1;
		sect3.neigFace[0] = face;
		return 1;
	}

	if ( (v2==face->v1) && (v3==face->v2) ) {
		sect3.neigBool[1] = 1;
		sect3.neigFace[1] = face;
		return 1;
	}
	if ( (v2==face->v2) && (v3==face->v1) )
		return 1;
	return 0;
} /*isCommonEdge*/


static int isSamePoint (PStrucFace3 face, int v1, int v2, int v3) {
	if ( (v1==face->v1) || (v1==face->v2) )
		return v1;
	if ( (v2==face->v1) || (v2==face->v2) )
		return v2;
	if ( (v3==face->v1) || (v3==face->v2) )
		return v3;
	return -1;
} /*isSamePoint*/


int faceIntersect32( PStrucFace3  face, int  v1, int  v2, int  v3 ) {
	/* proection on perpendicular to surfNormal */
	int i, bool=0;
	int nSameVert=-1; /* num.  of  same  vert  in  v[0-5]  in  sect3 */
	double p, dx, dy, dz, sin=0., cos=1.;

	sameVert = -1;
	/*********************/
	if (isWorkFace(face))
		return 0;
	/*********************/
	if (isCommonEdge(face, v1, v2, v3))
		return 0;
	/*********************/
	if ( (sameVert=isSamePoint(face, v3, v1, v2)) > -1 ) {
		if (sameVert == v1)
			nSameVert = 3;
		if (sameVert == v2)
			nSameVert = 4;
		if (sameVert == v3)
			nSameVert = 5;
	}

	/*********************/
	sect3.v[1].x = mesh3.vert[face->v1].x;
	sect3.v[1].y = mesh3.vert[face->v1].y;
	sect3.v[1].z = mesh3.vert[face->v1].z;
	sect3.v[2].x = mesh3.vert[face->v2].x;
	sect3.v[2].y = mesh3.vert[face->v2].y;
	sect3.v[2].z = mesh3.vert[face->v2].z;

	sect3.v[3].x = mesh3.vert[v1].x;
	sect3.v[3].y = mesh3.vert[v1].y;
	sect3.v[3].z = mesh3.vert[v1].z;
	sect3.v[4].x = mesh3.vert[v2].x;
	sect3.v[4].y = mesh3.vert[v2].y;
	sect3.v[4].z = mesh3.vert[v2].z;
	sect3.v[5].x = mesh3.vert[v3].x;
	sect3.v[5].y = mesh3.vert[v3].y;
	sect3.v[5].z = mesh3.vert[v3].z;

	dx = sect3.v[5].x;
	dy = sect3.v[5].y;
	dz = sect3.v[5].z;
	for (i=1; i<6; i++)
		shift(&sect3.v[i], dx, dy, dz);

	p = sqrt(sect3.v[0].z*sect3.v[0].z+sect3.v[0].y*sect3.v[0].y);
	if (p != 0) {
		cos = sect3.v[0].z/p; sin = -sect3.v[0].y/p;
		for (i=0; i<5; i++)
			rotate(&sect3.v[i], 'x', cos, sin);
	} else
		bool = 1;

	p = sqrt(sect3.v[0].x*sect3.v[0].x+sect3.v[0].z*sect3.v[0].z);
	if (p != 0) {
		cos = sect3.v[0].z/p; sin = -sect3.v[0].x/p;
	} else
		if (bool)
			errorExit3(2, "p==0  in  faceIntersect");
	for (i=0; i<5; i++)
		rotate(&sect3.v[i], 'y', cos, sin);

	p = sqrt(sect3.v[4].x*sect3.v[4].x+sect3.v[4].y*sect3.v[4].y);
	if (p != 0) {
		cos = sect3.v[4].x/p; sin = sect3.v[4].y/p;
		for (i=1; i<5; i++)
			rotate(&sect3.v[i], 'z', cos, sin);
	} else
		errorExit3(2, "p==0  in  faceIntersect");

	/***************************/
	/* calc proection */
	sect3.v[1].z = 0.;
	sect3.v[2].z = 0.;
	sect3.v[3].z = 0.;
	sect3.v[4].z = 0.;

	/***************************/
	if ( sameVert > -1) {/* case  sameVert */
		sect3.xSame = sect3.v[nSameVert].x;
		sect3.ySame = sect3.v[nSameVert].y;
	}

	if (intersectOne())
		return  1;
	else
		return  0;
} /*faceIntersect32*/


