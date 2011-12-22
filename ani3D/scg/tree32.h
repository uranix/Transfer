#ifndef H_TREE32_MESH3D
#define H_TREE32_MESH3D


#include"struct3.h"


/* exported  function */
PStrucFace3 addFace32(int v1, int v2 );
void remFace32( PStrucFace3  face );
double nearest32( int *vert, double x, double y, double z );
void vicinityFaces32( double x, double y, double z, double size );


#endif


