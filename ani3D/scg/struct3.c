#include"error3.h"
#include"memory3.h"
#include"tree3.h"
#include"user3.h"


/* extern  variables */
extern StrucTree3 tree3;
extern int boolRegul;
extern int nVertBS;


/* global  variables */
StrucMesh3 mesh3;



void init (void) {
	mesh3.nPoint = 0;
	mesh3.nTetra = 0;
	mesh3.surf = 0;
	tree3.nFace = 0;
	initMemory();
	tree3.root->flag = EMPTY;
	return;
} /*init*/


void addPoint (double x, double y, double z) {
	if (mesh3.nPoint >= mesh3.maxPoint)
		errorExit3(2, "mesh3.nPoint >= mesh3.maxPoint");
	mesh3.vert[mesh3.nPoint].x = x;
	mesh3.vert[mesh3.nPoint].y = y;
	mesh3.vert[mesh3.nPoint].z = z;
	mesh3.nPoint++;
	return;
} /*addPoint*/


void addTetra (int v1, int v2, int v3, int v4) {
	if (mesh3.nTetra >= mesh3.maxTetra)
		errorExit3(2, "mesh3.nTetra >= mesh3.maxTetra");
	mesh3.tetra[mesh3.nTetra].v1 = v1;
	mesh3.tetra[mesh3.nTetra].v2 = v2;
	mesh3.tetra[mesh3.nTetra].v3 = v3;
	mesh3.tetra[mesh3.nTetra].v4 = v4;
	mesh3.nTetra++;
	return;
} /*addTetra*/


