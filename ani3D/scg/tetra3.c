#include "tetra3.h"
#include "error3.h"

/* global  variables */
double minNear=0.85;
int nBadVert, badVert[10000], nSphereVert, sphereVert[10000];

void  addBadVert(int v) {
	if (nBadVert >= 10000)
		errorExit3(2, "nBadVert >= 10000");
	badVert[nBadVert++] = v;
	return;
} /*addBadVert*/

int isBadVert(int v) {
	int  i;
	for (i=0; i<nBadVert; i++)
		if (badVert[i] == v) return 1;
	return 0;
} /*isBadVert*/

void addSphereVert(int v) {
	if( nSphereVert >= 10000 )
		errorExit3(2, "nSphereVert >= 10000");
	sphereVert[nSphereVert++] = v;
	return;
} /*addSphereVert*/

int isSphereVert(int v) {
	int  i;
	for (i=0; i<nSphereVert; i++)
		if (sphereVert[i] == v) return  1;
	return  0;
} /*isSphereVert*/

