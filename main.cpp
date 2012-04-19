#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "AngularData.h"
#include "MeshData.h"
#include "kernels.h"

int main() {
	AngularData ad(1); /* maxk = 1, maxl = 2*/
	MeshData md("mesh.vol");

	DeviceAngularData dad(ad);
	DeviceMeshData dmd(md);

	REAL *f = (REAL *)deviceAlloc(dad.aslm * dmd.nP * sizeof(REAL));
	REAL *Af = (REAL *)deviceAlloc(dad.aslm * dmd.nP * sizeof(REAL));
	REAL *b = (REAL *)deviceAlloc(dad.aslm * dmd.nP * sizeof(REAL));

	REAL *_f = new REAL[dad.aslm * dmd.nP];
	REAL *_Af = new REAL[dad.aslm * dmd.nP];
	REAL *_b = new REAL[dad.aslm * dmd.nP];

	for (int i = 0; i < dad.aslm * dmd.nP; i++)
		_f[i] = 1;

	copyToDev(f, _f, dad.aslm * dmd.nP * sizeof(REAL));
	computeRhs(&dmd, &dad, f, Af, b);
	copyToHost(_b, b, dad.aslm * dmd.nP * sizeof(REAL));
	copyToHost(_Af, Af, dad.aslm * dmd.nP * sizeof(REAL));

	return 0;
}
