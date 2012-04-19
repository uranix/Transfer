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

	REAL *f = deviceAlloc(dmd.aslm * dad.nP * sizeof(REAL));
	REAL *Af = deviceAlloc(dmd.aslm * dad.nP * sizeof(REAL));
	REAL *b = deviceAlloc(dmd.aslm * dad.nP * sizeof(REAL));

	REAL *_f = new REAL[dmd.aslm * dad.nP];
	REAL *_Af = new REAL[dmd.aslm * dad.nP];
	REAL *_b = new REAL[dmd.aslm * dad.nP];

	for (int i = 0; i < dmd.aslm * dad.nP; i++)
		_f[i] = 1;

	copyToDev(f, _f, dmd.aslm * dad.nP * sizeof(REAL));
	computeRhs(dmd, dad, f, Af, b);
	copyFromDev(_b, b, dmd.aslm * dad.nP * sizeof(REAL));
	copyFromDev(_Af, Af, dmd.aslm * dad.nP * sizeof(REAL));

	return 0;
}
