#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "AngularData.h"
#include "MeshData.h"
#include "kernels.h"

int main(int argc, char **argv) {
	CudaContext::setDevice(0);
	AngularData ad(0); /* maxk = 1, maxl = 2*/
	MeshData md("mesh.vol");

	DeviceAngularData dad(ad);
	DeviceMeshData dmd(md);

	CudaContext *ctx = new CudaContext (7, &dmd, &dad);

	REAL *f = (REAL *)ctx->deviceAlloc(dad.aslm * dmd.nP * sizeof(REAL));
	REAL *b = (REAL *)ctx->deviceAlloc(dad.aslm * dmd.nP * sizeof(REAL));
	REAL *p = (REAL *)ctx->deviceAlloc(dad.aslm * dmd.nP * sizeof(REAL));
	REAL *Ap = (REAL *)ctx->deviceAlloc(dad.aslm * dmd.nP * sizeof(REAL));
	REAL *z = (REAL *)ctx->deviceAlloc(dad.aslm * dmd.nP * sizeof(REAL));
	REAL *r = (REAL *)ctx->deviceAlloc(dad.aslm * dmd.nP * sizeof(REAL));

	REAL *_f = new REAL[dad.aslm * dmd.nP];
	REAL *_Af = new REAL[dad.aslm * dmd.nP];
	REAL *_b = new REAL[dad.aslm * dmd.nP];

	idx N = dad.aslm * dmd.nP;

	double *Z = new double[N*N];

	FILE *file = fopen("matrix.txt", "w");
	ctx->computeRhs(b);
	ctx->copyToHost(_b, b, dad.aslm * dmd.nP * sizeof(REAL));

	for (int k = 0; k < dad.aslm * dmd.nP; k++) {

		for (int i = 0; i < dad.aslm * dmd.nP; i++)
			_f[i] = i==k;
		ctx->copyToDev(f, _f, dad.aslm * dmd.nP * sizeof(REAL));
		ctx->computeLhs(f, Ap);
		ctx->copyToHost(_Af, Ap, dad.aslm * dmd.nP * sizeof(REAL));
		for (int i = 0; i < dad.aslm * dmd.nP; i++)
			Z[k*N+i] = _Af[i];
		int j = k % dad.aslm;
		if (j < dad.slm) {
			for (int i = 0; i < dmd.nP; i++) 
				for (j = 0; j < dad.slm; j++) 
			{
				fprintf(file, "% 2.10e ", _Af[i*dad.aslm+j]);
			}
			fprintf(file, "% 2.10e ", _b[k]);
			fprintf(file, "\n");
		}
	}
	fclose(file);
	printf("Matrix dumped\n");

	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++) {
			if (fabs(Z[N*i+j]-Z[N*j+i]) > 1e-10)
				printf("Z[%d,%d] < %2.2e > Z[%d,%d]\n", i, j, fabs(Z[N*i+j]-Z[N*j+i]), j, i);
		}
	
	/*--------*/

	ctx->computeRhs(b);
	ctx->mulAdd(r, 0, b);
	ctx->computeLhs(f, Ap);
	ctx->addProd(r, Ap, -1);
	ctx->mulAdd(z, 0, r); /* z = invprec(r); */
	ctx->mulAdd(p, 0, z);
	int k = 0;
	double nrz = ctx->dot(r, z);
	while (1) {
		ctx->computeLhs(p, Ap);
		double alpha = nrz/ctx->dot(p, Ap);
		ctx->addProd(f, p, alpha);
		ctx->addProd(r, Ap, -alpha);
		double nr = ctx->norm(r);
		printf("norm r = %e\n", nr);
		if (nr < 1e-16)
			break;
		ctx->mulAdd(z, 0, r); /* z = invprec(r); */
		double nrz2 = ctx->dot(z, r);
		double beta = nrz2/nrz;
		nrz = nrz2;
		ctx->mulAdd(p, beta, z);
		k++;
	}

	/*--------*/
	
	ctx->copyToHost(_f, f, dad.aslm * dmd.nP * sizeof(REAL));

	return 0;
}
