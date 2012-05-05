#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "meshProcessor/mesh.h"

#include "AngularData.h"
#include "MeshData.h"
#include "CudaContext.h"
#include "Config.h"

REAL uround() {
	REAL x = 1.0;
	while ((REAL)1.0 + x > 1.0)
		x = 0.5 * x;
	printf("uround = %2.10e\n", x);
	return x;
}

void pbar(const char *msg, double val) {
	static double last;
	int width = 100;
	double delta = 1.0 / (double)width;
	if (val > last + delta || val < 1e-12 || val > 1 - 1e-12) {
		int ticks = (double)width * val + 0.4999999;
		printf("\r%s: %6.2f%% [", msg, 100. * val);
		for (int i = 0; i < ticks; i++)
			printf("#");
		for (int i = ticks; i < width; i++)
			printf(" ");
		printf("]");
		last = val;
	}
	if (val > 1 - 1e-12)
		printf("\n");
}

int main(int argc, char **argv) {
	if (argc != 2) {
		fprintf(stderr, "USAGE: %s <config-file>\n", argv[0]);
		return 1;
	}
	Config cfg(argv[1]);
	REAL eps = uround();
	AngularData ad(cfg.getMaxK()); 
	MeshData md(cfg);

	CudaContext *ctx = new CudaContext (cfg.getDevice(), md, ad);


	REAL *f = ctx->allocVector();
	REAL *b = ctx->allocVector();
	REAL *p = ctx->allocVector();
	REAL *Ap= ctx->allocVector();
	REAL *z = ctx->allocVector();
	REAL *r = ctx->allocVector();

	REAL *_f = new REAL[ctx->angdata->aslm * ctx->meshdata->nP];

	idx N = ctx->N();

	if (cfg.doDump()) {
		REAL *_Af = new REAL[N];
		REAL * _b = new REAL[N];

		double *Z = new double[N*N];

		FILE *file = fopen("slae.dat", "w");
		ctx->computeRhs(b);
		ctx->copyToHost(_b, b, N * sizeof(REAL));

		double q = 1.0 / (N-1), qq = 0;

		printf("N = %d\n", N);
		pbar("Matrix dump", 0);
		for (idx k = 0; k < N; k++) {
			idx j = k % ctx->angdata->aslm;
			pbar("Matrix dump", qq += q);
			if (j >= ad.slm)
				continue;
			for (idx i = 0; i < N; i++)
				_f[i] = i==k;
			ctx->copyToDev(f, _f, N * sizeof(REAL));
			ctx->computeLhs(f, Ap);
			ctx->copyToHost(_Af, Ap, N * sizeof(REAL));
			for (idx i = 0; i < N; i++) 
				Z[k * N + i] = _Af[i];
			for (idx i = 0; i < md.nP; i++) 
				for (idx jj = 0; jj < ad.slm; jj++) 
					fprintf(file, "% 2.10e ", _Af[i*ctx->angdata->aslm+jj]);
			fprintf(file, "% 2.10e ", _b[k]);
			fprintf(file, "\n");
		}
		fclose(file);
		printf("System dumped\n");
	/* Symmetry check */
		for (idx i=0; i<N; i++)
			for (idx j=0; j<N; j++) {
				REAL diff =	fabs(Z[N*i+j]-Z[N*j+i]);
				REAL allow = eps * (fabs(Z[N*i+j]) + fabs(Z[N*j+i]));
				if ( diff > 100 * allow )
					printf("Z[%d,%d] = %2.10e < %2.2e (%2.2f x uround) > Z[%d,%d] = %2.10e \n", i, j, Z[i*N+j], diff, diff/allow, j, i, Z[j*N+i]);
			}
	}

	/* ----- cgs ----- */

	for (idx i = 0; i < ctx->N(); i++)
		_f[i] = 0;
	ctx->copyToDev(f, _f, ctx->N() * sizeof(REAL));
	ctx->computeRhs(b);
	ctx->computeLhs(f, Ap);
	ctx->mulAdd(r, 0, b);
	ctx->addProd(r, Ap, -1);
	ctx->mulAdd(z, 0, r); /* z = invprec(r); */
	ctx->mulAdd(p, 0, z);
	int k = 0;
	double nrz = ctx->dot(r, z);
	while (k < 1000) {
		ctx->computeLhs(p, Ap);
		double alpha = nrz/ctx->dot(p, Ap);
		ctx->addProd(f, p, alpha);
		ctx->addProd(r, Ap, -alpha);
		double nr = ctx->norm(r);
		printf("k = %d norm r = %e\n", k, nr);
		if (nr < 1e-10)
			break;
		ctx->mulAdd(z, 0, r); /* z = invprec(r); */
		double nrz2 = ctx->dot(z, r);
		double beta = nrz2/nrz;
		nrz = nrz2;
		ctx->mulAdd(p, beta, z);
		k++;
	}

	/*------ cgs end -----*/

	ctx->copyToHost(_f, f, ctx->N() * sizeof(REAL));
	REAL *u[ad.slm];
	REAL *Wx, *Wy, *Wz;
	for (idx k=0; k < ad.slm; k++) {
		u[k] = new REAL[md.nP];
		for (idx i = 0; i < md.nP; i++)
			u[k][i] = _f[i*ctx->angdata->aslm + k];
	}
	Wx = new REAL[md.nT];
	Wy = new REAL[md.nT];
	Wz = new REAL[md.nT];

	md.ComputeFlux(u[0], Wx, Wy, Wz);

	md._m->saveVtk(cfg.getOutFilename(), sizeof(REAL), 3, 1, Wx, Wy, Wz, u[0]);

	delete ctx;

	return 0;
}
