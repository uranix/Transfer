#include "Eigenvalue.h"
#include <stdio.h>
#include <math.h>

double trieig(REAL *a, REAL *b, idx N) {
	if (N == 1)
		return a[0];
	if (N == 2) {
		double ap = 0.5 * (a[0] + a[1]);
		double am = 0.5 * (a[0] - a[1]);
		return ap + sqrt(am * am + b[0] * b[0]);
	}
	double q, lmax = a[0] + fabs(b[0]);
	for (idx i = 1; i < N-1; i++) {
		q = a[i] + fabs(b[i]) + fabs(b[i-1]);
		if (q > lmax)
			lmax = q;
	}
	q = a[N-1] + fabs(b[N-2]);
	if (q > lmax)
		lmax = q;

	int iters = 0;
	double dl = 1;
	double th = 100 * uround();
	while (iters < 100 && dl > th) {
		double x = lmax;
		double Po = 1;
		double P = a[0] - x;
		double Qo = 0;
		double Q = -1;
		for (idx j = 1; j < N; j++) {
			double Qs = Q;
			Q = -P + (a[j] - x) * Q - b[j-1]*b[j-1] * Qo;
			Qo = Qs;
			double Ps = P;
			P = (a[j] - x) * P - b[j-1]*b[j-1] * Po;
			Po = Ps;
		}
		dl = P/Q;
		lmax -= dl;
		iters ++;
	}

	if (dl > th)
		printf("Newton process did not converge!. dl = %e, lmax = %e\n", dl, lmax);
	else
		printf("Newton process converged to %2.10e in %d iters (N = %d)\n", lmax, iters, N);

	return lmax;
}

double eigest(CudaContext *ctx, REAL *v, REAL *v0, REAL *w, double *extra) {
	idx m = 30;
	REAL a[m];
	REAL b[m+1];
	REAL *s;
	ctx->scale(v, (REAL)1/ctx->norm(v));
	ctx->mulAdd(v0, 0, v);
	b[0] = 0;
	double lmax = 0, lold = 0, loold;
	idx j;
	printf("Estimating first eigenvalue\n");
	for (j = 0; j < m; j++) {
		ctx->computeLhs(v, w);
		ctx->mulAdd(v0, -b[j], w);
		s = v0;
		v0 = v;
		v = s;
		a[j] = ctx->dot(v, v0);
		ctx->addProd(v, v0, -a[j]);
		b[j+1] = ctx->norm(v);
		ctx->scale(v, (REAL)1/b[j+1]);
		loold = lold;
		lold = lmax;
		lmax = trieig(a, b, j+1);
		if (fabs(lmax - lold) < 1e-6 * lmax) {
			j++;
			break;
		}
	}
	double q = (lmax - lold)/(lold - loold);
	printf("Done %d Arnoldi iterations. Approximate lmax = %e (diff = %e, extrapolated diff = %e)\n", 
			j, lmax, lmax - lold, (lmax - lold)/(1-q));

	*extra = lmax + (lmax - lold) / (1-q);
	return lmax;
}
