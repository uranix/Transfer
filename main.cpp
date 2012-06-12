#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "meshProcessor/mesh.h"

#include "AngularData.h"
#include "MeshData.h"
#include "CudaContext.h"
#include "Config.h"
#include "Eigenvalue.h"

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

int preck = 0;

void invprec(CudaContext *ctx, REAL *x, REAL *b, REAL *Ap, REAL *p, REAL *r);

double U1(double r);
double W1(double r);
double T1(double r);
double U2(double r);
double W2(double r);
double T2(double r);

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

	REAL *precp = ctx->allocVector();
	REAL *precAp= ctx->allocVector();
	REAL *precr = ctx->allocVector();

	REAL *_f = new REAL[ctx->angdata->aslm * ctx->meshdata->nP];

	idx N = ctx->N();
	printf("N = %d\n", N);

	if (cfg.doDump()) {
		REAL *_Af = new REAL[N];
		REAL * _b = new REAL[N];

		double *Z = new double[N*N];

		FILE *file = fopen("slae.dat", "w");
		ctx->computeRhs(b);
		ctx->copyToHost(_b, b, N * sizeof(REAL));

		double q = 1.0 / N, qq = 0;

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
					fprintf(file, "% 2.14g ", _Af[i*ctx->angdata->aslm+jj]);
			fprintf(file, "% 2.14g ", _b[k]);
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


#if 0 /* ----- eig ----- */
	ctx->computeRhs(b);
	double lmax;
	eigest(ctx, b, f, Ap, &lmax);  
	printf("Using %2.10e as Lmax estimate\n", lmax); 

#endif

	/* ----- cgs ----- */

	for (idx i = 0; i < ctx->N(); i++)
		_f[i] = 0;
	ctx->copyToDev(f, _f, ctx->N() * sizeof(REAL));
	ctx->copyToDev(z, _f, ctx->N() * sizeof(REAL));
	ctx->computeRhs(b);
	ctx->computeLhs(f, Ap);
	ctx->mulAdd(r, 0, b);
	ctx->addProd(r, Ap, -1);
	if (cfg.usePrec())
		invprec(ctx, z, r, precAp, precp, precr);
	else
		ctx->mulAdd(z, 0, r);
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
		if (nr < 1e-8)
			break;
		if (cfg.usePrec())
			invprec(ctx, z, r, precAp, precp, precr);
		else
			ctx->mulAdd(z, 0, r);
		double nrz2 = ctx->dot(z, r);
		double beta = nrz2/nrz;
		nrz = nrz2;
		ctx->mulAdd(p, beta, z);
		k++;
	}
	
	printf("preck = %d\n", preck);

	/*------ cgs end -----*/

	ctx->copyToHost(_f, f, ctx->N() * sizeof(REAL));
	REAL *I[ad.slm];
	for (idx k=0; k < ad.slm; k++) {
		I[k] = new REAL[md.nP];
		for (idx i = 0; i < md.nP; i++)
			I[k][i] = _f[i*ctx->angdata->aslm + k];
	}

	REAL *U = new REAL[md.nP];
	REAL *W[3];
	REAL *T[3][3];

	for (int i = 0; i < 3; i++)
		W[i] = new REAL[md.nT];
	
	for (int i = 0; i < 3; i++) 
		for (int j = 0; j < 3; j++)
			T[i][j] = new REAL[md.nP];

	md.ComputeMoments(ad, I, U, T);
	md.ComputeFlux(T, W);

	REAL *aU = new REAL[md.nP];
	REAL *aW = new REAL[md.nP];
	REAL *aT = new REAL[md.nP];

	for (idx i=0; i < md.nP; i++) {
		if (cfg.solType() == 1) {
			Vector r = md._m->vertices[i]->r;
			double z = r.norm();
			aU[i] = U1(z);
			aW[i] = W1(z);
			aT[i] = T1(z);
		}
		if (cfg.solType() == 2) {
			Vector r = md._m->vertices[i]->r;
			double z = r.x;
			aU[i] = U2(z);
			aW[i] = W2(z);
			aT[i] = T2(z);
		}
	}

	md._m->saveVtk(cfg.getOutFilename(), sizeof(REAL), "W%v", 
			"U%sT%taU%saW%saT%s", 
			W[0], W[1], W[2],
			U, 
			T[0][0], T[0][1], T[0][2], 
				     T[1][1], T[1][2], 
				              T[2][2],
			aU, aW, aT);

	delete ctx;

	return 0;
}

void invprec(CudaContext *ctx, REAL *f, REAL *b, REAL *Ap, REAL *p, REAL *r) {
	ctx->mulAddProd(f, 0, b, 0);
	ctx->computeLhsDiag(f, Ap);
	ctx->mulAdd(r, 0, b);
	ctx->addProd(r, Ap, -1);
	ctx->mulAdd(p, 0, r);
	int k = 0;
	double nrz = ctx->norm(r);
	nrz *= nrz;
	while (k < 1000) {
		ctx->computeLhsDiag(p, Ap);
		double alpha = nrz/ctx->dot(p, Ap);
		ctx->addProd(f, p, alpha);
		ctx->addProd(r, Ap, -alpha);
		double nr = ctx->norm(r);
		printf("[prec] k = %d norm r = %e\n", k, nr);
		if (nr < 1e-8)
			break;
		double nrz2 = nr*nr; 
		double beta = nrz2/nrz;
		nrz = nrz2;
		ctx->mulAdd(p, beta, r);
		k++;
	}
	preck += k;
}

double U1(double r) {
	double et;
	double s;
	if (r < 0.2769230769230769230769) {
		et = r;
		s = 0;
		s = et * s + (-2.11880084041759575135e+010);
		s = et * s + (3.43393746825775788375e+010);
		s = et * s + (-2.46363480185502942141e+010);
		s = et * s + (1.02883536610813649387e+010);
		s = et * s + (-2.76824271261706411026e+009);
		s = et * s + (5.01849769148549257802e+008);
		s = et * s + (-6.22990397777505843927e+007);
		s = et * s + (5.27483038420710597643e+006);
		s = et * s + (-2.98283361773050237296e+005);
		s = et * s + (1.05962730403749002818e+004);
		s = et * s + (-2.30453411343367735557e+002);
		s = et * s + (-1.24666203694156616331e+001);
		s = et * s + (-1.09666605825782615283e-002);
		s = et * s + (1.21028908279436689570e+001);		
		return s;
	}
	if (r < 0.3) {
		et = 0.09 - r*r;
		et = et > 0? sqrt(et) : 0;
		s = 0;
		s = et * s + (-6.85773771975857260490e+014);
		s = et * s + (5.58519227599132096854e+014);
		s = et * s + (-2.03513365245128352228e+014);
		s = et * s + (4.37959693369701296837e+013);
		s = et * s + (-6.18970106013871872712e+012);
		s = et * s + (6.04614258589792795095e+011);
		s = et * s + (-4.18678249794868131335e+010);
		s = et * s + (2.07847863241272682988e+009);
		s = et * s + (-7.43405199711879191052e+007);
		s = et * s + (1.94387111288820816042e+006);
		s = et * s + (-4.08764309826525618172e+004);
		s = et * s + (9.88940366132652978487e+002);
		s = et * s + (5.37939946086077840228e-001);
		s = et * s + (5.33236102854060962220e+000);
		return s;
	}
	if (r < 0.34510867853474794933) {
		et = (r*r-0.09)/(3-0.09);
		et = et > 0? sqrt(et) : 0;
		s = 0;
		s = et * s + (1.79954478964540546536e+015);
		s = et * s + (-1.26999567405664114254e+015);
		s = et * s + (4.00973362277578285280e+014);
		s = et * s + (-7.47615493043627453883e+013);
		s = et * s + (9.15323574065913694227e+012);
		s = et * s + (-7.74347500675467556730e+011);
		s = et * s + (4.64168284791431416959e+010);
		s = et * s + (-1.99242070551878420453e+009);
		s = et * s + (6.14131176173976509910e+007);
		s = et * s + (-1.36566674796372251004e+006);
		s = et * s + (2.27468966487741300497e+004);
		s = et * s + (-3.49555801097923724719e+002);
		s = et * s + (-2.54513977718950306759e-001);
		s = et * s + (1.67386675966713077971e+000);
		return exp(s);
	}
	et = (r*r-0.09)/(3-0.09);
	et = et > 0? sqrt(et) : 0;
	s = 0;
	s = et * s + (-1.83672819173466640958e+002);
	s = et * s + (1.32783976645613893652e+003);
	s = et * s + (-4.23896487615352271604e+003);
	s = et * s + (7.78279861516373977843e+003);
	s = et * s + (-8.86784384396269128545e+003);
	s = et * s + (6.07270363746468148573e+003);
	s = et * s + (-1.72649452604486017227e+003);
	s = et * s + (-1.01430416195816327120e+003);
	s = et * s + (1.44935298413292195060e+003);
	s = et * s + (-8.34560715014520449044e+002);
	s = et * s + (2.83703619350034748739e+002);
	s = et * s + (-5.25611818440145473229e+001);
	s = et * s + (-3.60754314920493664929e+000);
	s = et * s + (1.69738499278806078916e+000);
	return exp(s);
}

double W1(double r) {
	double et;
	double s;
	if (r < 0.2769230769230769230769) {
		et = r;
		s = 0;
		s = et * s + (1.97740840297776218112e+009);
		s = et * s + (-3.18018297055580323840e+009);
		s = et * s + (2.26801732688725966848e+009);
		s = et * s + (-9.42522542328579489792e+008);
		s = et * s + (2.52589985087748407296e+008);
		s = et * s + (-4.56381193984082837504e+007);
		s = et * s + (5.65144884349181296640e+006);
		s = et * s + (-4.77407032550683705344e+005);
		s = et * s + (2.72220470575009103872e+004);
		s = et * s + (-9.71843268777485139968e+002);
		s = et * s + (5.37888246597859016704e+001);
		s = et * s + (-2.29150134574863908864e-001);
		s = et * s + (1.70044131552613531648e+000);
		s = et * s + (-6.98732380843836178432e-007);
		return s;
	}
	if (r < 0.3) {
		et = 0.09 - r*r;
		et = et > 0? sqrt(et) : 0;
		s = 0;
		s = et * s + (-4.71116314394159284224e+011);
		s = et * s + (3.90903031077199413248e+011);
		s = et * s + (-1.45934497357725122560e+011);
		s = et * s + (3.24433622412406816768e+010);
		s = et * s + (-4.79775790275191177216e+009);
		s = et * s + (5.00880550792231714816e+008);
		s = et * s + (-3.85410878886819004416e+007);
		s = et * s + (2.30699640029769203712e+006);
		s = et * s + (-1.20224324760459821056e+005);
		s = et * s + (6.15209993179449524224e+003);
		s = et * s + (2.25798077839204155392e+001);
		s = et * s + (-9.93663995788049514496e+001);
		s = et * s + (1.08253838861038780416e-004);
		s = et * s + (2.99884182278352601088e+000);
		return s;
	}
	if (r < 0.34510867853474794933) {
		et = (r*r-0.09)/(3-0.09);
		et = et > 0? sqrt(et) : 0;
		s = 0;
		s = et * s + (4.15948157989178572800e+011);
		s = et * s + (-2.81180642823378698240e+011);
		s = et * s + (8.38194737243012988928e+010);
		s = et * s + (-1.44218181168364519424e+010);
		s = et * s + (1.56603976159004065792e+009);
		s = et * s + (-1.08721185808400744448e+008);
		s = et * s + (4.44536574832441950208e+006);
		s = et * s + (-5.67840072770654765056e+004);
		s = et * s + (-6.80523373388829032448e+003);
		s = et * s + (1.00104437334385115136e+003);
		s = et * s + (9.27968474618780450816e+000);
		s = et * s + (-4.09977789007185444864e+001);
		s = et * s + (9.38935465590873391104e-005);
		s = et * s + (1.09822614644879015936e+000);
		return exp(s);
	}
	et = (r*r-0.09)/(3-0.09);
	et = et > 0? sqrt(et) : 0;
	s = 0;
	s = et * s + (-6.93285710426655162368e+002);
	s = et * s + (5.24549216049522671616e+003);
	s = et * s + (-1.78431474777974571008e+004);
	s = et * s + (3.59819734764197773312e+004);
	s = et * s + (-4.76983233318159777792e+004);
	s = et * s + (4.34514443514294566912e+004);
	s = et * s + (-2.74632843494618431488e+004);
	s = et * s + (1.16942466783187337216e+004);
	s = et * s + (-2.95794188106609459200e+003);
	s = et * s + (1.67973608841035317248e+002);
	s = et * s + (1.68129588323524050944e+002);
	s = et * s + (-5.93762432627748110336e+001);
	s = et * s + (1.10748429636623335424e+000);
	s = et * s + (1.07091913709027606528e+000);	
	return exp(s);
}

double T1(double r) {
	double et;
	double s;
	if (r < 0.2769230769230769230769) {
		et = r;
		s = 0;
		s = et * s + (-7.51559930774475112448e+008);
		s = et * s + (1.20628045139377078272e+009);
		s = et * s + (-8.59026433853490069504e+008);
		s = et * s + (3.56561417902624735232e+008);
		s = et * s + (-9.54651121290292363264e+007);
		s = et * s + (1.72333358798490238976e+007);
		s = et * s + (-2.13235115864461508608e+006);
		s = et * s + (1.79592125407696224256e+005);
		s = et * s + (-1.01585240161849229312e+004);
		s = et * s + (2.88837386659760668672e+002);
		s = et * s + (-7.82472308140753027072e+000);
		s = et * s + (-6.65482446294618210304e+000);
		s = et * s + (-3.71857943026299109376e-004);
		s = et * s + (4.03429461771170414592e+000);
		return s;
	}
	if (r < 0.3) {
		et = 0.09 - r*r;
		et = et > 0? sqrt(et) : 0;
		s = 0;
		s = et * s + (1.02719119372039962624e+011);
		s = et * s + (-8.49820920010218143744e+010);
		s = et * s + (3.17418244010389536768e+010);
		s = et * s + (-7.09532100036383997952e+009);
		s = et * s + (1.06237557097176317952e+009);
		s = et * s + (-1.13344766435211395072e+008);
		s = et * s + (9.01839944596069285888e+006);
		s = et * s + (-5.68070424387097985024e+005);
		s = et * s + (3.31048291034671546368e+004);
		s = et * s + (-2.25133638978220982272e+003);
		s = et * s + (-7.73709176640455442432e+000);
		s = et * s + (5.95867008723164594176e+001);
		s = et * s + (-4.61735886596592631808e-005);
		s = et * s + (2.05243236609067253760e+000);		
		return s;
	}
	if (r < 0.34510867853474794933) {
		et = (r*r-0.09)/(3-0.09);
		et = et > 0? sqrt(et) : 0;
		s = 0;
		s = et * s + (2.12125759591631159296e+012);
		s = et * s + (-1.51580831159405182976e+012);
		s = et * s + (4.86576731461709398016e+011);
		s = et * s + (-9.27841533116361080832e+010);
		s = et * s + (1.17203720919570956288e+010);
		s = et * s + (-1.03749846574271004672e+009);
		s = et * s + (6.68583108413558620160e+007);
		s = et * s + (-3.29330619153092640768e+006);
		s = et * s + (1.34681488895962529792e+005);
		s = et * s + (-4.13925200323301801984e+003);
		s = et * s + (-2.41683551364721803264e+001);
		s = et * s + (-1.35231020909353713664e+001);
		s = et * s + (-7.63025094969727123456e-005);
		s = et * s + (7.19025621861797855232e-001);		
		return exp(s);
	}
	et = (r*r-0.09)/(3-0.09);
	et = et > 0? sqrt(et) : 0;
	s = 0;
	s = et * s + (-1.51044850944340426752e+003);
	s = et * s + (1.15914338241695563776e+004);
	s = et * s + (-4.01579399609017630720e+004);
	s = et * s + (8.29752091793223450624e+004);
	s = et * s + (-1.13746604128831373312e+005);
	s = et * s + (1.08788250308062986240e+005);
	s = et * s + (-7.41857426352939991040e+004);
	s = et * s + (3.60857182450387517440e+004);
	s = et * s + (-1.22080628158646960128e+004);
	s = et * s + (2.65578220477532405760e+003);
	s = et * s + (-2.78046694923063885824e+002);
	s = et * s + (-1.49286693375251415040e+001);
	s = et * s + (7.59567228205288128512e-001);
	s = et * s + (6.92310233006297841664e-001);	
	return exp(s);
}

double xE1(double x);

double U2(double x) {
	double r = 2 * M_PI * (exp(-fabs(x)) - xE1(fabs(x)));
	if (-x > 0)
		return 4 * M_PI - r;
	return r;
}

double W2(double x) {
	double r = M_PI * ((1-fabs(x)) * exp(-fabs(x)) + fabs(x) * xE1(fabs(x)));
	return r;
}

double T2(double x) {
	double r = M_PI / 3.0 * ((2 - fabs(x) + x*x) * exp(-fabs(x)) - x * x * xE1(fabs(x)));
	if (-x > 0)
		return 4 * M_PI / 3.0 - r;
	return r;
}
