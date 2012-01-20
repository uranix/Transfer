#include "Spherical.h"

#include <math.h>
#include <string.h>

void legendre(int d, double *v) {
	for (int i=0; i<d*d; i++)
		v[i] = 0;
	v[0] = 1;
	v[d+1] = 1;
	for (int i=2; i<d; i++) {
		for (int j=0; j<i; j++) {
			int s = i*d+j;
			v[s] -= v[s-d-d]*(i-1)/i;
			v[s+1] += v[s-d]*(2*i-1)/i;
		}
	}
}

void diff(int d, double *v) {
	for (int i=0; i<d; i++) {
		int s = i*d;
		for (int j=1; j<=i; j++) 
			v[s+j-1] = j*v[s+j];
		v[s+i] = 0;
	}
}

void legendre_all(int d, double *v) {
	double *w = v;
	legendre(d, w);
	for (int i=0; i<d-1; i++) {
		memcpy(w + d*d, w, d*d*sizeof(double));
		w += d*d;
		diff(d, w);
	}
}

void coeff_all(int d, double *c) {
	for (int l=0; l<d; l++)
		for (int m=0; m<=l; m++) {
			double C = (4.*l+2.)/(m==0?2.:1.);
			for (int i = l + 1 - m; i<= l+m; i++)
				C /= i;
			c[d*l+m] = sqrt(C);
		}
}

Spherical::Spherical(int topow) : d(topow+1) {
	v = new double [d*d*d];
	c = new double [d*d];

	legendre_all(d, v);
	coeff_all(d, c);
}

Spherical::~Spherical() {
	delete[] v;
	delete[] c;
}

double Spherical::value(int l, int m, double x, double y, double z) {
	return value(l, m, acos(z), atan2(y,x));
}

double Spherical::value(int l, int m, double theta, double phi) {
	double r;
	if (m < 0) {
		m = -m;
		r = sin(m*phi);
	}
	else
		r = cos(m*phi);

	r *= pow(sin(theta), m);
	double P = 0;
	int s = d*d*m+d*l;
	double z = cos(theta);
	for (int k = d-1; k >=0; k--)
		P = z*P + v[s + k];
	return r*P*c[d*l+m];
}
