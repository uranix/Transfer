#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "AngularData.h"

void torat(double x, int *rp, int *rq) {
	int p = floor(x), p_ = 1;
	int q = 1, q_ = 0;
	double z = 1./(x - p);

	while (fabs(x - (1.0*p)/q) > 1e-10) {
		double a = floor(z), t;
		z = 1./(z-a);
		t = p;
		p = a*p + p_;
		p_= t;
		t = q;
		q = a*q + q_;
		q_ = t;
	}
	rp[0] = p;
	rq[0] = q;
}

char *frac(double x) {
	int p, q;
	torat(x, &p, &q);
	char *ret = (char *)malloc(1024);
	if (q == 1)
		sprintf(ret, "%d", p);
	else
		sprintf(ret, "%d/%d", p, q);

	return ret;
}

int main() {
	AngularData dat(1);
	return 0;
}
