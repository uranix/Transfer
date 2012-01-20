#include "AngularData.h"

#include "Spherical.h"
#include "LebedevQuad.h"

#include <math.h>

AngularData::AngularData(int maxk) {
	Spherical s(2*maxk);
	LebedevQuad q(4*maxk+2);

	slm = (1 + maxk)*(1 + 2*maxk);
	int aslm = slm; /* Not aligning in host version */
	omega = new REAL [3*aslm*aslm];
	omega_pos = new idx [3*aslm*aslm];
	for (int i=0; i<3*aslm*aslm; i++) {
		omega[i] = 0.;
		omega_pos[i] = -1;
	}

	for (int l1 = 0, lm1 = 0; l1 <=2; l1+=2)
		for (int m1 = -l1; m1 <= l1; m1++, lm1++)
			for (int l2 = 0, lm2 = 0; l2 <=2; l2+=2)
				for (int m2 = -l2; m2 <= l2; m2++, lm2++)
					for (int i = 0; i<3; i++) {
						for (int j = 0; j<3; j++) {
							REAL *r1 = (i==0)?q.x:(i==1)?q.y:q.z;
							REAL *r2 = (j==0)?q.x:(j==1)?q.y:q.z;
							REAL sum = 0;
							for (int p = 0; p < q.order; p++) {
								REAL t1 = s.value(l1, m1, q.x[p], q.y[p], q.z[p]);
								REAL t2 = s.value(l2, m2, q.x[p], q.y[p], q.z[p]);
								sum += q.w[p] * t1 * t2 * r1[p] * r2[p];
							}
							if ((fabs(sum) > 1e-10) || j == 2) {
								omega[i*aslm*aslm + aslm * lm1 + lm2] = sum;
								omega_pos[i*aslm*aslm + aslm * lm1 + lm2] = j;
								break;
							}
						}
					}
	Ox = new REAL[3*aslm*aslm];
	Oy = new REAL[3*aslm*aslm];
	Oz = new REAL[3*aslm*aslm];
}

AngularData::~AngularData() {
	delete[] omega;
	delete[] omega_pos;
	delete[] Ox;
	delete[] Oy;
	delete[] Oz;
}
