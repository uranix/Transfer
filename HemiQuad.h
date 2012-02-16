#ifndef __HEMIQUAD_H__
#define __HEMIQUAD_H__

struct HemiQuad {
	int order;
	double *x;
	double *y;
	double *z;
	double *w;
	HemiQuad(int mindegree);
	~HemiQuad();
};

#endif 
