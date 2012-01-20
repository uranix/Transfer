#ifndef __SPHERICAL_H__
#define __SPHERICAL_H__

class Spherical {
	int d;
	double *v;
	double *c;
public:	
	Spherical(int topow);
	~Spherical();
	double value(int l, int m, double x, double y, double z);
	double value(int l, int m, double theta, double phi);
};

#endif
