#ifndef __MESH3D__VECTOR_H__
#define __MESH3D__VECTOR_H__

#include <stdio.h>

struct Vector {
	double x, y, z;
	Vector(): x(0), y(0), z(0) {}
	Vector(const Vector &r): x(r.x), y(r.y), z(r.z) {}
	Vector(double _x, double _y, double _z): x(_x), y(_y), z(_z) {}
	void scale(double a) {
		x *= a;
		y *= a;
		z *= a;
	}
	void add(const Vector &r) {
		x += r.x;
		y += r.y;
		z += r.z;
	}
	void sub(const Vector &r) {
		x -= r.x;
		y -= r.y;
		z -= r.z;
	}
	void madd(double a, const Vector &r) {
		x += a*r.x;
		y += a*r.y;
		z += a*r.z;
	}
	double dot(const Vector &r) const {
		return x*r.x + y*r.y + z*r.z;
	}
	void cross(const Vector &b, Vector &c) const {
		c.x = y*b.z - z*b.y;
		c.y = z*b.x - x*b.z;
		c.z = x*b.y - y*b.x;
	}
	void dotscale(const Vector &r) {
		x *= r.x;
		y *= r.y;
		z *= r.z;
	}
	double norm() const {
		return sqrt(x*x+y*y+z*z);
	}
	void serialize(FILE *f) const {
		int foo;
		foo = fwrite(&x, sizeof(x), 1, f);
		foo = fwrite(&y, sizeof(y), 1, f);
		foo = fwrite(&z, sizeof(z), 1, f);
		foo ++;
	}
	Vector(FILE *f) {
		int foo;
		foo = fread(&x, sizeof(x), 1, f);
		foo = fread(&y, sizeof(y), 1, f);
		foo = fread(&z, sizeof(z), 1, f);
		foo ++;
	}
};

#endif 
