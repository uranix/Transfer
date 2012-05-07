#include "util.h"
#include <stdio.h>

REAL uround() {
	static REAL _uround = 0;
	if (_uround > 0)
		return _uround;
	printf("REAL = %s\n", sizeof(REAL) == 4 ? "float" : sizeof(REAL) == 8 ? "double" : "unknown!");
	REAL x = 1.0;
	while ((REAL)1.0 + x > (REAL)1.0)
		x = (REAL)0.5 * x;
	printf("uround = %2.10e\n", x);
	return _uround = x;
}
