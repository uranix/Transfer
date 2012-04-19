#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "AngularData.h"
#include "MeshData.h"

int main() {
	AngularData ad(1); /* maxk = 1, maxl = 2*/
	MeshData md("mesh.vol");

	DeviceAngularData dad(ad);
	DeviceMeshData dmd(md);

	return 0;
}
