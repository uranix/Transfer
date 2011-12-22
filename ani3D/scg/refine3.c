#include"struct3.h"
#include"memory3.h"
#include"tree3.h"
#include"tria32.h"
#include"user3.h"
#include"error3.h"
#include"region3.h"
#include"refine3.h"


/* extern  variables */
extern StrucMesh3  mesh3;
extern StrucRegion3  reg3;

static void calcNeigTetra2 (void) {
	int i, j, v1;

	for (i=0; i<mesh3.nPoint; i++)
		mesh3.neigTetra[i].n = 0;
	for (i=0; i<reg3.nTria; i++) {
		v1 = reg3.v1[i];
		j = mesh3.neigTetra[v1].n;
		if (j >= MAX_NEIGBOR)
			errorExit3(3, "j >= NAX_NEIGBOR");
		mesh3.neigTetra[v1].n++;
      mesh3.neigTetra[v1].neig[j] = i;

      v1 = reg3.v2[i];
      j = mesh3.neigTetra[v1].n;
      if (j >= MAX_NEIGBOR)
        errorExit3(3,"j >= NAX_NEIGBOR");
      mesh3.neigTetra[v1].n++;
      mesh3.neigTetra[v1].neig[j] = i;

      v1 = reg3.v3[i];
      j = mesh3.neigTetra[v1].n;
      if (j >= MAX_NEIGBOR)
        errorExit3(3,"j >= NAX_NEIGBOR");
      mesh3.neigTetra[v1].n++;
      mesh3.neigTetra[v1].neig[j] = i;
   }
   return;
} /*calcNeigTetra2*/


static void calcNeigbor2 (void) {
	int i, j, k, n, iTria, vert[3*MAX_NEIGBOR];

	for (i=0; i<mesh3.nPoint; i++) {
		if (i < mesh3.nLinePoint)
			mesh3.neigbor[i].n = -1;
		else
			mesh3.neigbor[i].n = 0;
	}
	for (i=0; i<mesh3.nPoint; i++) {
		n = 0;
		for (j=0; j<mesh3.neigTetra[i].n; j++) {
			iTria = mesh3.neigTetra[i].neig[j];
			vert[3*j+0] = reg3.v1[iTria];
			vert[3*j+1] = reg3.v2[iTria];
			vert[3*j+2] = reg3.v3[iTria];
		}
		for (j=0; j<3*mesh3.neigTetra[i].n; j++) {
			if (vert[j] == i)
				continue;
			iTria = 0;
			for (k=0; k<n; k++)
				if (vert[j] == mesh3.neigbor[i].neig[k]) {
					iTria = 1;
					break;
				}
			if (iTria) continue;
			mesh3.neigbor[i].neig[n] = vert[j];
			n++;
		}
		if (mesh3.neigbor[i].n < 0)
			mesh3.neigbor[i].n = -n;
		else
			mesh3.neigbor[i].n = n;
	}
   return;
} /*calcNeigbor2*/


void smoothingSurf (void) {
	int    i, j, k, n, nn;
	double x0, y0, z0, xx, yy, zz, x, y, z, u, v;
	double f1, du, dv, du1, dv1, du2, dv2, maxdu, maxdv, d1u, d1v, size;

	/* !!! alloc only for  2*mesh3.nPoint */
	mesh3.neigTetra = (PStrucNeigbor)myAlloc(2*mesh3.nPoint*S_StrucNeigbor);   // cbD
	mesh3.neigbor   = (PStrucNeigbor)myAlloc(2*mesh3.nPoint*S_StrucNeigbor);   // cbD
	calcNeigTetra2();
	calcNeigbor2();

	for (i=0; i<mesh3.nPoint; i++) {
		n = mesh3.neigbor[i].n;
		if (n > 0) {
			xx = yy = zz = 0.0;
			for (j=0; j<n; j++) {
				k = mesh3.neigbor[i].neig[j];
				nn = mesh3.neigbor[k].n;
				x = mesh3.vert[k].x;
				y = mesh3.vert[k].y;
				z = mesh3.vert[k].z;
				xx += x;  yy += y;  zz += z;
			}
			xx /= n;  yy /= n;  zz /= n;

			u = mesh3.vert[i].u;
			v = mesh3.vert[i].v;
			bounSurf(reg3.iSurf, u, v, &x, &y, &z);
			size = sizeFace(x, y, z);
			f1 = distance(xx, yy, zz, x, y, z);
			for (j=0; j<15; j++) { /* NEW METH */
				du = 0.2*(reg3.uMax - reg3.uMin);
				do {
					du /= 2.0;
					bounSurf(reg3.iSurf, u+du, v, &x0, &y0, &z0);
				} while (distance(x, y, z, x0, y0, z0) > 0.01*size);
				d1u = (distance(xx, yy, zz, x0, y0, z0) - f1)/du;

				dv = 0.2*(reg3.vMax - reg3.vMin);
				do {
					dv /= 2.0;
					bounSurf(reg3.iSurf, u, v+dv, &x0, &y0, &z0);
				} while (distance(x, y, z, x0, y0, z0) > 0.01*size);
				d1v = (distance(xx, yy, zz, x0, y0, z0) - f1)/dv;

				du = 0.6*(reg3.uMax - reg3.uMin);
				do {
					du /= 2.0;
					bounSurf(reg3.iSurf, u+du, v, &x0, &y0, &z0);
				} while (distance(x, y, z, x0, y0, z0) > 0.1*f1);

				dv = 0.6*(reg3.vMax - reg3.vMin);
				do {
					dv /= 2.0;
					bounSurf(reg3.iSurf, u, v+dv, &x0, &y0, &z0);
				} while (distance(x, y, z, x0, y0, z0) > 0.1*f1);

				maxdu = du;
				maxdv = dv;
				if (d1u == 0. && d1v == 0)
					errorExit3(3, " der == 0.   in  smooth ");
				if (fabs(d1u) > fabs(d1v)) {
					dv1 = maxdv;
					du1 = (-f1 - d1v*dv)/d1u;
					dv2 = -maxdv;
					du2 = (-f1 - d1v*dv)/d1u;
					if (fabs(du1) < fabs(du2)) {
						du = du1;
						dv = dv1;
					} else {
						du = du2;
						dv = dv2;
					}
					if (du > maxdu) {
						dv *= (maxdu/du);
						du = maxdu;
					}
					if (du < -maxdu) {
						dv *= (-maxdu/du);
						du = -maxdu;
					}
				} else {
					du1 = maxdu;
					dv1 = (-f1 - d1u*du)/d1v;
					du2 = maxdu;
					dv2 = (-f1 - d1u*du)/d1v;
					if (fabs(dv1) < fabs(dv2)) {
						du = du1;
						dv = dv1;
					} else {
						du = du2;
						dv = dv2;
					}
					if (dv > maxdv) {
						du *= (maxdv/dv);
						dv = maxdv;
					}
					if (dv < -maxdv) {
						du *= (-maxdv/dv);
						dv = -maxdv;
					}
				}
				u += du;
				v += dv;
				bounSurf(reg3.iSurf, u, v, &x, &y, &z);
				f1 = distance(xx, yy, zz, x, y, z);
			}/*NEW METH*/
			mesh3.vert[i].u = u;
			mesh3.vert[i].v = v;
			mesh3.vert[i].x = x;
			mesh3.vert[i].y = y;
			mesh3.vert[i].z = z;
		}
	}/* for i */

	free(mesh3.neigbor);
	free(mesh3.neigTetra);

	return;
} /*smoothingSurf*/


