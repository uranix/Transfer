int aft3dboundary (
		int nVVert, double *VVertxyz,
		int nLine, int *LineD, int *LineP, double *LineT,
		int nSurface, int *SurfL, int *SurfI, double *SurfT,
		int *pnVout, double *vertexout,
		int *pnFout, int *faceout, int *facecolor,
		int maxnV, int maxnF, int indexshift);
int aft3dsurfmeshrefiner(int *nV, double *vertex, int *nF, int *face, int *facematerial, int maxnV, int maxnF, int indexshift);

