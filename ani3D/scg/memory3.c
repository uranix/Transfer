#include"struct3.h"
#include"error3.h"
#include"memory3.h"
#include"tree3.h"


/* extern  variables */
extern  StrucMesh3  mesh3;
extern  StrucTree3  tree3;
extern  StrucTree3  tree32;


void *myAlloc (size_t n) {
    void  *p;

    p = (void *)malloc(n);
    if (p == NULL)
	errorExit3(2,"myAlloc");

    return  p;
} /*myAlloc*/


void initMemory (void) {
    int        i,allSize=0,size[100];
    size_t     maxMemory=0,allMemory=0,sMemory[100];
    char  *    pMemory=NULL;
    int maxTetra=0;

    /*   maxTetra = findMaxTetra();   */ //by Danilov only for tests
    maxTetra = 6400000;

    //	printf("maxTetra = %7d\n", maxTetra);

    maxMemory = 48*maxTetra;


    if (maxMemory == 0)
	errorExit3(2, "Memory");

    /* init  distribution  of  memory */
    size[0] = 1 * S_StrucVert3;    /*VERT*/
    size[1] = 5 * S_StrucTetra3;   /*TETRA*/

    allMemory = 0;
    for (i=0; i<2; i++) allSize += size[i];
    for (i=0; i<2; i++) {
	sMemory[i] = maxMemory*size[i] / allSize;
	allMemory += sMemory[i];
    }
    /* end  distribution  of  memory */

    /* init  memory  ptr */
    pMemory = (char*)myAlloc(allMemory);

    mesh3.maxPoint = sMemory[0] / S_StrucVert3;
    mesh3.vert = (PStrucVert3)pMemory;
    pMemory += mesh3.maxPoint*S_StrucVert3;
    mesh3.maxTetra = sMemory[1] / S_StrucTetra3;
    mesh3.tetra = (PStrucTetra3)pMemory;
    pMemory += mesh3.maxTetra*S_StrucTetra3;
    /* end  init  memory  ptr */

    tree3.root = myAlloc(S_StrucNode3);
    tree3.maxFace = 2 * mesh3.maxPoint;
    if (tree3.maxFace < 3000) tree3.maxFace = 3000;
    tree3.face = myAlloc(tree3.maxFace * sizeof(PStrucFace3));
    tree3.maxVicinityFace = 109990;
    tree3.vicinityFace = myAlloc(tree3.maxVicinityFace * sizeof(StrucFace4));

    tree32.root = myAlloc(S_StrucNode3);
    tree32.maxFace = mesh3.maxPoint/5;
    if (tree32.maxFace < 2000) tree32.maxFace = 2000;
    tree32.face = myAlloc(tree32.maxFace * sizeof(PStrucFace3));
    tree32.maxVicinityFace = 4409;
    tree32.vicinityFace = myAlloc(tree32.maxVicinityFace * sizeof(StrucFace4));

    return;
} /*initMemory*/

void freeMemory (void) {
    int i;
    free((void*)mesh3.vert);
    for (i=tree3.nFace-1; i>=0; i--)  remFace(tree3.face[i]);
    free(tree3.root);
    free(tree3.face);
    free(tree3.vicinityFace);
    for (i=tree32.nFace-1; i>=0; i--)  remFace(tree32.face[i]);
    free(tree32.root);
    free(tree32.face);
    free(tree32.vicinityFace);
    return;
} /*freeMemory*/


