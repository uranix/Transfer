#include"struct3.h"
#include"error3.h"
#ifdef DEBUG
#include"debug3.h"
#endif


/* extern  variables */
extern StrucMesh3  mesh3;
#ifdef DEBUG
extern  StrucDebug3  debug3;
#endif


void errorExit3 (int group, char *number) {
	switch (group) {
		case  GE_Memory:
			printf("\nOut  of  memory!\n");
			break;
		case  GE_Critical:
			printf("\nCritical  error  in  program!\n");
			break;
		case  GE_User:
			printf("\nError  in  function  of  user!\n");
			break;
		case  GE_Temp:
			printf("\nTempoparal  error!\n");
			break;
	}
	printf("\n\n\n%s\n",number);

	printf("\nRESULT : %10d %10d\n", mesh3.nPoint, mesh3.nTetra);
	/* outMesh(); */
	exit(1);
	return;
} /* errorExit3 */

