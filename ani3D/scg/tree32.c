#include"error3.h"
#include"memory3.h"
#include"struct3.h"
#include"section3.h"
#include"tree3.h"
#include"tree32.h"
#ifdef DEBUG
#include"debug3.h"
#endif


/* extern  variables */
extern StrucMesh3  mesh3;
extern StrucSect3  sect3;
#ifdef DEBUG
extern  StrucDebug3  debug3;
#endif


/* global  variables */
StrucTree3  tree32;
extern StrucTree3 tree3;


/* static  function */
static void addFaceArray( PStrucFace3  face );
static void remFaceArray( int  face );
static void remFaceRec( PStrucNode3  node, PStrucFace3  face );
static int sectCube( void );
static void vicinityFacesRec( PStrucNode3  node );


static void addFaceArray( PStrucFace3  face )
{
	int son,fath,i;
	PStrucFace3  cface,*tface=tree32.face;
	if( tree32.nFace >= tree32.maxFace )
		errorExit3(1,"tree32.nFace");

	tface[tree32.nFace] = face;  face->f = tree32.nFace;
	son=++tree32.nFace;
aaa:fath=son/2;

#ifdef DEBUG
	 if( mesh3.debug && debug3.bool.sortZ ){
		 if( (fath>0) && (tface[son-1]->z < tface[fath-1]->z)  ){
			 cface=tface[son-1];
			 tface[son-1]=tface[fath-1];
			 tface[fath-1]=cface;
			 i=cface->f;  cface->f=tface[son-1]->f;  tface[son-1]->f=i;
			 son=fath; if(son>1) goto aaa;
		 }
	 }
	 else{
		 if( (fath>0) && (tface[son-1]->s < tface[fath-1]->s)  ){
			 cface=tface[son-1];
			 tface[son-1]=tface[fath-1];
			 tface[fath-1]=cface;
			 i=cface->f;  cface->f=tface[son-1]->f;  tface[son-1]->f=i;
			 son=fath; if(son>1) goto aaa;
		 }
	 }
#else
	 if( (fath>0) && (tface[son-1]->s < tface[fath-1]->s)  ){
		 cface=tface[son-1];
		 tface[son-1]=tface[fath-1];
		 tface[fath-1]=cface;
		 i=cface->f;  cface->f=tface[son-1]->f;  tface[son-1]->f=i;
		 son=fath; if(son>1) goto aaa;
	 }
#endif

	 return;
};/*addFaceArray*/


static void remFaceArray( int  fath )
{
	int son1,son2,i;
	PStrucFace3 face,*tface;
	tface = tree32.face;
	if( (fath >= tree32.nFace) || (fath<0) )  errorExit3(2," tree32.nFace < 0   in  rem_face ");

	tree32.nFace--;
	tface[fath]=tface[tree32.nFace];
	tface[fath]->f=fath;
aaa:i = 0;
#ifdef DEBUG
	 if( mesh3.debug && debug3.bool.sortZ ){
		 if( 2*fath+1 < (tree32.nFace-1) ){
			 son1=2*fath+1;  son2=son1++;
			 if((tface[fath]->z<=tface[son1]->z)&&(tface[fath]->z<=tface[son2]->z)) i=0;
			 if((tface[fath]->z>=tface[son1]->z)&&(tface[son1]->z>=tface[son2]->z)) i=son2;
			 if((tface[fath]->z>=tface[son2]->z)&&(tface[son2]->z>tface[son1]->z)) i=son1;
			 if((tface[son1]->z>=tface[fath]->z)&&(tface[fath]->z>tface[son2]->z)) i=son2;
			 if((tface[son2]->z>=tface[fath]->z)&&(tface[fath]->z>tface[son1]->z)) i=son1;
		 }
		 if( 2*fath+1 == (tree32.nFace-1) ){
			 son1=2*fath+1;
			 if((tface[fath]->z>tface[son1]->z)) i=son1;else i=0;
		 }
	 }
	 else{
		 if( 2*fath+1 < (tree32.nFace-1) ){
			 son1=2*fath+1;  son2=son1++;
			 if((tface[fath]->s<=tface[son1]->s)&&(tface[fath]->s<=tface[son2]->s)) i=0;
			 if((tface[fath]->s>=tface[son1]->s)&&(tface[son1]->s>=tface[son2]->s)) i=son2;
			 if((tface[fath]->s>=tface[son2]->s)&&(tface[son2]->s>tface[son1]->s)) i=son1;
			 if((tface[son1]->s>=tface[fath]->s)&&(tface[fath]->s>tface[son2]->s)) i=son2;
			 if((tface[son2]->s>=tface[fath]->s)&&(tface[fath]->s>tface[son1]->s)) i=son1;
		 }
		 if( 2*fath+1 == (tree32.nFace-1) ){
			 son1=2*fath+1;
			 if((tface[fath]->s>tface[son1]->s)) i=son1;else i=0;
		 }
	 }
#else
	 if( 2*fath+1 < (tree32.nFace-1) ){
		 son1=2*fath+1;  son2=son1++;
		 if((tface[fath]->s<=tface[son1]->s)&&(tface[fath]->s<=tface[son2]->s)) i=0;
		 if((tface[fath]->s>=tface[son1]->s)&&(tface[son1]->s>=tface[son2]->s)) i=son2;
		 if((tface[fath]->s>=tface[son2]->s)&&(tface[son2]->s>tface[son1]->s)) i=son1;
		 if((tface[son1]->s>=tface[fath]->s)&&(tface[fath]->s>tface[son2]->s)) i=son2;
		 if((tface[son2]->s>=tface[fath]->s)&&(tface[fath]->s>tface[son1]->s)) i=son1;
	 }
	 if( 2*fath+1 == (tree32.nFace-1) ){
		 son1=2*fath+1;
		 if((tface[fath]->s>tface[son1]->s)) i=son1;else i=0;
	 }
#endif

	 if(i)
	 {  face=tface[i];  tface[i]=tface[fath];   tface[fath]=face;
		 son1=face->f;  face->f=tface[i]->f;  tface[i]->f=son1;
		 fath=i;   goto aaa;
	 }

	 return;
}; /*remFaceArray*/


#define EPS 1e-8
PStrucFace3 addFace32(int v1, int v2 )
{
	int          flag,j,d,d1,bool;
	double       x,y,z,x1=-1.0,y1=-1.0,z1=-1.0;
	double       xc=tree3.boxcx, yc=tree3.boxcy, zc=tree3.boxcz, side=tree3.boxsize;
	PStrucFace3  old_face,face;
	PStrucNode3  node;
	PStrucList3  list;

	node=tree32.root;   flag=node->flag;
	list = ( PStrucList3 ) node->faceNode;
	face = myAlloc( S_StrucFace3 );
	face->fail = 0;
	face->v1 = v1;   face->v2 = v2;  
	x = (mesh3.vert[v1].x+mesh3.vert[v2].x)/2.;
	y = (mesh3.vert[v1].y+mesh3.vert[v2].y)/2.;
	z = (mesh3.vert[v1].z+mesh3.vert[v2].z)/2.;
	face->x=x;   face->y=y;   face->z=z;
	face->s = distance(mesh3.vert[v1].x,mesh3.vert[v1].y,mesh3.vert[v1].z,mesh3.vert[v2].x,mesh3.vert[v2].y,mesh3.vert[v2].z);

	if ((fabs(x-tree3.boxcx)>tree3.boxsize+EPS) || (fabs(y-tree3.boxcy)>tree3.boxsize+EPS) || (fabs(z-tree3.boxcz)>tree3.boxsize+EPS)) {
		printf("\nbox: x: %le, y: %le, z: %le, size: %le\n", tree3.boxcx, tree3.boxcy, tree3.boxcz, tree3.boxsize);
		printf("x: %le, y: %le, z: %le\n", x, y, z);
		errorExit3(3,"x_y_z_insert");
	}
	addFaceArray( face );
	return face;

	while(flag==TREE){
		d = direction(x,y,z,xc,yc,zc);  side *= 0.5;
		center(&xc,&yc,&zc,side,d);
		node = (*list)[d];
		flag = node->flag;  list = ( PStrucList3 ) node->faceNode;
	}
	if (flag==EMPTY){
		node->flag=FACE; node->faceNode=face;
	}
	if (flag==FACE){
		bool = 1;  old_face = node->faceNode;
		x1 = old_face->x;  y1 = old_face->y;  z1 = old_face->z;
		if( (x1==x) && (y1==y) && (z1==z) )
			errorExit3(3," coincident   in  insert tree32");

		while(bool){
			node->flag = TREE;
			node->faceNode = myAlloc(S_StrucList3);;
			list = ( PStrucList3 ) node->faceNode;
			d = direction(x,y,z,xc,yc,zc);
			d1 = direction(x1,y1,z1,xc,yc,zc);
			side *= 0.5;

			if(d != d1){
				bool=0;
				for(j=0;j<8;j++){
					(*list)[j] = node = myAlloc(S_StrucNode3);
					if(d == j){  node->flag = FACE; node->faceNode = face;  }
					if(d1 == j){ node->flag = FACE; node->faceNode = old_face;}
					if( (d1 != j) && (d != j) )  node->flag = EMPTY;
				};/* for list */
			};/*( d != d1 )*/

			if(d == d1){
				for(j=0;j<8;j++){
					(*list)[j] = node = myAlloc(S_StrucNode3);
					if( d!=j )  node->flag = EMPTY;
				};/* for list */
				node = (*list)[d];  center(&xc,&yc,&zc,side,d);
			};/*( d == d1 )*/
		};/* while  bool */
	};/* if (flag==FACE) */
	addFaceArray( face );

	return( face );
};/*addFace32*/
#undef EPS

static void remFaceRec( PStrucNode3  node, PStrucFace3  face )
	/*  NEED    tree32.xc=0.5;  tree32.yc=0.5;  tree32.zc=0.5;  tree32.side=0.5;
		 tree32.x=face->x;  tree32.y=face->y;  tree32.z=face->z;
		 remFaceRec(tree32.root,face);
		 */
{
	int          flag,i;
	PStrucNode3  sub_node;
	PStrucList3  list;

	flag = node->flag;  list = ( PStrucList3 ) node->faceNode;
	switch (flag){
		case FACE:
			if( node->faceNode == face ){
				node->flag = EMPTY;   free(face);
			}
			else  errorExit3(3,"no  coincidence  in  removing");
			break;

		case TREE:
			i = direction(tree32.x,tree32.y,tree32.z,tree32.xc,tree32.yc,tree32.zc);
			tree32.side *= 0.5;   center(&tree32.xc,&tree32.yc,&tree32.zc,tree32.side,i);
			sub_node = (*list)[i];
			remFaceRec( sub_node, face );

			tree32.empty = 0;  tree32.fill = 0;
			for(i=0;i<8;i++){
				sub_node = (*list)[i];
				if( sub_node->flag == EMPTY)  tree32.empty++;
				if( sub_node->flag == FACE ){
					tree32.fill++;
					face = sub_node->faceNode;
				};
			};/* for  list */

			if ( (tree32.empty == 7) && (tree32.fill == 1) ){
				for(i=0;i<8;i++)  free( (*list)[i] );
				free( list );
				node->flag = FACE;  node->faceNode = face;
			};
			break;
	};/* switch */

	return;
}; /* remFaceRec */


void remFace32( PStrucFace3  face )
{
	remFaceArray(face->f);
	free(face);
	return;
	tree32.xc=tree3.boxcx;
	tree32.yc=tree3.boxcy;
	tree32.zc=tree3.boxcz;
	tree32.side=tree3.boxsize;
	tree32.x = face->x;  tree32.y = face->y;  tree32.z = face->z;
	remFaceRec( tree32.root, face );

	return;
}; /* remFace32 */


double nearest32( int *vert, double x, double y, double z )
{
	int          i,j,vn=0,v[2];
	double       p=0,dist=0.0;
	PStrucFace3  face;

	/*   vicinityFaces( x, y, z, size );
	*/
	for(i=0;i<tree32.nVicinityFace;i++){
		face = tree32.vicinityFace[i].face;
		v[0] = face->v1;  v[1] = face->v2; 
		for(j=0;j<2;j++){
			/*        
						 if( v[j] == sect3.work->v1 || v[j] == sect3.work->v2 )
						 continue;
						 */
			p = distance(mesh3.vert[v[j]].x,mesh3.vert[v[j]].y,mesh3.vert[v[j]].z,x,y,z);
			if((p < dist) || ((i==0)&&(j==0))){
				dist = p;
				vn = v[j];
			}
		}/* for  3  vert  of  face */
	}/* for  i < tree32.nVicinityFace */
	vert[0] = vn;

	return( dist );
}; /* nearest32 */


static int sectCube( void )
	/*  need  Cube     : tree32.xc;  tree32.yc;
		 tree32.zc;  tree32.side;
		 Vicinity : tree32.xVicinity;  tree32.yVicinity;
		 tree32.zVicinity;  tree32.sVicinity;
		 */
{
	if( tree32.xc+tree32.side < tree32.xVicinity-tree32.sVicinity )  return(0);
	if( tree32.xc-tree32.side > tree32.xVicinity+tree32.sVicinity )  return(0);
	if( tree32.yc+tree32.side < tree32.yVicinity-tree32.sVicinity )  return(0);
	if( tree32.yc-tree32.side > tree32.yVicinity+tree32.sVicinity )  return(0);
	if( tree32.zc+tree32.side < tree32.zVicinity-tree32.sVicinity )  return(0);
	if( tree32.zc-tree32.side > tree32.zVicinity+tree32.sVicinity )  return(0);

	return(1);

};/*sectCube*/


static void vicinityFacesRec( PStrucNode3  node )
	/*  NEED    tree32.xc=0.5;tree32.yc=0.5;tree32.zc=0.5;tree32.side=0.5;
		 tree32.nVicinityFace=0;
		 tree32.sVicinity , tree32.xVicinity , tree32.yVicinity ,
		 tree32.zVicinity
		 */
{
	int flag,i;
	PStrucList3 list;
	double x,y,z,s;

	flag = node->flag;   list = ( PStrucList3 ) node->faceNode;
	switch (flag){
		case FACE:
			x=node->faceNode->x;  y=node->faceNode->y;  z=node->faceNode->z;
			s = distanceS(x,y,z,tree32.xVicinity,tree32.yVicinity,tree32.zVicinity);
			if( s < tree32.ssVicinity ){
				tree32.vicinityFace[tree32.nVicinityFace].face = node->faceNode;
				tree32.vicinityFace[tree32.nVicinityFace].dist = s;
				tree32.nVicinityFace++;
				if( tree32.nVicinityFace > tree32.maxVicinityFace )
					errorExit3(3,"tree32.nVicinityFace");
			}
			break;

		case TREE:
			tree32.side *= 0.5;
			x = tree32.xc;  y = tree32.yc;  z = tree32.zc;
			s = tree32.side; center(&tree32.xc,&tree32.yc,&tree32.zc,tree32.side,0);
			if( sectCube() )  vicinityFacesRec( (*list)[0] );
			for(i=1;i<8;i++)
			{  tree32.side = s;   tree32.xc = x;   tree32.yc = y;  tree32.zc = z;
				center(&tree32.xc,&tree32.yc,&tree32.zc,tree32.side,i);
				if( sectCube() )  vicinityFacesRec( (*list)[i] );
			}
			break;
	}; /* switch  flag */

	return;
}; /* vicinityFacesRec */

void vicinityFacesRec_bf (double x, double y, double z, double size ) {
	int i;
	double x1, y1, z1, s;
	tree32.nVicinityFace = 0;
	for (i=0; i<tree32.nFace; i++) {
		x1 = tree32.face[i]->x;
		y1 = tree32.face[i]->y;
		z1 = tree32.face[i]->z;
		s = sqrt(distanceS(x, y, z,  x1, y1, z1));
		if ( s < size + tree32.face[i]->s/2.0 ) {
			tree32.vicinityFace[tree32.nVicinityFace].face = tree32.face[i];
			tree32.vicinityFace[tree32.nVicinityFace].dist = s*s;
			tree32.nVicinityFace++;
			if (tree32.nVicinityFace>tree32.maxVicinityFace)
				errorExit3(3,"tree32.nVicinityFace");
		}
	}
	return;
}

void vicinityFaces32( double x, double y, double z, double size )
{
	vicinityFacesRec_bf(x, y, z, size);
	return;
	tree32.xc=tree3.boxcx;
	tree32.yc=tree3.boxcy;
	tree32.zc=tree3.boxcz;
	tree32.side=tree3.boxsize;
	tree32.nVicinityFace = 0;
	tree32.ssVicinity= size*size;
	tree32.sVicinity = size;
	tree32.xVicinity = x;  tree32.yVicinity = y;  tree32.zVicinity = z;
	vicinityFacesRec( tree32.root );
	/*
		qsort(tree32.vicinityFace,tree32.nVicinityFace,sizeof(StrucFace4),compare);
		*/
	return;
}; /* vicinityFaces32 */


