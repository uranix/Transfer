#include "region3.h"
#include "error3.h"
#include "memory3.h"
#include "tree3.h"
#ifdef DEBUG
#include "debug3.h"
#endif

/* extern  variables */
extern StrucMesh3 mesh3;
#ifdef DEBUG
extern StrucDebug3 debug3;
#endif

/* global  variables */
StrucTree3	tree3;

/* static  function */
static void addFaceArray(PStrucFace3  face);
static void remFaceArray(int  face);
static void remFaceRec(PStrucNode3  node, PStrucFace3  face);
static int sectCube(void);
static void vicinityFacesRec(PStrucNode3  node);

int direction(double x,double y,double z, double xc,double yc,double zc) {
	int dx=0, dy=0, dz=0;
	if (x>xc) dx = 1;
	if (y>yc) dy = 2;
	if (z>zc) dz = 4;
	return dx+dy+dz;
}

void center(double *x,double *y,double *z, double side,int d) {
	int dx=0, dy=0, dz=0;
	if ((d<0) || (d>7)) errorExit3(2, "direction");
	dz = d/4; d = d%4;
	dy = d/2;
	dx = d%2;
	if (dx) x[0]+=side; else x[0]-=side;
	if (dy) y[0]+=side; else y[0]-=side;
	if (dz) z[0]+=side; else z[0]-=side;
}

double distance(double x,double y,double z, double xc,double yc,double zc) {
	return sqrt( (x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc) );
} /*distance*/

double distanceS(double x,double y,double z, double xc,double yc,double zc) {
	return (x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc);
} /*distanceS*/

double distance9(double x,double y,double z, double xc,double yc,double zc) {
	double p1, p2, p3;
	p1 = fabs(x-xc);
	p2 = fabs(y-yc);
	p3 = fabs(z-zc);
	if (p1<p2) p1 = p2;
	if (p1<p3) p1 = p3;
	return  p1;
}/*distance9*/

static void addFaceArray(PStrucFace3  face) {
	int son, fath, i;
	PStrucFace3 cface, *tface=tree3.face;
	if (tree3.nFace>=tree3.maxFace) errorExit3(1, "tree3.nFace");
	tface[tree3.nFace] = face;
	face->f = tree3.nFace;
	son=++tree3.nFace;
aaa:fath=son/2;
#ifdef DEBUG
	 if (mesh3.debug && debug3.bool.sortZ) {
		 if ((fath>0) && (tface[son-1]->z<tface[fath-1]->z)) {
			 cface=tface[son-1];
			 tface[son-1]=tface[fath-1];
			 tface[fath-1]=cface;
			 i=cface->f;
			 cface->f=tface[son-1]->f;
			 tface[son-1]->f=i;
			 son=fath;
			 if (son>1) goto aaa;
		 }
	 } else {
		 if ((fath>0) && (tface[son-1]->s<tface[fath-1]->s)) {
			 cface=tface[son-1];
			 tface[son-1]=tface[fath-1];
			 tface[fath-1]=cface;
			 i=cface->f;
			 cface->f=tface[son-1]->f;
			 tface[son-1]->f=i;
			 son=fath;
			 if (son>1) goto aaa;
		 }
	 }
#else
	 if ((fath>0) && (tface[son-1]->s<tface[fath-1]->s)) {
		 cface=tface[son-1];
		 tface[son-1]=tface[fath-1];
		 tface[fath-1]=cface;
		 i=cface->f;
		 cface->f=tface[son-1]->f;
		 tface[son-1]->f=i;
		 son=fath;
		 if (son>1) goto aaa;
	 }
#endif
}/*addFaceArray*/

static void remFaceArray(int  fath) {
	int son1, son2, i;
	PStrucFace3 face, *tface;
	tface = tree3.face;

	if ((fath>=tree3.nFace) || (fath<0)) errorExit3(2," tree3.nFace < 0   in   rem_face ");

	tree3.nFace--;
	tface[fath]=tface[tree3.nFace];
	tface[fath]->f=fath;
aaa:i = 0;
#ifdef DEBUG
	 if (mesh3.debug && debug3.bool.sortZ) {
		 if (2*fath+1<(tree3.nFace-1)) {
			 son1=2*fath+1;
			 son2=son1++;
			 if((tface[fath]->z<=tface[son1]->z)&&(tface[fath]->z<=tface[son2]->z)) i=0;
			 if((tface[fath]->z>=tface[son1]->z)&&(tface[son1]->z>=tface[son2]->z)) i=son2;
			 if((tface[fath]->z>=tface[son2]->z)&&(tface[son2]->z>tface[son1]->z)) i=son1;
			 if((tface[son1]->z>=tface[fath]->z)&&(tface[fath]->z>tface[son2]->z)) i=son2;
			 if((tface[son2]->z>=tface[fath]->z)&&(tface[fath]->z>tface[son1]->z)) i=son1;
		 }
		 if (2*fath+1==(tree3.nFace-1)) {
			 son1=2*fath+1;
			 if((tface[fath]->z>tface[son1]->z)) i=son1;
			 else i=0;
		 }
	 } else {
		 if (2*fath+1<(tree3.nFace-1)) {
			 son1=2*fath+1;
			 son2=son1++;
			 if((tface[fath]->s<=tface[son1]->s)&&(tface[fath]->s<=tface[son2]->s)) i=0;
			 if((tface[fath]->s>=tface[son1]->s)&&(tface[son1]->s>=tface[son2]->s)) i=son2;
			 if((tface[fath]->s>=tface[son2]->s)&&(tface[son2]->s>tface[son1]->s)) i=son1;
			 if((tface[son1]->s>=tface[fath]->s)&&(tface[fath]->s>tface[son2]->s)) i=son2;
			 if((tface[son2]->s>=tface[fath]->s)&&(tface[fath]->s>tface[son1]->s)) i=son1;
		 }
		 if (2*fath+1==(tree3.nFace-1)) {
			 son1=2*fath+1;
			 if((tface[fath]->s>tface[son1]->s)) i=son1;
			 else i=0;
		 }
	 }
#else
	 if (2*fath+1<(tree3.nFace-1)) {
		 son1=2*fath+1;
		 son2=son1++;
		 if((tface[fath]->s<=tface[son1]->s)&&(tface[fath]->s<=tface[son2]->s)) i=0;
		 if((tface[fath]->s>=tface[son1]->s)&&(tface[son1]->s>=tface[son2]->s)) i=son2;
		 if((tface[fath]->s>=tface[son2]->s)&&(tface[son2]->s>tface[son1]->s)) i=son1;
		 if((tface[son1]->s>=tface[fath]->s)&&(tface[fath]->s>tface[son2]->s)) i=son2;
		 if((tface[son2]->s>=tface[fath]->s)&&(tface[fath]->s>tface[son1]->s)) i=son1;
	 }
	 if (2*fath+1==(tree3.nFace-1)) {
		 son1=2*fath+1;
		 if((tface[fath]->s>tface[son1]->s)) i=son1;
		 else i=0;
	 }
#endif

	 if (i) {
		 face=tface[i];
		 tface[i]=tface[fath];
		 tface[fath]=face;
		 son1=face->f;
		 face->f=tface[i]->f;
		 tface[i]->f=son1;
		 fath=i;
		 goto aaa;
	 }
} /*remFaceArray*/


PStrucFace3 findFace(int v1, int v2, int v3) {
	int i, w1, w2, w3;
	for (i=0; i<tree3.nFace; i++) {
		w1 = tree3.face[i]->v1;
		w2 = tree3.face[i]->v2;
		w3 = tree3.face[i]->v3;
		if (  ((w3==v1) && (w2==v2) && (w1==v3)) ||
				((w3==v2) && (w2==v3) && (w1==v1)) ||
				((w3==v3) && (w2==v1) && (w1==v2))  )
			return tree3.face[i];
	}
	return 0;
}

#define EPS 1e-8
PStrucFace3	addFace(int v1, int v2, int v3, int twin, int color )
{
	int	flag, j, d, d1, bool;
	double	x, y, z, x1, y1, z1;
	double	xc=tree3.boxcx, yc=tree3.boxcy, zc=tree3.boxcz, side=tree3.boxsize;
	PStrucFace3	old_face, face;
	PStrucNode3	node;
	PStrucList3	list;

	node=tree3.root;
	flag=node->flag;
	list = ( PStrucList3 ) node->faceNode;
	face = myAlloc( S_StrucFace3 );
	face->color = color;
	face->v1 = v1;
	face->v2 = v2;
	face->v3 = v3;
	if (twin==0) {
		x = (mesh3.vert[v1].x+mesh3.vert[v2].x+mesh3.vert[v3].x)/3.;
		y = (mesh3.vert[v1].y+mesh3.vert[v2].y+mesh3.vert[v3].y)/3.;
		z = (mesh3.vert[v1].z+mesh3.vert[v2].z+mesh3.vert[v3].z)/3.;
		x = 0.27*mesh3.vert[v1].x+0.3*mesh3.vert[v2].x+0.43*mesh3.vert[v3].x;
		y = 0.27*mesh3.vert[v1].y+0.3*mesh3.vert[v2].y+0.43*mesh3.vert[v3].y;
		z = 0.27*mesh3.vert[v1].z+0.3*mesh3.vert[v2].z+0.43*mesh3.vert[v3].z;
	} else {
		x = 0.25*mesh3.vert[v1].x+0.3*mesh3.vert[v2].x+0.45*mesh3.vert[v3].x;
		y = 0.25*mesh3.vert[v1].y+0.3*mesh3.vert[v2].y+0.45*mesh3.vert[v3].y;
		z = 0.25*mesh3.vert[v1].z+0.3*mesh3.vert[v2].z+0.45*mesh3.vert[v3].z;
	}
	face->x=x;
	face->y=y;
	face->z=z;
	face->s = distance(mesh3.vert[v1].x,mesh3.vert[v1].y,mesh3.vert[v1].z,mesh3.vert[v2].x,mesh3.vert[v2].y,mesh3.vert[v2].z)
		+ distance(mesh3.vert[v1].x,mesh3.vert[v1].y,mesh3.vert[v1].z,mesh3.vert[v3].x,mesh3.vert[v3].y,mesh3.vert[v3].z)
		+ distance(mesh3.vert[v3].x,mesh3.vert[v3].y,mesh3.vert[v3].z,mesh3.vert[v2].x,mesh3.vert[v2].y,mesh3.vert[v2].z);
#ifdef  BADSIT
	if (twin==33) {
		x = (mesh3.vert[v1].x+mesh3.vert[v2].x+mesh3.vert[v3].x)/3.;
		y = (mesh3.vert[v1].y+mesh3.vert[v2].y+mesh3.vert[v3].y)/3.;
		z = (mesh3.vert[v1].z+mesh3.vert[v2].z+mesh3.vert[v3].z)/3.;
		face->x=x;
		face->y=y;
		face->z=z;
		face->s += 100.;  
	}
#endif

	if ((fabs(x-tree3.boxcx)>tree3.boxsize+EPS) || (fabs(y-tree3.boxcy)>tree3.boxsize+EPS) || (fabs(z-tree3.boxcz)>tree3.boxsize+EPS)) {
		printf("\nbox: x: %le, y: %le, z: %le, size: %le\n", tree3.boxcx, tree3.boxcy, tree3.boxcz, tree3.boxsize);
		printf("x: %le, y: %le, z: %le\n", x, y, z);
		errorExit3(3,"x_y_z_insert");
	}

	while (flag==TREE) {
		d = direction(x,y,z,xc,yc,zc);
		side *= 0.5;
		center(&xc,&yc,&zc,side,d);
		node = (*list)[d];
		flag = node->flag;
		list = ( PStrucList3 ) node->faceNode;
	}
	if (flag==EMPTY) {
		node->flag=FACE;
		node->faceNode=face;
	}
	if (flag==FACE) {
		bool = 1;
		old_face = node->faceNode;
		x1 = old_face->x;
		y1 = old_face->y;
		z1 = old_face->z;
		if ((x1==x) && (y1==y) && (z1==z)) {
		    if ( ((v1==old_face->v1) || (v1==old_face->v2) || (v1==old_face->v3)) &&
			 ((v2==old_face->v1) || (v2==old_face->v2) || (v2==old_face->v3)) &&
			 ((v3==old_face->v1) || (v3==old_face->v2) || (v3==old_face->v3)) )
			errorExit3(3," coincident   in  insert tree3");
		    else {
		    }
		}

		while (bool) {
			node->flag = TREE;
			node->faceNode = myAlloc(S_StrucList3);
			list = ( PStrucList3 ) node->faceNode;
			d = direction(x,y,z,xc,yc,zc);
			d1 = direction(x1,y1,z1,xc,yc,zc);
			side *= 0.5;

			if (d!=d1) {
				bool=0;
				for (j=0; j<8; j++) {
					(*list)[j] = node = myAlloc(S_StrucNode3);
					if (d==j) {node->flag = FACE; node->faceNode = face;}
					if (d1==j) {node->flag = FACE; node->faceNode = old_face;}
					if ((d1!=j)&&(d!=j)) node->flag = EMPTY;
				};/* for list */
			};/*( d != d1 )*/

			if (d==d1) {
				for (j=0; j<8; j++) {
					(*list)[j] = node = myAlloc(S_StrucNode3);
					if (d!=j) node->flag = EMPTY;
				};/* for list */
				node = (*list)[d];
				center(&xc,&yc,&zc,side,d);
			};/*( d == d1 )*/
		};/* while  bool */
	};/* if (flag==FACE) */
	addFaceArray(face);

	return(face);
} /*addFace*/
#undef EPS

static	void	remFaceRec( PStrucNode3  node, PStrucFace3  face )
	/*  NEED    tree3.xc=0.5;  tree3.yc=0.5;  tree3.zc=0.5;  tree3.side=0.5;
	 *  tree3.x=face->x;  tree3.y=face->y;  tree3.z=face->z;
	 *  remFaceRec(tree3.root,face);
	 */
{
	int	flag,i;
	PStrucNode3	sub_node;
	PStrucList3	list;

	flag = node->flag;
	list = ( PStrucList3 ) node->faceNode;
	switch (flag) {
		case FACE:
			if (node->faceNode==face) {
				node->flag = EMPTY;
				free(face);
			} else 
				errorExit3(3,"no  coincidence  in  removing");
			break;
		case TREE:
			i = direction(tree3.x,tree3.y,tree3.z,tree3.xc,tree3.yc,tree3.zc);
			tree3.side *= 0.5;
			center(&tree3.xc,&tree3.yc,&tree3.zc,tree3.side,i);
			sub_node = (*list)[i];
			remFaceRec( sub_node, face );

			tree3.empty = 0;
			tree3.fill = 0;
			for (i=0; i<8; i++) {
				sub_node = (*list)[i];
				if (sub_node->flag == EMPTY) tree3.empty++;
				if (sub_node->flag == FACE) {
					tree3.fill++;
					face = sub_node->faceNode;
				};
			};/* for  list */

			if ((tree3.empty==7)&&(tree3.fill==1)) {
				for (i=0; i<8; i++) free( (*list)[i] );
				free( list );
				node->flag = FACE;
				node->faceNode = face;
			};
			break;
	};/* switch */

	return;
} /* remFaceRec */

void	remFace( PStrucFace3  face )
{
	remFaceArray(face->f);
	tree3.xc=tree3.boxcx;
	tree3.yc=tree3.boxcy;
	tree3.zc=tree3.boxcz;
	tree3.side=tree3.boxsize;
	tree3.x = face->x;
	tree3.y = face->y;
	tree3.z = face->z;
	remFaceRec( tree3.root, face );

	return;
} /* remFace */

double	nearest3( int *vert, double x, double y, double z )
{
	int	i, j, vn=0, v[3];
	double	p=0, dist=0.0;
	PStrucFace3  face;
	/*   vicinityFaces( x, y, z, size );
	*/

	for (i=0; i<tree3.nVicinityFace; i++) {
		face = tree3.vicinityFace[i].face;
		v[0] = face->v1;
		v[1] = face->v2;
		v[2] = face->v3;
		for (j=0; j<3; j++) {
			p = distance(mesh3.vert[v[j]].x,mesh3.vert[v[j]].y,mesh3.vert[v[j]].z,x,y,z);
			if ((p<dist) || ((i==0)&&(j==0))) {
				dist = p;
				vn = v[j];
			}
		}/* for  3  vert  of  face */
	}/* for  i < tree3.nVicinityFace */
	vert[0] = vn;

	return( dist );
} /* nearest3 */

static	int	sectCube( void )
	/*  need  Cube     : tree3.xc;  tree3.yc;
	 *  tree3.zc;  tree3.side;
	 *  Vicinity : tree3.xVicinity;  tree3.yVicinity;
	 *  tree3.zVicinity;  tree3.sVicinity;
	 */
{
	if( tree3.xc+tree3.side < tree3.xVicinity-tree3.sVicinity )  return(0);
	if( tree3.xc-tree3.side > tree3.xVicinity+tree3.sVicinity )  return(0);
	if( tree3.yc+tree3.side < tree3.yVicinity-tree3.sVicinity )  return(0);
	if( tree3.yc-tree3.side > tree3.yVicinity+tree3.sVicinity )  return(0);
	if( tree3.zc+tree3.side < tree3.zVicinity-tree3.sVicinity )  return(0);
	if( tree3.zc-tree3.side > tree3.zVicinity+tree3.sVicinity )  return(0);

	return(1);
};/*sectCube*/

static	void	vicinityFacesRec( PStrucNode3  node )
	/*  NEED    tree3.xc=0.5;tree3.yc=0.5;tree3.zc=0.5;tree3.side=0.5;
	 *  tree3.nVicinityFace=0;
	 *  tree3.sVicinity , tree3.xVicinity , tree3.yVicinity ,
	 *  tree3.zVicinity
	 */
{
	int flag, i;
	PStrucList3 list;
	double x, y, z, s;

	flag = node->flag;
	list = ( PStrucList3 ) node->faceNode;
	switch (flag) {
		case FACE:
			x=node->faceNode->x;
			y=node->faceNode->y;
			z=node->faceNode->z;
			s = distanceS(x,y,z,tree3.xVicinity,tree3.yVicinity,tree3.zVicinity);
			if (s<tree3.ssVicinity) {
				tree3.vicinityFace[tree3.nVicinityFace].face = node->faceNode;
				tree3.vicinityFace[tree3.nVicinityFace].dist = s;
				tree3.nVicinityFace++;
				if (tree3.nVicinityFace>tree3.maxVicinityFace)
					errorExit3(3,"tree3.nVicinityFace");
			}
			break;
		case TREE:
			tree3.side *= 0.5;
			x = tree3.xc;
			y = tree3.yc;
			z = tree3.zc;
			s = tree3.side;
			center(&tree3.xc,&tree3.yc,&tree3.zc,tree3.side,0);
			if (sectCube()) vicinityFacesRec( (*list)[0] );
			for (i=1; i<8; i++) {
				tree3.side = s;
				tree3.xc = x;
				tree3.yc = y;
				tree3.zc = z;
				center(&tree3.xc,&tree3.yc,&tree3.zc,tree3.side,i);
				if (sectCube()) vicinityFacesRec( (*list)[i] );
			}
			break;
	}; /* switch  flag */

	return;
} /* vicinityFacesRec */

static	void	vicinityFacesRecFAF( PStrucNode3  node )
	/*  NEED    tree3.xc=0.5;tree3.yc=0.5;tree3.zc=0.5;tree3.side=0.5;
	 *  tree3.nVicinityFace=0;
	 *  tree3.sVicinity , tree3.xVicinity , tree3.yVicinity ,
	 *  tree3.zVicinity
	 */
{
	int	flag, i;
	PStrucList3	list;
	double	x, y, z, s;

	flag = node->flag;
	list = ( PStrucList3 ) node->faceNode;
	switch (flag) {
		case FACE:
			x=node->faceNode->x;
			y=node->faceNode->y;
			z=node->faceNode->z;
			s = distanceS(x,y,z,tree3.xVicinity,tree3.yVicinity,tree3.zVicinity);
			if (( s < tree3.ssVicinity) || (s-1.*distance(x,y,z,
							mesh3.vert[node->faceNode->v1].x,
							mesh3.vert[node->faceNode->v1].y,
							mesh3.vert[node->faceNode->v1].z)<
						tree3.ssVicinity/1.)){
				tree3.vicinityFace[tree3.nVicinityFace].face = node->faceNode;
				tree3.vicinityFace[tree3.nVicinityFace].dist = s;
				tree3.nVicinityFace++;
				if (tree3.nVicinityFace>tree3.maxVicinityFace)
					errorExit3(3,"tree3.nVicinityFace");
			}
			break;
		case TREE:
			tree3.side *= 0.5;
			x = tree3.xc;
			y = tree3.yc;
			z = tree3.zc;
			s = tree3.side;
			center(&tree3.xc,&tree3.yc,&tree3.zc,tree3.side,0);
			if (sectCube()) vicinityFacesRecFAF( (*list)[0] );
			for (i=1; i<8; i++) {
				tree3.side = s;
				tree3.xc = x;
				tree3.yc = y;
				tree3.zc = z;
				center(&tree3.xc,&tree3.yc,&tree3.zc,tree3.side,i);
				if (sectCube()) vicinityFacesRecFAF( (*list)[i] );
			}
			break;
	}; /* switch  flag */

	return;
} /* vicinityFacesRec */


int	compare( const void *p1, const void *p2 )
{
	StrucFace4	*f1, *f2;

	f1 = (StrucFace4 *)p1;
	f2 = (StrucFace4 *)p2;
	if (f1->dist<f2->dist)
		return  -1;
	else
		return  1;
}/*compare*/

void vicinityFacesRecFAF_bf (double x, double y, double z, double size ) {
	int i;
	double x1, y1, z1, s, ss;
	ss=size*size;
	tree3.nVicinityFace = 0;
	for (i=0; i<tree3.nFace; i++) {
		x1 = tree3.face[i]->x;
		y1 = tree3.face[i]->y;
		z1 = tree3.face[i]->z;
		s = distanceS(x, y, z,  x1, y1, z1);
		if ( s < ss ) {
			tree3.vicinityFace[tree3.nVicinityFace].face = tree3.face[i];
			tree3.vicinityFace[tree3.nVicinityFace].dist = s;
			tree3.nVicinityFace++;
			if (tree3.nVicinityFace>tree3.maxVicinityFace)
				errorExit3(3,"tree3.nVicinityFace");
		}
	}
	return;
}

void	vicinityFaces( double x, double y, double z, double size )
{
	vicinityFacesRecFAF_bf(x, y, z, size);
	return;
	tree3.xc=tree3.boxcx;
	tree3.yc=tree3.boxcy;
	tree3.zc=tree3.boxcz;
	tree3.side=tree3.boxsize;
	tree3.nVicinityFace = 0;
	tree3.ssVicinity= size*size;
	tree3.sVicinity = size;
	tree3.xVicinity = x;
	tree3.yVicinity = y;
	tree3.zVicinity = z;
	if (mesh3.surf)
		vicinityFacesRecFAF( tree3.root );
	else
		vicinityFacesRec( tree3.root );
	//	qsort(tree3.vicinityFace,tree3.nVicinityFace,sizeof(StrucFace4),compare);

	return;
}/* vicinityFaces */


static	void	findMaxTetraRec( void )
	/*  NEED    tree3.xc=0.5;tree3.yc=0.5;tree3.zc=0.5;tree3.side=0.5;
	 *  tree3.sVicinity = 0.
	 */
{
	int	i;
	double	x, y, z, s;

	if (mesh3.surf) s = mesh3.sminsize;
	else s = sizeFace(tree3.xc,tree3.yc,tree3.zc);

	if (s>5*tree3.side) {
		tree3.sVicinity += ((8.*tree3.side*tree3.side*tree3.side)/(0.117*s*s*s));
	} else {
		tree3.side *= 0.5;
		x = tree3.xc;
		y = tree3.yc;
		z = tree3.zc;
		s = tree3.side; 
		for (i=0; i<8; i++) {
			tree3.side = s;
			tree3.xc = x;
			tree3.yc = y;
			tree3.zc = z;
			center(&tree3.xc,&tree3.yc,&tree3.zc,tree3.side,i);
			findMaxTetraRec();
		}
	}         

	return;
}/*findMaxTetraRec*/


int	findMaxTetra( void )
{
	int	n;

	tree3.xc=tree3.boxcx;
	tree3.yc=tree3.boxcy;
	tree3.zc=tree3.boxcz;
	tree3.side=tree3.boxsize;
	tree3.sVicinity = 0;
	findMaxTetraRec();

	n = (int)tree3.sVicinity;

	return  n; 
}/*findMaxTetra*/

