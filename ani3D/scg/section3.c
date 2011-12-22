#include"struct3.h"
#include"error3.h"
#include"section3.h"
#include"tetra3.h"
#include"tree3.h"
#ifdef DEBUG
#include"debug3.h"
#endif


/* extern  variables */
extern StrucMesh3  mesh3;
#ifdef DEBUG
extern  StrucDebug3  debug3;
#endif


/* global  variables */
StrucSect3  sect3;
int     sameVert=-1;/* same  vert */


/* static  function */
static int pointInTriangle (void);
static int faceIntersect2 (PStrucVert3 v1, PStrucVert3 v2, PStrucVert3 v3, PStrucVert3 v4);

double determ3 (double a1, double a2, double a3, double b1, double b2, double b3, double c1, double c2, double c3);


void shift (PStrucVert3 v, double dx, double dy, double dz) {
	v->x -= dx;
	v->y -= dy;
	v->z -= dz;
	return;
} /*shift*/


void rotate (PStrucVert3 v, unsigned char exix, double cos, double sin) {
	double x, y, z;
	/*  around  Z
    X  ===  x*cos + y*sin
    Y  === -x*sin + y*cos
	 */
	switch (exix) {
		case 'x':
			z = v->z;   y = v->y;
			v->y =  y*cos + z*sin;
			v->z = -y*sin + z*cos;
			break;
		case 'y':
			x = v->x;   z = v->z;
			v->x =  x*cos + z*sin;
			v->z = -x*sin + z*cos;
			break;
		case 'z':
			x = v->x;   y = v->y;
			v->x =  x*cos + y*sin;
			v->y = -x*sin + y*cos;
			break;
	}

	return;
} /*rotate*/


void makeVector (PStrucVert3 v, int v1, int v2) {
	v->x = mesh3.vert[v1].x-mesh3.vert[v2].x;
	v->y = mesh3.vert[v1].y-mesh3.vert[v2].y;
   v->z = mesh3.vert[v1].z-mesh3.vert[v2].z;
   return;
} /*makeVector*/


int makeNormal (PStrucVert3 v, PStrucVert3 v1, PStrucVert3 v2) {
   double p;

   v->x = v1->y*v2->z-v2->y*v1->z;
   v->y = v1->z*v2->x-v2->z*v1->x;
   v->z = v1->x*v2->y-v2->x*v1->y;
   p = sqrt(v->x*v->x + v->y*v->y + v->z*v->z);
   if (p != 0.0) {
      v->x /= p; v->y /= p; v->z /= p;
      return  1;
   }
	return  0;
} /*makeNormal*/





int  pointInTriangle( void )
{
   double  x1,x2,y1,y2;
   double  p,x,y;
   double  t1,t2;

   x1 = sect3.v[4].x-sect3.v[3].x;   x2 = sect3.v[5].x-sect3.v[3].x;
   y1 = sect3.v[4].y-sect3.v[3].y;   y2 = sect3.v[5].y-sect3.v[3].y;
/* I forgot  to  write  equation */
   x = sect3.v[1].x-sect3.v[3].x;
   y = sect3.v[1].y-sect3.v[3].y;
   p = x1*y2-y1*x2;
   if( p == 0.0 )
     errorExit3(2,"imposible  p==0  in  pointInTriangle");
   t1 = (y2*x-x2*y)/p;
   t2 =-(y1*x-x1*y)/p;
   if( (t1 >= 0.0) && (t2 >= 0.0) && (t1+t2 <= 1.0) )
     return  1;

   return  0;
}/*pointInTriangle*/


int  faceIntersect2( PStrucVert3 v1, PStrucVert3 v2, PStrucVert3 v3, PStrucVert3 v4 )
{
   int     i;
   double  xInt=0.;
   double  dx,dy,dz,cos,sin,p;
   double  x1,y1;

   sect3.vv[0].x = v1->x;
   sect3.vv[0].y = v1->y;
   sect3.vv[1].x = v2->x;
   sect3.vv[1].y = v2->y;
   sect3.vv[2].x = v3->x;
   sect3.vv[2].y = v3->y;
   sect3.vv[3].x = v4->x;
   sect3.vv[3].y = v4->y;
   sect3.vv[4].x = sect3.xSame;
   sect3.vv[4].y = sect3.ySame;

   dx = sect3.vv[3].x;
   dy = sect3.vv[3].y;
   dz = 0.;
   for(i=0;i<5;i++)
     shift(&sect3.vv[i],dx,dy,dz);

   p = sqrt(sect3.vv[2].x*sect3.vv[2].x+sect3.vv[2].y*sect3.vv[2].y);
   if( p != 0. ){
      cos = sect3.vv[2].x/p;   sin = sect3.vv[2].y/p;
      for(i=0;i<3;i++)
	rotate(&sect3.vv[i],'z',cos,sin);
      if( (sect3.vv[4].x != 0.) || (sect3.vv[4].y != 0.) )
	rotate(&sect3.vv[4],'z',cos,sin);
   }
   else
     errorExit3(2,"p==0  in  faceIntersect2");

/***************************/
   if( (sect3.vv[0].y > 0.0) && (sect3.vv[1].y > 0.0) )
     return  0;
   if( (sect3.vv[0].y < 0.0) && (sect3.vv[1].y < 0.0) )
     return  0;

/***************************/
   x1 = sect3.vv[1].x-sect3.vv[0].x;
   y1 = sect3.vv[1].y-sect3.vv[0].y;

/***************************/
   if( fabs(y1) < 1.e-10 ){/*  case  one  line */

      if( (sect3.vv[0].x < sect3.vv[3].x) && (sect3.vv[1].x < sect3.vv[3].x) )
	return  0;
      if( (sect3.vv[0].x > sect3.vv[2].x) && (sect3.vv[1].x > sect3.vv[2].x) )
	return  0;
      if( sameVert > -1 ){/* case  sameVert */
	 if( (sect3.vv[0].x >= sect3.vv[3].x) && (sect3.vv[0].x <= sect3.vv[2].x) )
	   ;
	 else{
	    p = sect3.vv[0].x;
	    sect3.vv[0].x = sect3.vv[1].x;
	    sect3.vv[1].x = p;
	 }
	 if( (sect3.vv[1].x >= sect3.vv[3].x) && (sect3.vv[1].x <= sect3.vv[2].x) )
	   xInt = 0.5*(sect3.vv[0].x+sect3.vv[1].x);
	 if( sect3.vv[1].x < sect3.vv[3].x )
	   xInt = 0.5*(sect3.vv[0].x+sect3.vv[3].x);
	 if( sect3.vv[1].x > sect3.vv[2].x )
	   xInt = 0.5*(sect3.vv[0].x+sect3.vv[2].x);
	 p = distance9(xInt,0.,0.,sect3.vv[4].x,sect3.vv[4].y,0.);
   if( p < M_ERROR )
	   return  0;
      }

      return  1;
   }

/***************************/

/*  equation  of  segment:
                               x1             t1 >= 0.
               sect3.vv[0]  +  y1 * t1        t1 <= 1.

                              ... :
                       sect3.vv[0].y  +  y1 * t1  ===  0.0
     t1  ===  - sect3.vv[0].y / y1

    coordinates  of  point  of  intersection  with  exixX:
	      x  ===  sect3.vv[0].x  +   x1 * t1
	      y  ===  0.
    intersection  is  when:
	    sect3.vv[3].x  <=  x  <=  sect3.vv[2].x


*/
   xInt = sect3.vv[0].x+x1*(-sect3.vv[0].y/y1);

   if( (sect3.vv[3].x <= xInt) && (xInt <= sect3.vv[2].x) ){
      if( sameVert > -1 ){/* case  sameVert */
	 p = distance9(xInt,0.,0.,sect3.vv[4].x,sect3.vv[4].y,0.);
   if( p < M_ERROR )
	   return  0;
      }

      return  1;
   }

   return  0;
}/*faceIntersect2*/


int  intersectOne( void )
{
   double  p,sx,sy;

   if( (sect3.v[1].y < 0.) && (sect3.v[2].y < 0.) )
     return  0;
   if( fabs(sect3.v[1].x-sect3.v[2].x) < M_ERROR && fabs(sect3.v[1].y-sect3.v[2].y) < M_ERROR ){
      if( pointInTriangle() ){
	 if( sameVert > -1 ){/* case  sameVert */
	    p = distance9(sect3.v[1].x,sect3.v[1].y,0.,sect3.xSame,sect3.ySame,0.);
      if( p < M_ERROR )
	      return  0;
	 }
	 return  1;
      }
      return  0;
   }

   sx = sect3.v[1].x;
   sy = sect3.v[1].y;
   sect3.v[1].x = 0.5*(sect3.v[1].x+sect3.v[2].x);
   sect3.v[1].y = 0.5*(sect3.v[1].y+sect3.v[2].y);
   if( pointInTriangle() )
     return  1;
   sect3.v[1].x = sx;
   sect3.v[1].y = sy;

   if( faceIntersect2(&sect3.v[1],&sect3.v[2],&sect3.v[4],&sect3.v[5]) )
     return  1;
   if( faceIntersect2(&sect3.v[1],&sect3.v[2],&sect3.v[4],&sect3.v[3]) )
     return  1;
   if( faceIntersect2(&sect3.v[1],&sect3.v[2],&sect3.v[5],&sect3.v[3]) )
     return  1;

   return  0;
}/*intersectOne*/

