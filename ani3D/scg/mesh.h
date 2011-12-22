#ifndef MESH_H
#define MESH_H

#define PI 3.14159265358979f  
#define EPS 10e-6

class mesh
{
public:
// Coordinates of the mesh:
    double x, y, z;
    
// Angles:
    double pitch, yaw, roll;
    
    int nV, nF;
    double *Vertex;
    int *Index;
    
    mesh ();
    mesh (int, int);
    ~mesh ();

    bool hide;
    int in;

    int mesh_smv (char *name, int color);
    int mesh_front (int nV_, int nF_, double *Vertex_, int *Index_);
    int mesh_write_smv (char *name);

    int primitive_paral (double x, double y, double z, double MS, int nnV, int nnF);
    int primitive_sphere (double x, double MS, int nnV, int nnF);
    int primitive_cylinder (double r, double h, double MS, int nnV, int nnF);
    int move (double x, double y, double z, double rot_x, double rot_y, double rot_z);

    int intersect (mesh *msh);
    int triInt (mesh *msh, int numF, int *ret, double pt1[3], double pt2[3]);
    int isInside (double x, double y, double z);
    int insideTri (double *u, double *v, double *x, double *y);

    void printPoint (double *x);
    double dist_point_curve (double *point, double *curve, int len, int *ret);
};

#endif // MESH_H
