#ifndef TRITRI_H
#define TRITRI_H

int tri_tri_intersect(double V0[3],double V1[3],double V2[3],
                      double U0[3],double U1[3],double U2[3]) ;


int triangle_ray_intersect(double orig_pt[3], double dir[3], double U0[3], 
                           double U1[3],double U2[3], double int_pt[3]) ;

/* some macros */
#define CROSS(dest,v1,v2)                      \
              dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
              dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
              dest[2]=v1[0]*v2[1]-v1[1]*v2[0];

#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
#define VLEN(v)    (sqrt(DOT(v,v)))
#define VZERO(v)   (FZERO(DOT(v,v)))
#define SUB(dest,v1,v2)          \
            dest[0]=v1[0]-v2[0]; \
            dest[1]=v1[1]-v2[1]; \
            dest[2]=v1[2]-v2[2]; 

#define ADD(dest,v1,v2)          \
            dest[0]=v1[0]+v2[0]; \
            dest[1]=v1[1]+v2[1]; \
            dest[2]=v1[2]+v2[2]; 

#define SCALAR_MUL(dest, s, v) \
            dest[0]=v[0]*s; \
            dest[1]=v[1]*s; \
            dest[2]=v[2]*s; 

#endif
