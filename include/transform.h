#ifndef TRANSFORM_H
#define TRANSFORM_H

#include "matrix.h"

typedef struct
{
  float      x0 ;            /* center of transform */
  float      y0 ;
  float      z0 ;
  float      sigma ;         /* spread of transform */
  MATRIX     *m_L ;          /* transform matrix */
  MATRIX     *m_dL ;         /* gradient of fuctional wrt transform matrix */
  MATRIX     *m_last_dL ;    /* last time step for momentum */
} LINEAR_TRANSFORM, LT ;

typedef struct
{
  int               num_xforms ;      /* number linear transforms */
  LINEAR_TRANSFORM  *xforms ;   /* transforms */
} LINEAR_TRANSFORM_ARRAY, LTA ;

typedef struct
{
  int        type ;
  void       *xform ;
} TRANSFORM ;

int      LTAfree(LTA **plta) ;
LTA      *LTAread(char *fname) ;
int      LTAwrite(LTA *lta, char *fname) ;
LTA      *LTAalloc(int nxforms, MRI *mri) ;
int      LTAdivide(LTA *lta, MRI *mri) ;
MRI      *LTAtransform(MRI *mri_src, MRI *mri_dst, LTA *lta) ;
MATRIX   *LTAtransformAtPoint(LTA *lta, float x, float y,float z,MATRIX *m_L);
int      LTAworldToWorld(LTA *lta, float x, float y, float z,
                         float *px, float *py, float *pz);
int      LTAinverseWorldToWorld(LTA *lta, float x, float y, float z,
                         float *px, float *py, float *pz);
VECTOR   *LTAtransformPoint(LTA *lta, VECTOR *v_X, VECTOR *v_Y) ;
VECTOR   *LTAinverseTransformPoint(LTA *lta, VECTOR *v_X, VECTOR *v_Y) ;
double   LTAtransformPointAndGetWtotal(LTA *lta, VECTOR *v_X, VECTOR *v_Y) ;
MATRIX   *LTAinverseTransformAtPoint(LTA *lta, float x, float y, float z, 
                                     MATRIX *m_L) ;
MATRIX   *LTAworldTransformAtPoint(LTA *lta, float x, float y,float z,
                                   MATRIX *m_L);
int      LTAtoVoxelCoords(LTA *lta, MRI *mri) ;

#define TRANSFORM_ARRAY_TYPE    0
#define MORPH_3D_TYPE           1
#define MNI_TRANSFORM_TYPE      2

int      TransformFileNameType(char *fname) ;


#endif
