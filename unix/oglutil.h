#ifndef OGLUTIL_H
#define OGLUTIL_H

#include "mrisurf.h"

int    OGLUinit(MRI_SURFACE *mris, long frame_xsize, long frame_ysize) ;
int    OGLUcompile(MRI_SURFACE *mris, int *marked_vertices, int flags, 
                   float cslope) ;
void   OGLUsetLightingModel(float lite0, float lite1, float lite2, 
                            float lite3, float newoffset) ;
int    OGLUnoscale(void) ;
int    OGLUsetFOV(int fov) ;
int    OGLUsetCoordParms(double coord_thickness, double coord_spacing) ;

#define SCALE_FACTOR       0.55f

extern double oglu_fov ;

#define PATCH_FLAG   0x0001
#define TP_FLAG      0x0002   /* show tangent plane and principal directions */
                              /*  of marked vertices */
#define MESH_FLAG     0x0004
#define COORD_FLAG    0x0008   /* draw canonical coordinate system */
#define BW_FLAG       0x0010   /* only use black and white */
#define NEG_FLAG      0x0020   /* show negative vertices (flat maps only) */
#define NOBORDER_FLAG 0x0040   /* don't draw a border */

#define LIGHT_OFFSET  0.25

#endif
