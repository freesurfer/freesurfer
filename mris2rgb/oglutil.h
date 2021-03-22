/*
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


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
int    OGLUsetCoordParms(double coord_thickness, int num_lines) ;

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
#define VAL_FLAG      0x0080   /* paint values on surface */
#define TIME_FLAG     0x0100   /* paint latencies on surface */
#define STAN_FLAG     0x0200   /* paint spatio-temporal meg loc results */
#define LIN_FLAG      0x0400   /* paint linear analysis results */
#define CSCALE_FLAG   0x0800   /* draw the color scale */
#define LEG_FLAG      0x1000   /* Write a legend (min and max near the color scale) */

#define LIGHT_OFFSET  0.25

#define MARK_WHITE          1
#define MARK_YELLOW         2
#define MARK_BLUE           3
#define MARK_PURPLE         4
#define MARK_LIGHTGREEN     5
#define MARK_ORANGE         6
#define MARK_RED            7
#define MARK_GREEN          8
#define MARK_CYAN           9

#endif
