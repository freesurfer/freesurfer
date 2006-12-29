/**
 * @file  tritri.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:00 $
 *    $Revision: 1.6 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


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
