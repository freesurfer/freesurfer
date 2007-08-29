/**
 * @file  icosahedron.h
 * @brief Icosahedron utils
 *
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2007/08/29 16:16:11 $
 *    $Revision: 1.8 $
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


#ifndef ICOSAHEDRON_H
#define ICOSAHEDRON_H

#include "mrisurf.h"

#define ICO_VERTICES_PER_FACE  3

typedef struct
{
  float x, y, z ;
}
IC_VERTEX ;

typedef struct
{
  int vno[ICO_VERTICES_PER_FACE] ;
}
IC_FACE ;

typedef struct
{
  int       nvertices ;
  int       nfaces ;
  IC_VERTEX *vertices ;
  IC_FACE   *faces ;
}
ICOSAHEDRON ;

MRI_SURFACE *ic642_make_surface(int max_vertices, int max_faces) ;
MRI_SURFACE *ic42_make_surface(int max_vertices, int max_faces) ;
MRI_SURFACE *ic40962_make_surface(int max_vertices, int max_faces) ;
MRI_SURFACE *ic10242_make_surface(int max_vertices, int max_faces) ;
MRI_SURFACE *ic2562_make_surface(int max_vertices, int max_faces) ;
MRI_SURFACE *ic162_make_surface(int max_vertices, int max_faces) ;
MRI_SURFACE *ic163842_make_surface(int max_vertices, int max_faces) ;
MRI_SURFACE *ic12_make_surface(int max_vertices, int max_faces) ;
MRI_SURFACE *ICOread(char *fname) ;
MRI_SURFACE *ICOreadOverAlloc(char *fname, double pct_over) ;
int          ICOreadVertexPositions(MRI_SURFACE *mris, 
                                    char *fname, 
                                    int which) ;
MRI_SURFACE *ReadIcoByOrder(int IcoOrder, float RescaleFactor);
MRI_SURFACE *ReadIcoByNVtxs(int nIcoVtxs, float RescaleFactor);
int          IcoOrderFromNVtxs(int nIcoVtxs);
int          IcoNVtxsFromOrder(int IcoOrder);

#define ICO4_NVERTICES    2562
#define ICO4_NFACES       5120

extern IC_VERTEX ic2562_vertices[] ;
extern IC_FACE   ic2562_faces[] ;


// version of ic2562 used for testing mrishash.c, contributed by G. Wideman
MRI_SURFACE *ic2562_make_two_icos(float x1, float y1, float z1, float r1,
                                  float x2, float y2, float z2, float r2 );

ICOSAHEDRON *read_icosahedron(char *fname) ;
ICOSAHEDRON *read_icosahedron_by_order(int order) ;


#endif

