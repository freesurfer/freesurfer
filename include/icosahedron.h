/**
 * @file  icosahedron.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:08:59 $
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


#ifndef ICOSOHEDRON_H
#define ICOSOHEDRON_H

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
ICOSOHEDRON ;

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
int          ICOreadVertexPositions(MRI_SURFACE *mris, char *fname, int which) ;
MRI_SURFACE *ReadIcoByOrder(int IcoOrder, float RescaleFactor);
MRI_SURFACE *ReadIcoByNVtxs(int nIcoVtxs, float RescaleFactor);
int          IcoOrderFromNVtxs(int nIcoVtxs);
int          IcoNVtxsFromOrder(int IcoOrder);

#define ICO4_NVERTICES    2562
#define ICO4_NFACES       5120

extern IC_VERTEX ic2562_vertices[] ;
extern IC_FACE   ic2562_faces[] ;

#endif

