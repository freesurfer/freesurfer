#ifndef ICOSOHEDRON_H
#define ICOSOHEDRON_H

#include "mrisurf.h"

#define ICO_VERTICES_PER_FACE  3

typedef struct
{
  float x, y, z ;
} IC_VERTEX ;

typedef struct
{
  int vno[ICO_VERTICES_PER_FACE] ;
} IC_FACE ;

typedef struct 
{
  int       nvertices ;
  int       nfaces ;
  IC_VERTEX *vertices ;
  IC_FACE   *faces ;
} ICOSOHEDRON ;

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
int         ICOreadVertexPositions(MRI_SURFACE *mris, char *fname, int which) ;

#define ICO4_NVERTICES    2562
#define ICO4_NFACES       5120

extern IC_VERTEX ic2562_vertices[] ;
extern IC_FACE   ic2562_faces[] ;

#endif
