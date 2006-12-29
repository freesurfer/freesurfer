/**
 * @file  vertex.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:19 $
 *    $Revision: 1.2 $
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


#include "topology/vertex.h"

Vertex::Vertex(void) {
  fnum=maxfnum=0;
  f=n=0;
  vnum=maxvnum=0;
  v=e=0;
  marked = 0;
}

Vertex::~Vertex(void) {
  Clear();
}

void Vertex::Clear() {
  if (f) delete [] f;
  if (n) delete [] n;
  if (v) delete [] v;
  if (e) delete [] e;
  f=n=v=e=0;
  fnum = vnum=0;
  maxfnum = maxvnum=0;
}

const Vertex& Vertex::operator=(const Vertex &vertex) {
  // copying faces
  AllocateFaces(vertex.fnum);
  fnum=vertex.fnum;
  for (int p = 0 ; p < fnum ; p++) {
    f[p] = vertex.f[p];
    n[p] = vertex.n[p];
  }

  // copying vertices
  AllocateVertices(vertex.vnum);
  vnum=vertex.vnum;
  for (int p = 0 ; p < vnum ; p++) {
    v[p] = vertex.v[p];
    e[p] = vertex.e[p];
  }

  // copying other variables
  x=vertex.x;
  y=vertex.y;
  z=vertex.z;
  tx=vertex.tx;
  ty=vertex.ty;
  tz=vertex.tz;
  xorig=vertex.xorig;
  yorig=vertex.yorig;
  zorig=vertex.zorig;
  sx=vertex.sx;
  sy=vertex.sy;
  sz=vertex.sz;
  nx=vertex.nx;
  ny=vertex.ny;
  nz=vertex.nz;
  curv=vertex.curv;
  marked=vertex.marked;

  return vertex;
}










