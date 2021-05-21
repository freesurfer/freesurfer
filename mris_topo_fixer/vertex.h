/**
 * @brief topology fixer fastloop worker routines
 *
 */
/*
 * Original Author: F. Segonne
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


#ifndef TOPOLOGY_VERTEX_H
#define TOPOLOGY_VERTEX_H

#include "globals.h"

class Vertex
{
public:
  double x,y,z; // location
  double tx,ty,tz; // temporary coordinates
  double xorig,yorig,zorig; // original coordinates
  double sx,sy,sz; //spherical coordinates

  double nx,ny,nz; // normal (defined by the set of normals)

  int fnum,maxfnum; // number of adjacent faces
  int *f;   // list of adjacent faces
  int *n;   // list of vertex positions in adjacent faces
  int vnum,maxvnum; // number of adjacent vertices (we must have vnum = fnum for a valid closed surface)
  int *v;   // list of adjacent vertices (represent an edge)
  int *e;   // state of the adjacent vertices or of the corresponding edges

  double curv; // the curvature of the vertex

  int marked; // for computational purposes

  //constructor/destructor
  Vertex(void);
  ~Vertex(void);
  void Clear();
  const Vertex &operator=(const Vertex &v);

  inline int AllocateFaces(int mf)
  {
    if (maxfnum < mf)
    {
      maxfnum = mf;
      if (f) delete [] f;
      if (n) delete [] n;
      f = new int[maxfnum];
      n = new int[maxfnum];
      ASSERT(f != 0 && n != 0);
    }
    fnum=0;
    return 0;
  }
  inline int ExpandFaces(int nf)
  {
    if ((unsigned int)fnum + (unsigned int)nf 
	<= (unsigned int)maxfnum) 
      return 1; // enough free faces are available

    int *f_tmp = f, *n_tmp = n;
    maxfnum = fnum + nf ;

    f = new int[maxfnum];
    n = new int[maxfnum];
    ASSERT(f != 0 && n != 0);
    if (!f || !n) return -1;

    for (int k = 0 ; k < fnum; k++)
    {
      f[k]=f_tmp[k];
      n[k]=n_tmp[k];
    }
    delete [] f_tmp;
    delete [] n_tmp;
    return 0;
  }
  inline int AllocateVertices(int mv)
  {
    if (maxvnum < mv)
    {
      maxvnum = mv;
      if (v) delete [] v;
      if (e) delete [] e;
      v = new int[maxvnum];
      e = new int[maxvnum];
      ASSERT(v != 0 && e != 0);
    }
    vnum=0;
    return 0;
  }
  inline int ExpandVertices(int nv)
  {
    if (vnum + nv <= maxvnum) return 1; // enough free vertices are available

    int *v_tmp = v, *e_tmp = e;
    maxvnum = vnum + nv ;

    v = new int[maxvnum];
    e = new int[maxvnum];
    ASSERT(v != 0 && e != 0);
    if (!v || !e) return -1;

    for (int k = 0 ; k < vnum; k++)
    {
      v[k]=v_tmp[k];
      e[k]=e_tmp[k];
    }
    delete [] v_tmp;
    delete [] e_tmp;
    return 0;
  }
  inline void AddFace(int fn, int _n)
  {
    if (fnum == maxfnum) ExpandFaces(2);
    f[fnum]=fn;
    n[fnum++]=_n;
  }
  inline void AddEdge(int nv)
  {
    if (vnum==maxvnum) ExpandVertices(1);
    v[vnum]=nv;
    e[vnum++]=0;
  }
};

#endif
