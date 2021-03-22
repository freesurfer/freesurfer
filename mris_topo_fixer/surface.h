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


#ifndef TOPOLOGY_SURFACE_H
#define TOPOLOGY_SURFACE_H

#include "globals.h"
#include "vertex.h"
#include "face.h"
#include "loop.h"

#define UNKNOWN_TYPE_OF_SURFACE 0
#define CLOSED_SURFACE 1
#define OPEN_SURFACE 2
#define PATCH 3

class PatchDisk;

class Surface
{
private:
  int _InitSurfaceConnectivity(void);
  bool _InitFaceConnectivity(void);
  int _InitFaceCoordinates(void);
  int _FindFace(int n1, int n2, int fn) const;
  int _Allocate(int nv, int nf);
  void _OverAlloc(int nextrav,int nextraf);
  double _FaceDistance(int fdst, int fsrc);

public:
  // number of vertices, edges, and faces, and Euler characteristic
  int nvertices, nedges, nfaces, euler;
  // list of vertices
  int maxvertices;
  Vertex *vertices;
  // list of faces
  int maxfaces;
  Face *faces;
  //type of surface
  int type_of_surface;


  //constructor/destructor
  Surface(void);
  Surface(int nv, int nf);
  Surface(const string s);
  ~Surface(void);

  Surface *Clone() const;
  void Expand(int nextrav,int nextraf);
  void Center();
  void scale(double scaling_factor);
  int OpenFile(const string s,int verbose=0);
  int WriteFile(const string s, int verbose = -1) const;
  int GetDefectLabels(const string s);
  int OpenCurvatureFile(const string s);
  bool IsSurfaceValid(int verbose = 0);
  void PrintDefectInfo(int ndefect=-1);
  int InitSurface(void);
  int GetEuler(int &nv, int &ne, int &nf, int mark = -1);
  int GetEuler(int mark = -1);
  int GetEuler(const int *list_of_faces, int nfs);
  Surface *ExtractPatch(int mark,int nextravertices=0, int nextrafaces=0);
  void SetMarks(int mark);
  void SetMarks(const int *v, int nv, int mark);
  void Smooth(int niters);
  void Smooth(int niters, const int *v,int nv);
  void SmoothMarked(int niters,int mark);
  void ExpandMarks(int niters,int mark);

  ///////////////////////////////////////////////////////////////
  //
  //       For The Topology Correction
  //
  ///////////////////////////////////////////////////////////////

  //when the surface is a patch extracted from another surface
  int *vtrans_to,*ftrans_to,*vtrans_from,*ftrans_from;
  Surface *surface_source;
  PatchDisk *disk;

  double GetLoopLength(Loop &loop);
  void CutLoop(Loop &loop,int very_small_patch = 0);
  bool LoopValid(Loop &loop);
  void KnitPatch(Loop &loop, PatchDisk *pdisk);
  void IncreaseEuler(int nattempts,int maxinitface = -1);
  void CorrectTopology();
  int CutPatch(int seed=-1, int maxinitface = -1, int nattempts = 10, int very_small_patch = 0);

#if 0
  int computeVertexNormals(void);
  int computeFaceNormals(void);
  int computeFaceNormal(int n);

  int Center(void);
  int saveVertices(int type);
  int restoreVertices(int src_type);
  int computeCurvature(void);
  int scale(double scaling_factor);
  double projectOntoSphere(int niter);
  int edgesIntersection(int vno1, int vno2, int vno3, int vno4);
  int computeIntersection(void);
  int smoothOnSphere(int niter);
  int smoothOnlyMarkedVertices(int niter);
  int sphericalProjection(void);
  int initSHT(SHT& sht);
  int initHT(HT &ht);
#endif
};

#endif
