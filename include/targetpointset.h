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
 * targetpointset.h - this is a class that manages target point sets
 * for surfaces. A target point set allows the user to place points where 
 * the surface should be so that the surface placement optimization is
 * encouraged to place the surface there.
 */


// fsglm.h - include file for fsglm.c

#ifndef TARGETPOINTSET_H
#define TARGETPOINTSET_H
#undef X
#include <stdio.h>
#include "mrisurf.h"
#include "mrisutils.h"
#include "dmatrix.h"
#include "matrix.h"
#include "diag.h"
#include "volcluster.h"
#include "surfcluster.h"
#include "pointset.h"
#include <vector>
#include <array>
#include "json.h"
using json = nlohmann::json;
#include "dtk.fs.h"

class SurfacePointSet {
public:
  MRIS *surf;
  MRI *mri = NULL; // template from surf
  json *pPointSet;
  std::vector<SURFHOPLIST *> m_shl;
  std::vector<int> m_shl_vtxlist;
  int m_nhops=2;
  int m_fill_holes = 1;
  std::vector<int> psvtxlist; // vertices that map from the point set
  std::vector<int> vtxlist; // all vertices
  std::vector<int> npervtxlist; // number of target points that map to each vertex
  std::vector<std::vector<double>> txyzlist; // target coords for all vertices
  std::vector<std::vector<double>> dxyz; // delta from vertex to target
  std::vector<double> angle; // angle between surface normal and vector to target
  std::vector<double> dist; // dist from vertex to its target
  double AngleDegThresh = 60; // prune points where vector relative to normal is more than this (abs)
  int m_prune_by_angle = 1;
  int m_debug = 0;
  FILE *m_debug_fp = stdout;
  int MapPointSet(void);
  fsPointSet ConvertToPointSet(void); // for debugging, per-vertex, tkreg
  fsPointSet VerticesToPointSet(void); // for debugging, per-vertex, tkreg
  int Print(FILE *fp);
  MRI *MakeMask(void);
  double CostAndGrad(double weight, int ComputeGradient);
  DTK_TRACK_SET *ConvertToTrack(int nsteps);
  int WriteAsPatch(const char *fname,int ndil);
  int PruneByAngle(void);
  ~SurfacePointSet(){
    printf("SurfacePointSet::Destructor\n");
    if(mri) MRIfree(&mri);
    if(m_debug_fp != stdout && m_debug_fp != stderr) fclose(m_debug_fp);
    for(int n=0; n < m_shl.size(); n++) SurfHopListFree(&m_shl[n]);
  }
};


#endif



