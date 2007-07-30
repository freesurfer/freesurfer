/**
 * @file  surfcluster.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2007/07/30 23:10:29 $
 *    $Revision: 1.9 $
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



#ifndef _SURFCLUSTER_H
#define _SURFCLUSTER_H

#include "mri.h"
#include "mrisurf.h"
#include "label.h"

#undef SIGN
#define SIGN(x) (((x)>0)? 1.0 : -1.0 )

/* Surface Cluster Summary */
typedef struct
{
  int   clusterno;
  int   nmembers;
  float area;
  float maxval;
  int   vtxmaxval;
  float x,y,z;
  float xxfm,yxfm,zxfm;
  double pval_clusterwise; // from cluster simulation
  double pval_clusterwise_low; // from cluster simulation
  double pval_clusterwise_hi; // from cluster simulation
}
SURFCLUSTERSUM, SCS;

SCS *sclustMapSurfClusters(MRI_SURFACE *Surf, float thmin, float thmax,
                           int thsign, float minarea, int *nClusters,
                           MATRIX *XFM);
int sclustGrowSurfCluster(int ClustNo, int SeedVtx, MRI_SURFACE *Surf,
                          float thmin, float thmax, int thsign);
float sclustSurfaceArea(int ClusterNo, MRI_SURFACE *Surf, int *nvtxs) ;
float sclustSurfaceMax(int ClusterNo, MRI_SURFACE *Surf, int *vtxmax) ;
float sclustZeroSurfaceClusterNo(int ClusterNo, MRI_SURFACE *Surf);
float sclustZeroSurfaceNonClusters(MRI_SURFACE *Surf);
float sclustSetSurfaceValToClusterNo(MRI_SURFACE *Surf);
float sclustSetSurfaceValToCWP(MRI_SURFACE *Surf, SCS *scs);
float sclustCountClusters(MRI_SURFACE *Surf);
SCS *SurfClusterSummary(MRI_SURFACE *Surf, MATRIX *T, int *nClusters);
int DumpSurfClusterSum(FILE *fp, SCS *scs, int nClusters);
SCS *SortSurfClusterSum(SCS *scs, int nClusters);
int sclustReMap(MRI_SURFACE *Surf, int nClusters, SCS *scs_sorted);
double sclustMaxClusterArea(SURFCLUSTERSUM *scs, int nClusters);
SCS *sclustPruneByCWPval(SCS *ClusterList, int nclusters, 
			 double cwpvalthresh, int *nPruned);

#endif
