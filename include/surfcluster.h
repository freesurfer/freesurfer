/**
 * @brief routines for growing clusters on the surface
 *
 * routines for growing clusters on the surface
 * based on intensity thresholds and area threshold. Note: this
 * makes use of the undefval in the MRI_SURFACE structure.
 */
/*
 * Original Author: Doug Greve
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

#ifndef _SURFCLUSTER_H
#define _SURFCLUSTER_H

#include "mri.h"
#include "mrisurf.h"
#include "label.h"

#undef SIGN
#define SIGN(x) (((x)>0)? 1.0 : -1.0 )


#ifdef SURFCLUSTER_SRC
int FixSurfClusterArea = 1;
#else
extern int FixSurfClusterArea;
#endif

/* Surface Cluster Summary */
typedef struct
{
  int   clusterno;
  int   nmembers; //  number of vertices;
  float area;
  float weightvtx;  // vertex weighted weight
  float weightarea; // area   weighted weight
  float maxval;
  int   vtxmaxval;
  float x,y,z;
  float xxfm,yxfm,zxfm;
  float cx,cy,cz; // centroid
  float cxxfm,cyxfm,czxfm; // centroid
  double pval_clusterwise; // from cluster simulation
  double pval_clusterwise_low; // from cluster simulation
  double pval_clusterwise_hi; // from cluster simulation
}
SURFCLUSTERSUM, SCS;

SCS *sclustMapSurfClusters(MRI_SURFACE *Surf, float thmin, float thmax,
                           int thsign, float minarea, int *nClusters,
                           MATRIX *XFM, MRI *fwhmmap);
int sclustGrowSurfCluster(int ClustNo, int SeedVtx, MRI_SURFACE *Surf,
                          float thmin, float thmax, int thsign);
float sclustSurfaceArea(int ClusterNo, MRI_SURFACE *Surf, int *nvtxs) ;
float sclustWeight(int ClusterNo, MRI_SURFACE *Surf, MRI *mri, int UseArea);
float sclustSurfaceMax(int ClusterNo, MRI_SURFACE *Surf, int *vtxmax) ;
int sclustSurfaceCentroid(const int ClusterNo, const MRI_SURFACE *Surf, double *xyz);
float sclustZeroSurfaceClusterNo(int ClusterNo, MRI_SURFACE *Surf);
float sclustZeroSurfaceNonClusters(MRI_SURFACE *Surf);
float sclustSetSurfaceValToClusterNo(MRI_SURFACE *Surf);
float sclustSetSurfaceValToCWP(MRI_SURFACE *Surf, SCS *scs);
float sclustCountClusters(MRI_SURFACE *Surf);
SCS *SurfClusterSummary(MRI_SURFACE *Surf, MATRIX *T, int *nClusters, MRI *fwhmmap);
SCS *SurfClusterSummaryOld(MRI_SURFACE *Surf, MATRIX *T, int *nClusters);
int DumpSurfClusterSum(FILE *fp, SCS *scs, int nClusters);
SCS *SortSurfClusterSum(SCS *scs, int nClusters);
int sclustReMap(MRI_SURFACE *Surf, int nClusters, SCS *scs_sorted);
double sclustMaxClusterArea(SURFCLUSTERSUM *scs, int nClusters);
int sclustMaxClusterCount(SURFCLUSTERSUM *scs, int nClusters);
float sclustMaxClusterWeightVtx(SURFCLUSTERSUM *scs, int nClusters, int thsign);
SCS *sclustPruneByCWPval(SCS *ClusterList, int nclusters, 
			 double cwpvalthresh,int *nPruned, 
			 MRIS *surf);
int sclustAnnot(MRIS *surf, int NClusters);
int sclustGrowByDist(MRIS *surf, int seedvtxno, double dthresh, 
		     int shape, int vtxno, int *vtxlist);
int sclustSaveAsPointSet(char *fname, SCS *scslist, int NClusters, MRIS *surf);

#endif
