
#ifndef _SURFCLUSTER_H
#define _SURFCLUSTER_H

#include "mri.h"
#include "mrisurf.h"
#include "label.h"

#undef SIGN
#define SIGN(x) (((x)>0)? 1.0 : -1.0 )

/* Surface Cluster Summary */
typedef struct {
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
} SURFCLUSTERSUM, SCS;

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

#endif
