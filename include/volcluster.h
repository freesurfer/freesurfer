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



#ifndef _VOLCLUSTER_H
#define _VOLCLUSTER_H

#include "mri.h"
#include "label.h"

#undef SIGN
#define SIGN(x) (((x)>0)? 1.0 : -1.0 )

typedef struct
{
  int nmembers;
  int *col;
  int *row;
  int *slc;
  float *x;
  float *y;
  float *z;
  int maxmember;
  float maxval;
  float voxsize;
  double pval_clusterwise;
  double pval_clusterwise_low;
  double pval_clusterwise_hi;
}
VOLCLUSTER;

#ifdef VOLCLUSTER_SRC
int VolClustFixMNI = 0;
#else
extern int VolClustFixMNI;
#endif

VOLCLUSTER *clustAllocCluster(int nmembers);
int clustFreeCluster(VOLCLUSTER **ppvc);
VOLCLUSTER *clustCopyCluster(VOLCLUSTER *vc);

VOLCLUSTER **clustAllocClusterList(int nlist);
int clustFreeClusterList(VOLCLUSTER ***pppvclist, int nlist);
VOLCLUSTER **clustCopyClusterList(VOLCLUSTER **vclist, int nlist,
                                  VOLCLUSTER **vclist2);

int clustValueInRange(float val, float thmin, float thmax, int thsign);
MRI *clustInitHitMap(MRI *vol, int frame,
                     float thmin, float thmax, int thsign,
                     int *nhits, int **hitcol, int **hitrow, int **hitslc,
                     MRI *binmask, int maskframe);
int clustDumpCluster(FILE *fp, VOLCLUSTER *vc, MRI *vol, int frame);
int clustAddMember(VOLCLUSTER *vc, int col, int row, int slc);

VOLCLUSTER *clustGrow(int col0, int row0, int slc0,
                      MRI *HitMap, int AllowDiag);
int clustGrowOneVoxel(VOLCLUSTER *vc, int col0, int row0, int slc0,
                      MRI *HitMap, int AllowDiag);

int clustMaxMember(VOLCLUSTER *vc, MRI *vol, int frame, int thsign);


VOLCLUSTER **clustPruneBySize(VOLCLUSTER **vclist, int nlist,
                              float voxsize, float sizethresh,
                              int *nkeep);

VOLCLUSTER **clustPruneByDistance(VOLCLUSTER **vclist, int nlist,
                                  float distthresh, int *nkeep);

VOLCLUSTER **clustPruneByCWPval(VOLCLUSTER **vclist, int nlist,
				double cwpvalthresh, int *nkeep);

int clustCompareCluster(const void *a, const void *b);
VOLCLUSTER **clustSortClusterList(VOLCLUSTER **vclist, int nlist,
                                  VOLCLUSTER **vcsorted);

int clustComputeXYZ(VOLCLUSTER *vc, MATRIX *CRS2XYZ);
int clustComputeTal(VOLCLUSTER *vc, MATRIX *CRS2MNI);


MRI * clustClusterList2Vol(VOLCLUSTER **vclist, int nlist, MRI *tvol,
                           int frame, int ValOption);

LABEL *clustCluster2Label(VOLCLUSTER *vc, MRI *vol, int frame,
                          float colres, float rowres, float sliceres,
                          MATRIX *FSA2Func);

int clustDumpClusterList(FILE *fp, VOLCLUSTER **vclist, int nlist,
                         MRI *vol, int frame);
VOLCLUSTER **clustGetClusters(MRI *vol, int frame,
                              float threshmin, float threshmax,
                              int threshsign, float minclustsizemm3,
                              MRI *binmask, int *nClusters,
                              MATRIX *XFM);
int clustMaxClusterCount(VOLCLUSTER **VolClustList, int nClusters);
int clustDumpSummary(FILE *fp,VOLCLUSTER **VolClustList, int nClusters);

/*----------------------------------------------------------*/
typedef struct
{
  char simtype[100];  // perm, null-full, null-z
  char anattype[100]; // surface or volume
  char subject[100];  // when anattype==surf
  char hemi[10];      // when anattype==surf
  long seed;          // used for simulation
  char contrast[100]; // contrast name
  double thresh;
  double threshsign;  //0=abs,+1,-1
  double nullfwhm;    // smoothing of null simulation
  double varfwhm;     // amount of variance smoothing
  double searchspace; // in mm^2 for surf or mm^3 for vol
  int nreps;          // number of repetitions
  int *nClusters;
  double *MaxClusterSize;
  double *MaxClusterSizeVtx;
  double *MaxClusterWeightVtx;
  double *MaxClusterWeightArea;
  double *MaxSig;
  double *MaxStat;
  int mergedflag;     // Flag to indicate that two or more merged
  HISTOGRAM *mcs_pdf, *mcs_cdf; // max cluster size
  HISTOGRAM *ms_pdf, *ms_cdf;   // max sig
  double *grf_cdf; // for Gauss Rand Fields
  int FixGroupSubjectArea; // flag for keeping track
}
CLUSTER_SIM_DATA, CSD;

CLUSTER_SIM_DATA *CSDalloc(void);
int CSDallocData(CLUSTER_SIM_DATA *csd);
int CSDfreeData(CLUSTER_SIM_DATA *csd);
CSD *CSDcopy(CSD *csd, CSD *csdcopy);
CSD *CSDread(char *csdfile);
CSD *CSDmerge(CSD *csd1, CSD *csd2);
CSD *CSDreadMerge(char *csdfile, CSD *csd);
int CSDprintHeader(FILE *fp, CLUSTER_SIM_DATA *csd);
int CSDprint(FILE *fp, CSD *csd);
int CSDprintWeight(FILE *fp, CLUSTER_SIM_DATA *csd);
double CSDpvalMaxSig(double val, CSD *csd);
MRI *CSDpvalMaxSigMap(MRI *sig, CSD *csd, MRI *mask, MRI *vwsig, double *maxmaxsig, int Bonf);
double CSDpvalClustSize(CLUSTER_SIM_DATA *csd, double ClusterSize,
                        double ciPct, double *pvalLow, double *pvalHi);

int CSDcheckSimType(char *simtype);
int CSDpdf(CSD *csd, int nbins);
int CSDprintPDF(FILE *fp, CSD *csd);
int CSDwritePDF(char *fname, CSD *csd);

/*----------------------------------------------------------*/
typedef struct
{
  int nsim; /* number of simulation runs to generate table */
  long int seed; /* seed for random number generator */
  int nvox; /* number of voxels/vertices in search area */
  double totsize; /* total volume (mm^3) or area (mm^2) in search*/
  double fwhm;   /* fwhm in mm */
  int nsmooth;   /* number of smooth steps, surf only */
  double   ithr_lo, ithr_hi; /* intensity threshold range */
  int    n_ithr; /* Number ithreshs bet lo and hi*/
  double  *ithr; /* intensity thresholds*/
  char     ithr_sign[50]; /* abs, pos, neg*/
  int      ithr_signid; /* 0=abs, 1=pos, -1=neg*/
  double   sthr_lo, sthr_hi; /* cluster size threshold range */
  int    n_sthr; /* Number sthr's bet lo and hi*/
  double  *sthr; /* list of size thresholds*/
  int **hits;  /* hit[ithr][sthr] */

}
CLUSTER_HIT_TABLE, CHT;

CHT *CHTalloc(int n_ithr, double ithr_lo, double ithr_hi,
              int n_sthr, double sthr_lo, double sthr_hi);
int CHTfree(CHT **ppcht);
int CHTprint(FILE *fp, CHT *cht);
int CHTwrite(char *fname, CHT *cht);
CHT *CHTread(char *fname);
int CHTcompare(CHT *src, CHT *targ);
int CHTsetSignString(CHT *cht, char *ithr_sign);
int CHTsignId(char *ithr_sign);
int CSDwrite(char *fname, CSD *csd);

#endif
