
#ifndef _VOLCLUSTER_H
#define _VOLCLUSTER_H

#include "mri.h"
#include "label.h"

typedef struct {
  int nmembers;
  int *col;
  int *row;
  int *slc;
  float *x;
  float *y;
  float *z;
  int maxmember;
  float maxval;
} VOLCLUSTER;

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

/*----------------------------------------------------------*/
typedef struct {
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

} CLUSTER_HIT_TABLE, CHT;

CHT *CHTalloc(int n_ithr, double ithr_lo, double ithr_hi,
	      int n_sthr, double sthr_lo, double sthr_hi);
int CHTfree(CHT **ppcht);
int CHTprint(FILE *fp, CHT *cht);
int CHTwrite(char *fname, CHT *cht);
CHT *CHTread(char *fname);
int CHTcompare(CHT *src, CHT *targ);
int CHTsetSignString(CHT *cht, char *ithr_sign);
int CHTsignId(char *ithr_sign);

#endif
