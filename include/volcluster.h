
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

#endif
