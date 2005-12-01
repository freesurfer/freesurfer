#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_randist.h>
#include "mri.h"
#include "resample.h"
#include "transform.h"
#include "matrix.h"

#include "volcluster.h"

static int ConvertCRS2XYZ(int col, int row, int slc, MATRIX *CRS2XYZ,
        float *x, float *y, float *z);

/*----------------------------------------------------------------*/
VOLCLUSTER *clustAllocCluster(int nmembers)
{
  VOLCLUSTER *vc;

  vc = (VOLCLUSTER *) calloc(1, sizeof(VOLCLUSTER));

  if(nmembers == 0) return(vc);

  vc->nmembers = nmembers ;
  vc->col = (int *) calloc(nmembers, sizeof(int));
  vc->row = (int *) calloc(nmembers, sizeof(int));
  vc->slc = (int *) calloc(nmembers, sizeof(int));
  vc->x   = (float *) calloc(nmembers, sizeof(float));
  vc->y   = (float *) calloc(nmembers, sizeof(float));
  vc->z   = (float *) calloc(nmembers, sizeof(float));

  return(vc);
}
/*----------------------------------------------------------------*/
int clustFreeCluster(VOLCLUSTER **ppvc)
{
  VOLCLUSTER *vc;

  vc = *ppvc;

  if(vc->col != NULL) free(&vc->col);
  if(vc->row != NULL) free(&vc->row);
  if(vc->slc != NULL) free(&vc->slc);

  if(vc->x != NULL) free(&vc->x);
  if(vc->y != NULL) free(&vc->y);
  if(vc->z != NULL) free(&vc->z);

  free(ppvc);

  return(0);
}
/*----------------------------------------------------------------*/
VOLCLUSTER **clustAllocClusterList(int nlist)
{
  VOLCLUSTER **vclist;
  vclist = (VOLCLUSTER **) calloc( nlist , sizeof(VOLCLUSTER *));
  if(vclist == NULL){
    fprintf(stderr,"ERROR: clustAllocClusterList: could not alloc %d\n",nlist);
    return(NULL);
  }
  return(vclist);
}
/*------------------------------------------------------------------------*/
int clustFreeClusterList(VOLCLUSTER ***pppvclist, int nlist)
{
  int n;
  VOLCLUSTER **vclist;
  vclist = *pppvclist;
  
  for(n=0;n<nlist;n++)
    if(vclist[n] != NULL) clustFreeCluster(&vclist[n]);

  free(pppvclist);
  
  return(0);
}
/*------------------------------------------------------------------------*/
int clustDumpCluster(FILE *fp, VOLCLUSTER *vc, MRI *vol, int frame)
{
  int n;
  float val;

  for( n = 0; n < vc->nmembers; n++){
    fprintf(fp,"%4d  %3d %3d %3d ", n, vc->col[n], vc->row[n], vc->slc[n]);
    fprintf(fp,"%7.2f %7.2f %7.2f ", vc->x[n], vc->y[n], vc->z[n]);
    if(vol != NULL){
      val = MRIFseq_vox(vol,vc->col[n], vc->row[n], vc->slc[n], frame);
      fprintf(fp,"%12.5f\n",val);
    }
    else fprintf(fp,"\n");
  }
  return(0);
}

/*------------------------------------------------------------------------*/
int clustDumpClusterList(FILE *fp, VOLCLUSTER **vclist, int nlist, 
       MRI *vol, int frame)
{
  int n;

  for( n = 0; n < nlist; n++){
    fprintf(fp,"%3d %5d %5d %g-------------------------------\n",
      n,vclist[n]->nmembers,vclist[n]->maxmember, vclist[n]->maxval);
    clustDumpCluster(fp, vclist[n], vol, frame);
  }

  return(0);
}

/*------------------------------------------------------------------------*/
int clustValueInRange(float val, float thmin, float thmax, int thsign)
{

  if(thsign ==  0) val = fabs(val);
  if(thsign == -1) val = -val;

  if(thmax > 0){
    if(val >= thmin && val <= thmax) return 1;
  }
  else{
    if(val >= thmin ) return 1;
  }

  return(0);
}
/*------------------------------------------------------------------------*/
MRI *clustInitHitMap(MRI *vol, int frame, 
         float thmin, float thmax, int thsign, 
         int *nhits, int **hitcol, int **hitrow, int **hitslc,
         MRI *binmask, int maskframe)
{
  MRI *HitMap;
  int row, col, slc;
  float val;
  int nh, *hrow, *hcol, *hslc;
  int maskval;

  printf("maskframe = %d\n",maskframe);

  /* count the number of hits */
  nh = 0;
  for(col = 0; col < vol->width; col ++){
    for(row = 0; row < vol->height; row ++){
      for(slc = 0; slc < vol->depth; slc ++){
  if(binmask != NULL){
    maskval = MRIIseq_vox(binmask,col,row,slc,maskframe);
    if(maskval == 0) continue;
  }
  val = MRIFseq_vox(vol,col,row,slc,frame);
  if(clustValueInRange(val,thmin,thmax,thsign)) nh++;
      }
    }
  }
  printf("INFO: clustInitHitMap: found %d hits\n", nh );

  /* check that there are hits */
  if( nh == 0 ){
    fprintf(stderr,"ERROR: clustInitHitMap: no hits found\n");
    return(NULL);
  }
  
  /* allocate the rows cols and slices */
  hcol = (int *) calloc( nh, sizeof(int));
  if(hcol == NULL){
    fprintf(stderr,"ERROR: clustInitHitMap: could not alloc hit cols\n");
    MRIfree(&HitMap);
    return(NULL);
  }
  hrow = (int *) calloc( nh, sizeof(int));
  if(hrow == NULL){
    fprintf(stderr,"ERROR: clustInitHitMap: could not alloc hit rows\n");
    free(hcol);
    MRIfree(&HitMap);
    return(NULL);
  }
  hslc = (int *) calloc( nh, sizeof(int));
  if(hslc == NULL){
    fprintf(stderr,"ERROR: clustInitHitMap: could not alloc hit slices\n");
    free(hcol);
    free(hrow);
    MRIfree(&HitMap);
    return(NULL);
  }

  /* allocate the hit map */
  HitMap = MRIalloc(vol->width, vol->height, vol->depth, MRI_INT) ;
  if(HitMap == NULL){
    free(hcol);
    free(hrow);
    free(hslc);
    fprintf(stderr,"ERROR: clustInitHitMap: could not alloc HitMap\n");
    return(NULL);
  }

  /* Now go back through the volume to assign values */
  nh = 0;
  for(col = 0; col < vol->width; col ++){
    for(row = 0; row < vol->height; row ++){
      for(slc = 0; slc < vol->depth; slc ++){
  val = MRIFseq_vox(vol,col,row,slc,frame);
  if(binmask != NULL){
    maskval = MRIIseq_vox(binmask,col,row,slc,maskframe);
    //printf("%2d %2d %2d  %d %g\n",col,row,slc,maskval,val);
    if(maskval == 0){
      MRIIseq_vox(HitMap,col,row,slc,0) = 1;
      continue;
    }
  }
  if(clustValueInRange(val,thmin,thmax,thsign)){
    hcol[nh] = col;
    hrow[nh] = row;
    hslc[nh] = slc;
    MRIIseq_vox(HitMap,col,row,slc,0) = 0;
    nh ++;
  }
  else MRIIseq_vox(HitMap,col,row,slc,0) = 1;
      }
    }
  }
  printf("INFO: clustInitHitMap: found %d hits\n", nh );

  *hitcol = hcol;
  *hitrow = hrow;
  *hitslc = hslc;
  *nhits  = nh;

  return(HitMap);
}
/*------------------------------------------------------------------------*/
int clustAddMember(VOLCLUSTER *vc, int col, int row, int slc)
{
  int nmemb;
  int *tmp;
  float *ftmp;

  nmemb = vc->nmembers;

  tmp = (int *) realloc(vc->col, (nmemb+1) * sizeof(int));
  if(tmp == NULL){
    fprintf(stderr,"ERROR: clustAddMember: could not alloc %d\n",nmemb+1);
    return(1);
  }
  vc->col = tmp;

  tmp = (int *) realloc(vc->row, (nmemb+1) * sizeof(int));
  if(tmp == NULL){
    fprintf(stderr,"ERROR: clustAddMember: could not alloc %d\n",nmemb+1);
    return(1);
  }
  vc->row = tmp;

  tmp = (int *) realloc(vc->slc, (nmemb+1) * sizeof(int));
  if(tmp == NULL){
    fprintf(stderr,"ERROR: clustAddMember: could not alloc %d\n",nmemb+1);
    return(1);
  }
  vc->slc = tmp;

  ftmp = (float *) realloc(vc->x, (nmemb+1) * sizeof(float));
  if(ftmp == NULL){
    fprintf(stderr,"ERROR: clustAddMember: could not alloc %d\n",nmemb+1);
    return(1);
  }
  vc->x = ftmp;

  ftmp = (float *) realloc(vc->y, (nmemb+1) * sizeof(float));
  if(ftmp == NULL){
    fprintf(stderr,"ERROR: clustAddMember: could not alloc %d\n",nmemb+1);
    return(1);
  }
  vc->y = ftmp;

  ftmp = (float *) realloc(vc->z, (nmemb+1) * sizeof(float));
  if(ftmp == NULL){
    fprintf(stderr,"ERROR: clustAddMember: could not alloc %d\n",nmemb+1);
    return(1);
  }
  vc->z = ftmp;

  vc->col[nmemb] = col;
  vc->row[nmemb] = row;
  vc->slc[nmemb] = slc;

  vc->nmembers = nmemb+1;

  return(0);
}
/*------------------------------------------------------------------------*/
int clustGrowOneVoxel(VOLCLUSTER *vc, int col0, int row0, int slc0, 
          MRI *HitMap, int AllowDiag)
{
  int col, row, slc;
  int dcol, drow, dslc;
  int nadded;

  nadded = 0;
  for( dcol = -1; dcol <= +1; dcol++ ){
    for( drow = -1; drow <= +1; drow++ ){
      for( dslc = -1; dslc <= +1; dslc++ ){

  col = col0 + dcol;
  if(col < 0 || col >= HitMap->width) continue;
  
  row = row0 + drow;
  if(row < 0 || row >= HitMap->height) continue;

  slc = slc0 + dslc;
  if(slc < 0 || slc >= HitMap->depth) continue;

  if(!AllowDiag && dcol != 0 && drow != 0 && dslc != 0) continue;

  if(MRIIseq_vox(HitMap,col,row,slc,0)) continue;
  //printf("Adding %3d %3d %3d\n",col,row,slc);

  clustAddMember(vc,col,row,slc);
  MRIIseq_vox(HitMap,col,row,slc,0) = 1;
  nadded ++;

      }
    }
  }
  //printf("GrowOneFrom: %3d %3d %3d  %d\n",col0,row0,slc0,nadded);

  return(nadded);
}
/*------------------------------------------------------------------------*/
VOLCLUSTER *clustGrow(int col0, int row0, int slc0, MRI *HitMap, int AllowDiag)
{
  VOLCLUSTER *vc;
  int nthmember, nmembers_now;
  int col, row, slc;
  int nadded, nthpass, n;

  vc = (VOLCLUSTER *) calloc(1, sizeof(VOLCLUSTER));

  /* put the seed point in the cluster */
  clustAddMember(vc,col0,row0,slc0);
  MRIIseq_vox(HitMap,col0,row0,slc0,0) = 1;

  nthpass = 0;
  nadded = 1;
  while(nadded > 0){
    //printf("%4d  %5d  %d\n",nthpass,vc->nmembers,nadded);

    nadded = 0;
    nmembers_now = vc->nmembers;
    for(nthmember = 0; nthmember < nmembers_now; nthmember ++){
      col = vc->col[nthmember];
      row = vc->row[nthmember];
      slc = vc->slc[nthmember];
      n = clustGrowOneVoxel(vc, col, row, slc, HitMap, AllowDiag);
      //printf("Grown: %3d %3d %3d  %d\n",col,row,slc,n);      
      nadded += n;
    }

    nthpass ++;
  }

  return(vc);
} 
/*-------------------------------------------------------------------*/
int clustMaxMember(VOLCLUSTER *vc, MRI *vol, int frame, int thsign)
{
  int n;
  float val=0.0, val0;

  vc->maxval = 0.0;
  for( n = 0; n < vc->nmembers; n++){
    val0 = MRIFseq_vox(vol,vc->col[n], vc->row[n], vc->slc[n], frame);
    if(thsign ==  1) val = val0;
    if(thsign ==  0) val = fabs(val0);
    if(thsign == -1) val = -val0;
    if(vc->maxval < val){
      vc->maxval    = val0; /* keep orginal signed value */
      vc->maxmember = n;
    }
  }
  
  return(0);
}
/*------------------------------------------------------------------------*/
VOLCLUSTER **clustPruneBySize(VOLCLUSTER **vclist, int nlist, 
            float voxsize, float sizethresh, 
            int *nkeep)
{
  VOLCLUSTER **vcprune;
  int n;
  float clustersize;

  /* count the number to keep */
  *nkeep = 0;
  for(n=0; n < nlist ; n++){
    clustersize = vclist[n]->nmembers * voxsize;
    if( clustersize >= sizethresh) (*nkeep) ++;
  }

  vcprune = (VOLCLUSTER **) calloc( *nkeep , sizeof(VOLCLUSTER *));

  *nkeep = 0;
  for(n=0; n < nlist ; n++){
    clustersize = vclist[n]->nmembers * voxsize;
    if( clustersize >= sizethresh){
      vcprune[(*nkeep)] = clustCopyCluster(vclist[n]);
      (*nkeep) ++;
    }
  }

  return(vcprune);
}
/*------------------------------------------------------------------------*/
VOLCLUSTER **clustPruneByDistance(VOLCLUSTER **vclist, int nlist, 
          float distthresh, int *nkeep)
{
  VOLCLUSTER **vcprune;
  int n1, n2, nmax1, nmax2;
  float max1, max2;
  int keep;
  float x1, y1, z1;
  float x2, y2, z2;
  int pass;
  float d;

  vcprune = NULL;

  /* Two passes: (1) counts the number for alloc, (2) copies clusters */
  for(pass = 1; pass <= 2; pass ++){

    if(pass == 2){
      if( *nkeep == 0) return(NULL);
      vcprune = (VOLCLUSTER **) calloc( *nkeep , sizeof(VOLCLUSTER *));
    }

    /*-- Go through each cluster -- */
    *nkeep = 0;
    for(n1=0; n1 < nlist ; n1++){

      nmax1 = vclist[n1]->maxmember;
      max1 = vclist[n1]->maxval;
      x1   = vclist[n1]->x[nmax1];
      y1   = vclist[n1]->y[nmax1];
      z1   = vclist[n1]->z[nmax1];

      /*-- Compare to every other cluster -- */
      keep = 1;
      for(n2=0; n2 < nlist ; n2++){
  
  if(n1 == n2) continue; /* dont compare to self */
  
  nmax2 = vclist[n2]->maxmember;
  max2 = vclist[n2]->maxval;
  x2   = vclist[n2]->x[nmax2];
  y2   = vclist[n2]->y[nmax2];
  z2   = vclist[n2]->z[nmax2];
  
  /* Compute the distance from the max of one to the max of the
     other */
  d = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) );

  /* If the distance is less than threshold and the max of the
     first is less than the max of the second, throw out the
     first (dont worry about the second here */
  if(d < distthresh && max1 < max2) {
    //printf("Pruning %d: (%5.2f %5.2f %5.2f) (%5.2f %5.2f %5.2f) %g\n",
    // n1,x1,y1,z1,x2,y2,z2,d);
    keep = 0;
    break;
  }

      }/* end n2 loop */

      if(keep){
  if(pass == 2) vcprune[(*nkeep)] = clustCopyCluster(vclist[n1]);
  (*nkeep) ++;
      }

    }/* end n1 loop */
    
  }/* end pass loop */

  return(vcprune);
}
/*----------------------------------------------------------------*/
VOLCLUSTER *clustCopyCluster(VOLCLUSTER *vc)
{
  VOLCLUSTER *vc2;
  int ncopy;

  vc2 = clustAllocCluster(vc->nmembers);

  ncopy = vc->nmembers * sizeof(int);
  memcpy(vc2->col,vc->col,ncopy);
  memcpy(vc2->row,vc->row,ncopy);
  memcpy(vc2->slc,vc->slc,ncopy);

  ncopy = vc->nmembers * sizeof(float);
  memcpy(vc2->x,vc->x,ncopy);
  memcpy(vc2->y,vc->y,ncopy);
  memcpy(vc2->z,vc->z,ncopy);

  vc2->maxmember = vc->maxmember;
  vc2->maxval = vc->maxval;
  
  return(vc2);
}
/*----------------------------------------------------------------*/
VOLCLUSTER **clustCopyClusterList(VOLCLUSTER **vclist, int nlist,
          VOLCLUSTER **vclist2)
{
  int n;

  if(vclist2 == NULL) 
    vclist2 = clustAllocClusterList(nlist);

  for(n=0; n<nlist; n++){
    vclist2[n] = clustCopyCluster(vclist[n]);
  }
  
  return(vclist2);
}

/*----------------------------------------------------------------*/
int clustCompareCluster(const void *a, const void *b)
{
  VOLCLUSTER *vc1, *vc2;

  vc1 = *((VOLCLUSTER **)a);
  vc2 = *((VOLCLUSTER **)b);

  if(fabs(vc1->maxval) > fabs(vc2->maxval) ) return(-1);
  if(fabs(vc1->maxval) < fabs(vc2->maxval) ) return(+1);

  if(vc1->nmembers > vc2->nmembers) return(-1);
  if(vc1->nmembers < vc2->nmembers) return(+1);

  return(0);
}
/*----------------------------------------------------------------*/
VOLCLUSTER **clustSortClusterList(VOLCLUSTER **vclist, int nlist,
          VOLCLUSTER **vcsorted)
{

  if(vcsorted == NULL) 
    vcsorted = clustAllocClusterList(nlist);

  if(vclist != vcsorted)
    clustCopyClusterList(vclist, nlist, vcsorted);

  qsort((void *) vcsorted, nlist, sizeof(VOLCLUSTER **), clustCompareCluster);
  
  return(vcsorted);
}
/*----------------------------------------------------------------
  clustComputeXYZ() - computes the xyz coordinate of each member
  of a cluster given the 4x4 matrix that transforms the col, row, 
  and slice into x, y, and z.
  ----------------------------------------------------------------*/
int clustComputeXYZ(VOLCLUSTER *vc, MATRIX *CRS2XYZ)
{
  int n;

  for(n=0; n < vc->nmembers; n++){
    ConvertCRS2XYZ(vc->col[n],vc->row[n],vc->slc[n], CRS2XYZ,
       &(vc->x[n]), &(vc->y[n]), &(vc->z[n]) );

  }

  return(0);
}
/*----------------------------------------------------------------
  clustComputeTal() - computes the talairach xyz coordinate of each 
  member of a cluster given the 4x4 matrix that transforms the col, 
  row, and slice into MNI coorinates. The MNI coordinates are 
  transformed into talairach coordinates using a piece-wise linear
  transformation. See FixMNITal in transforms.c.
  ----------------------------------------------------------------*/
int clustComputeTal(VOLCLUSTER *vc, MATRIX *CRS2MNI)
{
  int n;

  for(n=0; n < vc->nmembers; n++){
    ConvertCRS2XYZ(vc->col[n],vc->row[n],vc->slc[n], CRS2MNI,
       &(vc->x[n]), &(vc->y[n]), &(vc->z[n]) );
    FixMNITal(  vc->x[n],    vc->y[n],   vc->z[n],
       &(vc->x[n]), &(vc->y[n]), &(vc->z[n]) );
  }

  return(0);
}
/*----------------------------------------------------------------*/
MRI * clustClusterList2Vol(VOLCLUSTER **vclist, int nlist, MRI *tvol,
         int frame, int ValOption)
{
  MRI *vol;
  int nthvc, n;
  VOLCLUSTER *vc;
  float val;

  vol = MRIallocSequence(tvol->width, tvol->height, tvol->depth, tvol->type,1);
  MRIcopyHeader(tvol, vol);

  for(nthvc = 0; nthvc < nlist; nthvc++){
    vc = vclist[nthvc];
    for(n = 0; n < vc->nmembers; n++){
  val = MRIFseq_vox(tvol,vc->col[n],vc->row[n],vc->slc[n],frame);
      if(ValOption == 1)
  MRIFseq_vox(vol,vc->col[n],vc->row[n],vc->slc[n],frame) = val;
      else
  MRIFseq_vox(vol,vc->col[n],vc->row[n],vc->slc[n],frame) = nthvc + 1;
    }
  }
  return(vol);
}
/*----------------------------------------------------------------*/
LABEL *clustCluster2Label(VOLCLUSTER *vc, MRI *vol, int frame,
        float colres, float rowres, float sliceres, 
        MATRIX *FSA2Func)
{
  LABEL *label;
  MATRIX *CRS2Func, *Func2CRS, *Func2FSA;
  MATRIX *xyzFunc, *xyzFSA;
  int n,nlabel;
  float xc, yc, zc, x, y, z, val;

  /* Compute the Cluster XYZ in Functional FOV space */
  Func2CRS = FOVQuantMatrix(vol->width, vol->height, vol->depth, 
          colres, rowres, sliceres); 
  CRS2Func = MatrixInverse(Func2CRS,NULL);
  clustComputeXYZ(vc,CRS2Func);

  /* Matrix to convert from volume to anatomical */
  Func2FSA = MatrixInverse(FSA2Func,NULL);

  /* Preallocate vectors */
  xyzFunc = MatrixAlloc(4, 1, MATRIX_REAL);
  xyzFunc->rptr[4][1] = 1;
  xyzFSA = MatrixAlloc(4, 1, MATRIX_REAL);

  /* First pass: count the number in the label */
  nlabel = 0;
  for(n = 0; n < vc->nmembers; n++){
    xc = vc->x[n];
    yc = vc->y[n];
    zc = vc->z[n];
    /* go through functional space in incr of 1 mm */
    for(x = xc - colres/2;   x <= xc + colres/2;  x += 1.0){
      for(y = yc - sliceres/2; y <= yc + sliceres/2; y += 1.0){
  for(z = zc - rowres/2;   z <= zc + rowres/2;   z += 1.0){
    nlabel ++;
  }
      }
    }
  }
  fprintf(stderr,"INFO: nlabel = %d\n",nlabel);

  /* Alloc the label */
  label = LabelAlloc(nlabel,NULL,NULL);

  /* Second pass: assign label values */
  nlabel = 0;
  for(n = 0; n < vc->nmembers; n++){
    xc = vc->x[n];
    yc = vc->y[n];
    zc = vc->z[n];
    val = MRIFseq_vox(vol,vc->col[n],vc->row[n],vc->slc[n],frame);

    /* go through functional space in incr of 1 mm */
    for(x = xc - colres/2; x <= xc + colres/2; x += 1.0){
      for(y = yc - sliceres/2; y <= yc + sliceres/2; y += 1.0){
  for(z = zc - rowres/2; z <= zc + rowres/2; z += 1.0){

    /* convert Functional XYZ FSA XYZ */
    xyzFunc->rptr[1][1] = x;
    xyzFunc->rptr[2][1] = y;
    xyzFunc->rptr[3][1] = z;
    MatrixMultiply(Func2FSA,xyzFunc,xyzFSA);
    
    /* assign fields to label */
    label->lv[nlabel].x = rint(xyzFSA->rptr[1][1]);
    label->lv[nlabel].y = rint(xyzFSA->rptr[2][1]);
    label->lv[nlabel].z = rint(xyzFSA->rptr[3][1]);
    label->lv[nlabel].stat   = val;
    nlabel ++;
  }
      }
    }
  }
  label->n_points = nlabel;

  MatrixFree(&xyzFunc);
  MatrixFree(&xyzFSA);
  MatrixFree(&Func2FSA);
  MatrixFree(&Func2CRS);
  MatrixFree(&CRS2Func);

  return(label);
}
/*----------------------------------------------------------------*/
/*--------------- STATIC FUNCTIONS BELOW HERE --------------------*/
/*----------------------------------------------------------------*/

/*----------------------------------------------------------------
  ConvertCRS2XYZ() - computes the xyz coordinate given the CRS and
  the transform matrix. This function just hides the matrix 
  operations.
  ----------------------------------------------------------------*/
static int ConvertCRS2XYZ(int col, int row, int slc, MATRIX *CRS2XYZ,
        float *x, float *y, float *z)
{
  MATRIX *crs, *xyz;

  crs = MatrixAlloc(4,1,MATRIX_REAL);
  crs->rptr[1][1] = (float) col;
  crs->rptr[2][1] = (float) row;
  crs->rptr[3][1] = (float) slc;
  crs->rptr[4][1] = 1;

  xyz = MatrixMultiply(CRS2XYZ,crs,NULL);

  *x = xyz->rptr[1][1];
  *y = xyz->rptr[2][1];
  *z = xyz->rptr[3][1];

  MatrixFree(&crs);
  MatrixFree(&xyz);

  return(0);
}

/*----------------------------------------------------*/
CHT *CHTalloc(int n_ithr, double ithr_lo, double ithr_hi,
	      int n_sthr, double sthr_lo, double sthr_hi)
{
  CHT *cht;
  double dithr, dsthr;
  int i,v;
  
  cht = (CHT *) calloc(1,sizeof(CHT));

  cht->n_ithr = n_ithr;
  cht->ithr = (double *) calloc(n_ithr,sizeof(double));
  if(n_ithr != 1) dithr = (ithr_hi - ithr_lo)/(n_ithr-1);
  else            dithr = 0;
  for(i=0; i < n_ithr; i++) cht->ithr[i] = ithr_lo + dithr*i;
  cht->ithr_lo = ithr_lo;
  cht->ithr_hi = ithr_hi;

  cht->n_sthr = n_sthr;
  cht->sthr = (double *) calloc(n_sthr,sizeof(double));
  if(n_sthr != 1) dsthr = (sthr_hi - sthr_lo)/(n_sthr-1);
  else            dsthr = 0;
  for(v=0; v < n_sthr; v++) cht->sthr[v] = sthr_lo + dsthr*v;
  cht->sthr_lo = sthr_lo;
  cht->sthr_hi = sthr_hi;

  cht->hits = (int **) calloc(n_ithr,sizeof(int*));
  for(i=0; i < n_ithr; i++)    
    cht->hits[i] = (int *) calloc(n_sthr,sizeof(int));

  return(cht);
}

/*----------------------------------------------------*/
int CHTfree(CHT **ppcht)
{
  CHT *cht;
  int i;

  cht = *ppcht;
  for(i=0; i < cht->n_ithr; i++) free(cht->hits[i]);
  free(cht->hits);
  free(cht->ithr);
  free(cht->sthr);
  free(*ppcht);
  *ppcht = NULL;

  return(0);
}
/*----------------------------------------------------*/
int CHTprint(FILE *fp, CHT *cht)
{
  int i,v;

  fprintf(fp,"# CHT 1\n");
  fprintf(fp,"# nsim        %d\n",cht->nsim);
  fprintf(fp,"# seed        %ld\n",cht->seed);
  fprintf(fp,"# nvox        %d\n",cht->nvox);
  fprintf(fp,"# totsize     %lf\n",cht->totsize);
  fprintf(fp,"# fwhm        %lf\n",cht->fwhm);
  fprintf(fp,"# nsmooth     %d\n",cht->nsmooth);
  fprintf(fp,"# n_ithr      %d\n",cht->n_ithr);
  fprintf(fp,"# ithr_lo     %lf\n",cht->ithr_lo);
  fprintf(fp,"# ithr_hi     %lf\n",cht->ithr_hi);
  fprintf(fp,"# ithr_sign   %s\n",cht->ithr_sign);
  fprintf(fp,"# n_sthr      %d\n",cht->n_sthr);
  fprintf(fp,"# sthr_lo     %lf\n",cht->sthr_lo);
  fprintf(fp,"# sthr_hi     %lf\n",cht->sthr_hi);
  //fprintf(fp,"# STARTDATA\n");

  fprintf(fp,"     ");
  for(v=0; v < cht->n_sthr; v++) 
    fprintf(fp,"%4.2lf ",cht->sthr[v]);
  fprintf(fp,"\n");

  for(i=0; i < cht->n_ithr; i++) {
    fprintf(fp,"%4.2lf ",cht->ithr[i]);
    for(v=0; v < cht->n_sthr; v++) 
      fprintf(fp,"%5d ",cht->hits[i][v]);
    fprintf(fp,"\n");
  }

  return(0);
}
/*----------------------------------------------------*/
int CHTwrite(char *fname, CHT *cht)
{
  FILE *fp;

  fp = fopen(fname,"w");
  if(fp == NULL){
    printf("ERROR: could not open %s for writing\n",fname);
    return(1);
  }

  CHTprint(fp,cht);

  fclose(fp);

  return(0);
}
/*----------------------------------------------------*/
CHT *CHTread(char *fname)
{
  CHT *cht;
  FILE *fp;
  char tag[1000];
  char tmpstr[1000];
  int nsim; /* number of simulation runs to generate table */
  long int seed; // Seed for random number generator
  int nvox; /* number of voxels/vertices in search area */
  double totsize; /* total volume (mm^3) or area (mm^2) in search*/
  double fwhm;   /* fwhm in mm */
  int nsmooth;   /* number of smooth steps, surf only */
  double   ithr_lo, ithr_hi; /* intensity threshold range */
  int    n_ithr; /* Number ithreshs bet lo and hi*/
  char     ithr_sign[50]; /* abs, pos, neg*/
  double   sthr_lo, sthr_hi; /* volume threshold range */
  int    n_sthr; /* Number sthreshs bet lo and hi*/
  int i,v;
  double dummy;

  fp = fopen(fname,"r");
  if(fp == NULL){
    printf("ERROR: could not open %s for reading\n",fname);
    return(NULL);
  }

  fgets(tmpstr,1000,fp);
  if(tmpstr[0] != '#'){
    printf("ERROR: %s is not formatted correctly (#)\n",fname);
    return(NULL);
  }

  sscanf(tmpstr,"%*s %s",tag);
  if(strcmp(tag,"CHT")){
    printf("ERROR: %s is not formatted correctly (CHT)\n",fname);
    return(NULL);
  }

  while(1){

    // Grab a line
    fgets(tmpstr,1000,fp);

    // Test whether still in the TAG section
    if(tmpstr[0] != '#') break;

    // Scan the tag
    sscanf(tmpstr,"%*s %s",tag);
    //printf("%s \n",tag);

    if(!strcmp(tag,"nsim"))    sscanf(tmpstr,"%*s %*s %d", &nsim);
    if(!strcmp(tag,"seed"))    sscanf(tmpstr,"%*s %*s %ld",&seed);
    if(!strcmp(tag,"nvox"))    sscanf(tmpstr,"%*s %*s %d", &nvox);
    if(!strcmp(tag,"totsize"))  sscanf(tmpstr,"%*s %*s %lf",&totsize);
    if(!strcmp(tag,"fwhm"))    sscanf(tmpstr,"%*s %*s %lf",&fwhm);
    if(!strcmp(tag,"nsmooth")) sscanf(tmpstr,"%*s %*s %d", &nsmooth);
    if(!strcmp(tag,"n_ithr"))  sscanf(tmpstr,"%*s %*s %d", &n_ithr);
    if(!strcmp(tag,"ithr_lo")) sscanf(tmpstr,"%*s %*s %lf",&ithr_lo);
    if(!strcmp(tag,"ithr_hi")) sscanf(tmpstr,"%*s %*s %lf",&ithr_hi);
    if(!strcmp(tag,"ithr_sign")) sscanf(tmpstr,"%*s %*s %s",ithr_sign);
    if(!strcmp(tag,"n_sthr"))  sscanf(tmpstr,"%*s %*s %d", &n_sthr);
    if(!strcmp(tag,"sthr_lo")) sscanf(tmpstr,"%*s %*s %lf",&sthr_lo);
    if(!strcmp(tag,"sthr_hi")) sscanf(tmpstr,"%*s %*s %lf",&sthr_hi);
  }
  
  cht = CHTalloc(n_ithr, ithr_lo, ithr_hi, n_sthr, sthr_lo, sthr_hi);
  cht->nsim    = nsim;
  cht->seed    = seed;
  cht->nvox    = nvox;
  cht->totsize = totsize;
  cht->fwhm    = fwhm;
  cht->nsmooth = nsmooth;
  strcpy(cht->ithr_sign,ithr_sign);
  if(!strcasecmp(ithr_sign,"pos"))      cht->ithr_signid = +1;
  else if(!strcasecmp(ithr_sign,"abs")) cht->ithr_signid =  0;
  else if(!strcasecmp(ithr_sign,"neg")) cht->ithr_signid = -1;
  else{
    printf("ERROR: ithr_sign = %s, unrecognized.\n",ithr_sign);
    printf("       ithr_sign must be pos, neg, or abs.\n");
    return(NULL);
  }

  // The first data line has already been skipped by the fgets above
  
  // Read in the hit table
  for(i=0; i < cht->n_ithr; i++) {
    // Skip first column as it is the ithr for that row
    fscanf(fp,"%lf ",&dummy); 
    for(v=0; v < cht->n_sthr; v++) 
      fscanf(fp,"%d ",&(cht->hits[i][v]));
  }

  fclose(fp);

  return(cht);
}

/*-----------------------------------------------------------------
  CHTcompare() - compares two CHTs. Returns 1 if they are different
  and 0 if they are the same. If an item in the targ is less than 0,
  then it's value is replaced with the value of the source. This 
  does not trigger a 1 return and allows for values to be filled in.
  Does not compare seeds.
  ----------------------------------------------------------------*/
int CHTcompare(CHT *src, CHT *targ)
{
  if(targ->nvox < 0) targ->nvox = src->nvox;
  else if(targ->nvox != src->nvox) return(1);

  if(targ->totsize < 0) targ->totsize = src->totsize;
  else if(targ->totsize != src->totsize) return(1);

  if(targ->fwhm < 0) targ->fwhm = src->fwhm;
  else if(targ->fwhm != src->fwhm) return(1);

  if(targ->nsmooth < 0) targ->nsmooth = src->nsmooth;
  else if(targ->nsmooth != src->nsmooth) return(1);

  if(targ->ithr_lo < 0) targ->ithr_lo = src->ithr_lo;
  else if(targ->ithr_lo != src->ithr_lo) return(1);

  if(targ->ithr_hi < 0) targ->ithr_hi = src->ithr_hi;
  else if(targ->ithr_hi != src->ithr_hi) return(1);

  if(targ->n_ithr < 0) targ->n_ithr = src->n_ithr;
  else if(targ->n_ithr != src->n_ithr) return(1);

  if(strlen(targ->ithr_sign)==0) 
    CHTsetSignString(targ, src->ithr_sign);
  else if(strcmp(targ->ithr_sign,src->ithr_sign)) return(1);

  if(targ->sthr_lo < 0) targ->sthr_lo = src->sthr_lo;
  else if(targ->sthr_lo != src->sthr_lo) return(1);

  if(targ->sthr_hi < 0) targ->sthr_hi = src->sthr_hi;
  else if(targ->sthr_hi != src->sthr_hi) return(1);

  if(targ->n_sthr < 0) targ->n_sthr = src->n_sthr;
  else if(targ->n_sthr != src->n_sthr) return(1);

  return(0);
}

/*--------------------------------------------------------------
  CHTsetSignString() - sets the sign string. If the string is
  unrecognized, returns 1, otherwise 0. If the string is NULL,
  abs is used.
  --------------------------------------------------------------*/
int CHTsetSignString(CHT *cht, char *ithr_sign)
{
  int ithr_signid;

  ithr_signid = CHTsignId(ithr_sign);
  if(ithr_signid == -100) return(1);

  cht->ithr_signid = ithr_signid;

  if(ithr_sign == NULL) strcpy(cht->ithr_sign,"abs");
  else                  strcpy(cht->ithr_sign,ithr_sign);

  return(0);
}

/*--------------------------------------------------------------
  CHTsignId() - converts the sign string to a numeric code. The
  code is set in the CHT structure and returned.
  --------------------------------------------------------------*/
int CHTsignId(char *ithr_sign)
{

  if(ithr_sign == NULL) return(0); // abs
  if(!strcasecmp(ithr_sign,"pos"))      return(+1);
  else if(!strcasecmp(ithr_sign,"abs")) return(0);
  else if(!strcasecmp(ithr_sign,"neg")) return(-1);

  printf("ERROR: ithr_sign = %s, unrecognized.\n",ithr_sign);
  printf("       ithr_sign must be pos, neg, or abs.\n");
  return(-100);
}

/*-------###########################################--------*/
/*-------###########################################--------*/
/*-------###########################################--------*/

/*--------------------------------------------------------------
  CSDread() - reads a cluster simulation data file. The format
  of this file is currently defined by mri_glmfit.
  --------------------------------------------------------------*/
CLUSTER_SIM_DATA *CSDread(char *csdfile)
{
  FILE *fp;
  CLUSTER_SIM_DATA *csd;
  char tag[1000], tmpstr[1000];
  int r,nthrep;

  fp = fopen(csdfile,"r");
  if(fp == NULL){
    printf("ERROR: CSDread(): could not open %s\n",csdfile);
    return(NULL);
  }

  csd = (CLUSTER_SIM_DATA *) calloc(sizeof(CLUSTER_SIM_DATA),1);
  csd->mergedflag = 0; // not a merged data 

  // Go through each input line
  nthrep = 0;
  while(1){
    r = fscanf(fp,"%s",tag);
    if(r==EOF) break;

    if(!strcmp(tag,"#")){
      // ----------- Header ------------
      fscanf(fp,"%s",tag);
      if(!strcmp(tag,"simtype")) fscanf(fp,"%s",csd->simtype);
      else if(!strcmp(tag,"anattype")){
	fscanf(fp,"%s",csd->anattype);
	fscanf(fp,"%s",csd->subject);
	fscanf(fp,"%s",csd->hemi);
      }
      else if(!strcmp(tag,"thresh")) fscanf(fp,"%lf",&(csd->thresh));
      else if(!strcmp(tag,"seed"))   fscanf(fp,"%ld",&(csd->seed));
      if(!strcmp(tag,"contrast"))    fscanf(fp,"%s",csd->contrast);
      else if(!strcmp(tag,"nsim")){
	fscanf(fp,"%d",&(csd->nreps));
	CSDallocData(csd);
      }
      else fgets(tmpstr,1000,fp); // not an interesting line, so get past it
    }
    else{
      // ----------- Data ------------
      //printf("%s \n",tag);
      fscanf(fp,"%d", &(csd->nClusters[nthrep]));
      fscanf(fp,"%lf",&(csd->MaxClusterSize[nthrep]));
      fscanf(fp,"%lf",&(csd->MaxSig[nthrep]));
      nthrep++;
    }
  }
  return(csd);
}
/*--------------------------------------------------------------
  CSDreadMerge() - reads in a CSD file and merges it with 
  another CSD. If the input csd is NULL, then it is the 
  same as CSDread(). The purpose is to be able to do somthing
  like this: csd = CSDreadMerge(csdfile, csd);
  --------------------------------------------------------------*/
CLUSTER_SIM_DATA *CSDreadMerge(char *csdfile, CSD *csd)
{
  CLUSTER_SIM_DATA *csd2=NULL, *csdmerged=NULL;

  if(csd == NULL){
    csd = CSDread(csdfile);
    return(csd);
  }
  
  csd2 = CSDread(csdfile);
  if(csd2 == NULL) return(NULL);

  csdmerged = CSDmerge(csd,csd2);
  if(csdmerged == NULL){
    CSDfreeData(csd2);
    return(NULL);
  }

  CSDcopy(csdmerged,csd);

  CSDfreeData(csd2);
  CSDfreeData(csdmerged);
  return(csd);
}

/*--------------------------------------------------------------
  CSDallocData() - allocates the arrays for a CSD (not the
  structure iteself).
  --------------------------------------------------------------*/
int CSDallocData(CLUSTER_SIM_DATA *csd)
{
  csd->nClusters = (int *) calloc(csd->nreps, sizeof(int));
  csd->MaxClusterSize = (double *) calloc(csd->nreps, sizeof(double));
  csd->MaxSig = (double *) calloc(csd->nreps, sizeof(double));
  return(0);
}
/*--------------------------------------------------------------
  CSDfreeData() - frees the data arrays and sets their pointers
  to NULL. Does not try to free the structure.
  --------------------------------------------------------------*/
int CSDfreeData(CLUSTER_SIM_DATA *csd)
{
  if(csd->nClusters){
    free(csd->nClusters); 
    csd->nClusters=NULL;
  }
  if(csd->MaxClusterSize){
    free(csd->MaxClusterSize); 
    csd->MaxClusterSize=NULL;
  }
  if(csd->MaxSig){
    free(csd->MaxSig); 
    csd->MaxSig = NULL;
  }
  return(0);
}
/*--------------------------------------------------------------
  CSDcopy() - copies csd into csdcopy. If csdcopy is NULL, it is
  allocated. If csdcopy is non-null, the data are freed and 
  then re-allocated.
  --------------------------------------------------------------*/
CSD *CSDcopy(CSD *csd, CSD *csdcopy)
{
  int nthrep;
  
  if(csdcopy == NULL)
    csdcopy = (CLUSTER_SIM_DATA *) calloc(sizeof(CLUSTER_SIM_DATA),1);
  else CSDfreeData(csdcopy);

  strcpy(csdcopy->simtype, csd->simtype);
  strcpy(csdcopy->anattype,csd->anattype);
  strcpy(csdcopy->subject, csd->subject);
  strcpy(csdcopy->hemi,    csd->hemi);
  csdcopy->thresh = csd->thresh;

  csdcopy->nreps = csd->nreps;
  CSDallocData(csdcopy);

  for(nthrep = 0; nthrep < csd->nreps; nthrep++){
    csdcopy->nClusters[nthrep]      = csd->nClusters[nthrep];
    csdcopy->MaxClusterSize[nthrep] = csd->MaxClusterSize[nthrep];
    csdcopy->MaxSig[nthrep]         = csd->MaxSig[nthrep];
  }
  return(csdcopy);
}
/*--------------------------------------------------------------
  CSDmerge() - merge two CSDs into one. Requires that: 
    (1) simtypes be the same
    (2) anattypes be the same
    (3) contrasts be the same
    (4) thresholds be the same
    (5) seeds be different
  The seed from the first is copied into the merge, and the
  mergeflag is set to 1.
  --------------------------------------------------------------*/
CSD *CSDmerge(CSD *csd1, CSD *csd2)
{
  int nthrep1, nthrep2, nthrep;
  CSD *csd;

  if(strcmp(csd1->simtype,csd2->simtype)){
    printf("ERROR: CSDmerge: CSDs have different sim types\n");
    return(NULL);
  }
  if(strcmp(csd1->anattype,csd2->anattype)){
    printf("ERROR: CSDmerge: CSDs have different anat types\n");
    return(NULL);
  }
  if(strcmp(csd1->contrast,csd2->contrast)){
    printf("ERROR: CSDmerge: CSDs have different contrasts\n");
    return(NULL);
  }
  if(csd1->thresh != csd2->thresh){
    printf("ERROR: CSDmerge: CSDs have different thresholds\n");
    return(NULL);
  }
  if(csd1->seed == csd2->seed){
    printf("ERROR: CSDmerge: CSDs have same seed\n");
    return(NULL);
  }
  
  csd = (CLUSTER_SIM_DATA *) calloc(sizeof(CLUSTER_SIM_DATA),1);
  strcpy(csd->simtype,  csd1->simtype);
  strcpy(csd->anattype, csd1->anattype);
  strcpy(csd->contrast, csd1->contrast);
  strcpy(csd->subject,  csd1->subject);
  strcpy(csd->hemi,     csd1->hemi);
  csd->thresh = csd1->thresh;
  csd->seed   = csd1->seed;
  csd->mergedflag = 1;

  csd->nreps = csd1->nreps + csd2->nreps;
  CSDallocData(csd);

  nthrep = 0;
  for(nthrep1 = 0; nthrep1 < csd1->nreps; nthrep1++){
    csd->nClusters[nthrep]      = csd1->nClusters[nthrep1];
    csd->MaxClusterSize[nthrep] = csd1->MaxClusterSize[nthrep1];
    csd->MaxSig[nthrep]         = csd1->MaxSig[nthrep1];
    nthrep++;
  }
  for(nthrep2 = 0; nthrep2 < csd2->nreps; nthrep2++){
    csd->nClusters[nthrep]      = csd2->nClusters[nthrep2];
    csd->MaxClusterSize[nthrep] = csd2->MaxClusterSize[nthrep2];
    csd->MaxSig[nthrep]         = csd2->MaxSig[nthrep2];
    nthrep++;
  }

  return(csd);
}
/*--------------------------------------------------------------
  CSDprint() - prints a CSD to the given stream.
  --------------------------------------------------------------*/
int CSDprint(FILE *fp, CLUSTER_SIM_DATA *csd)
{
  int nthrep;

  fprintf(fp,"# simtype %s\n",csd->simtype);
  if(!strcmp(csd->anattype,"surface"))
    fprintf(fp,"# anattype %s  %s %s\n",csd->anattype,csd->subject,csd->hemi);
  else
    fprintf(fp,"# anattype %s \n",csd->anattype);
  fprintf(fp,"# merged %d\n",csd->mergedflag);
  fprintf(fp,"# contrast %s\n",csd->contrast);
  fprintf(fp,"# seed %ld\n",csd->seed);
  fprintf(fp,"# thresh %lf\n",csd->thresh);
  fprintf(fp,"# nreps %d\n",csd->nreps);
  
  for(nthrep = 0; nthrep < csd->nreps; nthrep++){
    fprintf(fp,"%3d %3d %g %g \n",nthrep,csd->nClusters[nthrep],
	    csd->MaxClusterSize[nthrep],csd->MaxSig[nthrep]);
  }
  return(0);
}
/*--------------------------------------------------------------
  CSDpvalClustSize() - computes the emperical pvalue for a given
  cluster size, as well as the confidence interval. ciPct is the conf
  interval given in percent (eg, 90 for 90%). This means that there
  will be a 90% chance that the "true" pvalue for the cluster will lie
  between pvalLow and pvalHi (based on binomial distribution). If
  no item from the simulation is larger than ClusterSize, then
  it is assumed that 1 item is so that things dont break.
  --------------------------------------------------------------*/
double CSDpvalClustSize(CLUSTER_SIM_DATA *csd, double ClusterSize,
			double ciPct, double *pvalLow, double *pvalHi)
{
  int nthrep,nover, k, nlow, nhi;
  double pval,psum, pcilow, pcihi;

  // First, count the number of MaxClusters whose size is greater than
  // the one under test
  nover = 0;
  for(nthrep = 0; nthrep < csd->nreps; nthrep++)
    if(csd->MaxClusterSize[nthrep] > ClusterSize) nover++;

  // If none is over, then set nover = 1 so that things don't beak
  if(nover == 0) nover = 1;

  // Compute the nomial pvalue
  pval = (double)nover/csd->nreps;

  // Ranges for confidence interval
  pcihi  = ciPct/100;
  pcilow = 1-pcihi;

  // Loop thru all possible outcomes (k), computing the CDF of the 
  // binomial distribution, until the upper conf interval is reached.
  nlow = -1;
  k = 0;
  psum = gsl_ran_binomial_pdf(k,pval,csd->nreps);
  while(psum < pcihi && k <= csd->nreps){
    //printf("%3d %lf\n",k,psum);
    if(nlow < 0 && psum > pcilow) nlow = k;
    k++;
    psum += gsl_ran_binomial_pdf(k,pval,csd->nreps);
  }
  nhi = k;

  // Compute the pvalues at the lower and upper confidence intervals
  *pvalLow = (double)nlow/csd->nreps;
  *pvalHi  = (double)nhi/csd->nreps;

  //printf("csize=%lf  p=%lf  ci=%lf  nLow=%d pLow=%lf nHi=%d pHi=%lf\n",
  //	 ClusterSize,pval,ciPct,nlow,*pvalLow,nhi,*pvalHi);

  return(pval);
}

