#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

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
	      int n_vthr, double vthr_lo, double vthr_hi)
{
  CHT *cht;
  double dithr, dvthr;
  int i,v;
  
  cht = (CHT *) calloc(1,sizeof(CHT));

  cht->n_ithr = n_ithr;
  cht->ithr = (double *) calloc(n_ithr,sizeof(double));
  if(n_ithr != 1) dithr = (ithr_hi - ithr_lo)/(n_ithr-1);
  else            dithr = 0;
  for(i=0; i < n_ithr; i++) cht->ithr[i] = ithr_lo + dithr*i;
  cht->ithr_lo = ithr_lo;
  cht->ithr_hi = ithr_hi;

  cht->n_vthr = n_vthr;
  cht->vthr = (double *) calloc(n_vthr,sizeof(double));
  if(n_vthr != 1) dvthr = (vthr_hi - vthr_lo)/(n_vthr-1);
  else            dvthr = 0;
  for(v=0; v < n_vthr; v++) cht->vthr[v] = vthr_lo + dvthr*v;
  cht->vthr_lo = vthr_lo;
  cht->vthr_hi = vthr_hi;

  cht->hits = (int **) calloc(n_ithr,sizeof(int*));
  for(i=0; i < n_ithr; i++)    
    cht->hits[i] = (int *) calloc(n_vthr,sizeof(int));

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
  free(cht->vthr);
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
  fprintf(fp,"# nvox        %d\n",cht->nvox);
  fprintf(fp,"# totvol      %lf\n",cht->totvol);
  fprintf(fp,"# fwhm        %lf\n",cht->fwhm);
  fprintf(fp,"# nsmooth     %d\n",cht->nsmooth);
  fprintf(fp,"# n_ithr      %d\n",cht->n_ithr);
  fprintf(fp,"# ithr_lo     %lf\n",cht->ithr_lo);
  fprintf(fp,"# ithr_hi     %lf\n",cht->ithr_hi);
  fprintf(fp,"# ithr_sign   %s\n",cht->ithr_sign);
  fprintf(fp,"# n_vthr      %d\n",cht->n_vthr);
  fprintf(fp,"# vthr_lo     %lf\n",cht->vthr_lo);
  fprintf(fp,"# vthr_hi     %lf\n",cht->vthr_hi);
  //fprintf(fp,"# STARTDATA\n");

  fprintf(fp,"     ");
  for(v=0; v < cht->n_vthr; v++) 
    fprintf(fp,"%4.2lf ",cht->vthr[v]);
  fprintf(fp,"\n");

  for(i=0; i < cht->n_ithr; i++) {
    fprintf(fp,"%4.2lf ",cht->ithr[i]);
    for(v=0; v < cht->n_vthr; v++) 
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
  int nvox; /* number of voxels/vertices in search area */
  double totvol; /* total volume (mm^3) or area (mm^2) in search*/
  double fwhm;   /* fwhm in mm */
  int nsmooth;   /* number of smooth steps, surf only */
  double   ithr_lo, ithr_hi; /* intensity threshold range */
  int    n_ithr; /* Number ithreshs bet lo and hi*/
  char     ithr_sign[50]; /* abs, pos, neg*/
  double   vthr_lo, vthr_hi; /* volume threshold range */
  int    n_vthr; /* Number vthreshs bet lo and hi*/
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
    if(!strcmp(tag,"nvox"))    sscanf(tmpstr,"%*s %*s %d", &nvox);
    if(!strcmp(tag,"totvol"))  sscanf(tmpstr,"%*s %*s %lf",&totvol);
    if(!strcmp(tag,"fwhm"))    sscanf(tmpstr,"%*s %*s %lf",&fwhm);
    if(!strcmp(tag,"nsmooth")) sscanf(tmpstr,"%*s %*s %d", &nsmooth);
    if(!strcmp(tag,"n_ithr"))  sscanf(tmpstr,"%*s %*s %d", &n_ithr);
    if(!strcmp(tag,"ithr_lo")) sscanf(tmpstr,"%*s %*s %lf",&ithr_lo);
    if(!strcmp(tag,"ithr_hi")) sscanf(tmpstr,"%*s %*s %lf",&ithr_hi);
    if(!strcmp(tag,"ithr_sign")) sscanf(tmpstr,"%*s %*s %s",ithr_sign);
    if(!strcmp(tag,"n_vthr"))  sscanf(tmpstr,"%*s %*s %d", &n_vthr);
    if(!strcmp(tag,"vthr_lo")) sscanf(tmpstr,"%*s %*s %lf",&vthr_lo);
    if(!strcmp(tag,"vthr_hi")) sscanf(tmpstr,"%*s %*s %lf",&vthr_hi);
  }
  
  cht = CHTalloc(n_ithr, ithr_lo, ithr_hi, n_vthr, vthr_lo, vthr_hi);
  cht->nsim    = nsim;
  cht->nvox    = nvox;
  cht->totvol  = totvol;
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
    for(v=0; v < cht->n_vthr; v++) 
      fscanf(fp,"%d ",&(cht->hits[i][v]));
  }

  fclose(fp);

  return(cht);
}
