#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "mri.h"
#include "resample.h"
#include "matrix.h"

#include "volcluster.h"

static int ConvertMNI2Tal(float  xmni, float  ymni, float  zmni,
			  float *xtal, float *ytal, float *ztal);
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
  transformation. See ConvertMNI2Tal.
  ----------------------------------------------------------------*/
int clustComputeTal(VOLCLUSTER *vc, MATRIX *CRS2MNI)
{
  int n;

  for(n=0; n < vc->nmembers; n++){
    ConvertCRS2XYZ(vc->col[n],vc->row[n],vc->slc[n], CRS2MNI,
		   &(vc->x[n]), &(vc->y[n]), &(vc->z[n]) );
    ConvertMNI2Tal(  vc->x[n],    vc->y[n],   vc->z[n],
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
/*-----------------------------------------------------------
  ConvertMNI2Tal() - function to compute the "real" talairach
  coordinates from the MNI talaiarch coordinates. Has nothing
  to do with clustering.
  -----------------------------------------------------------*/
static int ConvertMNI2Tal(float  xmni, float  ymni, float  zmni,
			  float *xtal, float *ytal, float *ztal)
{
  MATRIX *T, *xyzMNI, *xyzTal;


  T = MatrixAlloc(4, 4, MATRIX_REAL);
  if(zmni >= 0.0){
    stuff_four_by_four(T, 
		       .9900,  .0000, .0000, 0,
		       .0000,  .9688, .0460, 0,
		       .0000, -.0485, .9189, 0,
		       .0000,  .0000, .0000, 1);
  }
  else {
    stuff_four_by_four(T, 
		       .9900,  .0000, .0000, 0,
		       .0000,  .9688, .0420, 0,
		       .0000, -.0485, .8390, 0,
		       .0000,  .0000, .0000, 1);
  }

  xyzMNI = MatrixAlloc(4, 1, MATRIX_REAL);
  xyzMNI->rptr[1][1] = xmni;
  xyzMNI->rptr[2][1] = ymni;
  xyzMNI->rptr[3][1] = zmni;
  xyzMNI->rptr[4][1] = 1.0;

  xyzTal = MatrixMultiply(T,xyzMNI,NULL);

  *xtal = xyzTal->rptr[1][1];
  *ytal = xyzTal->rptr[2][1];
  *ztal = xyzTal->rptr[3][1];

  MatrixFree(&T);
  MatrixFree(&xyzMNI);
  MatrixFree(&xyzTal);

  return(0);
}
