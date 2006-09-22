#include <stdio.h>
#include <stdlib.h>

#include "mri.h"
#include "diag.h"
#include "error.h"
#include "voxlist.h"

VOXEL_LIST  *
VLSTcreateInRegion(MRI *mri, float low_val, float hi_val , 
                   VOXEL_LIST *vl, int skip, int border_only, 
                   MRI_REGION *box)
{
  int   x, y, z, nvox, i, width, height, depth ;
  Real  val ;

  skip++ ;  /* next voxel + amount to skip */
  width = box->x+box->dx ; height = box->y+box->dy ; depth = box->z+box->dz ;
  for (nvox = 0, x = box->x ; x < width ; x+=skip)
	{
		for (y = box->y ; y < height ; y+=skip)
		{
			for (z = box->z ; z < depth ; z+=skip)
	    {
	      if (x == Gx && y == Gy && z == Gz)
					DiagBreak() ;
	      val = MRIgetVoxVal(mri, x, y, z, 0) ;
	      if (val > 0)
					DiagBreak() ;
	      if (x == Gx && y == Gy && z == Gz)
					DiagBreak() ;
	      if (val >= low_val && val <= hi_val)
					nvox++ ;
	    }
		}
	}

	if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
		printf("allocating %d voxel indices...\n", nvox) ;
  if (vl == NULL)
		vl = (VOXEL_LIST *)calloc(1, sizeof(VOXEL_LIST)) ;
  else if (vl->nvox < nvox)
  {
    free(vl->xi) ; free(vl->yi) ; free(vl->zi) ;
    vl->xi = vl->yi = zl->zi = NULL ;
  }
  if (vl->xi == NULL)
  {
		vl->xi = (int *)calloc(nvox, sizeof(int)) ;
		vl->yi = (int *)calloc(nvox, sizeof(int)) ;
		vl->zi = (int *)calloc(nvox, sizeof(int)) ;
  }
  if (!vl || !vl->xi || !vl->yi || !vl->zi)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d voxel list\n",
              Progname, nvox) ;
  vl->nvox = nvox ;
  for (nvox = 0, x = box->x ; x < width ; x+=skip)
	{
		for (y = box->y ; y < height ; y+=skip)
		{
			for (z = box->z ; z < depth ; z+=skip)
	    {
	      val = MRIgetVoxVal(mri, x, y, z, 0) ;
	      if (val >= low_val && val <= hi_val)
				{
					if (x == Gx && y == Gy && z == Gz)
						DiagBreak() ;
					if ((border_only == 0) || 
							(MRIneighbors(mri, x, y, z, 0) > 0))
					{
						i = nvox++ ;
						vl->xi[i] = x ; vl->yi[i] = y ; vl->zi[i] = z ;
					}
				}
	    }
		}
	}
  vl->mri = mri ;
  return(vl) ;
}
VOXEL_LIST *
VLSTcreate(MRI *mri, 
	   float low_val, 
	   float hi_val, 
	   VOXEL_LIST *vl, 
	   int skip, 
	   int border_only)
{
  int   x, y, z, nvox, i ;
  Real  val ;

  skip++ ;  /* next voxel + amount to skip */
  for (nvox = x = 0 ; x < mri->width ; x+=skip)
	{
		for (y = 0 ; y < mri->height ; y+=skip)
		{
			for (z = 0 ; z < mri->depth ; z+=skip)
	    {
	      if (x == Gx && y == Gy && z == Gz)
					DiagBreak() ;
	      val = MRIgetVoxVal(mri, x, y, z, 0) ;
	      if (val > 0)
					DiagBreak() ;
	      if (x == Gx && y == Gy && z == Gz)
					DiagBreak() ;
	      if (val >= low_val && val <= hi_val)
					nvox++ ;
	    }
		}
	}

	if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
		printf("allocating %d voxel indices...\n", nvox) ;
  if (vl == NULL)
		vl = (VOXEL_LIST *)calloc(1, sizeof(VOXEL_LIST)) ;
  else if (vl->nvox < nvox)
  {
    free(vl->xi) ; free(vl->yi) ; free(vl->zi) ;
    vl->xi = vl->yi = zl->zi = NULL ;
  }
  if (vl->xi == NULL)
  {
		vl->xi = (int *)calloc(nvox, sizeof(int)) ;
		vl->yi = (int *)calloc(nvox, sizeof(int)) ;
		vl->zi = (int *)calloc(nvox, sizeof(int)) ;
  }
  if (!vl || !vl->xi || !vl->yi || !vl->zi)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d voxel list\n",
              Progname, nvox) ;
  vl->nvox = nvox ;
  for (nvox = x = 0 ; x < mri->width ; x+=skip)
	{
		for (y = 0 ; y < mri->height ; y+=skip)
		{
			for (z = 0 ; z < mri->depth ; z+=skip)
	    {
	      val = MRIgetVoxVal(mri, x, y, z, 0) ;
	      if (val >= low_val && val <= hi_val)
				{
					if (x == Gx && y == Gy && z == Gz)
						DiagBreak() ;
					if ((border_only == 0) || 
							(MRIneighbors(mri, x, y, z, 0) > 0))
					{
						i = nvox++ ;
						vl->xi[i] = x ; vl->yi[i] = y ; vl->zi[i] = z ;
					}
				}
	    }
		}
	}
  vl->mri = mri ;
  return(vl) ;
}

int
VLSTfree(VOXEL_LIST **pvl)
{
  VOXEL_LIST *vl = *pvl ;
  *pvl = NULL ;

  if (!vl)
    return(ERROR_BADPARM) ;
  free(vl->xi) ;
  free(vl->yi) ;
  free(vl->zi) ;
  free(vl) ;
  return(NO_ERROR) ;
}
MRI *
VLSTcreateMri(VOXEL_LIST *vl, int val)
{
  int   i ;
  MRI *mri ;
	
  mri = MRIalloc(vl->mri->width, vl->mri->height, vl->mri->depth, MRI_UCHAR) ;
  MRIcopyHeader(vl->mri, mri) ;
  for (i = 0 ; i < vl->nvox ; i++)
    {
      MRIsetVoxVal(mri, vl->xi[i], vl->yi[i], vl->zi[i], 0, val) ;
    }
  return(mri) ;
}

MRI *
VLSTaddToMri(VOXEL_LIST *vl, MRI *mri, int val)
{
  int   i ;

  if (mri == NULL)
    {
      mri = MRIalloc(vl->mri->width, 
		     vl->mri->height, 
		     vl->mri->depth, 
		     MRI_UCHAR) ;
      MRIcopyHeader(vl->mri, mri) ;
    }
  for (i = 0 ; i < vl->nvox ; i++)
    {
      MRIsetVoxVal(mri, vl->xi[i], vl->yi[i], vl->zi[i], 0, val) ;
    }
  return(mri) ;
}
MRI *
VLSTtoMri(VOXEL_LIST *vl, MRI *mri)
{
  int   i ;
  Real  val ;

  if (mri == NULL)
    mri = MRIclone(vl->mri, NULL) ;

  for (i = 0 ; i < vl->nvox ; i++)
    {
      val = MRIgetVoxVal(vl->mri, vl->xi[i], vl->yi[i], vl->zi[i], 0) ;
      MRIsetVoxVal(mri, vl->xi[i], vl->yi[i], vl->zi[i], 0, val) ;
    }
  return(mri) ;
}

VOXEL_LIST *
VLSTdilate(VOXEL_LIST *vl, int mode, MRI *mri_exclude)
{
  MRI         *mri_current, *mri_new ; 
  int         i, xi, yi, zi, xk, yk, zk, nvox ;
  VOXEL_LIST  *vl_exp = NULL;

  // create volume to keep track of what voxels are in current list
  mri_current = VLSTcreateMri(vl, 1) ;
  mri_new = MRIclone(mri_current, NULL) ;

  // count how many vox will be in new list
  for (nvox = i = 0 ; i < vl->nvox ; i++)
    {
      for (xk = -1 ; xk <= 1 ; xk++)
	{
	  xi = vl->mri->xi[vl->xi[i]+xk] ;
	  for (yk = -1 ; yk <= 1 ; yk++)
	    {
	      yi = vl->mri->yi[vl->yi[i]+yk] ;
	      for (zk = -1 ; zk <= 1 ; zk++)
		{
		  zi = vl->mri->zi[vl->zi[i]+zk] ;
		  if (nint(MRIgetVoxVal(mri_current, xi, yi, zi,0)) == 0 && 
		      nint(MRIgetVoxVal(mri_new, xi, yi, zi,0)) == 0 &&
		      (mri_exclude == NULL || 
		       FZERO(MRIgetVoxVal(mri_exclude, xi, yi, zi, 0))))
		    {
		      MRIvox(mri_new, xi, yi, zi) = 1 ;
		      nvox++ ;
		    }
		}
	    }
	}
    }

  if (nvox == 0)
    return(NULL) ;

  if (mode == VL_DILATE_ADD)
    {
      vl_exp = (VOXEL_LIST *)calloc(1, sizeof(VOXEL_LIST)) ;
      vl_exp->nvox = vl->nvox + nvox ;
      vl_exp->mri = vl->mri ;
      vl_exp->xi = (int *)calloc(vl->nvox+nvox, sizeof(int)) ;
      vl_exp->yi = (int *)calloc(vl->nvox+nvox, sizeof(int)) ;
      vl_exp->zi = (int *)calloc(vl->nvox+nvox, sizeof(int)) ;
      if (!vl_exp || !vl_exp->xi || !vl_exp->yi || !vl_exp->zi)
	ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d voxel list\n",
		  Progname, nvox) ;
      for (i = 0 ; i < vl->nvox ; i++)
	{
	  vl_exp->xi[i] = vl->xi[i] ;
	  vl_exp->yi[i] = vl->yi[i] ;
	  vl_exp->zi[i] = vl->zi[i] ;
	}
      nvox = vl->nvox ;
    }
  else if (mode == VL_DILATE_REPLACE)
    {
      vl_exp = (VOXEL_LIST *)calloc(1, sizeof(VOXEL_LIST)) ;
      vl_exp->nvox = nvox ;
      vl_exp->mri = vl->mri ;
      vl_exp->xi = (int *)calloc(nvox, sizeof(int)) ;
      vl_exp->yi = (int *)calloc(nvox, sizeof(int)) ;
      vl_exp->zi = (int *)calloc(nvox, sizeof(int)) ;
      if (!vl_exp || !vl_exp->xi || !vl_exp->yi || !vl_exp->zi)
	ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d voxel list\n",
		  Progname, nvox) ;
      nvox = 0 ;
    }
  for (i = 0 ; i < vl->nvox ; i++)
    {
      for (xk = -1 ; xk <= 1 ; xk++)
	{
	  xi = vl->mri->xi[vl->xi[i]+xk] ;
	  for (yk = -1 ; yk <= 1 ; yk++)
	    {
	      yi = vl->mri->yi[vl->yi[i]+yk] ;
	      for (zk = -1 ; zk <= 1 ; zk++)
		{
		  zi = vl->mri->zi[vl->zi[i]+zk] ;
		  if (MRIvox(mri_new, xi, yi, zi) == 1) // add it 
		    {
		      MRIvox(mri_new, xi, yi, zi) = 0 ;  // only add it once
		      vl_exp->xi[nvox] = xi ; 
		      vl_exp->yi[nvox] = yi ; 
		      vl_exp->zi[nvox] = zi ;
		      nvox++ ;
		    }
		}
	    }
	}
    }

  MRIfree(&mri_current) ;
  MRIfree(&mri_new) ;

  return(vl_exp) ;
}



void VLSTcomputeStats(VOXEL_LIST *vl)
{
  double mean =0; 
  double std = 0;
  int i;
  double val;
  
  if(vl->mri == NULL || vl->nvox <= 0){
    vl->mean = 0;
    vl->std = 0;
    return;
  }

  for (i = 0 ; i < vl->nvox ; i++){
    val = MRIgetVoxVal(vl->mri, vl->xi[i], vl->yi[i], vl->zi[i], 0);
    mean += val;
    std += val*val;
  }
  
  mean /= (double)vl->nvox;
  std /= (double)vl->nvox;
  
  std = sqrt(std - mean*mean);

  vl->mean = mean;
  vl->std = std;
  
  return;
}
