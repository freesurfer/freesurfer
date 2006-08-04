#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "mrisurf.h"
#include "utils.h"
#include "error.h"
#include "proto.h"
#include "mri.h"
#include "macros.h"
#include "diag.h"
#include "volume_io.h"
#include "region.h"
#include "machine.h"
#include "fio.h"
#include "mri_identify.h" 
#include "mrisurf.h" 
#include "fmriutils.h" 
#include "gca.h"
#include "gcsa.h"
#include "fsgdf.h"
#include "icosahedron.h"
#include "gca.h"
#include "gcamorph.h"
#include "DICOMRead.h"
#include "cma.h"

/*
  mri_surfmask 
    --m mask.mgz --lh lh.mask.mgz --rh rh.mask.mgz --d dist.mgz
    --s subject 
    --t template.mgz  --r register.dat
    --thresh distthresh

  tkregister2 --reg register.dat --t template.mgz

*/

void convert_surf_to_vox(MRI_SURFACE* mris, MRI* vol);

char *Progname = "dngtester";
char *subject=NULL, *hemi=NULL, *surfname=NULL, *outsurfname=NULL;
char *SUBJECTS_DIR = NULL;
char tmpstr[2000];

int err,nl,n,c,r,s;
char *fsh;
MRI *TempVol=NULL,*TempVol2=NULL;
char *tempvolpath;
MRIS *lhwhite, *rhwhite;
MRIS *lhpial, *rhpial;
VERTEX vtx;
MRI_REGION *region;
int isleft,isright;
double dw,dp;
MRI *lhwhitedist,*lhpialdist,*rhwhitedist,*rhpialdist;
MRI *MaskVol,*LMaskVol,*RMaskVol;

/*----------------------------------------*/
int main(int argc, char **argv)
{
  tempvolpath = argv[1];
  subject = argv[2];
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");

  /* Read in the template volume header */
  // Use orig.mgz if not specified
  TempVol = MRIreadHeader(tempvolpath,MRI_VOLUME_TYPE_UNKNOWN);
  if(TempVol == NULL){
    printf("ERROR: reading %s header\n",tempvolpath);
    exit(1);
  }
  /* ------ Load subject's lh white surface ------ */
  sprintf(tmpstr,"%s/%s/surf/lh.white",SUBJECTS_DIR,subject);
  printf("\nReading lh white surface \n %s\n",tmpstr);
  lhwhite = MRISread(tmpstr);
  if(lhwhite == NULL){
    fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
    exit(1);
  }
  printf("\n");
  /* ------ Load subject's lh thickness ------ */
  sprintf(tmpstr,"%s/%s/surf/lh.thickness",SUBJECTS_DIR,subject);
  printf("Reading thickness %s\n",tmpstr);
  err = MRISreadCurvatureFile(lhwhite, tmpstr);
  if(err) exit(1);
  /* ------ Load subject's lh pial surface ------ */
  sprintf(tmpstr,"%s/%s/surf/lh.pial",SUBJECTS_DIR,subject);
  printf("\nReading lh pial surface \n %s\n",tmpstr);
  lhpial = MRISread(tmpstr);
  if(lhpial == NULL){
    fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
    exit(1);
  }
  printf("\n");
  /* Check that they have the same number of vertices */
  if(lhwhite->nvertices != lhpial->nvertices){
    printf("ERROR: lh white and pial have a different number of vertices (%d,%d)\n",
	   lhwhite->nvertices,lhpial->nvertices);
    exit(1);
  }

  /* ------ Load subject's rh white surface ------ */
  sprintf(tmpstr,"%s/%s/surf/rh.white",SUBJECTS_DIR,subject);
  printf("\nReading rh white surface \n %s\n",tmpstr);
  rhwhite = MRISread(tmpstr);
  if(rhwhite == NULL){
    fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
    exit(1);
  }
  printf("\n");
  /* ------ Load subject's rh thickness ------ */
  sprintf(tmpstr,"%s/%s/surf/rh.thickness",SUBJECTS_DIR,subject);
  printf("Reading thickness %s\n",tmpstr);
  err = MRISreadCurvatureFile(rhwhite, tmpstr);
  if(err) exit(1);
  /* ------ Load subject's rh pial surface ------ */
  sprintf(tmpstr,"%s/%s/surf/rh.pial",SUBJECTS_DIR,subject);
  printf("\nReading rh pial surface \n %s\n",tmpstr);
  rhpial = MRISread(tmpstr);
  if(rhpial == NULL){
    fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
    exit(1);
  }
  printf("\n");
  /* Check that they have the same number of vertices */
  if(rhwhite->nvertices != rhpial->nvertices){
    printf("ERROR: rh white and pial have a different number of vertices (%d,%d)\n",
	   rhwhite->nvertices,rhpial->nvertices);
    exit(1);
  }

  convert_surf_to_vox(lhpial, TempVol);
  convert_surf_to_vox(lhwhite,TempVol);
  MRIScomputeMetricProperties(lhpial);
  MRIScomputeMetricProperties(lhwhite);
  MRISsaveVertexPositions(lhwhite,ORIGINAL_VERTICES);
  MRISsaveVertexPositions(lhpial,ORIGINAL_VERTICES);

  convert_surf_to_vox(rhpial, TempVol);
  convert_surf_to_vox(rhwhite,TempVol);
  MRIScomputeMetricProperties(rhpial);
  MRIScomputeMetricProperties(rhwhite);
  MRISsaveVertexPositions(rhwhite,ORIGINAL_VERTICES);
  MRISsaveVertexPositions(rhpial,ORIGINAL_VERTICES);

  region = REGIONalloc();
  region->x = 0;
  region->y = 0;
  region->z = 0;
  region->dx = TempVol->width;
  region->dy = TempVol->height;
  region->dz = TempVol->depth;

  lhwhitedist = MRISbinarizeVolume(lhwhite,region,1.0,5.0);
  lhpialdist  = MRISbinarizeVolume(lhpial,region,1.0,5.0);
  rhwhitedist = MRISbinarizeVolume(rhwhite,region,1.0,5.0);
  rhpialdist  = MRISbinarizeVolume(rhpial,region,1.0,5.0);
  MRIwrite(lhwhitedist,"lh.dwhite.mgh");
  MRIwrite(lhpialdist,"lh.dpial.mgh");
  MRIwrite(rhwhitedist,"rh.dwhite.mgh");
  MRIwrite(rhpialdist,"rh.dpial.mgh");

  MaskVol = MRIallocSequence(TempVol->width,TempVol->height,TempVol->depth,MRI_INT,1);
  MRIcopyHeader(TempVol,MaskVol);
  LMaskVol = MRIallocSequence(TempVol->width,TempVol->height,TempVol->depth,MRI_INT,1);
  MRIcopyHeader(TempVol,LMaskVol);
  RMaskVol = MRIallocSequence(TempVol->width,TempVol->height,TempVol->depth,MRI_INT,1);
  MRIcopyHeader(TempVol,RMaskVol);

  for(c=0; c < TempVol->width; ++c){
    printf("%3d ",c); fflush(stdout);
    if(c%20 == 19) printf("\n");
    for(r=0; r < TempVol->height; ++r){
      for(s=0; s < TempVol->depth; ++s){

	dw = MRIgetVoxVal(lhwhitedist,c,r,s,0);
	dp = MRIgetVoxVal(lhpialdist, c,r,s,0);
	if(dw*dp < 0 && dw != 1000 && dp != 1000) isleft = 1;
	else                                      isleft = 0;

	dw = MRIgetVoxVal(rhwhitedist,c,r,s,0);
	dp = MRIgetVoxVal(rhpialdist, c,r,s,0);
	if(dw*dp < 0 && dw != 1000 && dp != 1000) isright = 1;
	else                                      isright = 0;

	if(isleft){
	  MRIsetVoxVal(MaskVol,c,r,s,0,1);
	  MRIsetVoxVal(LMaskVol,c,r,s,0,1);
	}
	if(isright){
	  MRIsetVoxVal(MaskVol,c,r,s,0,1);
	  MRIsetVoxVal(LMaskVol,c,r,s,0,1);
	}
      }
    }
  }
  printf("\n");

  MRIwrite(MaskVol,"mask.mgh");
  MRIwrite(LMaskVol,"lh.mask.mgh");
  MRIwrite(RMaskVol,"rh.mask.mgh");

  return(0);
  exit(0);
}

/*-----------------------------------------------------
  convert_surf_to_vox() - replace surface xyz with 
  values in volume index space.
  -----------------------------------------------------*/
void convert_surf_to_vox(MRI_SURFACE* mris, MRI* vol)
{
 double cx, cy, cz;
 Real vx, vy, vz;
 VERTEX* pvtx = &( mris->vertices[0] );
 unsigned int nvertices = (unsigned int)mris->nvertices;
 unsigned int ui;

 for(ui=0;ui < nvertices; ++ui, ++pvtx ) {
   cx = pvtx->x;
   cy = pvtx->y;
   cz = pvtx->z;
   MRIsurfaceRASToVoxel(vol, cx, cy, cz, &vx, &vy, &vz);
   pvtx->x = vx;
   pvtx->y = vy;
   pvtx->z = vz;
 } // next ui, pvtx
 return;
}

