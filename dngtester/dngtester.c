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

char *Progname = "dngtester";
char *subject=NULL, *hemi=NULL, *surfname=NULL, *outsurfname=NULL;
char *SUBJECTS_DIR = NULL;
char tmpstr[2000];
MRI *mri, *mri2, *template;
int err,nl,n,c,r,s;
char *fsh;
char *dcmfile, *dcmdir;
GCA *gca;
float hashres = 16;
MRI *TempVol=NULL;
char *tempvolpath;
MRIS *lhwhite, *rhwhite;
MRIS *lhpial, *rhpial;
MHT *lhwhite_hash, *rhwhite_hash;
MHT *lhpial_hash, *rhpial_hash;
MRI *DistVol=NULL,*MaskVol,*LMaskVol,*RMaskVol;
VERTEX vtx;
int  lhwvtx, lhpvtx, rhwvtx, rhpvtx;
int  lhwvtx2, lhpvtx2, rhwvtx2, rhpvtx2;
MATRIX *Vox2RAS, *CRS, *RAS;
float dlhw, dlhp, drhw, drhp;
float dlhw2, dlhp2, drhw2, drhp2;
float dthresh = 5.0;

/*
  mri_surfmask 
    --m mask.mgz --lh lh.mask.mgz --rh rh.mask.mgz --d dist.mgz
    --s subject 
    --t template.mgz  --r register.dat
    --thresh distthresh

  tkregister2 --reg register.dat --t template.mgz

*/



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
    printf("mri_surf2vol ERROR: reading %s header\n",tempvolpath);
    exit(1);
  }
  Vox2RAS = MRIxfmCRS2XYZtkreg(TempVol);
  printf("Vox2RAS (tkreg)\n");
  MatrixPrint(stdout,Vox2RAS);
  printf("\n");

  /* ------ Load subject's lh white surface ------ */
  sprintf(tmpstr,"%s/%s/surf/lh.white",SUBJECTS_DIR,subject);
  printf("\nReading lh white surface \n %s\n",tmpstr);
  lhwhite = MRISread(tmpstr);
  if(lhwhite == NULL){
    fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
    exit(1);
  }
  printf("\n");
  printf("Building hash of lh white\n");
  lhwhite_hash = MHTfillVertexTableRes(lhwhite, NULL,CURRENT_VERTICES,hashres);
  /* ------ Load subject's lh pial surface ------ */
  sprintf(tmpstr,"%s/%s/surf/lh.pial",SUBJECTS_DIR,subject);
  printf("\nReading lh pial surface \n %s\n",tmpstr);
  lhpial = MRISread(tmpstr);
  if(lhpial == NULL){
    fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
    exit(1);
  }
  printf("\n");
  printf("Building hash of lh pial\n");
  lhpial_hash = MHTfillVertexTableRes(lhpial, NULL,CURRENT_VERTICES,hashres);
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
  printf("Building hash of rh white\n");
  rhwhite_hash = MHTfillVertexTableRes(rhwhite, NULL,CURRENT_VERTICES,hashres);
  /* ------ Load subject's rh pial surface ------ */
  sprintf(tmpstr,"%s/%s/surf/rh.pial",SUBJECTS_DIR,subject);
  printf("\nReading rh pial surface \n %s\n",tmpstr);
  rhpial = MRISread(tmpstr);
  if(rhpial == NULL){
    fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
    exit(1);
  }
  printf("\n");
  printf("Building hash of rh pial\n");
  rhpial_hash = MHTfillVertexTableRes(rhpial, NULL,CURRENT_VERTICES,hashres);
  /* Check that they have the same number of vertices */
  if(rhwhite->nvertices != rhpial->nvertices){
    printf("ERROR: rh white and pial have a different number of vertices (%d,%d)\n",
	   rhwhite->nvertices,rhpial->nvertices);
    exit(1);
  }

  // Distance from each voxel to each of the four surfaces
  DistVol = MRIallocSequence(TempVol->width,TempVol->height,TempVol->depth,MRI_FLOAT,4);
  MRIcopyHeader(TempVol,DistVol);

  MaskVol = MRIallocSequence(TempVol->width,TempVol->height,TempVol->depth,MRI_INT,1);
  MRIcopyHeader(TempVol,MaskVol);

  LMaskVol = MRIallocSequence(TempVol->width,TempVol->height,TempVol->depth,MRI_INT,1);
  MRIcopyHeader(TempVol,LMaskVol);

  RMaskVol = MRIallocSequence(TempVol->width,TempVol->height,TempVol->depth,MRI_INT,1);
  MRIcopyHeader(TempVol,RMaskVol);

  CRS = MatrixAlloc(4,1,MATRIX_REAL);
  CRS->rptr[4][1] = 1;
  RAS = MatrixAlloc(4,1,MATRIX_REAL);
  RAS->rptr[4][1] = 1;

  for(c=0; c < TempVol->width; c++){
    printf("%3d ",c); fflush(stdout);
    if(c%20 == 19) printf("\n");
    for(r=0; r < TempVol->height; r++){
      for(s=0; s < TempVol->depth; s++){
	// Convert the CRS to RAS
	CRS->rptr[1][1] = c;
	CRS->rptr[2][1] = r;
	CRS->rptr[3][1] = s;
	RAS = MatrixMultiply(Vox2RAS,CRS,RAS);
	vtx.x = RAS->rptr[1][1];
	vtx.y = RAS->rptr[2][1];
	vtx.z = RAS->rptr[3][1];

	//MRIvoxelToSurfaceRAS(TempVol, c,r,s, &xs, &ys, &zs);
	//printf("[%g,%g,%g] [%g,%g,%g] \n",vtx.x,vtx.y,vtx.z,xs, ys, zs);

	// Get the index of the closest vertex in the 
	// lh.white, lh.pial, rh.white, rh.pial
	lhwvtx = MHTfindClosestVertexNo(lhwhite_hash,lhwhite,&vtx,&dlhw);
	//printf("%d %d %d lhwvtx %d\n",c,r,s,lhwvtx);
	lhpvtx = MHTfindClosestVertexNo(lhpial_hash, lhpial, &vtx,&dlhp);
	rhwvtx = MHTfindClosestVertexNo(rhwhite_hash,rhwhite,&vtx,&drhw);
	rhpvtx = MHTfindClosestVertexNo(rhpial_hash, rhpial, &vtx,&drhp);

	if(lhwvtx < 0) dlhw = -2*dthresh;
	if(lhpvtx < 0) dlhp = -2*dthresh;
	if(rhwvtx < 0) drhw = -2*dthresh;
	if(rhpvtx < 0) drhp = -2*dthresh;

	if(0 && lhwvtx >= 0){
	  vtx = lhwhite->vertices[lhwvtx];
	  printf("lhw:  %d  a=[%g %g %g]; b=[%g %g,%g]; %g [%d,%d,%d]\n",
		 lhwvtx,vtx.x,vtx.y,vtx.z,
		 RAS->rptr[1][1],RAS->rptr[2][1],RAS->rptr[3][1],
		 dlhw,c,r,s);
	  lhwvtx2 = MRISfindClosestVertex(lhwhite, vtx.x, vtx.y, vtx.z,&dlhw2);
	  vtx = lhwhite->vertices[lhwvtx2];
	  printf("lhw2: %d  a=[%g %g %g]; b=[%g %g,%g]; %g [%d,%d,%d]\n",
		 lhwvtx2,vtx.x,vtx.y,vtx.z,
		 RAS->rptr[1][1],RAS->rptr[2][1],RAS->rptr[3][1],
		 dlhw2,c,r,s);
	  if(lhwvtx != lhwvtx2) exit(1);
	  //exit(1);
	  fflush(stdout);
	  //lhpvtx2 = MRISfindClosestVertex(lhpial,  vtx.x, vtx.y, vtx.z,&dlhp2);
	  //rhwvtx2 = MRISfindClosestVertex(rhwhite, vtx.x, vtx.y, vtx.z,&drhw2);
	  //rhpvtx2 = MRISfindClosestVertex(rhpial,  vtx.x, vtx.y, vtx.z,&drhp2);
	}

	MRIsetVoxVal(DistVol,c,r,s,0,dlhw);
	MRIsetVoxVal(DistVol,c,r,s,1,dlhp);
	MRIsetVoxVal(DistVol,c,r,s,2,drhw);
	MRIsetVoxVal(DistVol,c,r,s,3,drhp);

	if(fabs(dlhw) < fabs(dthresh) && fabs(dlhp) < dthresh) 
	  MRIsetVoxVal(LMaskVol,c,r,s,0,1);
	if(fabs(drhw) < fabs(dthresh) && fabs(drhp) < dthresh) 
	  MRIsetVoxVal(RMaskVol,c,r,s,0,1);

	if((fabs(dlhw) < dthresh && fabs(dlhp) < dthresh) || 
	   (fabs(drhw) < dthresh && fabs(drhp) < dthresh))
	  MRIsetVoxVal(MaskVol,c,r,s,0,1);


      } // s
    } // r 
  } // c
  printf("\n");

  MRIwrite(DistVol,"dist.mgh");
  MRIwrite(MaskVol,"mask.mgh");
  MRIwrite(LMaskVol,"lh.mask.mgh");
  MRIwrite(RMaskVol,"rh.mask.mgh");

  return(0);
  exit(0);
}
