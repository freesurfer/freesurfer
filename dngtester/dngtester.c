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

int GCAMsampleInverseMorphRAS(GCA_MORPH *gcam, 
			      float  xAnat,  float  yAnat,  float  zAnat, 
			      float *xMorph, float *yMorph, float *zMorph,
			      MATRIX *vox2ras);
int GCAMsampleMorphRAS(GCA_MORPH *gcam, float xMorph, float yMorph, float zMorph,
		       float  *xAnat,  float  *yAnat,  float  *zAnat, 
		       MATRIX *vox2ras);

int GCAMmorphSurfToAtlas(MRIS *mris, GCA_MORPH *gcam);

char *Progname = "dngtester";
char *subject=NULL, *hemi=NULL, *surfname=NULL, *outsurfname=NULL;
char *SUBJECTS_DIR = NULL;
MRI *mri, *mri2, *template;
char tmpstr[2000];
int err;
GCA  *gca ;
MATRIX *V2Rsrc,*V2Rtemplate,*V2V;
MRIS *mris;
GCSA *gcsa;
char *fsh;
GCA_MORPH *gcam;
float c,r,s;

/*----------------------------------------*/
int main(int argc, char **argv)
{
  
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  subject     = argv[1];
  hemi        = argv[2];
  surfname    = argv[3];
  outsurfname = argv[4];
  subject = "fsr-tst";
  hemi = "lh";
  surfname = "white";
  outsurfname = "white.tal.nonlin2";


  sprintf(tmpstr,"%s/%s/mri/transforms/talairach.m3z",SUBJECTS_DIR,subject);
  printf("Reading %s\n",tmpstr);
  gcam = GCAMreadAndInvert(tmpstr);
  if(gcam == NULL) exit(1);
  Gdiag_no = 1;


  if(!strcmp(surfname,outsurfname)){
    printf("ERROR: input and output surfaces are the same\n");
    exit(1);
  }

  printf("subj=%s, hemi=%s, surf=%s, out=%s\n",subject,hemi,surfname,outsurfname);
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");

  sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,subject,hemi,surfname);
  printf("Reading %s\n",tmpstr);
  mris = MRISread(tmpstr);
  if(mris == NULL) exit(1);

  printf("Morphing\n");
  GCAMmorphSurfToAtlas(mris, gcam);

  sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,subject,hemi,outsurfname);
  printf("Writing to %s\n",tmpstr);
  MRISwrite(mris,tmpstr);

  exit(0);
  /*----------------------------------------------------------------*/

  fsh = getenv("FREESURFER_HOME");
  //set CPAtlas = $FREESURFER_HOME/average/$hemi.$GCS
  //set GCS     = curvature.buckner40.filled.desikan_killiany.gcs


  printf("Reading ico\n");
  mris = ReadIcoByOrder(7, 100);

  sprintf(tmpstr,"%s/average/%s.curvature.buckner40.filled.desikan_killiany.gcs",fsh,hemi);
  printf("Reading gcsa  %s\n",tmpstr);
  gcsa = GCSAread(tmpstr);
  mris->ct = gcsa->ct;

  printf("Building\n");
  GCSAbuildMostLikelyLabels(gcsa,mris);

  MRISmodeFilterAnnotations(mris, 2);
  sprintf(tmpstr,"./%s.aparc.annot",hemi);
  MRISwriteAnnotation(mris, tmpstr);



  exit(1);
  /*-------------------------------------------------*/

  //set GCA = RB40_talairach_2005_12_30.gca
  //$FREESURFER_HOME/average/$GCA
  printf("Reading %s\n",argv[2]) ;
  template = MRIreadHeader(argv[2],MRI_VOLUME_TYPE_UNKNOWN);
  if(!template) exit(1);

  printf("Reading %s\n",argv[1]) ;
  gca = GCAread(argv[1]) ;
  printf("Building\n");
  mri = GCAbuildMostLikelyVolume(gca, NULL) ;

  printf("Upsampling\n");
  V2Rsrc = MRIxfmCRS2XYZtkreg(mri);
  V2Rtemplate = MRIxfmCRS2XYZtkreg(template);
  V2V = MatrixMultiply(MatrixInverse(V2Rsrc,NULL),V2Rtemplate,NULL);
  MatrixPrint(stdout,V2V);

  mri2 = MRIallocSequence(template->width, template->height,template->depth,
                             MRI_UCHAR, 1);
  MRIcopyHeader(template, mri2) ;
  template->nframes = 1;
  MRIvol2Vol(mri, mri2, NULL, SAMPLE_NEAREST, 0);

  printf("Writings\n");
  MRIwrite(mri,"T1MLV.lowres.mgz");
  MRIwrite(mri2,"T1MLV.mgz");

  MRIfree(&mri);
  MRIfree(&mri2);

  mri = GCAbuildMostLikelyLabelVolume(gca);
  mri2 = MRIallocSequence(template->width, template->height,template->depth,
                             MRI_UCHAR, 1);
  MRIcopyHeader(template, mri2) ;
  template->nframes = 1;
  MRIvol2Vol(mri, mri2, NULL, SAMPLE_NEAREST, 0);

  MRIwrite(mri,"aseg.mlv.lowres.mgz");
  MRIwrite(mri2,"aseg.mlv.mgz");

  return(0);
}
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
int GCAMmorphSurfToAtlas(MRIS *mris, GCA_MORPH *gcam)
{
  int vtxno;
  VERTEX *v;
  float Mx, My, Mz;
  MATRIX *vox2ras;

  vox2ras = MRIxfmCRS2XYZtkreg(gcam->mri_xind);
  printf("Appling Inverse Morph \n");
  for(vtxno = 0; vtxno < mris->nvertices; vtxno++){
    v = &(mris->vertices[vtxno]);
    GCAMsampleInverseMorphRAS(gcam, v->x, v->y, v->z, &Mx, &My, &Mz, vox2ras);
    // pack it back into the vertex
    v->x = Mx;
    v->y = My;
    v->z = Mz;
  }
  return(0);
}
/*----------------------------------------------------------------------------
  GCAMsampleInverseMorphRAS() - 
  ----------------------------------------------------------------------------*/
int GCAMsampleInverseMorphRAS(GCA_MORPH *gcam, 
			      float  xAnat,  float  yAnat,  float  zAnat, 
			      float *xMorph, float *yMorph, float *zMorph,
			      MATRIX *vox2ras)
{
  static MATRIX *ras=NULL, *crs=NULL, *ras2vox=NULL;
  float  cMorph, rMorph, sMorph;
  float  cAnat, rAnat, sAnat;

  ras2vox = MatrixInverse(vox2ras,ras2vox);
  if(!ras){
    ras = MatrixAlloc(4,1,MATRIX_REAL);
    ras->rptr[4][1] = 1;
  }
  if(!crs){
    crs= MatrixAlloc(4,1,MATRIX_REAL);
    crs->rptr[4][1] = 1;
  }

  ras->rptr[1][1] = xAnat;
  ras->rptr[2][1] = yAnat;
  ras->rptr[3][1] = zAnat;
  crs = MatrixMultiply(ras2vox,ras,crs);
  cAnat = crs->rptr[1][1];
  rAnat = crs->rptr[2][1];
  sAnat = crs->rptr[3][1];
  
  err = GCAMsampleInverseMorph(gcam, cAnat, rAnat, sAnat, &cMorph, &rMorph, &sMorph);
  if(err) return(1);

  crs->rptr[1][1] = cMorph;
  crs->rptr[2][1] = rMorph;
  crs->rptr[3][1] = sMorph;
  ras = MatrixMultiply(vox2ras,crs,ras);

  *xMorph = ras->rptr[1][1];
  *yMorph = ras->rptr[2][1];
  *zMorph = ras->rptr[3][1];

  return(0);
}


/*----------------------------------------------------------------------------
  GCAMsampleMorphRAS() - given an RAS coord in the Morph space (xyzMorph),
  compute the RAS in the Anat space (xyzAnat). Has not been tested yet.
  ----------------------------------------------------------------------------*/
int GCAMsampleMorphRAS(GCA_MORPH *gcam, float xMorph, float yMorph, float zMorph,
		       float  *xAnat,  float  *yAnat,  float  *zAnat, 
		       MATRIX *vox2ras)
{
  static MATRIX *ras=NULL, *crs=NULL, *ras2vox=NULL;
  float  cMorph, rMorph, sMorph;
  float  cAnat, rAnat, sAnat;

  ras2vox = MatrixInverse(vox2ras,ras2vox);
  if(!ras){
    ras = MatrixAlloc(4,1,MATRIX_REAL);
    ras->rptr[4][1] = 1;
  }
  if(!crs){
    crs= MatrixAlloc(4,1,MATRIX_REAL);
    crs->rptr[4][1] = 1;
  }

  ras->rptr[1][1] = xMorph;
  ras->rptr[2][1] = yMorph;
  ras->rptr[3][1] = zMorph;
  crs = MatrixMultiply(ras2vox,ras,crs);
  cMorph = crs->rptr[1][1];
  rMorph = crs->rptr[2][1];
  sMorph = crs->rptr[3][1];
  
  err = GCAMsampleMorph(gcam, cMorph, rMorph, sMorph, &cAnat, &rAnat, &sAnat);
  if(err) return(1);

  crs->rptr[1][1] = cAnat;
  crs->rptr[2][1] = rAnat;
  crs->rptr[3][1] = sAnat;
  ras = MatrixMultiply(vox2ras,crs,ras);

  *xAnat = ras->rptr[1][1];
  *yAnat = ras->rptr[2][1];
  *zAnat = ras->rptr[3][1];

  return(0);
}



/*----------------------------------------------------------------------------*/


