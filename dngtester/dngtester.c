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

MATRIX *MRIcorVox2RAS(MATRIX *vox2ras);

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
  outsurfname = "white.tal.nonlin";


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
  GCAMmorphSurf(mris, gcam);

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
/*------------------------------------------------------------------------
  MRIcorVox2RAS() - computes vox2ras for the standard COR volume, ie,
  256^3, 1mm^3. The RAS is in "tkregister" space (also known as  "surface"
  space). 
  ------------------------------------------------------------------------*/
MATRIX *MRIcorVox2RAS(MATRIX *vox2ras)
{
  if(vox2ras==NULL) vox2ras = MatrixConstVal(0,4,4,NULL);
  vox2ras->rptr[1][1] = -1;
  vox2ras->rptr[1][4] = +128;
  vox2ras->rptr[2][3] = +1;
  vox2ras->rptr[2][4] = -128;
  vox2ras->rptr[3][2] = -1;
  vox2ras->rptr[3][4] = +128;
  vox2ras->rptr[4][4] = +1;
  return(vox2ras);
}

