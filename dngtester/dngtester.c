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

char *Progname = "dngtester";
char *subject=NULL, *hemi=NULL;
char *SUBJECTS_DIR = NULL;
MRI *mri, *mri2, *template;
char tmpstr[2000];
int err;
GCA  *gca ;
MATRIX *V2Rsrc,*V2Rtemplate,*V2V;
MRIS *mris;
GCSA *gcsa;
char *fsh;

/*----------------------------------------*/
int main(int argc, char **argv)
{

  fsh = getenv("FREESURFER_HOME");
  //set CPAtlas = $FREESURFER_HOME/average/$hemi.$GCS
  //set GCS     = curvature.buckner40.filled.desikan_killiany.gcs

  sprintf(tmpstr,"%s/average/lh.curvature.buckner40.filled.desikan_killiany.gcs",fsh);
  printf("Reading gcsa  %s\n",tmpstr);
  gcsa = GCSAread(tmpstr);

  printf("Reading ico\n");
  mris = ReadIcoByOrder(7, 100);
  printf("Building\n");
  GCSAbuildMostLikelyLabels(gcsa,mris);
  mris->ct = gcsa->ct;
  MRISwriteAnnotation(mris, "./lh.aparc.annot");
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
