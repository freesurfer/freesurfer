#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
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
#include "fsenv.h"
#include "qdecutils.h"
#include "dti.h"
#include "registerio.h"

// setenv SUBJECTS_DIR /space/greve/1/users/greve/subjects
// /autofs/space/greve_001/users/greve/dev/trunk/dngtester
// ./dngtester register.dat func.mgz func-in-m3z.mgh
// ./dngtester identity.dat ~/subjects/fbirn-anat-101/mri/orig.mgz orig-in-m3z.mgh
// tkregister2 --targ orig-in-m3z.mgh --mov func-in-m3z.mgh --regheader --reg tmp.reg
// Atlas: $SUBJECTS_DIR/avgtst/mri/T1MLV.mgz


LTA *TransformRegDat2LTA(MRI *targ, MRI *mov, MATRIX *R);

char *Progname = "dngtester";

char *outvolfile=NULL;
char *subject = NULL;
char *regfile = NULL;
TRANSFORM *Rtransform;  //types : M3D, M3Z, LTA, FSLMAT, DAT, OCT(TA), XFM
FSENV     *fsenv;
char gcamfile[1000];
char origfile[1000];
char *sourcefile;
GCAM      *gcam;
MRI       *mri_source;
MRI       *mri_dst; 
MRI       *mri_targ; 
MATRIX *V, *W, *m_tmp;
float ipr, bpr, intensity;
MATRIX *R;
int float2int;

/*----------------------------------------*/
int main(int argc, char **argv)
{
  int err;
    
  if(argc != 4){
    printf("USAGE: dngtester regfile sourcevol outvol\n");
    exit(1);
  }

  regfile = argv[1];
  sourcefile = argv[2];
  outvolfile = argv[3];
  
  // 2) LOAD IN TALAIRACH.M3Z AND READ IN AS A GCAM STRUCTURE 
  fsenv = FSENVgetenv(); 

  printf("Reading in source %s \n",sourcefile);
  mri_source = MRIread(sourcefile);
  if(!mri_source) exit(1);

  printf("Reading in reg \n");
  err = regio_read_register(regfile, &subject, &ipr, &bpr, 
				&intensity, &R, &float2int);
  if(err) exit(1);

  sprintf(gcamfile,"%s/%s/mri/transforms/talairach.m3z",fsenv->SUBJECTS_DIR,subject);
  printf("Reading GCAM %s\n",gcamfile);
  gcam = GCAMread(gcamfile);
 if(gcam == NULL) exit(1);
  printf("Done Reading GCAM\n");

  sprintf(origfile,"%s/%s/mri/orig.mgz",fsenv->SUBJECTS_DIR,subject);
  printf("Reading orig %s\n",origfile);
  mri_targ = MRIread(origfile);
  if(mri_targ == NULL) exit(1);
  printf("Done Reading orig\n");

  printf("Converting to LTA \n");
  Rtransform = (TRANSFORM *)calloc(sizeof(TRANSFORM),1);
  Rtransform->xform = (void *)TransformRegDat2LTA(mri_targ, mri_source, R);

  printf("Applying linear transform to gcam\n");
  GCAMapplyTransform(gcam, Rtransform);  //voxel2voxel
  printf("Done applying transform\n");
 
  // 4) RESAMPLING USING THE FUNCTIONAL IMAGE AS
  //    THE SOURCE OF THE TRANSFORMATION
  // classical functional - precise address ??
   //GCAMrasToVox(gcam, atlas);	// no more useful...
  // because GCAMmorphToAtlas assumes gcam is a Vox2Vox transformation :
  printf("Applying 3d transform to source\n");
  mri_dst = GCAMmorphToAtlas(mri_source, gcam, NULL, -1);
  printf("Done 3d transform to func\n");

  MRIwrite(mri_dst,outvolfile);
  
  
  printf("Finished\n");

  return(0);
  exit(0);
}

/*---------------------------------------------------------------*/
LTA *TransformRegDat2LTA(MRI *targ, MRI *mov, MATRIX *R)
{
  LTA *lta;
  MATRIX *vox2vox; // Targ->Mov
  MATRIX *Ttarg, *Tmov, *invTmov;

  Ttarg = MRIxfmCRS2XYZtkreg(targ);
  Tmov  = MRIxfmCRS2XYZtkreg(mov);
  invTmov = MatrixInverse(Tmov,NULL);

  // vox2vox = invTmov * R * Ttarg
  vox2vox = MatrixMultiply(invTmov,R,NULL);
  MatrixMultiply(vox2vox,Ttarg,vox2vox);

  lta = LTAalloc(1,NULL);
  lta->type = LINEAR_VOX_TO_VOX;
  lta->xforms[0].type = LINEAR_VOX_TO_VOX;
  getVolGeom(targ,&lta->xforms[0].src);
  getVolGeom(mov,&lta->xforms[0].dst);
  lta->xforms[0].m_L = MatrixCopy(vox2vox,NULL);

  MatrixFree(&Ttarg);
  MatrixFree(&Tmov);
  MatrixFree(&invTmov);
  MatrixFree(&vox2vox);

  return(lta);
}

