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
#include "timer.h"

// setenv SUBJECTS_DIR /space/greve/1/users/greve/subjects
// /autofs/space/greve_001/users/greve/dev/trunk/dngtester
// ./dngtester register.dat func.mgz func-in-m3z.mgh
// ./dngtester identity.dat ~/subjects/fbirn-anat-101/mri/orig.mgz orig-in-m3z.mgh
// tkregister2 --targ orig-in-m3z.mgh --mov func-in-m3z.mgh --regheader --reg tmp.reg
// Atlas: $SUBJECTS_DIR/avgtst/mri/T1MLV.mgz

MRI *MRIvote(MRI *in, MRI *mask, MRI *vote);


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
MRI *mri, *mri2, *mask=NULL;
MATRIX *V, *W, *m_tmp;
float ipr, bpr, intensity;
MATRIX *R;
int float2int;
int err;
struct timeb  mytimer;
int msecFitTime;
double a,b;

/*----------------------------------------*/
int main(int argc, char **argv)
{

  mri  = MRIread(argv[1]);
  mask = MRIread(argv[2]);
  mri2 = MRIvote(mri, mask, NULL);
  MRIwrite(mri2,"vote.mgh");

  return(0);
  exit(0);
}

/*---------------------------------------------------------------
  MRIvote() - select the most frequently occuring value measured
  across frames in each voxel. NOT TESTED YET!
  ---------------------------------------------------------------*/
MRI *MRIvote(MRI *in, MRI *mask, MRI *vote)
{
  int c, r, s, f, f0, ncols, nrows, nslices,nframes;
  float m;
  double vmax,v,v0;
  int runlen, runlenmax;
  MRI *sorted;

  if(0 && in->type != MRI_INT && in->type != MRI_SHORT && 
     in->type != MRI_LONG && in->type != MRI_UCHAR){
    printf("ERROR: MRIvote(): input is not of integer class\n");
    return(NULL);
  }

  sorted = MRIsort(in,mask,NULL);
  if(sorted == NULL) return(NULL);

  ncols   = in->width;
  nrows   = in->height;
  nslices = in->depth;
  nframes = in->nframes;

  if(vote==NULL){
    vote = MRIallocSequence(ncols, nrows, nslices, in->type, 1);
    if(vote==NULL){
      printf("ERROR: MRIvote: could not alloc\n");
      return(NULL);
    }
    MRIcopyHeader(in,vote);
    vote->nframes = 1;
  }
  if(in->type != vote->type){
    printf("ERROR: MRIvote: type mismatch\n");
    return(NULL);
  }
  if(vote->width != ncols   || vote->height  != nrows ||
     vote->depth != nslices || vote->nframes != 1){
    printf("ERROR: MRIvote: dimension mismatch\n");
    return(NULL);
  }

  vmax = 0;
  runlenmax = 0;
  for(s=0; s<nslices; s++){
    for(r=0; r<nrows; r++){
      for(c=0; c<ncols; c++) {
	if(mask){
	  m = MRIgetVoxVal(mask,c,r,s,0);
	  if(m < 0.5) continue;
	}
	v0 = MRIgetVoxVal(sorted,c,r,s,0); // value at start of run
	f0 = 0;                            // frame at start of run
	f = 1;
	while(f < nframes){
	  v = MRIgetVoxVal(sorted,c,r,s,f);
	  if(v0 != v){
	    // new value is different than that of run start
	    runlen = f - f0; // runlength for v0
	    if(runlenmax < runlen){
	      runlenmax = runlen;
	      vmax = v0;
	      v0 = v;
	      f0 = f;
	    }
	  }
	  f++;
	}
	// Need to do this one more time in case last value
	// has the longest run
	runlen = f - f0;
	if(runlenmax < runlen){
	  runlenmax = runlen;
	  vmax = v0;
	  v0 = v;
	  f0 = f;
	}
	MRIsetVoxVal(vote,c,r,s,0,vmax);
	// Should probably keep track of max run length
      } // cols
    } // rows
  } // slices

  MRIfree(&sorted);
  return(vote);
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

