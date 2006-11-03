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
#include "evschutils.h"
#include "matrix.h"
#include "matfile.h"
#include "randomfields.h"

// setenv SUBJECTS_DIR /space/greve/1/users/greve/subjects
// /autofs/space/greve_001/users/greve/dev/trunk/dngtester
// ./dngtester register.dat func.mgz func-in-m3z.mgh
// ./dngtester identity.dat ~/subjects/fbirn-anat-101/mri/orig.mgz orig-in-m3z.mgh
// tkregister2 --targ orig-in-m3z.mgh --mov func-in-m3z.mgh --regheader --reg tmp.reg
// Atlas: $SUBJECTS_DIR/avgtst/mri/T1MLV.mgz

MRI *MRIsetSliceNo(MRI *mri, MRI *out);
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
MRI *mri, *mri2, *mask=NULL;
MATRIX *V, *W, *m_tmp;
float ipr, bpr, intensity;
MATRIX *R;
int float2int;
int err,n;
struct timeb  mytimer;
int msecFitTime;
double a,b,zthresh,pthresh;

MRI *fMRIvariance2(MRI *fmri, float DOF, int RmMean, MRI *var);

/*----------------------------------------*/
int main(int argc, char **argv)
{

  printf("Reading\n");
  mri = MRIread(argv[1]);
  printf("chunck %d\n",mri->ischunked);

  mri2 = fMRIvariance(mri, -1, 1, NULL);

  printf("Var1 --------\n");
  TimerStart(&mytimer) ;
  for(n=0; n < 1; n++)  fMRIvariance(mri, -1, 1, mri2);
  msecFitTime = TimerStop(&mytimer) ;
  printf("t = %lf\n",(double)msecFitTime/n);
  MRIwrite(mri2,"var1.mgh");

  printf("Var2 --------\n");
  TimerStart(&mytimer) ;
  for(n=0; n < 1; n++)  fMRIvariance2(mri, -1, 1, mri2);
  msecFitTime = TimerStop(&mytimer) ;
  printf("t = %lf\n",(double)msecFitTime/n);
  MRIwrite(mri2,"var2.mgh");

  printf("Done\n");

  exit(0);


  mri2 = MRIvote(mri, mask, NULL);
  MRIwrite(mri2,"vote.mgh");

  return(0);
  exit(0);
}

/*--------------------------------------------------------*/
MRI *fMRIvariance2(MRI *fmri, float DOF, int RmMean, MRI *var)
{
  int c, r, s, f;
  double val,sumsqval, sumval;
  int nvox_per_row, nvox_per_slice, bytes_per_vol;
  void *p;
  val = 0;

  if(DOF < 0) DOF = fmri->nframes;

  if(var==NULL){
    var = MRIallocSequence(fmri->width, fmri->height, fmri->depth,
                           MRI_FLOAT, 1);
    if(var==NULL){
      printf("ERROR: fMRIvariance: could not alloc\n");
      return(NULL);
    }
    MRIcopyHeader(fmri,var);
  }
  else{
    if(var->width  != fmri->width ||
       var->height != fmri->height ||
       var->depth  != fmri->depth){
      printf("ERROR: fMRIvariance: output dimension mismatch\n");
      return(NULL);
    }
    if(var->type != MRI_FLOAT){
      printf("ERROR: fMRIvariance: structure passed is not MRI_FLOAT\n");
      return(NULL);
    }
  }

  nvox_per_row = fmri->width;
  nvox_per_slice = fmri->width * fmri->height;
  bytes_per_vol = fmri->width * fmri->height * fmri->depth * fmri->bytes_per_vox;
  for(c=0; c < fmri->width; c++){
    for(r=0; r < fmri->height; r++){
      for(s=0; s < fmri->depth; s++){
        sumval = 0;
        sumsqval = 0;
	p = fmri->chunk + (c + r*nvox_per_row + s*nvox_per_slice)*fmri->bytes_per_vox;
        for(f=0; f < fmri->nframes; f++){
	  if(fmri->ischunked){
	    switch(fmri->type){
	    case MRI_UCHAR: val = (double)(*((char *)p)); break;
	    case MRI_SHORT: val = (double)(*((short*)p)); break;
	    case MRI_INT:   val = (double)(*((int  *)p)); break;
	    case MRI_LONG:  val = (double)(*((long *)p)); break;
	    case MRI_FLOAT: val = (double)(*((float*)p)); break;
	    }
	  } 
	  else val = MRIgetVoxVal(fmri, c, r, s, f);
	  //printf("%d  %lu   %g  %g\n",f,(long int)p,val,MRIgetVoxVal(fmri,c,r,s,f));
          sumsqval += (val*val);
          if(RmMean) sumval += val;
	  p += bytes_per_vol;
        }
        MRIFseq_vox(var,c,r,s,0) = sumsqval/DOF;
        if(RmMean)
          MRIFseq_vox(var,c,r,s,0) -= ((sumval/DOF)*(sumval/DOF));
      }
    }
  }

  return(var);
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
/*
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
*/

/*---------------------------------------------------------------*/
MRI *MRIsetSliceNo(MRI *mri, MRI *out)
{
  int c, r, s, f, n, ncols, nrows, nslices,nframes;
  void   *pmri=NULL, *pout=NULL;
  int sz, szout;

  ncols   = mri->width;
  nrows   = mri->height;
  nslices = mri->depth;
  nframes = mri->nframes;

  if(out==NULL){
    out = MRIallocSequence(ncols, nrows, nslices, mri->type, nframes);
    if(out==NULL){
      printf("ERROR: MRIsetSliceNo: could not alloc output\n");
      return(NULL);
    }
    MRIcopyHeader(mri,out); // ordinarily would need to change nframes
  }
  if(out->width != ncols   || out->height != nrows ||
     out->depth != nslices || out->nframes != nframes){
    printf("ERROR: MRIsetSliceNo: dimension mismatch\n");
    return(NULL);
  }

  // Number of bytes in the mri data types
  sz   = MRIsizeof(mri->type);
  szout = MRIsizeof(out->type);

  n = 0;
  for(f=0; f<nframes; f++){
    for(s=0; s<nslices; s++){
      for(r=0; r<nrows; r++){
	// Pointers to the start of the column
	pmri  = (void *) mri->slices[n][r];
	pout  = (void *) out->slices[n][r];
        for(c=0; c<ncols; c++) {
	  MRIdbl2ptr(s, pout, out->type);
	  pmri += sz;
	  pout  += szout;
        } // cols
      } // rows
      n++;
    } // slices
  } // frames

  return(out);
}

/*--------------------------------------------------------*/
double pcluster(double clustersize, double vthresh, double fwhm, 
		double searchsize, int dim)
{
  double p,W,dLh,Em,beta,pvthresh,Pnk;

  W = fwhm/sqrt(4.0*log(2.0));
  dLh = pow(W,-dim);

  Em = searchsize * pow(2*M_PI,-(dim+1)/2.0) * dLh * 
    (pow(vthresh,dim-1.0) - 1) *  exp(-pow(vthresh,2)/2); 

  pvthresh = 1.0-sc_cdf_gaussian_Q(-vthresh,1);

  beta = pow( (gamma(dim/2.0+1)*Em) / (searchsize*pvthresh),2.0/dim );

  // Prob than n >= k, ie, the number of voxels in a cluster >= csize
  Pnk = exp(-beta * pow(clustersize,2.0/dim));
  p = 1 - exp(-Em*Pnk);

  #if 0
  printf("csize = %lf\n",clustersize);
  printf("vthresh = %lf\n",vthresh);
  printf("fwhm = %lf\n",fwhm);
  printf("searchsize = %lf\n",searchsize);
  printf("dim = %d\n",dim);
  printf("W = %lf dLh %lf\n",W,dLh);
  printf("Em = %lf\n",Em);
  printf("pvthresh = %lf\n",pvthresh);
  printf("beta = %lf\n",beta);
  printf("Pnk = %lf\n",Pnk);
  printf("p = %lf\n",p);
  #endif

  return(p);
}
