/**
 * @file  dngtester.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2007/04/06 06:17:04 $
 *    $Revision: 1.29 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


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
#include "mri2.h"
#include "annotation.h"

// setenv SUBJECTS_DIR /space/greve/1/users/greve/subjects
// /autofs/space/greve_001/users/greve/dev/trunk/dngtester
// ./dngtester register.dat func.mgz func-in-m3z.mgh
// ./dngtester identity.dat ~/subjects/fbirn-anat-101/mri/orig.mgz orig-in-m3z.mgh
// tkregister2 --targ orig-in-m3z.mgh --mov func-in-m3z.mgh --regheader --reg tmp.reg
// Atlas: $SUBJECTS_DIR/avgtst/mri/T1MLV.mgz

MRI *MRIsetSliceNo(MRI *mri, MRI *out);

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
char *SUBJECTS_DIR;
char tmpstr[2000];
MRIS *surf;

MRI *fMRIvariance2(MRI *fmri, float DOF, int RmMean, MRI *var);
void printrgb(void);
MRI *GetMyMask(MRIS *surf);
MRI *MRISdilateMask(MRIS *surf, MRI *mask, int annotidmask, int niters);


/*----------------------------------------*/
int main(int argc, char **argv) {
  int err;
  int *nunits;

  subject = "tl-wm";
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  sprintf(tmpstr,"%s/%s/surf/lh.white",SUBJECTS_DIR,subject);
  printf("\nReading lh white surface \n %s\n",tmpstr);
  surf = MRISread(tmpstr);
  if(!surf) exit(1);

  sprintf(tmpstr,"%s/%s/label/lh.aparc.annot",SUBJECTS_DIR,subject);
  printf("Loading annotations from %s\n",tmpstr);
  err = MRISreadAnnotation(surf, tmpstr);
  if(err) exit(1);
  mri = GetMyMask(surf);
  MRISdilateMask(surf, mri, 28, 20);

  nunits = (int *)calloc(surf->ct->nentries, sizeof(int)) ;
  nunits[28] = 5;
  MRISdivideAnnotation(surf, nunits) ;
  MRISwriteAnnotation(surf, "lh.fbirn.annot") ;

  return(0);

  printrgb();
  return(0);

  printf("Reading\n");
  mri = MRIread(argv[1]);
  printf("chunck %d\n",mri->ischunked);

  mri2 = fMRIvariance(mri, -1, 1, NULL);

  printf("Var1 --------\n");
  TimerStart(&mytimer) ;
  for (n=0; n < 1; n++)  fMRIvariance(mri, -1, 1, mri2);
  msecFitTime = TimerStop(&mytimer) ;
  printf("t = %lf\n",(double)msecFitTime/n);
  MRIwrite(mri2,"var1.mgh");

  printf("Var2 --------\n");
  TimerStart(&mytimer) ;
  for (n=0; n < 1; n++)  fMRIvariance2(mri, -1, 1, mri2);
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
MRI *fMRIvariance2(MRI *fmri, float DOF, int RmMean, MRI *var) {
  int c, r, s, f;
  double val,sumsqval, sumval;
  int nvox_per_row, nvox_per_slice, bytes_per_vol;
  void *p;
  val = 0;

  if (DOF < 0) DOF = fmri->nframes;

  if (var==NULL) {
    var = MRIallocSequence(fmri->width, fmri->height, fmri->depth,
                           MRI_FLOAT, 1);
    if (var==NULL) {
      printf("ERROR: fMRIvariance: could not alloc\n");
      return(NULL);
    }
    MRIcopyHeader(fmri,var);
  } else {
    if (var->width  != fmri->width ||
        var->height != fmri->height ||
        var->depth  != fmri->depth) {
      printf("ERROR: fMRIvariance: output dimension mismatch\n");
      return(NULL);
    }
    if (var->type != MRI_FLOAT) {
      printf("ERROR: fMRIvariance: structure passed is not MRI_FLOAT\n");
      return(NULL);
    }
  }

  nvox_per_row = fmri->width;
  nvox_per_slice = fmri->width * fmri->height;
  bytes_per_vol = fmri->width * fmri->height * fmri->depth * fmri->bytes_per_vox;
  for (c=0; c < fmri->width; c++) {
    for (r=0; r < fmri->height; r++) {
      for (s=0; s < fmri->depth; s++) {
        sumval = 0;
        sumsqval = 0;
        p = fmri->chunk + (c + r*nvox_per_row + s*nvox_per_slice)*fmri->bytes_per_vox;
        for (f=0; f < fmri->nframes; f++) {
          if (fmri->ischunked) {
            switch (fmri->type) {
            case MRI_UCHAR:
              val = (double)(*((char *)p));
              break;
            case MRI_SHORT:
              val = (double)(*((short*)p));
              break;
            case MRI_INT:
              val = (double)(*((int  *)p));
              break;
            case MRI_LONG:
              val = (double)(*((long *)p));
              break;
            case MRI_FLOAT:
              val = (double)(*((float*)p));
              break;
            }
          } else val = MRIgetVoxVal(fmri, c, r, s, f);
          //printf("%d  %lu   %g  %g\n",f,(long int)p,val,MRIgetVoxVal(fmri,c,r,s,f));
          sumsqval += (val*val);
          if (RmMean) sumval += val;
          p += bytes_per_vol;
        }
        MRIFseq_vox(var,c,r,s,0) = sumsqval/DOF;
        if (RmMean)
          MRIFseq_vox(var,c,r,s,0) -= ((sumval/DOF)*(sumval/DOF));
      }
    }
  }

  return(var);
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
MRI *MRIsetSliceNo(MRI *mri, MRI *out) {
  int c, r, s, f, n, ncols, nrows, nslices,nframes;
  void   *pmri=NULL, *pout=NULL;
  int sz, szout;

  ncols   = mri->width;
  nrows   = mri->height;
  nslices = mri->depth;
  nframes = mri->nframes;

  if (out==NULL) {
    out = MRIallocSequence(ncols, nrows, nslices, mri->type, nframes);
    if (out==NULL) {
      printf("ERROR: MRIsetSliceNo: could not alloc output\n");
      return(NULL);
    }
    MRIcopyHeader(mri,out); // ordinarily would need to change nframes
  }
  if (out->width != ncols   || out->height != nrows ||
      out->depth != nslices || out->nframes != nframes) {
    printf("ERROR: MRIsetSliceNo: dimension mismatch\n");
    return(NULL);
  }

  // Number of bytes in the mri data types
  sz   = MRIsizeof(mri->type);
  szout = MRIsizeof(out->type);

  n = 0;
  for (f=0; f<nframes; f++) {
    for (s=0; s<nslices; s++) {
      for (r=0; r<nrows; r++) {
        // Pointers to the start of the column
        pmri  = (void *) mri->slices[n][r];
        pout  = (void *) out->slices[n][r];
        for (c=0; c<ncols; c++) {
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
                double searchsize, int dim) {
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

void printrgb(void)
{
  double f,r,g,b;

  for(f=0; f<64; f++){
    if(f < 25)      r =  0.0;
    else if(f < 40) r = -1.5 + .065*f;
    else if(f < 57) r =  1.0;
    else            r =  4.5 - .065*f;

    if(f < 9)       g =  0.0;
    else if(f < 23) g = -0.5 + .065*f;
    else if(f < 41) g =  1.0;
    else if(f < 55) g =  3.5 - .065*f;
    else            g = 0;

    if(f < 7)       b =  0.5 + .065*f;
    else if(f < 25) b =  1.0;
    else if(f < 40) b =  2.5 - .065*f;
    else            b =  0.0;

    printf("%6.4f %6.4f %6.4f %6.4f\n",f,r,g,b);

  }
}
/*-------------------------------------------------------------
  MRI *GetMyMask(MRIS *surf) - creates a mask at the intersection
  of SFG and CAcing, PAcing, and RAcing.
  -------------------------------------------------------------*/
MRI *GetMyMask(MRIS *surf)
{
  int vtxno, annot, annotid, nnbrs, nbrvtx, nthnbr, nbrannotid; 
  VERTEX *vtx;
  int superiorfrontal, posteriorcingulate;
  int caudalanteriorcingulate, rostralanteriorcingulate;
  MRI *mri;

  superiorfrontal = 28;
  posteriorcingulate = 23;
  caudalanteriorcingulate = 2;
  rostralanteriorcingulate = 26;

  mri = MRIalloc(surf->nvertices, 1, 1, MRI_INT) ;

  for (vtxno = 0; vtxno < surf->nvertices; vtxno++) {

    vtx = &(surf->vertices[vtxno]);
    annot = surf->vertices[vtxno].annotation;
    CTABfindAnnotation(surf->ct, annot, &annotid);

    if(annotid != superiorfrontal && 
       annotid != posteriorcingulate &&
       annotid != caudalanteriorcingulate && 
       annotid != rostralanteriorcingulate) continue;

    nnbrs = surf->vertices[vtxno].vnum;
    for (nthnbr = 0; nthnbr < nnbrs; nthnbr++){
      nbrvtx = surf->vertices[vtxno].v[nthnbr];
      if (surf->vertices[nbrvtx].ripflag) continue; //skip ripped vtxs
      annot = surf->vertices[nbrvtx].annotation;
      CTABfindAnnotation(surf->ct, annot, &nbrannotid);
      if(annotid == superiorfrontal && 
	 (nbrannotid == posteriorcingulate ||
	  nbrannotid == caudalanteriorcingulate ||
	  nbrannotid == rostralanteriorcingulate) ){
	MRIsetVoxVal(mri,vtxno,0,0,0,1);
      }
    }
  }

  MRIwrite(mri,"mymask.mgh");

  return(mri);
}

MRI *MRISdilateMask(MRIS *surf, MRI *mask, int annotidmask, int niters)
{
  int vtxno, annot, annotid, nnbrs, nbrvtxno, nthnbr, nthiter;
  VERTEX *vtx;
  MRI *mri1, *mri2;
  float val;

  //superiorfrontal = 28;

  mri1 = MRIcopy(mask,NULL);
  mri2 = MRIcopy(mask,NULL);

  for(nthiter = 0; nthiter < niters; nthiter++){
    printf("iter %d\n",nthiter);

    for (vtxno = 0; vtxno < surf->nvertices; vtxno++) {
      vtx = &(surf->vertices[vtxno]);

      if(annotidmask > -1){
	annot = surf->vertices[vtxno].annotation;
	CTABfindAnnotation(surf->ct, annot, &annotid);
	if(annotid != annotidmask){
	  MRIsetVoxVal(mri2,vtxno,0,0,0,0);
	  continue;
	}
      }

      val = MRIgetVoxVal(mri1,vtxno,0,0,0);

      if(val){
	MRIsetVoxVal(mri2,vtxno,0,0,0,1);
	continue;
      }
      // If it gets here, the vtx is in the annot and has not been set 
      nnbrs = surf->vertices[vtxno].vnum;
      for (nthnbr = 0; nthnbr < nnbrs; nthnbr++){
	nbrvtxno = surf->vertices[vtxno].v[nthnbr];
	if (surf->vertices[nbrvtxno].ripflag) continue; //skip ripped vtxs
	val = MRIgetVoxVal(mri1,nbrvtxno,0,0,0);
	if(val){
	  MRIsetVoxVal(mri2,vtxno,0,0,0,1);
	  continue;
	}
      }
    }
    MRIcopy(mri2,mri1);
  
  }

  MRIfree(&mri2);

  MRIwrite(mri1,"mymask2.mgh");

  return(mri1);
}
