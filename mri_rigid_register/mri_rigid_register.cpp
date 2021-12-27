/*
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


//
// mri_rigid_register.c
//
// original author:
// original date  :
//
// Warning: Do not edit the following four lines.  CVS maintains them.
//
////////////////////////////////////////////////////////////////////

/*E*

  Now works reasonably well with non-coronal anisotropic volumes -
  with the following exception:

  BUG: Gak, why do

  $ ~/dev/mri_rigid_register/mri_rigid_register /tmp/ebeth/BERA13_recon/synth.mgh /tmp/ebeth/BERA13_recon/cor-rot30/ /tmp/ebeth/BERA13_recon/dmrr_synthmgh2corrot30cor-corras.lta

  and

  $ ~/dev/mri_rigid_register/mri_rigid_register /tmp/ebeth/BERA13_recon/cor-rot30/ /tmp/ebeth/BERA13_recon/synth.mgh /tmp/ebeth/BERA13_recon/dmrr_corrot30cor2synthmgh-corras.lta

  give completely different results??? (one finds the 30degree
  rotation, but the other is close to I4 + c_ras!!)

  Already M_reg from estimate_rigid_regmatrix() is wrong for the
  second case.

M_reg
 0.996   0.091  -0.003   1.081;
-0.091   0.996  -0.013   19.870;
 0.002   0.013   1.000   7.840;
 0.000   0.000   0.000   1.000;

  Luckily, folks seem to prefer to register flash to mprage than vice
  versa...

  I don't see the problem.  It appears the old (before-my-tampering)
  version fails similarly.  Which makes sense, since it starts getting
  it wrong before it gets to the parts I was tampering with.

  Complaints to ebeth@nmr.mgh.harvard.edu :(

  Ugh, if one way works, we could always register that way then return
  the inverse if nec., could have a thing that decides which is the
  "bigger" one - is size the problem, or format?


  About the estimate_rigid_regmatrix() incarnations below - because of
  this bugginess, I reverted to the oldest copy I could find.  Others,
  #if-0'd out, are mine that was written to match only over a set of
  brain points, and the same thing only with my changes commented out.
  But none of them works.

  *E*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrimorph.h"
#include "mri_conform.h"
#include "utils.h"
#include "timer.h"
#include "matrix.h"
#include "transform.h"
#include "version.h"

//#define LINEAR_CORONAL_RAS_TO_CORONAL_RAS       21

#include "region.h"

static int write_iterations = 0 ;

static MRI *MRR_MRIhybrid(MRI *mri_src, MRI *mri_target, MRI *mri_xformed);
/*static */
MATRIX *MRR_VoxelXformToCoronalRasXform(MRI *mri_src, MRI *mri_dst, MATRIX *m_voxel_xform, MATRIX *m_coronalras_xform);
//static int maxabs(int a, int b, int c); //E/ unquestionably the wrong tool

static int apply_transform(MRI *mri_in, MRI *mri_target, MATRIX *M_reg, MRI *mri_xformed) ;

static char *B1_fname = NULL ;
static int noskull = 0 ;
static double tol=1e-10 ;

#if 0
// someone else's functions, never called
MRI *MRIsadd(MRI *mri1, MRI *mri2, MRI *mri_dst) ;
MRI *MRIsscalarMul(MRI *mri_src, MRI *mri_dst, float scalar) ;
MRI *MRIssqrt(MRI *mri_src, MRI *mri_dst) ;
#endif

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static int window_flag = 0 ;

static void usage_exit(int code) ;

static int niter=10 ; //not used(?)
static int voxel_xform = 0;
static float thresh = 25 ;
//static int conform = 0 ; //not implemented
static int sinc_flag = 1;
static int sinchalfwindow = 3;
static int max_ndelta = 3;
static int noscale = 1 ;

static char *residual_name = NULL ;

#define MAX_IMAGES 100
#define MIN_T1   5  // avoid singularity at T1=0
#define MIN_PD   5
#define MIN_ITER 5
#define MAX_ITER 5000

#define MAX_WRITE 10000
static int nwrite = 0 ;
static char *write_volumes[MAX_WRITE] ;
static char *out_fnames[MAX_WRITE] ;
static int apply_xform = 0 ;
static int exit_if_writefile_not_found = 1 ;

static void estimate_rigid_regmatrix(MRI *mri_source, MRI *mri_target, MATRIX *M_reg, MRI *mri_B1);

int
main(int argc, char *argv[]) {
  char   **av ;
  int    ac, nargs, i ;
  MRI    *mri_src, *mri_target, *mri_B1 ;
  char   *src_fname, *out_fname, *target_fname ;
  int    msec, minutes, seconds;
  Timer start ;
  LTA    *lta ;
  MATRIX *M_reg, *vox_s2vox_t; // *m_coronalras_src2trg;

  nargs = handleVersionOption(argc, argv, "mri_rigid_register");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;

  ac = argc ; //not used
  av = argv ; //not used
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    usage_exit(1) ;

  src_fname = argv[1] ;
  target_fname = argv[2] ;
  out_fname = argv[3] ;

  printf("reading source from %s...\n", src_fname) ;
  mri_src = MRIread(src_fname) ;
  printf("reading target from %s...\n", target_fname) ;
  mri_target = MRIread(target_fname) ;
  if (!mri_src)
    ErrorExit(ERROR_NOFILE, "%s: could not read src MRI from %s",
              Progname, src_fname) ;
  if (!mri_target)
    ErrorExit(ERROR_NOFILE, "%s: could not read dst MRI from %s",
              Progname, target_fname) ;

  if (B1_fname) {
    mri_B1  = MRIread(B1_fname) ;
    if (!mri_B1)
      ErrorExit(ERROR_NOFILE, "%s: could not read B1 mapfrom %s", Progname, B1_fname) ;
  } else
    mri_B1 = NULL ;


  //E/ someone was working on this conform piece:
#if 0
  if (conform) {
    MRI *mri_tmp ;

    printf("embedding and interpolating volume\n") ;
    mri_tmp = MRIconform(mri_flash[nvolumes]) ;
    /*      MRIfree(&mri_src) ;*/
    mri_flash[nvolumes] = mri_tmp ;
  }
#endif

  M_reg = MatrixIdentity(4,(MATRIX *)NULL);

  estimate_rigid_regmatrix(mri_src,mri_target,M_reg, mri_B1);

  //E/ output of this is the tx matrix from src to target in RAS
  //space, i.e. mm - the actual ras_s2ras_t matrix - e.g. identity4
  //for "mrr vol1 vol1 out.lta"

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    printf("M_reg\n");
    MatrixPrint(stdout,M_reg);
  }

  printf("writing registration matrix to %s...\n", out_fname);
  lta = LTAalloc(1,NULL) ;
//#define _MRR_DEBUG
#ifdef _MRR_DEBUG
  fprintf(stderr, "mrr: vox_s2vox_t = \n");
  MatrixPrint(stderr, vox_s2vox_t);
#endif

  if (voxel_xform) {
    // M_reg is the RAS-to-RAS transform
    vox_s2vox_t = MRIrasXformToVoxelXform(mri_src, mri_target, M_reg, NULL);
    lta->type = LINEAR_VOXEL_TO_VOXEL ;
    MatrixCopy(vox_s2vox_t, lta->xforms[0].m_L) ;
    MatrixFree(&vox_s2vox_t);
  } else {
    // keep as it is if ras-to-ras
    lta->type = LINEAR_RAS_TO_RAS;
    MatrixCopy(M_reg, lta->xforms[0].m_L);
#if 0
    lta->type = LINEAR_CORONAL_RAS_TO_CORONAL_RAS ;
    //E/ LTAvoxelTransformToCoronalRasTransform(lta) ;
    // That didn't take voxel size into account.
    // Replacement fn turns vox2vox tx into cor-ras_to_cor-ras tx:
    m_coronalras_src2trg =
      MRR_VoxelXformToCoronalRasXform
      (mri_src, mri_target, vox_s2vox_t, NULL);
    MatrixCopy(m_coronalras_src2trg, lta->xforms[0].m_L) ;
#endif
  }
  // save src and target information in lta
  getVolGeom(mri_src, &lta->xforms[0].src);
  getVolGeom(mri_target, &lta->xforms[0].dst);
  LTAwriteEx(lta,out_fname) ;


  //E/ Maybe need cases here, too?  -out_like targ or templ or other?
  //("conform"?)

  for (i = 0 ; i < nwrite ; i++) {
    MRI *mri_in, *mri_xformed ;

    mri_in = MRIread(write_volumes[i]) ;

    if (!mri_in) {
      if (exit_if_writefile_not_found)
        ErrorExit(ERROR_NOFILE, "%s: could not read volume %s", Progname, write_volumes[i]) ;
      else {
        fprintf(stderr, "%s: could not read volume %s - skipping...\n", Progname, write_volumes[i]) ;
        continue;
      }
    }
#if 1
    mri_xformed = MRR_MRIhybrid(mri_in, mri_target, NULL);
    apply_transform(mri_in, mri_target, M_reg, mri_xformed) ;
#else
    //E/ Someone could work on this when things quiet down - to make
    //this option act more like mri_transform -
    if (apply_out_type == 1) // -out_like mri_target
    {
      mri_xformed = MRR_MRIhybrid(mri_in, mri_target, NULL);
      apply_transform(mri_in, mri_target, M_reg, mri_xformed) ;
    } else if (apply_out_type == 2) // -out_like tkmedit template - coronal 1mm^3 256^3 c_ras=0 but not UCHAR
      // for compatibility with mri_transform, which won't know about target unless you specify -out_like.
    {
      mri_xformed = MRIalloc(256,256,256,  mri_in->type);
      mri_xformed->x_r = -1.0;
      mri_xformed->x_a =  0.0;
      mri_xformed->x_s =  0.0;
      mri_xformed->y_r =  0.0;
      mri_xformed->y_a =  0.0;
      mri_xformed->y_s = -1.0;
      mri_xformed->z_r =  0.0;
      mri_xformed->z_a =  1.0;
      mri_xformed->z_s =  0.0;
      mri_xformed->c_r =  0.0;
      mri_xformed->c_a =  0.0;
      mri_xformed->c_s =  0.0; //?

      apply_transform(mri_in, mri_xformed, M_reg, mri_xformed) ;
    }
#endif

    printf("writing transformed volume to %s...\n", out_fnames[i]) ;
    MRIcopyPulseParameters(mri_in, mri_xformed) ;
    MRIwrite(mri_xformed, out_fnames[i]) ;
    MRIfree(&mri_xformed) ;
    MRIfree(&mri_in) ;
  }

  if (apply_xform) // this used to be in e_r_r() - why?

    // this could use a little tinkering to make it more like
    // mri_transform should be - e.g. maybe let user specify the
    // target space

  {
    MRI *mri_xformed ;

    //E/    mri_xformed = MRIclone(mri_src, NULL) ;
    mri_xformed = MRR_MRIhybrid(mri_src, mri_target, NULL);
    apply_transform(mri_src, mri_target, M_reg, mri_xformed) ;
    printf("writing transformed volume to xformed.mgh...\n") ;
    MRIwrite(mri_xformed, "xformed.mgh") ;
    MRIfree(&mri_xformed) ;
  }

  MRIfree(&mri_src) ;
  MRIfree(&mri_target) ;

  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("parameter estimation took %d minutes and %d seconds.\n", minutes, seconds) ;

  exit(0);
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "sinc")) {
    sinchalfwindow = atoi(argv[2]);
    sinc_flag = 1;
    nargs = 1;
    printf("using sinc interpolation with windowwidth of %d\n",
           2*sinchalfwindow);
  } else if (!stricmp(option, "trilinear")) {
    sinc_flag = 0;
    printf("using trilinear interpolation\n");
  } else if (!stricmp(option, "scale")) {
    noscale = 0 ;
    printf("computing intensity scaling...\n") ;
  } else if (!stricmp(option, "noskull")) {
    noskull = 1 ;
    printf("assuming skull stripped images: discounting 0 voxels in fine alignment\n") ;
  } else if (!stricmp(option, "write")) {
    write_iterations = atoi(argv[2]) ;
    nargs = 1 ;
    printf("writing snapshots every %d iterations\n", write_iterations) ;
  } else if (!stricmp(option, "window")) {
    window_flag = 1 ;
    printf("applying hanning window to volumes...\n") ;
  } else if (!stricmp(option, "B1")) {
    B1_fname = argv[2] ;
    nargs = 1 ;
    printf("using B1 map from %s to weight registration...\n", B1_fname) ;
  } else if (!stricmp(option, "tol")) {
    sscanf(argv[2], "%lf", &tol) ;
    nargs = 1 ;
    printf("using tol = %2.2e\n", tol) ;
  }
  /*  else if (!stricmp(option, "noconform")) // not implemented
  {
    conform = 0 ;
    printf("inhibiting isotropic volume interpolation\n") ;
    } */
  else switch (toupper(*option)) {
    case 'N':
      niter = atoi(argv[2]) ;
      printf("performing estimation/motion correction %d times...\n", niter) ;
      // not used
      nargs = 1 ;
      break ;
    case 'W':
      write_volumes[nwrite] = argv[2] ;
      out_fnames[nwrite] = argv[3] ;
      printf("transforming volume %s and writing to %s....\n",
             write_volumes[nwrite], out_fnames[nwrite]) ;
      nwrite++ ;
      nargs = 2 ;
      break ;
    case 'A':
      apply_xform = 1 ;
      break ;
    case 'F':
      //E/trying to emulate mri_transform behavior - you'd hate to have
      //to rerun all of mri_rigid_register if one of your files is missing
      exit_if_writefile_not_found = 0 ;
      break ;
    case 'V':
      voxel_xform = 1 ;
      fprintf(stderr, "CAVEAT: e.g. mri_transform may not handle type LINEAR_VOXEL_TO_VOXEL well - consider writing out a LINEAR_CORONAL_RAS_TO_CORONAL_RAS transform.  Invoke mri_transform -h for more info.\n");
      break ;
    case 'R':
      residual_name = argv[2] ;
      printf("writing out residuals to %s...\n", residual_name) ;
      nargs = 1 ;
      break ;
    case 'T':
      thresh = atof(argv[2]) ;
      nargs = 1 ;
      printf("ignoring locations in which all images are less than %2.2f\n",
             thresh) ;
      break ;
    case '?':
    case 'U':
      usage_exit(0) ;
      break ;
    default:
      printf("unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
usage_exit(int code) {
  fprintf(stderr, "usage: %s [options] <src volume> <target volume> <transform fname>\n\n", Progname) ;
  fprintf(stderr, "Consider _not running with -v for vox2vox transform.  The default, coronal-ras-to-coronal-ras is pretty useful...\n\n") ;
  exit(code) ;
}


#define MAX_VOX ((64*64*64)+1)

#define NSTEP   11

/*E* INCORRECT/BUG/FIXTHIS, but the idea was to use this to see if
  brain is more or less MRI_CORONAL or MRI_SAGITTAL, i.e. -y_s the
  biggest of the ._s's. or not, i.e. +-y (which?) is superior - really
  we just want max but need to know if we meant +y or -y. */
int
maxabs(int a, int b, int c) {
  int A,B,C;
  A=abs(a);
  B=abs(b);
  C=abs(c);
  if (A>B)
    if (A>C) return 0;
    else return 2;
  else
    if (B>C) return 1;
    else return 2;
}
/*E*/





static void
estimate_rigid_regmatrix(MRI *mri_source, MRI *mri_target, MATRIX *M_reg, MRI *mri_B1) {
  double   xf, yf, zf, tx, ty, tz, ax, ay, az, ca, sa, val1, val2, err, sse, best_sse, dt=0.1, da=RADIANS(0.005), scale_2_to_1, new_scale ;
  int      x, y, z, txi, tyi, tzi, axi, ayi, azi, indx, stepindx, changed, pass, num;
  int      width=mri_source->width, height=mri_source->height, depth=mri_source->depth, dx=10, dy=10, dz=10, nvalues, iter=0;
#if 1
  int      step[NSTEP]={1024,512,256,128,64,32,16,8,4,2,1}, scale, max_wt_indx ;
#else
  int      nstep=1, step[1]={1}, scale;
#endif
  MATRIX   *vox2ras_source, *ras2vox_source, *vox2ras_target, *ras2vox_target, *vox_s2vox_t;
  MATRIX   *M_reg_bak, *M_reg_opt, *M_tmp, *M_delta, *M_delta1, *M_delta2, *M_delta3, *M_delta4, *M_delta5, *M_delta6, *M_tmp2;
  MATRIX   *voxmat1, *voxmat2 ;
  double   voxval1[MAX_VOX], voxval2[MAX_VOX], sum1, sum2, wt[MAX_VOX],
  max_wt ;
  MRI      *mri_source_edge, *mri_target_edge ;

  if (noskull) {
    MRI *mri_tmp, *mri_tmp2 ;
    int i ;

#define VOXELS_FROM_EDGE 2
    mri_tmp = MRIbinarize(mri_source, NULL, 1, 128, 0) ;
    mri_tmp2 = NULL ;
    for (i = 0 ;i < VOXELS_FROM_EDGE ; i++) {
      mri_tmp2 = MRIdilate(mri_tmp, mri_tmp2) ;
      MRIcopy(mri_tmp2, mri_tmp) ;
    }
    mri_source_edge = mri_tmp ;
    MRIfree(&mri_tmp2) ;

    mri_tmp = MRIbinarize(mri_target, NULL, 1, 128, 0) ;
    mri_tmp2 = NULL ;
    for (i = 0 ;i < VOXELS_FROM_EDGE ; i++) {
      mri_tmp2 = MRIdilate(mri_tmp, mri_tmp2) ;
      MRIcopy(mri_tmp2, mri_tmp) ;
    }
    mri_target_edge = mri_tmp ;
    MRIfree(&mri_tmp2) ;
  }


  vox2ras_source = MRIgetVoxelToRasXform(mri_source) ;
  vox2ras_target = MRIgetVoxelToRasXform(mri_target) ;
  ras2vox_source = MatrixInverse(vox2ras_source, NULL);
  ras2vox_target = MatrixInverse(vox2ras_target, NULL);
  vox_s2vox_t = MatrixIdentity(4,NULL);

  nvalues = 0;
  scale_2_to_1 = 0.0 ;
  for (sum1 = sum2 = 0.0, z = 0 ; z < depth ; z++)
    for (y = 0 ; y < height ; y++)
      for (x = 0 ; x < width ; x++) {
        if ((x%dx==0)&&(y%dy==0)&&(z%dz==0)) {
          nvalues++;
        }
        MRIsampleVolume(mri_source, x, y, z, &val1) ;
        MRIsampleVolume(mri_target, x, y, z, &val2) ;
        sum1 += val1 ;
        sum2 += val2 ;
      }

  if (nvalues >= MAX_VOX)
    ErrorExit(ERROR_NOMEMORY, "%s: # of samples %d exceeds max %d",Progname,nvalues,MAX_VOX) ;

  if (noscale)
    scale_2_to_1 = 1 ;
  else
    scale_2_to_1 = sum1 / sum2 ;
  printf("nvalues = %d, initial scaling = %2.2f\n", nvalues, scale_2_to_1);

  voxmat1 = MatrixAlloc(4,nvalues,MATRIX_REAL);

  max_wt_indx = indx = 0;
  max_wt = 0 ;
  for (z = 0 ; z < depth ; z++)
    for (y = 0 ; y < height ; y++)
      for (x = 0 ; x < width ; x++) {
        if ((x%dx==0)&&(y%dy==0)&&(z%dz==0)) {
          indx++;
          voxmat1->rptr[1][indx] = x;
          voxmat1->rptr[2][indx] = y;
          voxmat1->rptr[3][indx] = z;
          voxmat1->rptr[4][indx] = 1;
          xf=x;
          yf=y;
          zf=z;
          MRIsampleVolume(mri_source, xf, yf, zf, &val1) ;
          voxval1[indx] = val1;
          /*
                voxval1[indx] = MRISvox(mri_source, x, y, z);
          */
          if (mri_B1 == NULL)
            wt[indx] = 1 ;
          else    /* use inverse of field-map (squared) to wt registration */
          {
            MRIsampleVolume(mri_B1, xf, yf, zf, &val1) ;
            if (FZERO(val1))  /* should never happen */
              val1 = 1e-3 ;
            wt[indx] = 1.0 / (val1*val1) ;
            if (wt[indx] > max_wt) {
              max_wt = wt[indx] ;
              max_wt_indx = indx ;
            }
          }
        }
      }

  if (Gdiag & DIAG_SHOW)
    printf("M_reg (initial)\n");
  MatrixPrint(stdout,M_reg);

  M_delta = MatrixIdentity(4,NULL);
  M_delta1 = MatrixIdentity(4,NULL);
  M_delta2 = MatrixIdentity(4,NULL);
  M_delta3 = MatrixIdentity(4,NULL);
  M_delta4 = MatrixIdentity(4,NULL);
  M_delta5 = MatrixIdentity(4,NULL);
  M_delta6 = MatrixIdentity(4,NULL);

  M_reg_opt = MatrixCopy(M_reg, NULL);
  M_reg_bak = MatrixCopy(M_reg_opt, NULL);
  M_tmp = MatrixCopy(M_reg, NULL);
  M_tmp2 = MatrixCopy(M_tmp, NULL) ;

  MatrixMultiply(M_reg,vox2ras_source,M_tmp);
  MatrixMultiply(ras2vox_target,M_tmp,vox_s2vox_t);
  if (write_iterations > 0) {
    char fname[STRLEN] ;
    MRI  *mri_aligned ;

    sprintf(fname, "step%04d", iter) ;
    mri_aligned = MRIlinearTransform(mri_source, NULL, vox_s2vox_t);
    MRIwriteImageViews(mri_target, "target", IMAGE_SIZE) ;
    MRIwriteImageViews(mri_aligned, "step0000", IMAGE_SIZE) ;
    MRIfree(&mri_aligned) ;
  }

  best_sse = 1e30;
  sse = 0;
  voxmat2 = MatrixMultiply(vox_s2vox_t,voxmat1,NULL);
  printf("initial voxel transform:\n") ;
  MatrixPrint(stdout, vox_s2vox_t) ;
  scale = step[0] ;
  for (indx=1; indx<=nvalues; indx++) {
    if (indx == Gdiag_no)
      DiagBreak() ;
    if ((voxmat1->rptr[1][indx] == Gx) &&
        (voxmat1->rptr[2][indx] == Gy) &&
        (voxmat1->rptr[3][indx] == Gz))
      DiagBreak() ;
    xf=voxmat2->rptr[1][indx];
    yf=voxmat2->rptr[2][indx];
    zf=voxmat2->rptr[3][indx];
    MRIsampleVolume(mri_target, xf, yf, zf, &val2) ;
    voxval2[indx] = scale_2_to_1 * val2;
    val1 = voxval1[indx];
    err = val1-val2;
#if 1
    /* ignore background voxels when doing fine alignment */
    if (dt*scale <= 2 && (FZERO(val1) || FZERO(val2)))
      err = 0 ;
#endif
    sse += err*err*wt[indx];
  }
  sse /= nvalues;
  best_sse = sse ;
  printf("initial sse = %2.2f (%2.1f)\n", best_sse, sqrt(best_sse)) ;

  /* first find best translation */
  axi = ayi = azi = 0 ;
  scale = 256 ;
  for (stepindx=0; stepindx<NSTEP; stepindx++) {
    scale = step[stepindx] ;
    do {
      changed = 0 ;
      for (txi = -3; txi <= 3; txi++)
        for (tyi = -3; tyi <= 3; tyi++)
          for (tzi = -3; tzi <= 3; tzi++) {
            tx = txi*dt*scale;
            ty = tyi*dt*scale;
            tz = tzi*dt*scale;
            //      printf("checking translation (%2.1f, %2.1f, %2.1f)\n", tx, ty, tz) ;
            M_delta1->rptr[1][4]=tx;
            M_delta1->rptr[2][4]=ty;
            M_delta1->rptr[3][4]=tz;
            M_delta2->rptr[2][2]=1;
            M_delta2->rptr[2][3]=0;
            M_delta2->rptr[3][2]=0;
            M_delta2->rptr[3][3]=1;
            MatrixMultiply(M_delta2,M_delta1,M_delta5);
            M_delta3->rptr[1][1]=1;
            M_delta3->rptr[1][3]=0;
            M_delta3->rptr[3][1]=0;
            M_delta3->rptr[3][3]=1;
            MatrixMultiply(M_delta3,M_delta5,M_delta6);
            M_delta4->rptr[1][1]=1;
            M_delta4->rptr[1][2]=0;
            M_delta4->rptr[2][1]=0;
            M_delta4->rptr[2][2]=1;
            MatrixMultiply(M_delta4,M_delta6,M_delta);
            MatrixMultiply(M_delta,M_reg_bak,M_reg);

            MatrixMultiply(M_reg,vox2ras_source,M_tmp);
            MatrixMultiply(ras2vox_target,M_tmp,vox_s2vox_t);

            MatrixMultiply(vox_s2vox_t,voxmat1,voxmat2);
            sse = 0;
            for (indx=1; indx<=nvalues; indx++) {
              if (indx == Gdiag_no)
                DiagBreak() ;
              if ((voxmat1->rptr[1][indx] == Gx) &&
                  (voxmat1->rptr[2][indx] == Gy) &&
                  (voxmat1->rptr[3][indx] == Gz))
                DiagBreak() ;
              xf=voxmat2->rptr[1][indx];
              yf=voxmat2->rptr[2][indx];
              zf=voxmat2->rptr[3][indx];
              MRIsampleVolume(mri_target, xf, yf, zf, &val2) ;
              voxval2[indx] = scale_2_to_1 * val2;
              val1 = voxval1[indx];
              err = val1-val2;
              sse += err*err*wt[indx];
            }
            sse /= nvalues;
            if (sse<best_sse-tol) {
              iter++ ;
              if (iter == Gdiag_no)
                DiagBreak() ;
              best_sse = sse;
              if (FZERO(best_sse))
                DiagBreak() ;
              MatrixCopy(M_reg, M_reg_opt);
              if ((write_iterations > 0) && (iter % write_iterations) == 0) {
                char fname[STRLEN] ;
                MRI *mri_aligned ;

                MatrixMultiply(M_reg_opt,vox2ras_source,M_tmp);
                MatrixMultiply(ras2vox_target,M_tmp,vox_s2vox_t);
                sprintf(fname, "step%04d", iter) ;
                mri_aligned = MRIlinearTransform(mri_source, NULL, vox_s2vox_t);
                printf("writing snapshots to %s...\n", fname) ;
                MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
                MRIfree(&mri_aligned) ;
              }

              /*
              printf("M_delta\n"); MatrixPrint(stdout,M_delta);
              */
              /*
              printf("M_delta1\n"); MatrixPrint(stdout,M_delta1);
              printf("M_delta2\n"); MatrixPrint(stdout,M_delta2);
              printf("M_delta3\n"); MatrixPrint(stdout,M_delta3);
              printf("M_delta4\n"); MatrixPrint(stdout,M_delta4);
              printf("M_delta5\n"); MatrixPrint(stdout,M_delta5);
              printf("M_delta6\n"); MatrixPrint(stdout,M_delta6);
              */
              if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
                printf("M_reg_opt\n");
                MatrixPrint(stdout,M_reg_opt);
                printf("vox_s2vox_t\n");
                MatrixPrint(stdout,vox_s2vox_t);
              }
              changed = 1;
            }
          }
    } while (changed) ;
  }
  printf("scale %d: after searching for optimal translation %2.1f, sse = %2.2f (%2.1f)\n", scale, dt*scale, best_sse, sqrt(best_sse)) ;
  MatrixPrint(stdout, M_reg_opt) ;


  if (noskull)
    printf("ignoring skull-stripped voxels...\n") ;
  for (stepindx=0; stepindx<NSTEP; stepindx++) {
    scale = step[stepindx];
    changed = 1;
    pass = 0;

    while (changed) {
      pass++;
      changed = 0;
      MatrixCopy(M_reg_opt, M_reg_bak);
      for (txi = -1; txi <= 1; txi++)
        for (tyi = -1; tyi <= 1; tyi++)
          for (tzi = -1; tzi <= 1; tzi++)
            for (axi = -1; axi <= 1; axi++)
              for (ayi = -1; ayi <= 1; ayi++)
                for (azi = -1; azi <= 1; azi++)
                  if (((txi!=0)+(tyi!=0)+(tzi!=0)+(axi!=0)+(ayi!=0)+(azi!=0))<=max_ndelta) {
                    tx = txi*dt*scale;
                    ty = tyi*dt*scale;
                    tz = tzi*dt*scale;
                    ax = axi*da*scale;
                    ay = ayi*da*scale;
                    az = azi*da*scale;
                    M_delta1->rptr[1][4]=tx;
                    M_delta1->rptr[2][4]=ty;
                    M_delta1->rptr[3][4]=tz;
                    ca = cos(ax);
                    sa = sin(ax);
                    M_delta2->rptr[2][2]=ca;
                    M_delta2->rptr[2][3]=-sa;
                    M_delta2->rptr[3][2]=sa;
                    M_delta2->rptr[3][3]=ca;
                    MatrixMultiply(M_delta2,M_delta1,M_delta5);
                    ca = cos(ay);
                    sa = sin(ay);
                    M_delta3->rptr[1][1]=ca;
                    M_delta3->rptr[1][3]=-sa;
                    M_delta3->rptr[3][1]=sa;
                    M_delta3->rptr[3][3]=ca;
                    MatrixMultiply(M_delta3,M_delta5,M_delta6);
                    ca = cos(az);
                    sa = sin(az);
                    M_delta4->rptr[1][1]=ca;
                    M_delta4->rptr[1][2]=-sa;
                    M_delta4->rptr[2][1]=sa;
                    M_delta4->rptr[2][2]=ca;
                    MatrixMultiply(M_delta4,M_delta6,M_delta);
                    MatrixMultiply(M_delta,M_reg_bak,M_reg);

                    MatrixMultiply(M_reg,vox2ras_source,M_tmp);
                    MatrixMultiply(ras2vox_target,M_tmp,vox_s2vox_t);

                    MatrixMultiply(vox_s2vox_t,voxmat1,voxmat2);
                    sse = 0;
                    for (indx=1; indx<=nvalues; indx++) {
                      if (indx == Gdiag_no)
                        DiagBreak() ;
                      if ((voxmat1->rptr[1][indx] == Gx) &&
                          (voxmat1->rptr[2][indx] == Gy) &&
                          (voxmat1->rptr[3][indx] == Gz))
                        DiagBreak() ;
                      xf=voxmat2->rptr[1][indx];
                      yf=voxmat2->rptr[2][indx];
                      zf=voxmat2->rptr[3][indx];
                      MRIsampleVolume(mri_target, xf, yf, zf, &val2) ;
                      voxval2[indx] = scale_2_to_1 * val2;
                      val1 = voxval1[indx];
                      err = val1-val2;

                      /* ignore background voxels when doing fine alignment */
                      if (noskull) {
                        double b1, b2, xf1, yf1, zf1 ;

                        MRIsampleVolume(mri_target_edge, xf, yf, zf, &b2) ;
                        xf1 = voxmat1->rptr[1][indx];
                        yf1 = voxmat1->rptr[2][indx];
                        zf1 = voxmat1->rptr[3][indx];
                        MRIsampleVolume(mri_source_edge, xf1, yf1, zf1, &b1) ;

                        if (b1 > 0.1 && b2 > 0.1 && (FZERO(val1) || FZERO(val2)))
                          err = 0 ;
                      }

                      sse += (err*err*wt[indx]);
                    }
                    sse /= nvalues;
                    if (sse<best_sse-tol) {
                      iter++ ;
                      if (iter == Gdiag_no)
                        DiagBreak() ;
                      best_sse = sse;
                      if (FZERO(best_sse))
                        DiagBreak() ;
                      MatrixCopy(M_reg, M_reg_opt);
                      if ((write_iterations > 0) && (iter % write_iterations) == 0) {
                        char fname[STRLEN] ;
                        MRI *mri_aligned ;

                        MatrixMultiply(M_reg_opt,vox2ras_source,M_tmp);
                        MatrixMultiply(ras2vox_target,M_tmp,vox_s2vox_t);
                        sprintf(fname, "step%04d", iter) ;
                        mri_aligned = MRIlinearTransform(mri_source, NULL, vox_s2vox_t);
                        printf("writing snapshots to %s...\n", fname) ;
                        MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
                        MRIfree(&mri_aligned) ;
                      }

                      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
                        printf("%d (%d) %f %f %f %f %f %f sse = %f (%f)\n",scale,pass,tx,ty,tz,ax,
                               ay,az,sse,sqrt(sse));
                      /*
                                printf("M_delta\n"); MatrixPrint(stdout,M_delta);
                      */
                      /*
                                printf("M_delta1\n"); MatrixPrint(stdout,M_delta1);
                                printf("M_delta2\n"); MatrixPrint(stdout,M_delta2);
                                printf("M_delta3\n"); MatrixPrint(stdout,M_delta3);
                                printf("M_delta4\n"); MatrixPrint(stdout,M_delta4);
                                printf("M_delta5\n"); MatrixPrint(stdout,M_delta5);
                                printf("M_delta6\n"); MatrixPrint(stdout,M_delta6);
                      */
                      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
                        printf("M_reg_opt\n");
                        MatrixPrint(stdout,M_reg_opt);
                        printf("vox_s2vox_t\n");
                        MatrixPrint(stdout,vox_s2vox_t);
                      }
                      changed = 1;
                    }
                  }

    }

    MatrixMultiply(M_reg,vox2ras_source,M_tmp);
    MatrixMultiply(ras2vox_target,M_tmp,vox_s2vox_t);

#if 1
    for (new_scale = 0.0, num = 0.0, z = 0 ; z < depth ; z++) {
      for (y = 0 ; y < height ; y++) {
        for (x = 0 ; x < width ; x++) {
          M_tmp->rptr[1][1] = x;
          M_tmp->rptr[2][1] = y;
          M_tmp->rptr[3][1] = z;
          M_tmp->rptr[4][1] = 1;
          MatrixMultiply(vox_s2vox_t,M_tmp,M_tmp2);
          xf = M_tmp2->rptr[1][1] ;
          yf = M_tmp2->rptr[2][1] ;
          zf = M_tmp2->rptr[3][1] ;

          MRIsampleVolume(mri_source, x, y, z, &val1) ;
          MRIsampleVolume(mri_target, xf, yf, zf, &val2) ;
          if (FZERO(val1)|| FZERO(val2))
            continue ;
          num++ ;
          new_scale += val1 / val2 ;
          if (!std::isfinite(new_scale))
            DiagBreak()  ;
        }
      }
    }

#else
    for (num = 0, new_scale = 0.0, indx=1; indx<=nvalues; indx++) {
      xf=voxmat2->rptr[1][indx];
      yf=voxmat2->rptr[2][indx];
      zf=voxmat2->rptr[3][indx];
      MRIsampleVolume(mri_target, xf, yf, zf, &val2) ;
      val1 = voxval1[indx];
      /* ignore background voxels when doing fine alignment */
      if (!FZERO(val1) && !FZERO(val2)) {
        num++ ;
        new_scale += val1 / val2 ;
      }
    }
#endif
#if 0
    if (!FZERO(num)) {
      new_scale /= (double)num ;
      if (noscale)
        scale_2_to_1 = 1 ;
      else
        scale_2_to_1 = new_scale ;
      if (!finite(scale_2_to_1))
        DiagBreak() ;
      printf("setting image scaling factor to %2.2f\n", scale_2_to_1) ;
      if (scale_2_to_1 > 100)
        DiagBreak() ;
      if (stepindx >= 7)
        DiagBreak() ;
    }
#endif

    printf("step %d: best_sse = %f (%f)\n",stepindx,best_sse,sqrt(best_sse));
    if (Gdiag & DIAG_SHOW) {
      printf("M_reg_opt\n");
      MatrixPrint(stdout,M_reg_opt);
      /*
            printf("vox_s2vox_t\n"); MatrixPrint(stdout,vox_s2vox_t);
      */
    }
  }

  if (write_iterations > 0) {
    char fname[STRLEN] ;
    MRI  *mri_aligned ;
    MATRIX *vox2ras_source, *ras2vox_source, *vox2ras_target, *ras2vox_target, *m_tmp ;

    vox2ras_source = MRIgetVoxelToRasXform(mri_source) ;
    vox2ras_target = MRIgetVoxelToRasXform(mri_target) ;
    ras2vox_source = MatrixInverse(vox2ras_source, NULL);
    ras2vox_target = MatrixInverse(vox2ras_target, NULL);
    m_tmp = MatrixMultiply(M_reg_opt,vox2ras_source, NULL);
    vox_s2vox_t = MatrixMultiply(ras2vox_target,m_tmp, NULL);

    sprintf(fname, "step%04d", iter+1) ;
    mri_aligned = MRIlinearTransform(mri_source, NULL, vox_s2vox_t);
    MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
    MRIfree(&mri_aligned) ;
    MatrixFree(&vox2ras_source) ;
    MatrixFree(&ras2vox_source) ;
    MatrixFree(&vox2ras_target) ;
    MatrixFree(&ras2vox_target) ;
    MatrixFree(&m_tmp) ;
  }

  /*E* I (re)moved this:
  if (apply_xform)
  {
    MRI *mri_xformed ;

    mri_xformed = MRIclone(mri_source, NULL) ;
    apply_transform(mri_source, mri_target, M_reg_opt, mri_xformed) ;
    MRIwrite(mri_xformed, "xformed.mgh") ;
    MRIfree(&mri_xformed) ;
  }
  *E*/

  if (noskull) {
    MRIfree(&mri_source_edge) ;
    MRIfree(&mri_target_edge) ;
  }
  MatrixCopy(M_reg_opt, M_reg);
}




#if 0
//E/ my "contributions" commented-out
static void
estimate_rigid_regmatrix(MRI *mri_source, MRI *mri_target, MATRIX *M_reg) {
  /*E* 14 aug 2002 ebeth

    Actually, all my contribution is commented out for the moment for
    debugging purposes - or, rather, bbox bounds are set to the whole
    volume instead of just the skull.

    For coronal volumes (and anything where superior is within
    45degrees of y -- or -y, which is a BUG), take sample points only
    in the superior .75 of
    MRIfindApproximateSkullBoundingBox(mri_source, thresh=50, &bbox);

    For anything else, it looks at the whole thing (or could make it
    use the whole MRIfindApproximateSkullBoundingBox).

    It performs microscopically better.

  */

  MRI_REGION bbox;
  int xplusdx,yplusdy,zplusdz;
  /*E* sorry, these are not the same d? as the variables d?  below -
    bbox range goes from ? to ?+d? - these are _those.  */

  double   xf, yf, zf, tx, ty, tz, ax, ay, az, ca, sa, val1, val2, err, sse, best_sse, dt=0.01, da=RADIANS(0.005), tol=0.00001;
  int      x, y, z, txi, tyi, tzi, axi, ayi, azi, indx, stepindx, changed, pass;
  int dx=10, dy=10, dz=10, nvalues;

#if 1
  /*
    int      nstep=8, step[8]={32,16,8,4,2,1}, scale;
  */
  int      step[NSTEP]={1024,512,256,128,64,32,16,8,4,2,1}, scale;
#else
  int      nstep=1, step[1]={1}, scale;
#endif
  MATRIX   *vox2ras_source, *ras2vox_source, *vox2ras_target, *ras2vox_target, *vox_s2vox_t;
  MATRIX   *M_reg_bak, *M_reg_opt, *M_tmp, *M_delta, *M_delta1, *M_delta2, *M_delta3, *M_delta4, *M_delta5, *M_delta6;
  MATRIX   *voxmat1, *voxmat2;
  double   voxval1[MAX_VOX], voxval2[MAX_VOX];

  /*E*
  fprintf(stderr,"TEST!!!!\n");
  //  M_reg = MatrixIdentity(4,NULL); // why did this line break it??
  / * * /
  M_reg->rptr[2][2] = 0.866;
  M_reg->rptr[3][3] = 0.866;
  M_reg->rptr[2][3] = 0.5;
  M_reg->rptr[3][2] = -0.5;
  / * * /
  MatrixPrint(stderr, M_reg);
  return;
  *E*/

  vox2ras_source = MRIgetVoxelToRasXform(mri_source) ;
  vox2ras_target = MRIgetVoxelToRasXform(mri_target) ;
  ras2vox_source = MatrixInverse(vox2ras_source, NULL);
  ras2vox_target = MatrixInverse(vox2ras_target, NULL);
  vox_s2vox_t = MatrixIdentity(4,NULL);

  bbox = *REGIONalloc();

  bbox.x = bbox.y = bbox.z = 0;
  bbox.dx = mri_source->width;
  bbox.dy = mri_source->height;
  bbox.dz = mri_source->depth;

#if 0
  if (MRIfindApproximateSkullBoundingBox(mri_source, /*thresh*/ 50, &bbox)
      != NO_ERROR)
    exit(1);

  fprintf(stderr, "Skull1's ApproximateSkullBoundingBox: bbox.x,y,z;dx,dy,dz = %d,%d,%d;%d,%d,%d\n", bbox.x,bbox.y,bbox.z, bbox.dx,bbox.dy,bbox.dz);

  //E/ If there's no xyzc_ras, it does _not fall-back on
  //slice_direction - but I _could write it that way.


  if (mri_source->ras_good_flag) {
    switch (maxabs(mri_source->x_s,mri_source->y_s,mri_source->z_s)) {
    case 1: //coronal or sagittal
      /*E*
      if (mri_source->y_s>0) bbox.y += bbox.dy/4.;  //probably not, for coronal
      bbox.dy *=.75;
      fprintf(stderr, "Volume is coronal or sagittal or within 45 degrees - BUG: or 180deg off of that - for this case, matching on grid points only in the superior 3/4 of skull bounding box is implemented.\n");
      break;
      *E*/
    case 0:
    case 2:
    default:
      fprintf(stderr, "Matching on grid points over entire image.\n");
      bbox.x = bbox.y = bbox.z = 0;
      bbox.dx = mri_source->width;
      bbox.dy = mri_source->height;
      bbox.dz = mri_source->depth;
      break;
    }
  }
  /*E*
  else
    {
      switch (mri_source->slice_direction)
  {
  case MRI_CORONAL:
   bbox.dy *=.75;
   fprintf(stderr, "MRI_CORONAL case - so do this matching on grid points only in the superior 3/4 of skull bounding box.\n");
   break;
  case MRI_SAGITTAL:
   //   bbox.dx *= .75;
   //   fprintf(stderr, "Untested case - MRI_SAGITTAL\n");
   fprintf(stderr, "MRI_SAGITTAL case - finding some brain not implemented, so matching on grid points over entire image.\n");
   bbox.x = bbox.y = bbox.z = 0;
   bbox.dx = mri_source->width;
   bbox.dy = mri_source->height;
   bbox.dz = mri_source->depth;
   break;
  case MRI_HORIZONTAL:
   //   bbox.dz *= .75;
   //   fprintf(stderr, "Untested case - MRI_HORIZONTAL\n");
   fprintf(stderr, "MRI_HORIZONTAL case - finding some brain not implemented, so matching on grid points over entire image.\n");
   bbox.x = bbox.y = bbox.z = 0;
   bbox.dx = mri_source->width;
   bbox.dy = mri_source->height;
   bbox.dz = mri_source->depth;
   break;
  default: //contract in every direction?
   fprintf(stderr, "Finding some brain not implemented for this case, so matching on grid points over entire image.\n");
   bbox.x = bbox.y = bbox.z = 0;
   bbox.dx = mri_source->width;
   bbox.dy = mri_source->height;
   bbox.dz = mri_source->depth;
   break;
  }
      bbox.x = bbox.y = bbox.z = 0;
      bbox.dx = mri_source->width;
      bbox.dy = mri_source->height;
      bbox.dz = mri_source->depth;
    }
    *E*/
#endif


  xplusdx = bbox.x + bbox.dx;
  yplusdy = bbox.y + bbox.dy;
  zplusdz = bbox.z + bbox.dz;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stderr, "x,y,z;xplusdx,yplusdy,zplusdz = %d,%d,%d;%d,%d,%d\n", bbox.x,bbox.y,bbox.z, xplusdx,yplusdy,zplusdz);
  }

  nvalues = 0;
  for (z = bbox.z ; z < zplusdz ; z+=dz)
    for (y = bbox.y ; y < yplusdy ; y+=dy)
      for (x = bbox.x ; x < xplusdx ; x+=dx)
        nvalues++;

  printf("nvalues = %d\n", nvalues);

  voxmat1 = MatrixAlloc(4,nvalues,MATRIX_REAL);
  voxmat2 = MatrixCopy(voxmat1, NULL);

  indx = 0;
  for (z = bbox.z ; z < zplusdz ; z+=dz)
    for (y = bbox.y ; y < yplusdy ; y+=dy)
      for (x = bbox.x ; x < xplusdx ; x+=dx) {
        indx++;
        voxmat1->rptr[1][indx] = x;
        voxmat1->rptr[2][indx] = y;
        voxmat1->rptr[3][indx] = z;
        voxmat1->rptr[4][indx] = 1;
        xf=x;
        yf=y;
        zf=z;
        MRIsampleVolume(mri_source, xf, yf, zf, &val1) ;
        voxval1[indx] = val1;
        // voxval1[indx] = MRISvox(mri_source, x, y, z);
      }

  if (Gdiag & DIAG_SHOW)
    printf("M_reg (initial)\n");
  MatrixPrint(stdout,M_reg);

  M_delta = MatrixIdentity(4,NULL);
  M_delta1 = MatrixIdentity(4,NULL);
  M_delta2 = MatrixIdentity(4,NULL);
  M_delta3 = MatrixIdentity(4,NULL);
  M_delta4 = MatrixIdentity(4,NULL);
  M_delta5 = MatrixIdentity(4,NULL);
  M_delta6 = MatrixIdentity(4,NULL);

  M_reg_opt = MatrixCopy(M_reg, NULL);
  M_reg_bak = MatrixCopy(M_reg_opt, NULL);
  M_tmp = MatrixCopy(M_reg, NULL);

  best_sse = 1e30;
  for (stepindx=0; stepindx<NSTEP; stepindx++) {
    scale = step[stepindx];
    changed = 1;
    pass = 0;
    while (changed) {
      pass++;
      changed = 0;
      MatrixCopy(M_reg_opt, M_reg_bak);
      for (txi = -1; txi <= 1; txi++)
        for (tyi = -1; tyi <= 1; tyi++)
          for (tzi = -1; tzi <= 1; tzi++)
            for (axi = -1; axi <= 1; axi++)
              for (ayi = -1; ayi <= 1; ayi++)
                for (azi = -1; azi <= 1; azi++)
                  if (((txi!=0)+(tyi!=0)+(tzi!=0)+(axi!=0)+(ayi!=0)+(azi!=0))<=max_ndelta) {
                    tx = txi*dt*scale;
                    ty = tyi*dt*scale;
                    tz = tzi*dt*scale;
                    ax = axi*da*scale;
                    ay = ayi*da*scale;
                    az = azi*da*scale;
                    M_delta1->rptr[1][4]=tx;
                    M_delta1->rptr[2][4]=ty;
                    M_delta1->rptr[3][4]=tz;
                    ca = cos(ax);
                    sa = sin(ax);
                    M_delta2->rptr[2][2]=ca;
                    M_delta2->rptr[2][3]=-sa;
                    M_delta2->rptr[3][2]=sa;
                    M_delta2->rptr[3][3]=ca;
                    MatrixMultiply(M_delta2,M_delta1,M_delta5);
                    ca = cos(ay);
                    sa = sin(ay);
                    M_delta3->rptr[1][1]=ca;
                    M_delta3->rptr[1][3]=-sa;
                    M_delta3->rptr[3][1]=sa;
                    M_delta3->rptr[3][3]=ca;
                    MatrixMultiply(M_delta3,M_delta5,M_delta6);
                    ca = cos(az);
                    sa = sin(az);
                    M_delta4->rptr[1][1]=ca;
                    M_delta4->rptr[1][2]=-sa;
                    M_delta4->rptr[2][1]=sa;
                    M_delta4->rptr[2][2]=ca;
                    MatrixMultiply(M_delta4,M_delta6,M_delta);
                    MatrixMultiply(M_delta,M_reg_bak,M_reg);

                    MatrixMultiply(M_reg,vox2ras_source,M_tmp);
                    MatrixMultiply(ras2vox_target,M_tmp,vox_s2vox_t);

                    MatrixMultiply(vox_s2vox_t,voxmat1,voxmat2);
                    sse = 0;
                    for (indx=1; indx<=nvalues; indx++) {
                      xf=voxmat2->rptr[1][indx];
                      yf=voxmat2->rptr[2][indx];
                      zf=voxmat2->rptr[3][indx];
                      MRIsampleVolume(mri_target, xf, yf, zf, &val2) ;
                      voxval2[indx] = val2;
                      val1 = voxval1[indx];
                      err = val1-val2;
                      sse += err*err;
                    }
                    sse /= nvalues;
                    if (sse<best_sse-tol) {
                      best_sse = sse;
                      MatrixCopy(M_reg, M_reg_opt);
                      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
                        printf("%d (%d) %f %f %f %f %f %f sse = %f (%f)\n",scale,pass,tx,ty,tz,ax,
                               ay,az,sse,sqrt(sse));
                      /*
                                printf("M_delta\n"); MatrixPrint(stdout,M_delta);
                      */
                      /*
                                printf("M_delta1\n"); MatrixPrint(stdout,M_delta1);
                                printf("M_delta2\n"); MatrixPrint(stdout,M_delta2);
                                printf("M_delta3\n"); MatrixPrint(stdout,M_delta3);
                                printf("M_delta4\n"); MatrixPrint(stdout,M_delta4);
                                printf("M_delta5\n"); MatrixPrint(stdout,M_delta5);
                                printf("M_delta6\n"); MatrixPrint(stdout,M_delta6);
                      */
                      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
                        printf("M_reg_opt\n");
                        MatrixPrint(stdout,M_reg_opt);
                        printf("vox_s2vox_t\n");
                        MatrixPrint(stdout,vox_s2vox_t);
                      }
                      changed = 1;
                    }
                  }

    }
    printf("step %d: best_sse = %f (%f)\n",stepindx,best_sse,sqrt(best_sse));
    if (Gdiag & DIAG_SHOW) {
      printf("M_reg_opt\n");
      MatrixPrint(stdout,M_reg_opt);
      //      printf("vox_s2vox_t\n"); MatrixPrint(stdout,vox_s2vox_t);
    }
  }

  MatrixCopy(M_reg_opt, M_reg);
}
#endif





#if 0
//E/ my stuff dumb but in there
static void
estimate_rigid_regmatrix(MRI *mri_source, MRI *mri_target, MATRIX *M_reg) {
  /*E* 14 aug 2002
    ebeth changing where sample points are chosen - using
    MRIfindApproximateSkullBoundingBox(mri_source, thresh=50, &bbox);

    bounding box addition is implemented and tested for things closest
    to this direction and neither, otherwise:

    x_ras -1.000000 0.000000 0.000000
    y_ras 0.000000 0.000000 -1.000000
    z_ras 0.000000 1.000000 0.000000

    1 aug 2002

    took one MRI_CORONAL test brain and zeroed a grid of pixels from the
    ApproximateSkullBoundingBox to see where that sits.  Looking at the
    superior 3/4 of that seems be pretty good.  See
    ~ebeth/build/mri_rigid_register.c.compare for the comparison between
    matches for mine and for the prev rev.

    I've mostly tagged these additions w/ *E*
  */
  /*E*/ MRI_REGION bbox;
  /*E*/
  int xplusdx,yplusdy,zplusdz;

  double   xf, yf, zf, tx, ty, tz, ax, ay, az, ca, sa, val1, val2, err, sse, best_sse, dt=0.01, da=RADIANS(0.005), tol=0.00001;
  int      x, y, z, txi, tyi, tzi, axi, ayi, azi, indx, stepindx, changed, pass;

  /*E* unused: int width=mri_source->width, height=mri_source->height,
    depth=mri_source->depth; */

  /*E* We're taking points over a smaller volume, so we might want to
    lower dx,dy,dz - but so far we're still getting slightly better
    results even with so few points.  */
  int dx=10, dy=10, dz=10, nvalues;

#if 1
  /*
    int      nstep=8, step[8]={32,16,8,4,2,1}, scale;
  */
  int      step[NSTEP]={1024,512,256,128,64,32,16,8,4,2,1}, scale;
#else
  int      nstep=1, step[1]={1}, scale;
#endif
  MATRIX   *vox2ras_source, *ras2vox_source, *vox2ras_target, *ras2vox_target, *vox_s2vox_t;
  MATRIX   *M_reg_bak, *M_reg_opt, *M_tmp, *M_delta, *M_delta1, *M_delta2, *M_delta3, *M_delta4, *M_delta5, *M_delta6;
  MATRIX   *voxmat1, *voxmat2;
  double   voxval1[MAX_VOX], voxval2[MAX_VOX];


  vox2ras_source = MRIgetVoxelToRasXform(mri_source) ;
  vox2ras_target = MRIgetVoxelToRasXform(mri_target) ;
  ras2vox_source = MatrixInverse(vox2ras_source, NULL);
  ras2vox_target = MatrixInverse(vox2ras_target, NULL);
  vox_s2vox_t = MatrixIdentity(4,NULL);

  /*E*/
  bbox = *REGIONalloc();
  if (MRIfindApproximateSkullBoundingBox(mri_source, /*thresh*/ 50, &bbox)
      != NO_ERROR)
    exit(1);
  /*E*/

  fprintf(stderr, "Skull1's ApproximateSkullBoundingBox: bbox.x,y,z;dx,dy,dz = %d,%d,%d;%d,%d,%d\n", bbox.x,bbox.y,bbox.z, bbox.dx,bbox.dy,bbox.dz);

  /*E* Eyeballing it, looks like almost the same LTA output!

    coords are: col/row/slice = width/height/depth.

    Everything in ma_two_point has:

    x_ras -1.000000 0.000000 0.000000
    y_ras 0.000000 0.000000 -1.000000
    z_ras 0.000000 1.000000 0.000000

    Ick, could check orientation (MRI_CORONAL vs. MRI_HORIZONTAL
    vs. MRI_SAGITTAL) to see...  Or ImageDirCos?

    We just care which is inferior->superior direction, and for
    y_ras=0,0,-1, we want bbox.dy*=.75 (for our test image,
    $MA2/new/003001-*, this is the top 3/4 of skull bbox - gets most
    of the top of the head - top corners and edges extend only a
    little out of the skull into black

    should check *_ras or mri_src->slice_direction = MRI_FOO to
    decide.

  */


  if (mri_source->ras_good_flag>0) {
    switch (maxabs(mri_source->x_s,mri_source->y_s,mri_source->z_s)) {
    case 0:
      //   if (mri_source->x_s>0) bbox.x += bbox.dx/4.;
      //   bbox.dx *= .75;
      fprintf(stderr, "Untested case - |mri_source->x_s| greatest of |*_s| - let's just do this the old way\n");
      bbox.x = bbox.y = bbox.z = 0;
      bbox.dx = mri_source->width;
      bbox.dy = mri_source->height;
      bbox.dz = mri_source->depth;
    case 1:
      if (mri_source->y_s>0) bbox.y += bbox.dy/4.;  //probably not, for coronal
      bbox.dy *=.75;
      fprintf(stderr, "Whew.  We're pretty well equipped to handle the case where |mri_source->y_s| is the greatest of the |*_s|\n");
      break;
    case 2:
      //   if (mri_source->x_s>0) bbox.z += bbox.dz/4.; //?
      //   bbox.dx *= .75;
      fprintf(stderr, "Untested case - |mri_source->z_s| greatest of |*_s| - let's just do this the old way\n");
      bbox.x = bbox.y = bbox.z = 0;
      bbox.dx = mri_source->width;
      bbox.dy = mri_source->height;
      bbox.dz = mri_source->depth;
      break;
    default:
      break;
    }
  } else {
    switch (mri_source->slice_direction) {
    case MRI_CORONAL:
      bbox.dy *=.75;
      fprintf(stderr, "Whew.  We're pretty well equipped to handle the MRI_CORONAL case\n");
      break;
    case MRI_SAGITTAL:
      //   bbox.dx *= .75;
      //   fprintf(stderr, "Untested case - MRI_SAGITTAL\n");
      fprintf(stderr, "Untested case - MRI_SAGITTAL - let's just do this the old way\n");
      bbox.x = bbox.y = bbox.z = 0;
      bbox.dx = mri_source->width;
      bbox.dy = mri_source->height;
      bbox.dz = mri_source->depth;
      break;
    case MRI_HORIZONTAL:
      //   bbox.dz *= .75;
      //   fprintf(stderr, "Untested case - MRI_HORIZONTAL\n");
      fprintf(stderr, "Untested case - MRI_HORIZONTAL - let's just do this the old way\n");
      bbox.x = bbox.y = bbox.z = 0;
      bbox.dx = mri_source->width;
      bbox.dy = mri_source->height;
      bbox.dz = mri_source->depth;
      break;
    default: //contract in every direction?
      fprintf(stderr, "Untested case - let's just do this the old way\n");
      bbox.x = bbox.y = bbox.z = 0;
      bbox.dx = mri_source->width;
      bbox.dy = mri_source->height;
      bbox.dz = mri_source->depth;
      break;
    }
  }

  xplusdx = bbox.x + bbox.dx;
  yplusdy = bbox.y + bbox.dy;
  zplusdz = bbox.z + bbox.dz;
  fprintf(stderr, "x,y,z;xplusdx,yplusdy,zplusdz = %d,%d,%d;%d,%d,%d\n", bbox.x,bbox.y,bbox.z, xplusdx,yplusdy,zplusdz);

  nvalues = 0;
  for (z = bbox.z ; z < zplusdz ; z+=dz)
    for (y = bbox.y ; y < yplusdy ; y+=dy)
      for (x = bbox.x ; x < xplusdx ; x+=dx)
        nvalues++;

  printf("nvalues = %d\n", nvalues);

  voxmat1 = MatrixAlloc(4,nvalues,MATRIX_REAL);
  voxmat2 = MatrixCopy(voxmat1, NULL);


  indx = 0;
  for (z = bbox.z ; z < zplusdz ; z+=dz)
    for (y = bbox.y ; y < yplusdy ; y+=dy)
      for (x = bbox.x ; x < xplusdx ; x+=dx) {
        indx++;
        voxmat1->rptr[1][indx] = x;
        voxmat1->rptr[2][indx] = y;
        voxmat1->rptr[3][indx] = z;
        voxmat1->rptr[4][indx] = 1;
        xf=x;
        yf=y;
        zf=z;
        MRIsampleVolume(mri_source, xf, yf, zf, &val1) ;
        voxval1[indx] = val1;
        // voxval1[indx] = MRISvox(mri_source, x, y, z);
      }

  if (Gdiag & DIAG_SHOW)
    printf("M_reg (initial)\n");
  MatrixPrint(stdout,M_reg);

  M_delta = MatrixIdentity(4,NULL);
  M_delta1 = MatrixIdentity(4,NULL);
  M_delta2 = MatrixIdentity(4,NULL);
  M_delta3 = MatrixIdentity(4,NULL);
  M_delta4 = MatrixIdentity(4,NULL);
  M_delta5 = MatrixIdentity(4,NULL);
  M_delta6 = MatrixIdentity(4,NULL);

  M_reg_opt = MatrixCopy(M_reg, NULL);
  M_reg_bak = MatrixCopy(M_reg_opt, NULL);
  M_tmp = MatrixCopy(M_reg, NULL);

  best_sse = 1e30;
  for (stepindx=0; stepindx<NSTEP; stepindx++) {
    scale = step[stepindx];
    changed = 1;
    pass = 0;
    while (changed) {
      pass++;
      changed = 0;
      MatrixCopy(M_reg_opt, M_reg_bak);
      for (txi = -1; txi <= 1; txi++)
        for (tyi = -1; tyi <= 1; tyi++)
          for (tzi = -1; tzi <= 1; tzi++)
            for (axi = -1; axi <= 1; axi++)
              for (ayi = -1; ayi <= 1; ayi++)
                for (azi = -1; azi <= 1; azi++)
                  if (((txi!=0)+(tyi!=0)+(tzi!=0)+(axi!=0)+(ayi!=0)+(azi!=0))<=max_ndelta) {
                    tx = txi*dt*scale;
                    ty = tyi*dt*scale;
                    tz = tzi*dt*scale;
                    ax = axi*da*scale;
                    ay = ayi*da*scale;
                    az = azi*da*scale;
                    M_delta1->rptr[1][4]=tx;
                    M_delta1->rptr[2][4]=ty;
                    M_delta1->rptr[3][4]=tz;
                    ca = cos(ax);
                    sa = sin(ax);
                    M_delta2->rptr[2][2]=ca;
                    M_delta2->rptr[2][3]=-sa;
                    M_delta2->rptr[3][2]=sa;
                    M_delta2->rptr[3][3]=ca;
                    MatrixMultiply(M_delta2,M_delta1,M_delta5);
                    ca = cos(ay);
                    sa = sin(ay);
                    M_delta3->rptr[1][1]=ca;
                    M_delta3->rptr[1][3]=-sa;
                    M_delta3->rptr[3][1]=sa;
                    M_delta3->rptr[3][3]=ca;
                    MatrixMultiply(M_delta3,M_delta5,M_delta6);
                    ca = cos(az);
                    sa = sin(az);
                    M_delta4->rptr[1][1]=ca;
                    M_delta4->rptr[1][2]=-sa;
                    M_delta4->rptr[2][1]=sa;
                    M_delta4->rptr[2][2]=ca;
                    MatrixMultiply(M_delta4,M_delta6,M_delta);
                    MatrixMultiply(M_delta,M_reg_bak,M_reg);

                    MatrixMultiply(M_reg,vox2ras_source,M_tmp);
                    MatrixMultiply(ras2vox_target,M_tmp,vox_s2vox_t);

                    MatrixMultiply(vox_s2vox_t,voxmat1,voxmat2);
                    sse = 0;
                    for (indx=1; indx<=nvalues; indx++)
                      /*E* building the table of grid points */
                    {
                      xf=voxmat2->rptr[1][indx];
                      yf=voxmat2->rptr[2][indx];
                      zf=voxmat2->rptr[3][indx];
                      MRIsampleVolume(mri_target, xf, yf, zf, &val2) ;
                      voxval2[indx] = val2;
                      val1 = voxval1[indx];
                      err = val1-val2;
                      sse += err*err;
                    }
                    sse /= nvalues;
                    if (sse<best_sse-tol) {
                      best_sse = sse;
                      MatrixCopy(M_reg, M_reg_opt);
                      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
                        printf("%d (%d) %f %f %f %f %f %f sse = %f (%f)\n",scale,pass,tx,ty,tz,ax,
                               ay,az,sse,sqrt(sse));
                      /*
                                printf("M_delta\n"); MatrixPrint(stdout,M_delta);
                      */
                      /*
                                printf("M_delta1\n"); MatrixPrint(stdout,M_delta1);
                                printf("M_delta2\n"); MatrixPrint(stdout,M_delta2);
                                printf("M_delta3\n"); MatrixPrint(stdout,M_delta3);
                                printf("M_delta4\n"); MatrixPrint(stdout,M_delta4);
                                printf("M_delta5\n"); MatrixPrint(stdout,M_delta5);
                                printf("M_delta6\n"); MatrixPrint(stdout,M_delta6);
                      */
                      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
                        printf("M_reg_opt\n");
                        MatrixPrint(stdout,M_reg_opt);
                        printf("vox_s2vox_t\n");
                        MatrixPrint(stdout,vox_s2vox_t);
                      }
                      changed = 1;
                    }
                  }

    }
    printf("step %d: best_sse = %f (%f)\n",stepindx,best_sse,sqrt(best_sse));
    if (Gdiag & DIAG_SHOW) {
      printf("M_reg_opt\n");
      MatrixPrint(stdout,M_reg_opt);
      //      printf("vox_s2vox_t\n"); MatrixPrint(stdout,vox_s2vox_t);
    }
  }
  if (apply_xform) {
    MRI *mri_xformed ;

    mri_xformed = MRIclone(mri_source, NULL) ;
    apply_transform(mri_source, mri_target, M_reg_opt, mri_xformed) ;
    MRIwrite(mri_xformed, "xformed.mgh") ;
    MRIfree(&mri_xformed) ;
  }

  MatrixCopy(M_reg_opt, M_reg);
}
#endif

static int
apply_transform(MRI *mri_in, MRI *mri_target, MATRIX *M_reg, MRI *mri_xformed) {
  MATRIX *vox2ras_source, *ras2vox_source, *vox2ras_target, *ras2vox_target,
  *vox_s2vox_t;
  MATRIX   *m, *m_tmp;
  VECTOR   *v_in, *v_target ;
  int      x, y, z, width, height, depth ;
  double   xs, ys, zs, val ;

  vox2ras_source = MRIgetVoxelToRasXform(mri_in) ;
  vox2ras_target = MRIgetVoxelToRasXform(mri_target) ;
  ras2vox_source = MatrixInverse(vox2ras_source, NULL);
  ras2vox_target = MatrixInverse(vox2ras_target, NULL);
  // M_reg is RAS-to-RAS
  m_tmp = MatrixMultiply(M_reg,vox2ras_source, NULL);
  vox_s2vox_t =MatrixMultiply(ras2vox_target,m_tmp, NULL);
  m = MatrixInverse(vox_s2vox_t, NULL) ;
  MatrixFree(&vox2ras_source) ;
  MatrixFree(&ras2vox_source) ;
  MatrixFree(&vox2ras_target) ;
  MatrixFree(&ras2vox_target) ;
  MatrixFree(&m_tmp) ;
  MatrixFree(&vox_s2vox_t) ;
  MRIcopyPulseParameters(mri_in, mri_xformed) ;

  v_in = MatrixAlloc(4, 1, MATRIX_REAL) ;
  v_target = MatrixAlloc(4, 1, MATRIX_REAL) ;
  v_in->rptr[4][1] = v_target->rptr[4][1] = 1.0 ;

  width = mri_xformed->width ;
  height = mri_xformed->height ;
  depth = mri_xformed->depth ; //E/ Changing this from mri_target->* was correct.

  for (x = 0 ; x < width ; x++) {
    V3_X(v_target) = x ;
    for (y = 0 ; y < height ; y++) {
      V3_Y(v_target) = y ;
      for (z = 0 ; z < depth ; z++) {
        V3_Z(v_target) = z ;
        MatrixMultiply(m, v_target, v_in) ;
        xs = (V3_X(v_in)) ;
        ys = (V3_Y(v_in)) ;
        zs = (V3_Z(v_in)) ;
        MRIsampleVolume(mri_in, xs, ys, zs, &val) ;
        //E/ MRIsV() safely returns 0 for coords out-of-volume
        switch (mri_xformed->type) {
          //E/ MRI?vox() segfaults on coords out-of-volume:
        case MRI_UCHAR:
          MRIvox(mri_xformed, x, y, z) = (unsigned char)nint(val) ;
          break ;
        case MRI_SHORT:
          MRISvox(mri_xformed, x, y, z) = (short)nint(val) ;
          break ;
        case MRI_USHRT:
          MRIUSvox(mri_xformed, x, y, z) = (unsigned short)nint(val) ;
          break ;
        case MRI_FLOAT:
          MRIFvox(mri_xformed, x, y, z) = (float)(val) ;
          break ;
        default:
          ErrorExit(ERROR_UNSUPPORTED, "apply_transform: unsupported dst type %d",
                    mri_xformed->type) ;
          break ;
        }
      }
    }
  }

  MatrixFree(&v_in) ;
  MatrixFree(&v_target) ;
  MatrixFree(&m) ;

  return(NO_ERROR) ;
}



MATRIX *MRR_VoxelXformToCoronalRasXform(MRI *mri_src, MRI *mri_dst, MATRIX *M_vox_src2trg, MATRIX *M_cor_src2trg)
//E/Now using this plus MRIrasXformToVoxelXform() instead of
//MRR_LTArasTransformToCoronalRasTransform() in case anyone wants a
//voxel xform.
{
  MATRIX *D_src, *D_trg, *D_src_inv, *m_tmp;

  //E/ D's are the vox2coronalras matrices generated by
  //MRIxfmCRS2XYZtkreg().  Thanks, Doug.

  //E/ M_cor_src2trg = D_trg * M_vox_src2trg * D_src_inv
  D_src = MRIxfmCRS2XYZtkreg(mri_src);
  D_trg = MRIxfmCRS2XYZtkreg(mri_dst);
  D_src_inv = MatrixInverse(D_src, NULL);

  m_tmp = MatrixMultiply(D_trg, M_vox_src2trg, NULL);
  M_cor_src2trg = MatrixMultiply(m_tmp, D_src_inv, NULL);

//#define _MRR_VX2CRX
#ifdef _MRR_VX2CRX
  fprintf(stderr, "Entering MRR_VoxelXformToCoronalRasXform:\n");

  fprintf(stderr, "M_vox_src2trg = \n");
  MatrixPrint(stderr, M_vox_src2trg);

  fprintf(stderr, "D_src = \n");
  MatrixPrint(stderr, D_src);
  fprintf(stderr, "D_src_inv = \n");
  MatrixPrint(stderr, D_src_inv);
  fprintf(stderr, "D_trg = \n");
  MatrixPrint(stderr, D_trg);

  fprintf(stderr, "M_cor_src2trg = D_trg * M_vox_src2trg * D_src_inv = \n");
  MatrixPrint(stderr, M_cor_src2trg);

  fprintf(stderr, "Exiting MRR_VoxelXformToCoronalRasXform.\n");
#endif

  return M_cor_src2trg;
}



MRI *
MRR_MRIhybrid(MRI *mri_src, MRI *mri_target, MRI *mri_xformed)
//E/ Using this (when transforming) instead of how it used to be,
//(MRIclone for -a) or (used to be MRIalloc for -w).  Differences are:
//MRIclone misses src->type and MRIalloc misses MRIcopyHeader().
{
  if (!mri_xformed)
    mri_xformed =
      MRIallocSequence(mri_target->width,
                       mri_target->height,
                       mri_target->depth,
                       mri_src->type, //E/ !!
                       mri_target->nframes);
  MRIcopyHeader(mri_target, mri_xformed) ; //E/ !!
  //  mri_xformed->c_r = 0.0;
  //  mri_xformed->c_a = 0.0;
  //  mri_xformed->c_s = 0.0;
  return(mri_xformed) ;
}


//E/ These three aren't used but they must be important to someone.

#if 0

MRI *
MRIsadd(MRI *mri1, MRI *mri2, MRI *mri_dst) {
  int     width, height, depth, x, y, z ;
  short   *p1, *p2, *pdst ;

  width = mri1->width ;
  height = mri1->height ;
  depth = mri1->depth ;

  if (!mri_dst) {
    mri_dst = MRIalloc(width, height, depth, mri1->type) ;
    MRIcopyHeader(mri1, mri_dst) ;
  }

  for (z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      p1 = &MRISvox(mri1, 0, y, z) ;
      p2 = &MRISvox(mri2, 0, y, z) ;
      pdst = &MRISvox(mri_dst, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
        *pdst++ = *p1++ + *p2++ ;
    }
  }
  return(mri_dst) ;
}

MRI *
MRIssqrt(MRI *mri_src, MRI *mri_dst) {
  int     width, height, depth, x, y, z, frame ;
  short   *psrc, *pdst ;
  unsigned short   *psrc2, *pdst2 ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  for (frame = 0 ; frame < mri_src->nframes ; frame++) {
    for (z = 0 ; z < depth ; z++) {
      for (y = 0 ; y < height ; y++) {
        switch (mri_src->type) {
        case MRI_SHORT:
          psrc = &MRISseq_vox(mri_src, 0, y, z, frame) ;
          pdst = &MRISseq_vox(mri_dst, 0, y, z, frame) ;
          for (x = 0 ; x < width ; x++)
            *pdst++ = sqrt(*psrc++) ;
          break ;
        case MRI_USHRT:
          psrc2 = &MRIUSseq_vox(mri_src, 0, y, z, frame) ;
          pdst2 = &MRIUSseq_vox(mri_dst, 0, y, z, frame) ;
          for (x = 0 ; x < width ; x++)
            *pdst2++ = sqrt(*psrc2++) ;
          break ;
        default:
          ErrorReturn(NULL,
                      (ERROR_UNSUPPORTED,
                       "MRIssqrt: unsupported type %d", mri_src->type)) ;
        }
      }
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIsscalarMul(MRI *mri_src, MRI *mri_dst, float scalar) {
  int     width, height, depth, x, y, z, frame ;
  short   *psrc, *pdst ;
  unsigned short   *psrc2, *pdst2 ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  for (frame = 0 ; frame < mri_src->nframes ; frame++) {
    for (z = 0 ; z < depth ; z++) {
      for (y = 0 ; y < height ; y++) {
        switch (mri_src->type) {
        case MRI_SHORT:
          psrc = &MRISseq_vox(mri_src, 0, y, z, frame) ;
          pdst = &MRISseq_vox(mri_dst, 0, y, z, frame) ;
          for (x = 0 ; x < width ; x++)
            *pdst++ = *psrc++ * scalar ;
          break ;
        case MRI_USHRT:
          psrc2 = &MRIUSseq_vox(mri_src, 0, y, z, frame) ;
          pdst2 = &MRIUSseq_vox(mri_dst, 0, y, z, frame) ;
          for (x = 0 ; x < width ; x++)
            *pdst2++ = *psrc2++ * scalar ;
          break ;
        default:
          ErrorReturn(NULL,
                      (ERROR_UNSUPPORTED,
                       "MRIsscalarMul: unsupported type %d", mri_src->type)) ;
        }
      }
    }
  }
  return(mri_dst) ;
}

#endif




#if 0
/*E* this below doesn't belong here but it's code to steal from or a
  function to call(?) - MRIsynthesizeWeightedVolume was a tip/
  guideline for a way to rewrite mri_rigid_register "if you're feeling
  bold" - above, we use MRIfindApproximateSkullBoundingBox() instead
  and take the superior 75% */

MRI *
MRIsynthesizeWeightedVolume(MRI *mri_T1, MRI *mri_PD, float w5, float TR5,
                            float w30, float TR30, float target_wm, float TE) {
  MRI *mri_dst ;
  MRI_REGION box ;
  float      x0, y0, z0, min_real_val ;
  HISTOGRAM *h_mri, *h_smooth ;
  int        mri_peak, n, min_real_bin, x, y, z, width, height, depth ;
  MRI       *mri30, *mri5 ;
  double    mean_PD ;
  double    val30, val5, val ;

  mean_PD = MRImeanFrame(mri_PD, 0) ;
  /*  MRIscalarMul(mri_PD, mri_PD, 1000.0f/mean_PD) ;*/
  mri30 = MRIsynthesize(mri_T1, mri_PD, NULL, TR30, RADIANS(30), TE) ;
  mri5 = MRIsynthesize(mri_T1, mri_PD, NULL, TR5, RADIANS(5), TE) ;
  width = mri30->width ;
  height = mri30->height ;
  depth = mri30->depth ;

  mri_dst = MRIalloc(width, height, depth, MRI_FLOAT) ;
  MRIcopyHeader(mri_T1, mri_dst) ;

  h_mri = MRIhistogram(mri30, 100) ;
  h_smooth = HISTOsmooth(h_mri, NULL, 2) ;
  mri_peak = HISTOfindHighestPeakInRegion(h_smooth, 0, h_smooth->nbins) ;
  min_real_bin = HISTOfindNextValley(h_smooth, mri_peak) ;
  min_real_val = h_smooth->bins[min_real_bin] ;

  MRIfindApproximateSkullBoundingBox(mri30, min_real_val, &box) ;
  x0 = box.x+box.dx/3 ;
  y0 = box.y+box.dy/3 ;
  z0 = box.z+box.dz/2 ;
  printf("using (%.0f, %.0f, %.0f) as brain centroid...\n",x0, y0, z0) ;
  box.dx /= 4 ;
  box.x = x0 - box.dx/2;
  box.dy /= 4 ;
  box.y = y0 - box.dy/2;
  box.dz /= 4 ;
  box.z = z0 - box.dz/2;


  printf("using box (%d,%d,%d) --> (%d, %d,%d) "
         "to find MRI wm\n", box.x, box.y, box.z,
         box.x+box.dx-1,box.y+box.dy-1, box.z+box.dz-1) ;

  h_mri = MRIhistogramRegion(mri30, 0, NULL, &box) ;
  for (n = 0 ; n < h_mri->nbins-1 ; n++)
    if (h_mri->bins[n+1] > min_real_val)
      break ;
  HISTOclearBins(h_mri, h_mri, 0, n) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    HISTOplot(h_mri, "mri.histo") ;
  mri_peak = HISTOfindLastPeak(h_mri, HISTO_WINDOW_SIZE,MIN_HISTO_PCT);
  mri_peak = h_mri->bins[mri_peak] ;
  printf("before smoothing, mri peak at %d\n", mri_peak) ;
  h_smooth = HISTOsmooth(h_mri, NULL, 2) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    HISTOplot(h_smooth, "mri_smooth.histo") ;
  mri_peak = HISTOfindLastPeak(h_smooth, HISTO_WINDOW_SIZE,MIN_HISTO_PCT);
  mri_peak = h_mri->bins[mri_peak] ;
  printf("after smoothing, mri peak at %d\n", mri_peak) ;
  HISTOfree(&h_smooth) ;
  HISTOfree(&h_mri) ;


  return(mri_dst) ;
}
#endif
