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

static int apply_transform(MRI *mri_in, MRI *mri_target, MATRIX *M_reg, MRI *mri_xformed) ;
MRI *MRIsadd(MRI *mri1, MRI *mri2, MRI *mri_dst) ;
MRI *MRIsscalarMul(MRI *mri_src, MRI *mri_dst, float scalar) ;
MRI *MRIssqrt(MRI *mri_src, MRI *mri_dst) ;

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;
static int align = 1 ;
static int window_flag = 0 ;
static MORPH_PARMS  parms ;

static void usage_exit(int code) ;

static int niter=10 ;
static float thresh = 25 ;
static int conform = 0 ;
static int sinc_flag = 1;
static int sinchalfwindow = 3;

static char *residual_name = NULL ;

#define MAX_IMAGES 100
#define MIN_T1   5  /* avoid singularity at T1=0 */
#define MIN_PD   5
#define MIN_ITER 5
#define MAX_ITER 5000

#define MAX_WRITE 10000
static int nwrite = 0 ;
static char *write_volumes[MAX_WRITE] ;
static char *out_fnames[MAX_WRITE] ;

static void estimate_rigid_regmatrix(MRI *mri_source, MRI *mri_traget, MATRIX *M_reg);

int
main(int argc, char *argv[])
{
  char   **av ;
  int    ac, nargs, i ;
  MRI    *mri_src, *mri_target ;
  char   *src_fname, *out_fname, *target_fname ;
  int    msec, minutes, seconds;
  struct timeb start ;
  MATRIX *M_reg ;
  LTA    *lta ;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  TimerStart(&start) ;
  
  parms.dt = 1e-6 ;
  parms.tol = 1e-5 ;
  parms.momentum = 0.0 ;
  parms.niterations = 20 ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    usage_exit(1) ;

  src_fname = argv[1] ;
  target_fname = argv[2] ;
  out_fname = argv[3] ;
  
  printf("reading source from %s...\n", src_fname) ;
  mri_src = MRIread(src_fname) ;
  printf("reading target from %s...\n", target_fname) ;
  mri_target = MRIread(target_fname) ;
  if (!mri_src)
    ErrorExit(ERROR_NOFILE, "%s: could not read src MRI from %s", src_fname) ;
  if (!mri_target)
    ErrorExit(ERROR_NOFILE, "%s: could not read dst MRI from %s", target_fname) ;

#if 0
  if (conform)
  {
    MRI *mri_tmp ;
    
    printf("embedding and interpolating volume\n") ;
    mri_tmp = MRIconform(mri_flash[nvolumes]) ;
    /*      MRIfree(&mri_src) ;*/
    mri_flash[nvolumes] = mri_tmp ;
  }
#endif

  M_reg = MatrixIdentity(4,(MATRIX *)NULL);
  
  estimate_rigid_regmatrix(mri_src,mri_target,M_reg);
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    printf("M_reg\n"); MatrixPrint(stdout,M_reg);
  }
  
  printf("writing registration matrix to %s...\n", out_fname);
  lta = LTAalloc(1,NULL) ;
  MatrixCopy(M_reg,lta->xforms[0].m_L) ;
  lta->type = LINEAR_RAS_TO_RAS ;
  LTAwrite(lta,out_fname) ;
  MRIfree(&mri_src) ;

  for (i = 0 ; i < nwrite ; i++)
  {
    MRI *mri_in, *mri_xformed ;

    mri_in = MRIread(write_volumes[i]) ;
    if (!mri_in)
      ErrorExit(ERROR_NOFILE, "%s: could not read volume %s", Progname, write_volumes[i]) ;
    mri_xformed = MRIalloc(mri_target->width, mri_target->height, mri_target->depth,
                           mri_in->type) ;
    apply_transform(mri_in, mri_target, M_reg, mri_xformed) ;
    printf("writing transformed volume to %s...\n", out_fnames[i]) ;
    MRIwrite(mri_xformed, out_fnames[i]) ;
    MRIfree(&mri_xformed) ; MRIfree(&mri_in) ;
  }
  
  msec = TimerStop(&start) ;
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
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;
  
  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "dt"))
  {
    parms.dt = atof(argv[2]) ;
    nargs = 1 ;
    printf("using dt = %2.3e\n", parms.dt) ;
  }
  else if (!stricmp(option, "tol"))
  {
    parms.tol = atof(argv[2]) ;
    nargs = 1 ;
    printf("using tol = %2.3e\n", parms.tol) ;
  }
  else if (!stricmp(option, "conform"))
  {
    conform = 1 ;
    printf("interpolating volume to be isotropic 1mm^3\n") ;
  }
  else if (!stricmp(option, "sinc"))
  {
    sinchalfwindow = atoi(argv[2]);
    sinc_flag = 1;
    nargs = 1;
    printf("using sinc interpolation with windowwidth of %d\n",
      2*sinchalfwindow);
  }
  else if (!stricmp(option, "trilinear"))
  {
    sinc_flag = 0;
    printf("using trilinear interpolation\n");
  }
  else if (!stricmp(option, "window"))
  {
    window_flag = 1 ;
    printf("applying hanning window to volumes...\n") ;
  }
  else if (!stricmp(option, "noconform"))
  {
    conform = 0 ;
    printf("inhibiting isotropic volume interpolation\n") ;
  }
  else switch (toupper(*option))
  {
  case 'N':
    niter = atoi(argv[2]) ;
    printf("performing estimation/motion correction %d times...\n", niter) ;
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
  case 'M':
    parms.momentum = atof(argv[2]) ;
    nargs = 1 ;
    printf("using momentum = %2.3f\n", parms.momentum) ;
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
  case 'A':
    align = 1 ;
    printf("aligning volumes before averaging...\n") ;
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
usage_exit(int code)
{
  printf("usage: %s [options] <src volume> <target volume> <transform fname>\n", Progname) ;
  exit(code) ;
}


MRI *
MRIsadd(MRI *mri1, MRI *mri2, MRI *mri_dst)
{
  int     width, height, depth, x, y, z ;
  short   *p1, *p2, *pdst ;

  width = mri1->width ;
  height = mri1->height ;
  depth = mri1->depth ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, mri1->type) ;
    MRIcopyHeader(mri1, mri_dst) ;
  }

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      p1 = &MRISvox(mri1, 0, y, z) ;
      p2 = &MRISvox(mri2, 0, y, z) ;
      pdst = &MRISvox(mri_dst, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
        *pdst++ = *p1++ + *p2++ ;
    }
  }
  return(mri_dst) ;
}

#define MAX_VOX 262145

#define NSTEP   11

static void
estimate_rigid_regmatrix(MRI *mri_source, MRI *mri_target, MATRIX *M_reg)
{
  double   xf, yf, zf, tx, ty, tz, ax, ay, az, ca, sa, val1, val2, err, sse, best_sse, dt=0.01, da=RADIANS(0.005), tol=0.00001;
  int      x, y, z, txi, tyi, tzi, axi, ayi, azi, indx, stepindx, changed, pass;
  int      width=mri_source->width, height=mri_source->height, depth=mri_source->depth, dx=10, dy=10, dz=10, nvalues;
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

  nvalues = 0;
  for (z = 0 ; z < depth ; z++)
  for (y = 0 ; y < height ; y++)
  for (x = 0 ; x < width ; x++)
  {
    if ((x%dx==0)&&(y%dy==0)&&(z%dz==0))
    {
      nvalues++;
    }
  }

  voxmat1 = MatrixAlloc(4,nvalues,MATRIX_REAL); 
  voxmat2 = MatrixCopy(voxmat1, NULL);

  indx = 0; 
  for (z = 0 ; z < depth ; z++)
  for (y = 0 ; y < height ; y++)
  for (x = 0 ; x < width ; x++)
  { 
    if ((x%dx==0)&&(y%dy==0)&&(z%dz==0))
    { 
      indx++;
      voxmat1->rptr[1][indx] = x;
      voxmat1->rptr[2][indx] = y;
      voxmat1->rptr[3][indx] = z;
      voxmat1->rptr[4][indx] = 1;
      voxval1[indx] = MRISvox(mri_source, x, y, z); 
    }
  }
  
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

  best_sse = 10000000;
  for (stepindx=0; stepindx<NSTEP; stepindx++) 
  {
    scale = step[stepindx];
    changed = 1;
    pass = 0;
    while (changed)
    {
      pass++;
      changed = 0;
      MatrixCopy(M_reg_opt, M_reg_bak);
      for (txi = -1; txi <= 1; txi++)
      for (tyi = -1; tyi <= 1; tyi++)
      for (tzi = -1; tzi <= 1; tzi++)
      for (axi = -1; axi <= 1; axi++)
      for (ayi = -1; ayi <= 1; ayi++)
      for (azi = -1; azi <= 1; azi++)
      {
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
        {
          xf=voxmat2->rptr[1][indx]; yf=voxmat2->rptr[2][indx]; zf=voxmat2->rptr[3][indx];
          MRIsampleVolume(mri_target, xf, yf, zf, &val2) ;
          voxval2[indx] = val2;
          val1 = voxval1[indx];
          err = val1-val2;
          sse += err*err;
        }
        sse /= nvalues;
        if (sse<best_sse-tol)
        {
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
          if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          {
            printf("M_reg\n"); MatrixPrint(stdout,M_reg);
            printf("vox_s2vox_t\n"); MatrixPrint(stdout,vox_s2vox_t);
          }
          changed = 1;
        }
      }
     
    }
    printf("step %d: sse = %f (%f)\n",stepindx,sse,sqrt(sse));
    if (Gdiag & DIAG_SHOW)
    {
      printf("M_reg\n"); MatrixPrint(stdout,M_reg);
      printf("vox_s2vox_t\n"); MatrixPrint(stdout,vox_s2vox_t);
    }
  }
  MatrixCopy(M_reg_opt, M_reg);
}

MRI *
MRIssqrt(MRI *mri_src, MRI *mri_dst)
{
  int     width, height, depth, x, y, z, frame ;
  short   *psrc, *pdst ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  for (frame = 0 ; frame < mri_src->nframes ; frame++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        switch (mri_src->type)
        {
        case MRI_SHORT:
          psrc = &MRISseq_vox(mri_src, 0, y, z, frame) ;
          pdst = &MRISseq_vox(mri_dst, 0, y, z, frame) ;
          for (x = 0 ; x < width ; x++)
            *pdst++ = sqrt(*psrc++) ;
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
MRIsscalarMul(MRI *mri_src, MRI *mri_dst, float scalar)
{
  int     width, height, depth, x, y, z, frame ;
  short   *psrc, *pdst ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  for (frame = 0 ; frame < mri_src->nframes ; frame++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        switch (mri_src->type)
        {
        case MRI_SHORT:
          psrc = &MRISseq_vox(mri_src, 0, y, z, frame) ;
          pdst = &MRISseq_vox(mri_dst, 0, y, z, frame) ;
          for (x = 0 ; x < width ; x++)
            *pdst++ = *psrc++ * scalar ;
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
static int
apply_transform(MRI *mri_in, MRI *mri_target, MATRIX *M_reg, MRI *mri_xformed)
{
  MATRIX   *vox2ras_source, *ras2vox_source, *vox2ras_target, *ras2vox_target, *vox_s2vox_t;
  MATRIX   *m, *m_tmp ;
  VECTOR   *v_in, *v_target ;
  int      x, y, z, width, height, depth ;
  Real     xs, ys, zs, val ;


  vox2ras_source = MRIgetVoxelToRasXform(mri_in) ;
  vox2ras_target = MRIgetVoxelToRasXform(mri_target) ;
  ras2vox_source = MatrixInverse(vox2ras_source, NULL);
  ras2vox_target = MatrixInverse(vox2ras_target, NULL);
  m_tmp = MatrixMultiply(M_reg,vox2ras_source, NULL);
  vox_s2vox_t =MatrixMultiply(ras2vox_target,m_tmp, NULL);
  m = MatrixInverse(vox_s2vox_t, NULL) ;
  MatrixFree(&m_tmp) ;

  v_in = MatrixAlloc(4, 1, MATRIX_REAL) ;
  v_target = MatrixAlloc(4, 1, MATRIX_REAL) ;
  v_in->rptr[4][1] = v_target->rptr[4][1] = 1.0 ;
  width = mri_target->width ; height = mri_target->height ; depth = mri_target->depth ;
  for (x = 0 ; x < width ; x++)
  {
    V3_X(v_target) = x ;
    for (y = 0 ; y < height ; y++)
    {
      V3_Y(v_target) = y ;
      for (z = 0 ; z < depth ; z++)
      {
        V3_Z(v_target) = z ;
        MatrixMultiply(m, v_target, v_in) ;
        xs = (V3_X(v_in)) ; ys = (V3_Y(v_in)) ; zs = (V3_Z(v_in)) ;
        MRIsampleVolume(mri_in, xs, ys, zs, &val) ;
        switch (mri_xformed->type)
        {
        case MRI_UCHAR:
          MRIvox(mri_xformed, x, y, z) = (unsigned char)nint(val) ;
          break ;
        case MRI_SHORT:
          MRISvox(mri_xformed, x, y, z) = (short)nint(val) ;
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

  MatrixFree(&v_in) ; MatrixFree(&v_target) ;
  MatrixFree(&vox2ras_source) ;
  MatrixFree(&vox2ras_target) ;
  MatrixFree(&ras2vox_source) ;
  MatrixFree(&ras2vox_target) ;
  MatrixFree(&vox_s2vox_t) ;
  MatrixFree(&m) ;
  return(NO_ERROR) ;
}
