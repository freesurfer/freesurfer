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


static double FLASHforwardModel(double flip_angle, double TR, double PD, 
                                double T1) ;

static double estimate_ms_params(MRI **mri_flash, MRI **mri_flash_synth, int nvolumes, MRI *mri_T1, MRI *mri_PD, MRI *mri_sse, MATRIX **M_reg);
static void estimate_rigid_regmatrix(MRI *mri_source, MRI *mri_traget, MATRIX *M_reg);

static float        tr = 0, te = 0, fa = 0 ;

static int write_iterations=0;

int
main(int argc, char *argv[])
{
  char   **av, fname[STRLEN] ;
  int    ac, nargs, i ;
  MRI    *mri_flash[MAX_IMAGES], *mri_flash_synth[MAX_IMAGES], *mri_T1, *mri_PD, *mri_sse ;
  char   *in_fname, *out_dir ;
  int    msec, minutes, seconds;
  int    nvolumes ;
  struct timeb start ;
  MATRIX *M_reg[MAX_IMAGES];
  double rms ;

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

  if (argc < 4)
    usage_exit(1) ;

  out_dir = argv[argc-1] ;
#if 0
  FileNameOnly(out_T1_fname, fname) ;
  FileNameRemoveExtension(fname, fname) ;
  strcpy(parms.base_name, fname) ;
#endif

  nvolumes = 0 ;
  for (i = 1 ; i < argc-1 ; i++)
  {
    if (argv[i][0] == '-')
    {
      if (!stricmp(argv[i]+1, "te"))
        te = atof(argv[i+1]) ;
      else if (!stricmp(argv[i]+1, "tr"))
        tr = atof(argv[i+1]) ;
      else if (!stricmp(argv[i]+1, "fa"))
        fa = atof(argv[i+1]) ;
      else
        ErrorExit(ERROR_BADPARM, "%s: unsupported MR parameter %s",
                  Progname, argv[i]+1) ;
      i++ ;  /* skip parameter */
      continue ;
    }

    in_fname = argv[i] ;
    printf("reading %s...", in_fname) ;

    mri_flash[nvolumes] = MRIread(in_fname) ;
    mri_flash[nvolumes]->register_mat = MRIgetVoxelToRasXform(mri_flash[nvolumes]);

    if (!mri_flash[nvolumes])
      ErrorExit(Gerror, "%s: MRIread(%s) failed", Progname, in_fname) ;
    if (tr > 0)
    {
      mri_flash[nvolumes]->tr = tr ;
    }
    if (te > 0)
    {
      mri_flash[nvolumes]->te = te ;
    }
    if (fa > 0)
    {
      mri_flash[nvolumes]->flip_angle = fa ;
    }
    printf("TE = %2.2f, TR = %2.2f, alpha = %2.2f\n", mri_flash[nvolumes]->te, 
           mri_flash[nvolumes]->tr, mri_flash[nvolumes]->flip_angle) ;
    mri_flash[nvolumes]->flip_angle = RADIANS(mri_flash[nvolumes]->flip_angle);
    if (conform)
    {
      MRI *mri_tmp ;

      printf("embedding and interpolating volume\n") ;
      mri_tmp = MRIconform(mri_flash[nvolumes]) ;
      /*      MRIfree(&mri_src) ;*/
      mri_flash[nvolumes] = mri_tmp ;
    }
    mri_flash_synth[nvolumes] = MRIclone(mri_flash[nvolumes], NULL) ;
    mri_flash_synth[nvolumes]->register_mat = MRIgetVoxelToRasXform(mri_flash[nvolumes]);
    if (FZERO(mri_flash[nvolumes]->tr) || 
        FZERO(mri_flash[nvolumes]->flip_angle))
      ErrorExit(ERROR_BADPARM, "%s: invalid TR or FA for image %d:%s",
                Progname, nvolumes, in_fname) ;
    nvolumes++ ;
  }
  printf("using %d FLASH volumes to estimate tissue parameters.\n", nvolumes) ;
  mri_T1 = MRIclone(mri_flash[0], NULL) ;
  mri_PD = MRIclone(mri_flash[0], NULL) ;
  mri_sse = MRIclone(mri_flash[0], NULL) ;

  {
    int j, iter ; /* should exit when M_reg's converge. */
    LTA *lta;

    Progname = argv[0] ;
    ErrorInit(NULL, NULL, NULL) ;
    DiagInit(NULL, NULL, NULL) ;

    for (j=0;j<nvolumes;j++) M_reg[j] = MatrixIdentity(4,(MATRIX *)NULL);
 
    for (iter=0; iter<niter; iter++) 
    {
      printf("parameter estimation/motion correction iteration %d of %d\n", iter+1, niter) ;
      rms = estimate_ms_params(mri_flash, mri_flash_synth, nvolumes, mri_T1, mri_PD, mri_sse, M_reg) ;
      printf("parameter rms = %2.3f\n", rms) ;
  
      for (j=0;j<nvolumes;j++)
      {
        printf("estimating rigid alignment for FLASH volume #%d...\n", j+1) ;
        estimate_rigid_regmatrix(mri_flash_synth[j],mri_flash[j],M_reg[j]);
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
        {
          printf("M_reg[%d]\n",j); MatrixPrint(stdout,M_reg[j]);
        }
  
      }
 
      if ((write_iterations>0) && (iter%write_iterations==0))
      {
        sprintf(fname,"%s/T1-%d.mgh",out_dir,iter);
        printf("writing T1 esimates to %s...\n", fname) ;
        MRIwrite(mri_T1, fname) ;
        sprintf(fname,"%s/PD-%d.mgh",out_dir,iter);
        printf("writing PD estimates to %s...\n", fname) ;
        MRIwrite(mri_PD, fname) ;
        sprintf(fname,"%s/sse-%d.mgh",out_dir,iter);
        printf("writing residual sse to %s...\n", fname) ;
        MRIwrite(mri_sse, fname) ;
 
        for (j=0;j<nvolumes;j++)
        {
          sprintf(fname,"%s/vol%d-%d.mgh",out_dir,j,iter);
          printf("writing synthetic images to %s...\n", fname);
          MRIwrite(mri_flash_synth[j], fname) ;
          sprintf(fname,"%s/vol%d-%d.lta",out_dir,j,iter); 
          printf("writing regisration matrix to %s...\n", fname);
          lta = LTAalloc(1,NULL) ;
          MatrixCopy(M_reg[j],lta->xforms[0].m_L) ;
          lta->type = LINEAR_RAS_TO_RAS ;
          LTAwrite(lta,fname) ;
        }
      }
    }
    sprintf(fname,"%s/T1.mgh",out_dir);
    printf("writing T1 esimates to %s...\n", fname) ;
    MRIwrite(mri_T1, fname) ;
    sprintf(fname,"%s/PD.mgh",out_dir);
    printf("writing PD estimates to %s...\n", fname) ;
    MRIwrite(mri_PD, fname) ;
    sprintf(fname,"%s/sse.mgh",out_dir);
    printf("writing residual sse to %s...\n", fname) ;
    MRIwrite(mri_sse, fname) ;

    for (j=0;j<nvolumes;j++)
    {
      sprintf(fname,"%s/vol%d.mgh",out_dir,j);
      printf("writing synthetic images to %s...\n", fname);
      MRIwrite(mri_flash_synth[j], fname) ;
      sprintf(fname,"%s/vol%d.lta",out_dir,j);
      printf("writing registration matrix to %s...\n", fname);
      lta = LTAalloc(1,NULL) ;
      MatrixCopy(M_reg[j],lta->xforms[0].m_L) ;
      lta->type = LINEAR_RAS_TO_RAS ;
      LTAwrite(lta,fname) ;
    }
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
  else if (!stricmp(option, "tr"))
  {
    tr = atof(argv[2]) ;
    nargs = 1;
  }
  else if (!stricmp(option, "te"))
  {
    te = atof(argv[2]) ;
    nargs = 1;
  }
  else if (!stricmp(option, "fa"))
  {
    fa = atof(argv[2]) ;
    nargs = 1;
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
    write_iterations = atoi(argv[2]) ;
    printf("writing out intermediate results every %d iterations...\n",
           write_iterations) ;
    nargs = 1 ;
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
  printf("usage: %s [options] <volume> ... <output T1 volume> <output PD volume>\n", Progname) ;
  printf("\t-a    rigid alignment of input volumes before averaging\n") ;
  exit(code) ;
}


static double
FLASHforwardModel(double flip_angle, double TR, double PD, double T1)
{
  double FLASH, E1 ;
  double  CFA, SFA ;


  CFA = cos(flip_angle) ; SFA = sin(flip_angle) ;
  E1 = exp(-TR/T1) ;
      
  FLASH = PD * SFA ;
  if (!DZERO(T1))
    FLASH *= (1-E1)/(1-CFA*E1);
  return(FLASH) ;
}

#define T1_MIN 10.0
#define T1_MAX 5000.0
#define MAX_NVALS 10000
#define MAX_NVOLS 100

static double
estimate_ms_params(MRI **mri_flash, MRI **mri_flash_synth, int nvolumes, MRI *mri_T1, MRI *mri_PD, MRI *mri_sse, MATRIX **M_reg)
{
  double   SignalTableValues[MAX_NVALS][MAX_NVOLS], SignalTableT1[MAX_NVALS], SignalTableNorm[MAX_NVALS], ImageValues[MAX_NVOLS], total_sse ;
  double   se, best_se, ss, sse, err, val, norm, T1, PD, xf, yf, zf ;
  int      i, j, x, y, z, indx, min_indx, max_indx, best_indx, center_indx, stepindx;
  int      width=mri_T1->width, height=mri_T1->height, depth=mri_T1->depth, nvalues=5000, nevals;
  int      nstep=11, step[11]={1024,512,256,128,64,32,16,8,4,2,1};
  MRI      *mri ;
  MATRIX   *vox2ras[MAX_NVOLS], *ras2vox[MAX_NVOLS], *voxvec1, *voxvec2, *rasvec1, *rasvec2;

  voxvec1 = MatrixAlloc(4,1,MATRIX_REAL); voxvec1->rptr[4][1] = 1.0;
  voxvec2 = MatrixCopy(voxvec1, NULL);
  rasvec1 = MatrixCopy(voxvec1, NULL);
  rasvec2 = MatrixCopy(voxvec1, NULL);
  for (j = 0 ; j < nvolumes ; j++)
  {
    vox2ras[j] = MatrixCopy(mri_flash[j]->register_mat, NULL);
    ras2vox[j] = MatrixInverse(vox2ras[j], NULL);
  }

  PD = 1;
  for (i=0; i<nvalues; i++)
  {
    T1 = SignalTableT1[i] = T1_MIN+i*(T1_MAX-T1_MIN)/(nvalues-1);  
    ss = 0;
    for (j = 0 ; j < nvolumes ; j++)
    {
      mri = mri_flash[j] ;
      val = SignalTableValues[i][j] = FLASHforwardModel(mri->flip_angle, mri->tr, PD, T1) ;
      ss += val*val;
    }
    norm = SignalTableNorm[i] = sqrt(ss);
    if (norm>0) for (j = 0 ; j < nvolumes ; j++) SignalTableValues[i][j] /= norm;
  }

  total_sse = 0 ;
  for (z = 0 ; z < depth ; z++)
  for (y = 0 ; y < height ; y++)
  for (x = 0 ; x < width ; x++)
  {
    ss = 0;
    for (j = 0 ; j < nvolumes ; j++)
    {
      mri = mri_flash[j] ;
#if 0
      val = ImageValues[j] = MRISvox(mri, x, y, z) ;
#endif
      voxvec1->data[0]=x; voxvec1->data[1]=y; voxvec1->data[2]=z;
      MatrixMultiply(vox2ras[j],voxvec1,rasvec1);
      MatrixMultiply(M_reg[j],rasvec1,rasvec2);
      MatrixMultiply(ras2vox[j],rasvec2,voxvec2);
      xf=voxvec2->data[0]; yf=voxvec2->data[1]; zf=voxvec2->data[2];
      MRIsampleVolume(mri, xf, yf, zf, &val) ;
      ImageValues[j] = val;
      ss += val*val;
    }
    norm = sqrt(ss);
    if (norm>0) for (j = 0 ; j < nvolumes ; j++) ImageValues[j] /= norm;
   
    min_indx = best_indx = 0;
    max_indx = nvalues-1; 
    best_indx = -1;
    center_indx = -1;
    best_se = 10000000;
    nevals = 0;
    for (stepindx=0; stepindx<nstep; stepindx++) 
    {
      for (indx=min_indx; indx<=max_indx; indx+=step[stepindx])
      if (indx!=center_indx)
      {
        se = 0;
        for (j = 0 ; j < nvolumes ; j++)
        {
          err = ImageValues[j]-SignalTableValues[indx][j];
          se += err*err; 
        }
        if (se<best_se)
        {
          best_se = se;
          best_indx = indx;
        }
        nevals++;
      }
      min_indx = MAX(best_indx-step[stepindx]/2,1);
      max_indx = MIN(best_indx+step[stepindx]/2,nvalues);
      center_indx = best_indx;
    }
    T1 = MRISvox(mri_T1, x, y, z) = (short)nint(SignalTableT1[best_indx]);
    PD = MRISvox(mri_PD, x, y, z) = (short)nint(norm/SignalTableNorm[best_indx]);
    for (j = 0 ; j < nvolumes ; j++)
    {
      mri = mri_flash_synth[j] ;
      MRISvox(mri, x, y, z) = (short)nint(PD*SignalTableNorm[best_indx]*SignalTableValues[best_indx][j]);
    }
    sse = 0;
    for (j = 0 ; j < nvolumes ; j++)
    {
      err = MRISvox(mri_flash_synth[j], x, y, z)-ImageValues[j]*norm;
      sse += err*err; 
    }
    total_sse += sse ;
    MRISvox(mri_sse, x, y, z) = (short)nint(sqrt(sse));
  }
  total_sse = sqrt(total_sse / (width*height*depth)) ;
  return(total_sse) ;
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

static void
estimate_rigid_regmatrix(MRI *mri_source, MRI *mri_target, MATRIX *M_reg)
{
  double   xf, yf, zf, tx, ty, tz, ax, ay, az, ca, sa, val1, val2, err, sse, best_sse, dt=0.01, da=RADIANS(0.005), tol=0.00001;
  int      x, y, z, txi, tyi, tzi, axi, ayi, azi, indx, stepindx, changed, pass;
  int      width=mri_source->width, height=mri_source->height, depth=mri_source->depth, dx=10, dy=10, dz=10, nvalues;
#if 1
/*
  int      nstep=6, step[6]={32,16,8,4,2,1}, scale;
*/
  int      nstep=5, step[5]={16,8,4,2,1}, scale;
#else
  int      nstep=1, step[1]={1}, scale;
#endif
  MATRIX   *vox2ras_source, *ras2vox_source, *vox2ras_target, *ras2vox_target, *vox_s2vox_t;
  MATRIX   *M_reg_bak, *M_reg_opt, *M_tmp, *M_delta, *M_delta1, *M_delta2, *M_delta3, *M_delta4, *M_delta5, *M_delta6;
  MATRIX   *voxmat1, *voxmat2;
  double   voxval1[MAX_VOX], voxval2[MAX_VOX];

  vox2ras_source = MatrixCopy(mri_source->register_mat, NULL);
  vox2ras_target = MatrixCopy(mri_target->register_mat, NULL);
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
  for (stepindx=0; stepindx<nstep; stepindx++) 
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
/*
          printf("M_reg\n"); MatrixPrint(stdout,M_reg);
          printf("vox_s2vox_t\n"); MatrixPrint(stdout,vox_s2vox_t);
*/
          changed = 1;
        }
      }
     
    }
    printf("step %d: sse = %f (%f)\n",stepindx,sse,sqrt(sse));
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
