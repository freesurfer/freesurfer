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

//E/ Maybe should update these three to not demand MRI_SHORT - but they're never called.
MRI *MRIsadd(MRI *mri1, MRI *mri2, MRI *mri_dst) ;
MRI *MRIsscalarMul(MRI *mri_src, MRI *mri_dst, float scalar) ;
MRI *MRIssqrt(MRI *mri_src, MRI *mri_dst) ;

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

static int kalpha = 0 ;
static int debug_slice = -1 ;
static int correct_PD = 0 ;

char *Progname ;
static int align = 1 ;

static void usage_exit(int code) ;

static int niter=10 ;
static float thresh = 25 ;
static int conform = 0 ;

static int InterpMethod = SAMPLE_TRILINEAR;  /*E* prev default behavior */
static int sinchalfwindow = 3;

static char *residual_name = NULL ;

#define MAX_IMAGES 100


static double FLASHforwardModel(double flip_angle, double TR, double PD, 
                                double T1) ;

static double estimate_ms_params(MRI **mri_flash, MRI **mri_flash_synth, int nvolumes, MRI *mri_T1, MRI *mri_PD, MRI *mri_sse, MATRIX **M_reg);
static double estimate_ms_params_with_kalpha(MRI **mri_flash, MRI **mri_flash_synth, int nvolumes, MRI *mri_T1, MRI *mri_PD, MRI *mri_fa, MRI *mri_sse, MATRIX **M_reg);
static void estimate_rigid_regmatrix(MRI *mri_source, MRI *mri_traget, MATRIX *M_reg);

static float        tr = 0, te = 0, fa = 0 ;

static int write_iterations=0;

static int average_volumes_with_different_echo_times(MRI **mri_flash, MRI **mri_all_flash, int nvolumes_total) ;

static MRI *estimate_T2star(MRI **mri_all_flash, int nvolumes, MRI *mri_T1, MRI *mri_PD) ;
static MRI *compute_T2star_map(MRI **mri_flash, int nvolumes, int *scan_types) ;

int
main(int argc, char *argv[])
{
  char   **av, fname[STRLEN] ;
  int    ac, nargs, i ;
  MRI    *mri_flash[MAX_IMAGES], *mri_flash_synth[MAX_IMAGES], *mri_all_flash[MAX_IMAGES],
		     *mri_T1, *mri_PD, *mri_sse, *mri_fa, *mri_T2star ;
  char   *in_fname, *out_dir ;
  int    msec, minutes, seconds, nvolumes, nvolumes_total ;
  struct timeb start ;
  MATRIX *M_reg[MAX_IMAGES];
  double rms ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_ms_fitparms.c,v 1.16 2003/10/14 21:39:00 fischl Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  TimerStart(&start) ;
  
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
        fa = RADIANS(atof(argv[i+1])) ;
      else
        ErrorExit(ERROR_BADPARM, "%s: unsupported MR parameter %s",
                  Progname, argv[i]+1) ;
      i++ ;  /* skip parameter */
      continue ;
    }

    in_fname = argv[i] ;
    printf("reading %s...", in_fname) ;

    mri_flash[nvolumes] = MRIread(in_fname) ;
    if (mri_flash[nvolumes] == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read volume %s",
                Progname, in_fname) ;

    mri_flash[nvolumes]->register_mat = MRIgetVoxelToRasXform(mri_flash[nvolumes]);

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
           mri_flash[nvolumes]->tr, DEGREES(mri_flash[nvolumes]->flip_angle)) ;
#if 0
    mri_flash[nvolumes]->flip_angle = RADIANS(mri_flash[nvolumes]->flip_angle);
#endif
    if (conform)
    {
      MRI *mri_tmp ;

      printf("embedding and interpolating volume\n") ;
      mri_tmp = MRIconform(mri_flash[nvolumes]) ;
      /*      MRIfree(&mri_src) ;*/
      mri_flash[nvolumes] = mri_tmp ;
    }
		mri_all_flash[nvolumes] = mri_flash[nvolumes] ;  /* for multi-echo, mri_flash will be avg across echo time */
    if (FZERO(mri_flash[nvolumes]->tr) || 
        FZERO(mri_flash[nvolumes]->flip_angle))
      ErrorExit(ERROR_BADPARM, "%s: invalid TR or FA for image %d:%s",
                Progname, nvolumes, in_fname) ;
    //E/    if (mri_flash[nvolumes]->type != MRI_SHORT)
    //E/      ErrorExit(ERROR_BADPARM, "%s: input volumes must be 16 bits/voxel", Progname) ;

    nvolumes++ ;
  }

	nvolumes_total = nvolumes ;   /* all volumes read in */
	nvolumes = average_volumes_with_different_echo_times(mri_flash, mri_all_flash, nvolumes_total) ;
	if (nvolumes == 2)
		niter = 1 ;  /* don't bother motion-correcting when we only have 2 volumes */
	for (i = 0 ; i < nvolumes ; i++)
	{
    mri_flash_synth[i] = MRIclone(mri_flash[i], NULL) ;
    mri_flash_synth[i]->register_mat = MRIgetVoxelToRasXform(mri_flash[i]);
	}
	
  {
    int i, j ;
		
    for (i = 0 ; i < nvolumes ; i++)
    {
      for (j = i+1 ; j < nvolumes ; j++)
      {
				if ((mri_flash[i]->width != mri_flash[j]->width) ||
						(mri_flash[i]->height != mri_flash[j]->height) ||
						(mri_flash[i]->depth != mri_flash[j]->depth) ||
						(mri_flash[i]->type != mri_flash[j]->type))
					ErrorExit(ERROR_BADPARM, "%s:\nvolumes %d (type %d) and %d (type %d) don't match (%d x %d x %d) vs (%d x %d x %d)\n",
										Progname, i, mri_flash[i]->type, j, mri_flash[j]->type, mri_flash[i]->width, 
										mri_flash[i]->height, mri_flash[i]->depth, 
										mri_flash[j]->width, mri_flash[j]->height, mri_flash[j]->depth) ;
      }
    }
  }
  
  printf("using %d FLASH volumes to estimate tissue parameters.\n", nvolumes) ;
  mri_T1 = MRIclone(mri_flash[0], NULL) ;
  mri_PD = MRIclone(mri_flash[0], NULL) ;
  mri_sse = MRIclone(mri_flash[0], NULL) ;
  if (kalpha)
  {
    mri_fa = MRIalloc(mri_T1->width, mri_T1->height, mri_T1->depth, MRI_FLOAT) ;
    MRIcopyHeader(mri_T1, mri_fa) ;
  }
  else
    mri_fa = NULL ;
  
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
      if (mri_fa)
				rms = estimate_ms_params_with_kalpha(mri_flash, mri_flash_synth, nvolumes, 
																						 mri_T1, mri_PD, mri_fa, mri_sse, M_reg) ;
      else
				rms = estimate_ms_params(mri_flash, mri_flash_synth, nvolumes, 
																 mri_T1, mri_PD, mri_sse, M_reg) ;
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
        printf("writing T1 estimates to %s...\n", fname) ;
        MRIwrite(mri_T1, fname) ;
        sprintf(fname,"%s/PD-%d.mgh",out_dir,iter);
        printf("writing PD estimates to %s...\n", fname) ;
        MRIwrite(mri_PD, fname) ;
        sprintf(fname,"%s/sse-%d.mgh",out_dir,iter);
        printf("writing residual sse to %s...\n", fname) ;
        MRIwrite(mri_sse, fname) ;
				if (mri_fa)
				{
					sprintf(fname,"%s/fa-%d.mgh",out_dir,iter);
					printf("writing fa map to %s...\n", fname) ;
					MRIwrite(mri_fa, fname) ;
				}
				
        for (j=0;j<nvolumes;j++)
        {
          sprintf(fname,"%s/vol%d-%d.mgh",out_dir,j,iter);
          printf("writing synthetic images to %s...\n", fname);
          MRIwrite(mri_flash_synth[j], fname) ;
          sprintf(fname,"%s/vol%d-%d.lta",out_dir,j,iter); 
          printf("writing regisration matrix to %s...\n", fname);
          lta = LTAalloc(1,NULL) ;
          MatrixCopy(M_reg[j],lta->xforms[0].m_L) ;
					// add src and dst information
          getVolGeom(mri_flash_synth[j], &lta->xforms[0].src);
					getVolGeom(mri_flash[j], &lta->xforms[0].dst);
					lta->type = LINEAR_RAS_TO_RAS ;
          LTAwriteEx(lta,fname) ;
        }
      }
    }
    sprintf(fname,"%s/PD.mgh",out_dir);
    printf("writing PD estimates to %s...\n", fname) ;
    MRIwrite(mri_PD, fname) ;
		mri_T2star = estimate_T2star(mri_all_flash, nvolumes_total, mri_T1, mri_PD) ;
		if (mri_T2star)
		{
			if  (correct_PD)
			{
				sprintf(fname,"%s/PDcorrected.mgh",out_dir);
				printf("writing corrected PD estimates to %s...\n", fname) ;
				MRIwrite(mri_PD, fname) ;
			}
			sprintf(fname,"%s/T2star.mgh",out_dir);
			printf("writing T2star estimates to %s...\n", fname) ;
			MRIwrite(mri_T2star, fname) ;
		}
    sprintf(fname,"%s/T1.mgh",out_dir);
    printf("writing T1 estimates to %s...\n", fname) ;
    MRIwrite(mri_T1, fname) ;
    sprintf(fname,"%s/sse.mgh",out_dir);
    printf("writing residual sse to %s...\n", fname) ;
    MRIwrite(mri_sse, fname) ;
    if (mri_fa)
    {
      sprintf(fname,"%s/fa.mgh",out_dir);
      printf("writing fa map to %s...\n", fname) ;
      MRIwrite(mri_fa, fname) ;
    }
    for (j=0;j<nvolumes;j++)
    {
      sprintf(fname,"%s/vol%d.mgh",out_dir,j);
      printf("writing synthetic images to %s...\n", fname);
      MRIwrite(mri_flash_synth[j], fname) ;
      sprintf(fname,"%s/vol%d.lta",out_dir,j);
      printf("writing registration matrix to %s...\n", fname);
      lta = LTAalloc(1,NULL) ;
      MatrixCopy(M_reg[j],lta->xforms[0].m_L) ;
      // add src and dst information
      getVolGeom(mri_flash_synth[j], &lta->xforms[0].src);
      getVolGeom(mri_flash[j], &lta->xforms[0].dst);
      lta->type = LINEAR_RAS_TO_RAS ;
      LTAwriteEx(lta,fname) ;
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
  if (!stricmp(option, "kalpha"))
  {
    kalpha = 1 ;
    printf("fitting for inhomogeneities in transmit field...\n") ;
  }
  else if (!stricmp(option, "debug_slice"))
  {
    debug_slice = atoi(argv[2]) ;
    nargs = 1 ;
    printf("debugging slice %d...\n", debug_slice) ;
  }
  else if (!stricmp(option, "debug_voxel"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)...\n", Gx, Gy, Gz) ;
  }
  else if (!stricmp(option, "correct"))
  {
		correct_PD = 1 ;
    printf("correcting PD by T2* estimates\n") ;
  }
  else if (!stricmp(option, "conform"))
  {
    conform = 1 ;
    printf("interpolating volume to be isotropic 1mm^3\n") ;
  }
  else if (!stricmp(option, "window"))
  {
    printf("window option not implemented\n");
    /*E* window_flag = 1 ; */
  }

  /*E* Interpolation method.  Default is trilinear, other options are
    nearest, cubic, sinc.  You can say -foo or -interp foo.  For sinc,
    you can say -interp sinc 3 or -interp sinc -hw 3 or -sinc 3 or
    -sinc -hw 3.  Maybe -hw 3 should imply sinc, but right now it
    doesn't.  */

  else if (!stricmp(option, "st") ||
	   !stricmp(option, "sample") ||
	   !stricmp(option, "sample_type") ||
	   !stricmp(option, "interp"))
  {
    InterpMethod = MRIinterpCode(argv[2]) ;
    nargs = 1;
    if (InterpMethod==SAMPLE_SINC)
    {
      if ((argc<4) || !strncmp(argv[3],"-",1)) /*E* i.e. no sinchalfwindow value supplied */
      {
	printf("using sinc interpolation (default windowwidth is 6)\n");
      }
      else
      {
	sinchalfwindow = atoi(argv[3]);
	nargs = 2;
	printf("using sinc interpolation with windowwidth of %d\n", 2*sinchalfwindow);
      }
    }
  }
  else if (!stricmp(option, "sinc"))
  {
    InterpMethod = SAMPLE_SINC;
    if ((argc<3) || !strncmp(argv[2],"-",1)) /*E* i.e. no sinchalfwindow value supplied */
    {
      printf("using sinc interpolation (default windowwidth is 6)\n");
    }
    else
    {
      sinchalfwindow = atoi(argv[2]);
      nargs = 1;
      printf("using sinc interpolation with windowwidth of %d\n", 2*sinchalfwindow);
    }
  }
  else if (!stricmp(option, "sinchalfwindow") ||
	   !stricmp(option, "hw"))
  {
    /*E* InterpMethod = SAMPLE_SINC; //? */
    sinchalfwindow = atoi(argv[2]);
    nargs = 1;
    printf("using sinc interpolation with windowwidth of %d\n", 2*sinchalfwindow);
  }
  else if (!stricmp(option, "trilinear"))
  {
    InterpMethod = SAMPLE_TRILINEAR;
    printf("using trilinear interpolation\n");
  }
  else if (!stricmp(option, "cubic"))
  {
    InterpMethod = SAMPLE_CUBIC;
    printf("using cubic interpolation\n");
  }
  else if (!stricmp(option, "nearest"))
  {
    InterpMethod = SAMPLE_NEAREST;
    printf("using nearest-neighbor interpolation\n");
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
    fa = RADIANS(atof(argv[2])) ;
    nargs = 1;
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
  printf("usage: %s [options] <volume> ... <output directory>\n", Progname) ;
  printf("this program takes an arbitrary # of FLASH images as input, and estimates\n"
	 "the T1 and PD values of the data for voxel, as well as a linear transform\n"
	 "aligning each of the images. The T1 and PD maps are written into <output directory>\n"
	 "together with synthetic volumes names vol?.mgh, one for each of the input\n"
	 "volumes. All the output volumes are generated in the common (motion-corrected) space.\n");
  printf("Note that TR, TE and the flip angle are read directly from the image header.\n"
	 "If this information is not available, it can be specified on the command line using\n"
	 "-tr <TR in msec> -te <TE in msec> -fa <flip angle in degrees> before each volume.\n");
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
#define T1_MAX 7000.0
#define MAX_NVALS 5000
#define MAX_NVOLS MAX_IMAGES
#define FAk_MAX    2.0
#define FAk_MIN    0.5

static double
estimate_ms_params(MRI **mri_flash, MRI **mri_flash_synth, int nvolumes, MRI *mri_T1, MRI *mri_PD, MRI *mri_sse, MATRIX **M_reg)
{
  double   SignalTableValues[MAX_NVALS][MAX_NVOLS], SignalTableT1[MAX_NVALS], SignalTableNorm[MAX_NVALS], ImageValues[MAX_NVOLS], total_sse ;
  double   se, best_se, ss, sse, err, val, norm, T1, PD, xf, yf, zf ;
  int      i, j, x, y, z, indx, min_indx, max_indx, best_indx, center_indx, stepindx;
  int      width=mri_T1->width, height=mri_T1->height, depth=mri_T1->depth, nvalues=MAX_NVALS, nevals;
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
		if (x == 0 && y == 0 && z == 4)
			DiagBreak() ;
    ss = 0;
    for (j = 0 ; j < nvolumes ; j++)
    {
      mri = mri_flash[j] ;
#if 0
      //E/ val = ImageValues[j] = MRISvox(mri, x, y, z) ;
      val = ImageValues[j] = MRIgetVoxVal(mri, x, y, z, 0);
#endif
      voxvec1->data[0]=x; voxvec1->data[1]=y; voxvec1->data[2]=z;
      MatrixMultiply(vox2ras[j],voxvec1,rasvec1);
      MatrixMultiply(M_reg[j],rasvec1,rasvec2);
      MatrixMultiply(ras2vox[j],rasvec2,voxvec2);
      xf=voxvec2->data[0]; yf=voxvec2->data[1]; zf=voxvec2->data[2];
      if(InterpMethod==SAMPLE_SINC)
	MRIsincSampleVolume(mri, xf, yf, zf, sinchalfwindow, &val) ;
      else
	MRIsampleVolumeType(mri, xf, yf, zf, &val, InterpMethod) ;
      ImageValues[j] = val;
      ss += val*val;
    }
		if (ImageValues[0] == 76 && ImageValues[1] == 113)
			DiagBreak() ;
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
      max_indx = MIN(best_indx+step[stepindx]/2,nvalues-1);
      center_indx = best_indx;
    }

    //E/ T1 = MRISvox(mri_T1, x, y, z) = (short)nint(SignalTableT1[best_indx]);
    T1 = SignalTableT1[best_indx];
    MRIsetVoxVal(mri_T1, x, y, z, 0, T1);

    PD = norm/SignalTableNorm[best_indx];
    if ((short)PD < 0)
      PD = (double)(0x7fff-1) ;
    //E/ MRISvox(mri_PD, x, y, z) = PD ;
    MRIsetVoxVal(mri_PD, x, y, z, 0, PD);
    for (j = 0 ; j < nvolumes ; j++)
    {
      mri = mri_flash_synth[j] ;
      //E/ MRISvox(mri, x, y, z) = (short)nint(PD*SignalTableNorm[best_indx]*SignalTableValues[best_indx][j]);
      MRIsetVoxVal(mri, x, y, z, 0, PD*SignalTableNorm[best_indx]*SignalTableValues[best_indx][j]);
    }
    sse = 0;
    for (j = 0 ; j < nvolumes ; j++)
    {
      //E/ err = MRISvox(mri_flash_synth[j], x, y, z)-ImageValues[j]*norm;
      err = MRIgetVoxVal(mri_flash_synth[j], x, y, z, 0)-ImageValues[j]*norm;
      sse += err*err; 
    }
    total_sse += sse ;
    //E/ MRISvox(mri_sse, x, y, z) = (short)nint(sqrt(sse));
    MRIsetVoxVal(mri_sse, x, y, z, 0, sqrt(sse));
		if (T1 >= 4999 && ImageValues[0] > 70 && ImageValues[1] > 70)
			DiagBreak() ;
  }
  total_sse = sqrt(total_sse / (width*height*depth)) ;
  return(total_sse) ;
}

#define UNITY_FAk_INDEX  nint((1 - FAk_MIN) * (nvalues-1) / (FAk_MAX - FAk_MIN))
static double
estimate_ms_params_with_kalpha(MRI **mri_flash, MRI **mri_flash_synth, int nvolumes, MRI *mri_T1, MRI *mri_PD, MRI *mri_fa, MRI *mri_sse, MATRIX **M_reg)
{
  double   ***SignalTableValues, SignalTableT1[MAX_NVALS], **SignalTableNorm, ImageValues[MAX_NVOLS], total_sse, SignalTableFAk[MAX_NVALS], FAk, last_sse, last_FAk ;
  double   se, best_se, ss, sse, err, val, norm, T1, PD, xf, yf, zf, last_T1, last_PD ;
  int      i, j, k, x, y, z, T1_indx, min_T1_indx, max_T1_indx, best_T1_indx, center_T1_indx, 
           FA_indx, min_FA_indx, max_FA_indx, best_FA_indx, center_FA_indx,
		       stepindx;
  int      width=mri_T1->width, height=mri_T1->height, depth=mri_T1->depth, nvalues=MAX_NVALS/3, 
		       nevals, niter;
  int      nstep=11, step[11]={1024,512,256,128,64,32,16,8,4,2,1};
  MRI      *mri ;
  MATRIX   *vox2ras[MAX_NVOLS], *ras2vox[MAX_NVOLS], *voxvec1, *voxvec2, *rasvec1, *rasvec2;

  SignalTableValues = (double ***)calloc(nvalues, sizeof(double **)) ;
  SignalTableNorm = (double **)calloc(nvalues, sizeof(double *)) ;
  for (i = 0 ; i < nvalues ; i++)
  {
    SignalTableValues[i] = (double **)calloc(nvalues, sizeof(double *)) ;
    SignalTableNorm[i] = (double *)calloc(nvalues, sizeof(double)) ;
    for (j = 0 ; j < nvalues ; j++)
    {
      SignalTableValues[i][j] = (double *)calloc(nvolumes, sizeof(double)) ;
    }
  }
  
  voxvec1 = MatrixAlloc(4,1,MATRIX_REAL); voxvec1->rptr[4][1] = 1.0;
  voxvec2 = MatrixCopy(voxvec1, NULL);
  rasvec1 = MatrixCopy(voxvec1, NULL);
  rasvec2 = MatrixCopy(voxvec1, NULL);
  for (k = 0 ; k < nvolumes ; k++)
  {
    vox2ras[k] = MatrixCopy(mri_flash[k]->register_mat, NULL);
    ras2vox[k] = MatrixInverse(vox2ras[k], NULL);
  }

  PD = 1;
  for (i=0; i<nvalues; i++)
  {
    T1 = SignalTableT1[i] = T1_MIN+i*(T1_MAX-T1_MIN)/(nvalues-1);  
    for (j = 0 ; j < nvalues ; j++)
    {
      FAk = SignalTableFAk[j] = FAk_MIN+(double)j*(FAk_MAX-FAk_MIN)/(double)(nvalues-1);  
      ss = 0;
      for (k = 0 ; k < nvolumes ; k++)
      {
	mri = mri_flash[k] ;
	val = SignalTableValues[i][j][k] = FLASHforwardModel(FAk*mri->flip_angle, mri->tr, PD, T1) ;
	ss += val*val;
      }
      norm = SignalTableNorm[i][j] = sqrt(ss);
      if (norm>0) 
	for (k = 0 ; k < nvolumes ; k++) 
	  SignalTableValues[i][j][k] /= norm;
    }
  }

  total_sse = 0 ;
  for (z = 0 ; z < depth ; z++)
  for (y = 0 ; y < height ; y++)
  for (x = 0 ; x < width ; x++)
  {
    if (x == Gx && y == Gy && z == Gz)
      DiagBreak() ;
    if (debug_slice >= 0 && z != debug_slice)
      continue ;
    if (x == 0 && y == 0 && z == 4)
      DiagBreak() ;
    ss = 0;
    for (k = 0 ; k < nvolumes ; k++)
    {
      mri = mri_flash[k] ;
#if 0
      //E/ val = ImageValues[k] = MRISvox(mri, x, y, z) ;
      val = ImageValues[k] = MRIgetVoxVal(mri, x, y, z, 0) ;
#endif
      voxvec1->data[0]=x; voxvec1->data[1]=y; voxvec1->data[2]=z;
      MatrixMultiply(vox2ras[k],voxvec1,rasvec1);
      MatrixMultiply(M_reg[k],rasvec1,rasvec2);
      MatrixMultiply(ras2vox[k],rasvec2,voxvec2);
      xf=voxvec2->data[0]; yf=voxvec2->data[1]; zf=voxvec2->data[2];
      if(InterpMethod==SAMPLE_SINC)
	MRIsincSampleVolume(mri, xf, yf, zf, sinchalfwindow, &val) ;
      else
	MRIsampleVolumeType(mri, xf, yf, zf, &val, InterpMethod) ;
      ImageValues[k] = val;
      ss += val*val;
    }
		if (ImageValues[0] == 76 && ImageValues[1] == 113)
			DiagBreak() ;
    norm = sqrt(ss);
		if (FZERO(norm))
			continue ;   /* no real data */

    if (norm>0) 
			for (k = 0 ; k < nvolumes ; k++) 
				ImageValues[k] /= norm;
   
    min_T1_indx = best_T1_indx = 0; max_T1_indx = nvalues-1; 
    best_T1_indx = -1; center_T1_indx = -1; 
    min_FA_indx = 0; max_FA_indx = nvalues-1; 
    best_FA_indx = -1; center_FA_indx = -1; best_se = 10000000;
		best_FA_indx = UNITY_FAk_INDEX  ;
    nevals = 0; 
#if 0
    best_T1_indx = nint((1000.0 - T1_MIN) * (nvalues-1) / (T1_MAX - T1_MIN)) ;
#endif

    sse = -1 ; niter = 0 ;
    last_T1 = last_FAk = last_PD = 0 ;  /* for compiler warning */
    do
    {
      for (stepindx=0; stepindx<nstep; stepindx++) 
      {
	for (T1_indx=min_T1_indx; T1_indx<=max_T1_indx; T1_indx+=step[stepindx])
	{
	  if (T1_indx!=center_T1_indx)
	  {
	    se = 0;
	    for (k = 0 ; k < nvolumes ; k++)
	    {
	      err = ImageValues[k]-SignalTableValues[T1_indx][best_FA_indx][k];
	      se += err*err; 
	    }
	    if (se<best_se)
	    {
	      best_se = se;
	      best_T1_indx = T1_indx;
	    }
	    nevals++;
	  }
	}
	min_T1_indx = MAX(best_T1_indx-step[stepindx]/2,1);
	max_T1_indx = MIN(best_T1_indx+step[stepindx]/2,nvalues-1);
	center_T1_indx = best_T1_indx;
      }
      for (stepindx=0; stepindx<nstep; stepindx++) 
      {
	for (FA_indx=min_FA_indx; FA_indx<=max_FA_indx; FA_indx+=step[stepindx])
	{
	  if (FA_indx!=center_FA_indx)
	  {
	    se = 0;
	    for (k = 0 ; k < nvolumes ; k++)
	    {
	      err = ImageValues[k]-SignalTableValues[best_T1_indx][FA_indx][k];
	      se += err*err; 
	    }
	    if (se<best_se)
	    {
	      best_se = se;
	      best_FA_indx = FA_indx;
	    }
	    nevals++;
	  }
	}
	min_FA_indx = MAX(best_FA_indx-step[stepindx]/2,1);
	max_FA_indx = MIN(best_FA_indx+step[stepindx]/2,nvalues-1);
	center_FA_indx = best_FA_indx;
      }
      
      //E/ T1 = MRISvox(mri_T1, x, y, z) = (short)nint(SignalTableT1[best_T1_indx]);
      T1 = SignalTableT1[best_T1_indx];
      MRIsetVoxVal(mri_T1, x, y, z, 0, T1);
      
      FAk = MRIFvox(mri_fa, x, y, z) = SignalTableFAk[best_FA_indx];
      PD = norm/SignalTableNorm[best_T1_indx][best_FA_indx];
      if ((short)PD < 0)
	PD = (double)(0x7fff-1) ;
      //E/ MRISvox(mri_PD, x, y, z) = PD ;
      MRIsetVoxVal(mri_PD, x, y, z, 0, PD);
      for (k = 0 ; k < nvolumes ; k++)
      {
	mri = mri_flash_synth[k] ;
	//E/ MRISvox(mri, x, y, z) = (short)nint(PD*SignalTableNorm[best_T1_indx][best_FA_indx] * SignalTableValues[best_T1_indx][best_FA_indx][k]);
	MRIsetVoxVal(mri, x, y, z, 0,
		     PD * SignalTableNorm[best_T1_indx][best_FA_indx]
		     * SignalTableValues[best_T1_indx][best_FA_indx][k]);
      }
      last_sse = sse ;
      sse = 0 ;
      for (k = 0 ; k < nvolumes ; k++)
      {
	//E/ err = MRISvox(mri_flash_synth[k], x, y, z)-ImageValues[k]*norm;
	err = MRIgetVoxVal(mri_flash_synth[k], x, y, z, 0)-ImageValues[k]*norm;
	sse += err*err; 
      }
      if (last_sse < 0)
	last_sse = sse+1 ;  /* first time */
      //E/ MRISvox(mri_sse, x, y, z) = (short)nint(sqrt(sse));
      MRIsetVoxVal(mri_sse, x, y, z, 0, sqrt(sse));
      if (sse > last_sse)  /* revert to old values */
      {
	T1 = last_T1 ; FAk = last_FAk ; PD = last_PD ;
	best_T1_indx = nint((T1 - T1_MIN) * (nvalues-1) / (T1_MAX - T1_MIN)) ;
	best_FA_indx = nint((FAk - FAk_MIN) * (nvalues-1) / (FAk_MAX - FAk_MIN)) ;
	//E/ MRISvox(mri_PD, x, y, z) = PD ;
	MRIsetVoxVal(mri_PD, x, y, z, 0, PD) ;
	for (k = 0 ; k < nvolumes ; k++)
	{
	  mri = mri_flash_synth[k] ;
	  //E/ MRISvox(mri, x, y, z) = (short)nint(PD*SignalTableNorm[best_T1_indx][best_FA_indx] * SignalTableValues[best_T1_indx][best_FA_indx][k]);
	  MRIsetVoxVal(mri, x, y, z, 0,
		       PD * SignalTableNorm[best_T1_indx][best_FA_indx]
		       * SignalTableValues[best_T1_indx][best_FA_indx][k]);
	}
	sse = last_sse ;
      }
      else
      {
	last_T1 = T1 ; last_FAk = FAk ; last_PD = PD ;
      }
      
    } while ((sse < last_sse) && (niter++ < 4));
    total_sse += sse ;
  }
  total_sse = sqrt(total_sse / (width*height*depth)) ;
  for (i = 0 ; i < nvalues ; i++)
  {
    for (j = 0 ; j < nvalues ; j++)
    {
      free(SignalTableValues[i][j]) ;
    }
    free(SignalTableValues[i]) ;
    free(SignalTableNorm[i]) ;
  }
  free(SignalTableValues) ;
  free(SignalTableNorm) ;
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
      //E/ voxval1[indx] = MRISvox(mri_source, x, y, z); 
      voxval1[indx] = MRIgetVoxVal(mri_source, x, y, z, 0); 
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
	  if(InterpMethod==SAMPLE_SINC)
	    MRIsincSampleVolume(mri_target, xf, yf, zf, sinchalfwindow, &val2) ;
	  else
	    MRIsampleVolumeType(mri_target, xf, yf, zf, &val2, InterpMethod) ;
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
double
dFLASH_dk(MRI *mri_T1, MRI *mri_PD, MRI *mri_fa, double TR, double flip_angle, int x, int y, int z)
{
  double T1, PD, dk, k, e_TR_T1_minus_1, cos_ka, numer, denom, e_TR_T1 ;
  
  //E/ T1 = MRISvox(mri_T1, x, y, z) ;
  T1 = MRIgetVoxVal(mri_T1, x, y, z, 0) ;
  //E/ PD = MRISvox(mri_PD, x, y, z) ;
  PD = MRIgetVoxVal(mri_PD, x, y, z, 0) ;
  k = MRIFvox(mri_fa, x, y, z) ;
  
  e_TR_T1 = exp(TR/T1) ;
  e_TR_T1_minus_1 = e_TR_T1 - 1 ;
  cos_ka = cos(k*flip_angle) ;
  
  numer = PD * flip_angle * e_TR_T1_minus_1 * (e_TR_T1*cos_ka-1) ;
  denom = (e_TR_T1-cos_ka) ; denom*= denom ;
  if (FZERO(denom))
    denom = 1 ;
  dk = numer /denom ;
  return(dk) ;
}
#define PARAMETERS_MATCH(mri1, mri2)  ((mri1->tr == mri2->tr) && (mri1->flip_angle == mri2->flip_angle))

static int
average_volumes_with_different_echo_times(MRI **mri_flash, MRI **mri_all_flash, int nvolumes_total)
{
	int i, j, nvolumes, averaged[MAX_IMAGES], navgs ;
	MRI    *mri_avg ;

	memset(averaged, 0, sizeof(averaged)) ;

	for (nvolumes = i = 0 ; i < nvolumes_total ; i++)
	{
		if (averaged[i])
			continue ;
		averaged[i] = 1 ;
		mri_avg = MRIcopy(mri_all_flash[i], NULL) ;
    mri_avg->register_mat = MRIgetVoxelToRasXform(mri_all_flash[nvolumes]);
		navgs = 1 ;
		for (j = i+1 ; j < nvolumes_total ; j++)
		{
			if (averaged[j])
				continue ;
			if (PARAMETERS_MATCH(mri_all_flash[j], mri_avg) == 0)
				continue ;
			MRIaverage(mri_all_flash[j], navgs, mri_avg) ;
			averaged[j] = 1 ;
			navgs++ ;
		}
		mri_flash[nvolumes] = mri_avg ;
		nvolumes++ ;
	}


	return(nvolumes) ;
}

static MRI *
estimate_T2star(MRI **mri_flash, int nvolumes, MRI *mri_T1, MRI *mri_PD)
{
	int    i, j, nechoes, processed[MAX_IMAGES], nprocessed, x, y, z, different_te, width, depth, height ;
	MRI    *mri_T2star = NULL ;
	double decay, T2star ;
	Real   PD ;

	/* first decide whether T2* can be estimated at all */
	different_te = 0 ;
	for (i = 0 ; different_te == 0 && i < nvolumes ; i++)
	{
		for (j = i+1 ; different_te == 0 && j < nvolumes ; j++)
		{
			if (mri_flash[i]->te != mri_flash[j]->te)
				different_te = 1 ;
		}
	}
	if (different_te == 0)
		return(NULL) ;  /* can't estimate T2* */

	memset(processed, 0, sizeof(processed)) ;

	for (nprocessed = nechoes = i = 0 ; i < nvolumes ; i++)
	{
		if (processed[i])
			continue ;
		processed[i] = nprocessed+1 ;
		nechoes = 1 ;

		for (j = i+1 ; j < nvolumes ; j++)
		{
			if (processed[j])
				continue ;
			if (PARAMETERS_MATCH(mri_flash[i], mri_flash[j]) == 0)
				continue ;
			processed[j] = nprocessed+1 ;
		}
		nprocessed++ ;
	}
	printf("estimating T2* with %d different acquisitions, each with %d echoes...\n",
				 nprocessed, nvolumes/nprocessed) ;
	mri_T2star = compute_T2star_map(mri_flash, nvolumes, processed) ;

	/* now update PD map to take out T2* component */
	if (correct_PD)
	{
		width = mri_T2star->width ; height = mri_T2star->height ; depth = mri_T2star->depth ;
		for (x = 0 ; x < width ; x++)
		{
			for (y = 0 ; y < height ; y++)
			{
				for (z = 0 ; z < depth ; z++)
				{
					PD = MRIgetVoxVal(mri_PD, x, y, z, 0);
					if (PD < 100)
						MRIsetVoxVal(mri_T2star, x, y, z, 0, 0) ;
					if (x == Gx && y == Gy && z == Gz)
						DiagBreak() ;
					T2star = MRIgetVoxVal(mri_T2star, x, y, z, 0) ;
					if (FZERO(T2star))
						continue ;
					for (i = 0, decay = 0.0 ; i < nvolumes ; i++)
						decay += exp(-mri_flash[i]->te / T2star) ;
					decay /= (float)nvolumes ;
					PD /= decay ;
					MRIsetVoxVal(mri_PD, x, y, z, 0, PD) ;
				}
			}
		}
	}

	return(mri_T2star) ;
}

static MRI *
compute_T2star_map(MRI **mri_flash, int nvolumes, int *scan_types)
{
	MATRIX *mX, *mXpinv = NULL ;
	VECTOR *vY, *vParms = NULL ;
	int    x, y, z, e, width, height, depth, nscans, i ;
	MRI    *mri_T2star ;
	float  T2star, cond ;
	Real   val ;

	for (i = nscans = 0 ; i < nvolumes ; i++)
	{
		if (scan_types[i] > nscans)
			nscans = scan_types[i] ;
	}

	width = mri_flash[0]->width ;
	height = mri_flash[0]->height ;
	depth = mri_flash[0]->depth ;
	mri_T2star = MRIalloc(width, height, depth, MRI_FLOAT) ;
	if (!mri_T2star)
		ErrorExit(ERROR_NOMEMORY, "%s: could not allocate T2* map", Progname) ;
	MRIcopyHeader(mri_flash[0], mri_T2star) ;

	mX = MatrixAlloc(nvolumes, nscans+1, MATRIX_REAL) ;
	vY = VectorAlloc(nvolumes, MATRIX_REAL) ;
	vParms = VectorAlloc(nscans+1, MATRIX_REAL) ;
	for (e = 0 ; e < nvolumes ; e++)
	{
		*MATRIX_RELT(mX, e+1, 1) = -mri_flash[e]->te ;
		for (i = 1 ; i <= nscans ; i++)  /* which multi-echo set does this volume belong to */
		{
			if (scan_types[e] == i)
				*MATRIX_RELT(mX, e+1, i+1) = 1 ;
			else
				*MATRIX_RELT(mX, e+1, i+1) = 0 ;
		}
	}
	mXpinv = MatrixPseudoInverse(mX, mXpinv) ;
	if (!mXpinv)
		ErrorReturn(NULL, (ERROR_BADPARM, "%s: could not invert matrix for T2* estimation", Progname)) ;

	cond = MatrixConditionNumber(mX) ;
	for (x = 0 ; x < width ; x++)
	{
		for (y = 0 ; y < height ; y++)
		{
			for (z = 0 ; z < depth ; z++)
			{
				if (x == Gx && y == Gy && z == Gz)
					DiagBreak() ;
				for (e = 0 ; e < nvolumes ; e++)
				{
					MRIsampleVolumeType(mri_flash[e], x, y, z, &val, SAMPLE_NEAREST) ;
					VECTOR_ELT(vY, e+1) = log(val) ;
				}
				vParms = MatrixMultiply(mXpinv, vY, vParms) ;
				if (!FZERO(*MATRIX_RELT(vParms, 1, 1)))
					T2star = 1 / *MATRIX_RELT(vParms, 1, 1) ;
				else
					T2star = 0 ;
				if (T2star > 10000 || T2star < -1000)
					DiagBreak() ;
				MRIsetVoxVal(mri_T2star, x, y, z, 0, T2star) ;
			}
		}
	}

	MatrixFree(&mX) ; VectorFree(&vY) ; MatrixFree(&mXpinv) ;
	VectorFree(&vParms) ;
	return(mri_T2star) ;
}
