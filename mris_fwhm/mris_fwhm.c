/**
 * @file  mris_fwhm.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: zkaufman $
 *    $Date: 2016/10/14 20:40:04 $
 *    $Revision: 1.40.2.2 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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


/*
BEGINHELP

This program has two functions:
  1. Apply surface-based smoothing so surface overlay data.
     This function overlaps with that in mri_surf2surf
  2. Estimate the smoothness of a surface-based data set.

--i input

Input data. Format must be something readable by mri_convert
(eg, mgh, mgz, img, nii). Alternately, one can synthesize
white gaussian noise with --synth and --synth-frames.

--subject subject (--s)

Subject whose surface the input is defined on. Can use --s instead of
--subject.

--surf surfname

Compute AR1 on surface surfname. Default is white.

--mask maskfile

Compute AR1 only over voxels in the given mask. Format can be anything
accepted by mri_convert. See also --label.

--mask-inv

Invert mask, ie, compute AR1 only over voxels outside the given mask.

--label label

Use label as a mask. Can be inverted with --mask-inv.

--cortex

Use hemi.cortex.label as a mask. Can be inverted with --mask-inv.

--hemi hemi (--h)

Hemifield that the input is defined on. Legal values are lh and rh.
Can use --h instead of --hemi.

--X x.mat

Detrend data with the matrix in x.mat. Ie, y = (I-inv(X'*X)*X')*y, where
y is the input. x.mat must be a matlab4 matrix.

--detrend order

Detrend data with polynomial regressors upto order. If no output is specified,
then order=0 by default. If an output is specified, then no detrending is done.

--sum sumfile

Prints ascii summary to sumfile.

--fwhm fwhm

Smooth input by fwhm mm.

--niters-only <nitersfile>

Only report the number of iterations needed to achieve the FWHM given
by fwhm. If nitersfile is specified, the number of iterations is 
written to the file.

--o outfile

Save (possibly synthesized and/or smoothed) data to outfile. Automatically
detects format. Format must be one accepted as by mri_convert. Note: do
not use analyze or nifit as these cannot store more than 32k in a dimension.
mri_surf2surf can store surface data in those formats. When an output
is specified, detrending is turned off. Normally the mean is removed.

--synth

Synthesize input with white gaussian noise. Ten frames are used by default,
but this can be changed with --synth-frames.

--synth-frames nframes

Synthesize input with white gaussian noise with the given number of frames.
Implies --synth.


ENDHELP
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
double round(double x);
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "macros.h"
#include "utils.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "fmriutils.h"
#include "mrisurf.h"
#include "mri2.h"
#include "fio.h"
#include "version.h"
#include "annotation.h"
#include "cmdargs.h"
#include "timer.h"
#include "matfile.h"
#include "randomfields.h"
#include "icosahedron.h"
#include "pdf.h"
#include "matfile.h"
#include "fsenv.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mris_fwhm.c,v 1.40.2.2 2016/10/14 20:40:04 zkaufman Exp $";
char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

char *subject=NULL, *hemi=NULL, *SUBJECTS_DIR=NULL;
char *surfname="white";
char *surfpath=NULL;
char *inpath=NULL;
char *outpath=NULL;
char *sumfile=NULL;
char *datfile=NULL;
char *ar1datfile=NULL;
char tmpstr[2000];
MRI *InVals=NULL, *mritmp;

char *maskpath=NULL;
char *labelpath=NULL;
MRI *mask=NULL;
LABEL *label=NULL;
int maskinv = 0;

MRI *mritmp=NULL;

MRIS *surf;
double infwhm = 0, ingstd = 0;
int synth = 0, nframes = 10;
int SynthSeed = -1;
int nitersonly=0;

char *Xfile=NULL;
MATRIX *X=NULL;
int DetrendOrder = -1;
int DoDetrend = 1;
int SmoothOnly = 0;
int DoSqr = 0; // take square of input before smoothing

char *ar1fname = NULL;
char *arNfname = NULL;

int FixGroupAreaTest(MRIS *surf, char *outfile);
char *GroupAreaTestFile = NULL;
char *nitersfile = NULL;
double DHiters2fwhm(MRIS *surf, int vtxno, int niters, char *outfile);
int DHvtxno=0, DHniters=0;
char *DHfile = NULL;
int UseCortexLabel = 0;
int   prunemask = 0;
float prune_thr = FLT_MIN;
char *outmaskpath=NULL;
int arNHops = 0;
int nthreads = 1;

typedef struct {
  int cvtx; // center vertex
  int nhops;
  char *hit; // vector to indicate whether vertex has been hit
  int *nperhop; // number of vertices per hop
  int **vtxlist; // list of vertices for each hop
  int *nperhop_alloced; // number of vertices alloced per hop
} SURFHOPLIST;
SURFHOPLIST *SetSurfHopListAlloc(MRI_SURFACE *Surf, int nHops);
SURFHOPLIST *SetSurfHopList(int CenterVtx, MRI_SURFACE *Surf, int nHops);
int SurfHopListFree(SURFHOPLIST **shl0);
MRI *MRISarN(MRIS *surf, MRI *src, MRI *mask, MRI *arN, int N);
MRI *MRISsmoothKernel(MRIS *surf, MRI *src, MRI *mask, MRI *mrikern, MATRIX *globkern, int SqrFlag, MRI *out);

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs, niters=0, Ntp, n, err;
  double fwhm = 0, ar1mn, ar1std, ar1max, avgvtxarea,ftmp, fwhmDH;
  double InterVertexDistAvg, InterVertexDistStdDev;
  MRI *ar1=NULL;
  FILE *fp;
  MATRIX *globkern;
  double gstd, gk0;

  if(0){
#ifdef _OPENMP
    omp_set_num_threads(7);
#endif
    InVals = MRIread(argv[1]);
    surf = MRISread(argv[2]);
    if(0){
      sscanf(argv[3],"%lf",&fwhm);
      gstd = fwhm/sqrt(log(256.0));
      printf("fwhm = %g, gstd = %g\n",fwhm,gstd);
      globkern = GaussianVector(15, 0, gstd, 1, NULL);
      gk0 = globkern->rptr[1][1];
      for(n = 0; n < 15; n++)  globkern->rptr[n+1][1] = globkern->rptr[n+1][1]/gk0;
      sprintf(tmpstr,"globkern.fwhm%02d.mtx",(int)fwhm);
      MatrixWriteTxt(tmpstr, globkern);
    }
    globkern = MatrixReadTxt(argv[3], NULL);
    mritmp = MRISsmoothKernel(surf, InVals, NULL, NULL, globkern, 1, NULL);
    MRIwrite(mritmp,argv[4]);
    exit(0);
  }

  nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd,2000);

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  if (argc == 0) usage_exit();
  parse_commandline(argc, argv);
  check_options();
  if (checkoptsonly) return(0);

  if (SynthSeed < 0) SynthSeed = PDFtodSeed();

  sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,subject,hemi,surfname);
  surfpath = strcpyalloc(tmpstr);

  if (debug) dump_options(stdout);

  surf = MRISread(surfpath);
  if (surf == NULL) {
    printf("ERROR: could not read %s\n",surfpath);
    exit(1);
  }
  MRIScomputeMetricProperties(surf);
  InterVertexDistAvg    = surf->avg_vertex_dist;
  InterVertexDistStdDev = surf->std_vertex_dist;
  avgvtxarea = surf->avg_vertex_area;

  printf("%s %s %s\n",subject,hemi,surfname);
  printf("Number of vertices %d\n",surf->nvertices);
  printf("Number of faces    %d\n",surf->nfaces);
  printf("Total area         %lf\n",surf->total_area);
  if (surf->group_avg_surface_area > 0) 
    printf("GroupSurface %f\n",surf->group_avg_surface_area);
  else                                 
    printf("GroupSurface 0\n");
  //if(getenv("FIX_VERTEX_AREA") != NULL) printf("FIX_VERTEX_AREA 1\n");
  //else                                  printf("FIX_VERTEX_AREA 0\n");
  printf("AvgVtxArea       %lf\n",avgvtxarea);
  printf("AvgVtxDist       %lf\n",InterVertexDistAvg);
  printf("StdVtxDist       %lf\n",InterVertexDistStdDev);

  if(GroupAreaTestFile != NULL){
    FixGroupAreaTest(surf,GroupAreaTestFile);
    exit(0);
  }

  if(DHniters > 0){
    DHiters2fwhm(surf, DHvtxno, DHniters, DHfile);
    exit(0);
  }

  if(nitersonly && infwhm > 0) {
    niters = MRISfwhm2niters(infwhm,surf);
    printf("niters %d\n",niters);
    if(nitersfile){
      fp = fopen(nitersfile,"w");
      fprintf(fp,"%d\n",niters);
      fclose(fp);
    }
    exit(0);
  }

  if (!synth) {
    InVals = MRIread(inpath);
    if (InVals == NULL) exit(1);
    if(InVals->type != MRI_FLOAT){
      printf("Changing input type to float\n");
      mritmp = MRISeqchangeType(InVals, MRI_FLOAT, 0, 0, 0);
      MRIfree(&InVals);
      InVals = mritmp;
    }
  } else {
    printf("Synthesizing %d frames, Seed = %d\n",nframes,SynthSeed);
    InVals = MRIrandn(surf->nvertices, 1, 1, nframes,0, 1, NULL);
  }

  if(DoSqr){
    printf("Computing square of input\n");
    MRIsquare(InVals,NULL,InVals);
  }

  if(labelpath) {
    label = LabelRead(subject, labelpath);
    if (label == NULL) {
      printf("ERROR reading %s\n",labelpath);
      exit(1);
    }
    mask = MRISlabel2Mask(surf, label, NULL);
    mritmp = mri_reshape(mask, InVals->width,InVals->height, InVals->depth, 1);
    MRIfree(&mask);
    mask = mritmp;
  }
  if (maskpath) {
    printf("Loading mask %s\n",maskpath);
    mask = MRIread(maskpath);
    if (mask==NULL) exit(1);
  }
  if (mask) {
    if (maskinv) MRImaskInvert(mask,mask);
    printf("Found %d voxels in mask\n",MRInMask(mask));
  }

  if(prunemask){
    printf("Pruning voxels by thr: %e\n", prune_thr);
    mask = MRIframeBinarize(InVals,FLT_MIN,mask);
  }
  if(outmaskpath){
    printf("Saving final mask to %s\n",outmaskpath);
    err = MRIwrite(mask,outmaskpath);
    if(err) exit(1);
  }


  if(DoDetrend) {
    Ntp = InVals->nframes;
    printf("Polynomial detrending, order = %d\n",DetrendOrder);
    X = MatrixAlloc(Ntp,DetrendOrder+1,MATRIX_REAL);
    for (n=0;n<Ntp;n++) X->rptr[n+1][1] = 1.0;
    ftmp = Ntp/2.0;
    if (DetrendOrder >= 1)
      for (n=0;n<Ntp;n++) X->rptr[n+1][2] = (n-ftmp)/ftmp;
    if (DetrendOrder >= 2)
      for (n=0;n<Ntp;n++) X->rptr[n+1][3] = pow((n-ftmp),2.0)/(ftmp*ftmp);
  } 
  else printf("Not Polynomial detrending\n");

  if(X) {
    if (X->rows != InVals->nframes) {
      printf("ERROR: dimension mismatch between X and input\n");
      exit(1);
    }
    mritmp = fMRIdetrend(InVals,X);
    if (mritmp == NULL) exit(1);
    MRIfree(&InVals);
    InVals = mritmp;
  }

  if (infwhm > 0) {
    niters = MRISfwhm2niters(infwhm,surf);
    printf("Smoothing input by fwhm=%lf, gstd=%lf, niters=%d \n",
           infwhm,ingstd,niters);
    InVals = MRISsmoothMRI(surf, InVals, niters, mask,InVals);
    if(InVals == NULL) exit(1);
    if(SmoothOnly) {
      printf("Only smoothing, so saving and exiting now\n");
      err = MRIwrite(InVals,outpath);
      exit(err);
    }
  }

  printf("Computing spatial AR1 \n");
  ar1 = MRISar1(surf, InVals, mask, NULL);
  if(ar1fname)  MRIwrite(ar1,ar1fname);
  if(arNfname)  {
    MRI *arN;
    printf("Computing ARN over %d hops\n",arNHops);
    arN = MRISarN(surf, InVals, mask, NULL, arNHops);
    MRIwrite(arN,arNfname);
    printf("done Computing ARN\n");
  }

  // Average AR1 over all vertices
  RFglobalStats(ar1, mask, &ar1mn, &ar1std, &ar1max);
  printf("ar1mn = %g, ar1std = %g, ar1max = %g\n",ar1mn, ar1std, ar1max);
  printf("avg vertex dist %g\n",surf->avg_vertex_dist);
  fwhm = MRISfwhmFromAR1(surf, ar1mn);
  printf("avg vertex dist %g\n",surf->avg_vertex_dist);

  printf("fwhm = %lf\n",fwhm);
  fflush(stdout);

  if (sumfile) {
    fp = fopen(sumfile,"w");
    if (fp == NULL) {
      printf("ERROR: opening %s\n",sumfile);
      exit(1);
    }
    dump_options(fp);
    if (infwhm>0) fprintf(fp,"inniters %d\n",niters);
    fprintf(fp,"Number of vertices %d\n",surf->nvertices);
    fprintf(fp,"Number of faces    %d\n",surf->nfaces);
    fprintf(fp,"Total area         %lf\n",surf->total_area);
    if (surf->group_avg_surface_area > 0) fprintf(fp,"GroupSurface %f\n",surf->group_avg_surface_area);
    else                                 fprintf(fp,"GroupSurface 0\n");
    //if (getenv("FIX_VERTEX_AREA") != NULL) fprintf(fp,"FIX_VERTEX_AREA 1\n");
    //else                                  fprintf(fp,"FIX_VERTEX_AREA 0\n");
    fprintf(fp,"AvgVtxArea       %lf\n",avgvtxarea);
    fprintf(fp,"AvgVtxDist       %lf\n",InterVertexDistAvg);
    fprintf(fp,"StdVtxDist       %lf\n",InterVertexDistStdDev);
    fprintf(fp,"ar1mn %lf\n",ar1mn);
    fprintf(fp,"ar1std %lf\n",ar1std);
    fprintf(fp,"fwhm %lf\n",fwhm);
    fclose(fp);
  }

  if(datfile) {
    fp = fopen(datfile,"w");
    if (fp == NULL) {
      printf("ERROR: opening %s\n",datfile);
      exit(1);
    }
    fprintf(fp,"%lf ",fwhm);
    if(infwhm>0) {
      fprintf(fp,"%lf %d ",infwhm,niters);
      fwhmDH = DHiters2fwhm(surf, 10000, niters, NULL);
      fprintf(fp,"%lf ",fwhmDH);
      fwhmDH = DHiters2fwhm(surf, 20000, niters, NULL);
      fprintf(fp,"%lf ",fwhmDH);
    }
    fprintf(fp,"\n");
    fclose(fp);
  }

  if(ar1datfile) {
    fp = fopen(ar1datfile,"w");
    if (fp == NULL) {
      printf("ERROR: opening %s\n",ar1datfile);
      exit(1);
    }
    fprintf(fp,"%lf %lf\n",ar1mn,ar1std);
    fclose(fp);
  }

  if(outpath) {
    err = MRIwrite(InVals,outpath);
    if(err) exit(1);
  }

  exit(0);
}
/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;

  if (argc < 1) usage_exit();

  nargc   = argc;
  pargv = argv;
  while (nargc > 0) {

    option = pargv[0];
    if (debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))  print_help() ;
    else if (!strcasecmp(option, "--version")) print_version() ;
    else if (!strcasecmp(option, "--debug"))   debug = 1;
    else if (!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if (!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;
    else if (!strcasecmp(option, "--synth")) synth = 1;
    else if (!strcasecmp(option, "--nosynth")) synth = 0;
    else if (!strcasecmp(option, "--no-detrend")) DoDetrend = 0;
    else if (!strcasecmp(option, "--prune"))    prunemask = 1;
    else if (!strcasecmp(option, "--no-prune")) prunemask = 0;
    else if (!strcasecmp(option, "--prune_thr")){
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%f",&prune_thr); 
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--sqr")) DoSqr = 1;
    else if (!strcasecmp(option, "--fast")) setenv("USE_FAST_SURF_SMOOTHER","1",1);
    else if (!strcasecmp(option, "--no-fast")) setenv("USE_FAST_SURF_SMOOTHER","0",1);
    else if (!strcasecmp(option, "--smooth-only") || !strcasecmp(option, "--so")) {
      DoDetrend = 0;
      SmoothOnly = 1;
    }
    else if (!strcasecmp(option, "--mask-inv")) maskinv = 1;
    else if (!strcasecmp(option, "--niters-only")){
      nitersonly = 1; synth = 1;
      if(nargc > 0 && !CMDisFlag(pargv[0])){
        nitersfile = pargv[0];
        nargsused++;
      }
    }
    else if (!strcasecmp(option, "--s") || !strcasecmp(option, "--subject")) {
      if (nargc < 1) CMDargNErr(option,1);
      subject = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--sd")) {
      if(nargc < 1) CMDargNErr(option,1);
      FSENVsetSUBJECTS_DIR(pargv[0]);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--dh")){
      if(nargc < 3) CMDargNErr(option,3);
      sscanf(pargv[0],"%d",&DHvtxno);
      sscanf(pargv[1],"%d",&DHniters);
      DHfile = pargv[2];
      nargsused = 3;
    } 
    else if (!strcasecmp(option, "--h") || !strcasecmp(option, "--hemi")) {
      if (nargc < 1) CMDargNErr(option,1);
      hemi = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--lh")) hemi = "lh";
    else if (!strcasecmp(option, "--rh")) hemi = "rh";

    else if (!strcasecmp(option, "--surf")) {
      if (nargc < 1) CMDargNErr(option,1);
      surfname = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--i")) {
      if (nargc < 1) CMDargNErr(option,1);
      inpath = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--mask")) {
      if (nargc < 1) CMDargNErr(option,1);
      maskpath = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--out-mask")) {
      if (nargc < 1) CMDargNErr(option,1);
      outmaskpath = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--label")) {
      if (nargc < 1) CMDargNErr(option,1);
      if(fio_FileExistsReadable(pargv[0]))
	labelpath = fio_fullpath(pargv[0]); // defeat LabelRead()
      else labelpath = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--cortex")) {
      UseCortexLabel = 1;
    } 
    else if (!strcasecmp(option, "--sum")) {
      if (nargc < 1) CMDargNErr(option,1);
      sumfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--dat")) {
      if (nargc < 1) CMDargNErr(option,1);
      datfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--ar1dat")) {
      if (nargc < 1) CMDargNErr(option,1);
      ar1datfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--ar1")) {
      if (nargc < 1) CMDargNErr(option,1);
      ar1fname = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--arN")) {
      if (nargc < 2) CMDargNErr(option,2);
      sscanf(pargv[0],"%d",&arNHops);
      arNfname = pargv[1];
      nargsused = 2;
    }
    else if (!strcasecmp(option, "--fwhm")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&infwhm);
      ingstd = infwhm/sqrt(log(256.0));
      nargsused = 1;
    } else if (!strcasecmp(option, "--synth-frames")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nframes);
      synth = 1;
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--threads") || !strcasecmp(option, "--nthreads") ){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nthreads);
      #ifdef _OPENMP
      omp_set_num_threads(nthreads);
      #endif
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--o")) {
      if (nargc < 1) CMDargNErr(option,1);
      outpath = pargv[0];
      DoDetrend = 0;
      nargsused = 1;
    } else if (!strcasecmp(option, "--group-area-test")) {
      if (nargc < 1) CMDargNErr(option,1);
      GroupAreaTestFile = pargv[0];
      synth = 1;
      nargsused = 1;
    } else if (!strcasecmp(option, "--X")) {
      if (nargc < 1) CMDargNErr(option,1);
      Xfile = pargv[0];
      //X = MatrixReadTxt(Xfile, NULL);
      X = MatlabRead(Xfile);
      DoDetrend = 0;
      nargsused = 1;
    } else if (!strcasecmp(option, "--detrend")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&DetrendOrder);
      if (DetrendOrder > 2) {
        printf("ERROR: cannot have detrending order > 2\n");
        exit(1);
      }
      if(DetrendOrder < 0) DoDetrend = 0;
      nargsused = 1;
    } else {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if (CMDsingleDash(option))
        fprintf(stderr,"       Did you really mean -%s ?\n",option);
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* ------------------------------------------------------ */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
static void print_usage(void) {
  printf("USAGE: %s\n",Progname) ;
  printf("Smooths surface data and/or estimates FWHM\n");
  printf("   --i input\n");
  printf("   --subject subject (--s)\n");
  printf("   --hemi hemi (--h), or --lh or --rh\n");
  printf("   --surf surf <white>\n");
  printf("   --label labelfile\n");
  printf("   --cortex : used hemi.cortex.label\n");
  printf("   --mask maskfile\n");
  printf("   --X x.mat : matlab4 detrending matrix\n");
  printf("   --detrend order : polynomial detrending (default 0, turned off with output)\n");
  printf("   --smooth-only : only smooth (implies --no-detrend)\n");
  printf("   --no-detrend : turn of poly detrending \n");
  printf("   --sqr : compute square of input before smoothing\n");
  printf("   --sum sumfile\n");
  printf("   --dat datfile (only contains fwhm)\n");
  printf("   --ar1dat ar1datfile (contains ar1mean ar1std)\n");
  printf("   --ar1 ar1vol : save spatial ar1 as an overlay\n");
  printf("   --prune - remove any voxel that is zero in any subject (after any inversion)\n");
  printf("   --no-prune - do not prune (default)\n");
  printf("   --out-mask outmask : save final mask\n");
  printf("   \n");
  printf("   --fwhm fwhm : apply before measuring\n");
  printf("   --niters-only <niters> : only report on niters for fwhm\n");
  printf("   --o output\n");
  printf("\n");
  printf("   --sd SUBJECTS_DIR \n");
  printf("   --synth \n");
  printf("   --synth-frames nframes : default is 10 \n");
  printf("   --threads nthreads\n");
  printf("\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  printf("%s\n", vcid) ;
  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
printf("\n");
printf("This program has two functions:\n");
printf("  1. Apply surface-based smoothing so surface overlay data.\n");
printf("     This function overlaps with that in mri_surf2surf\n");
printf("  2. Estimate the smoothness of a surface-based data set.\n");
printf("\n");
printf("--i input\n");
printf("\n");
printf("Input data. Format must be something readable by mri_convert\n");
printf("(eg, mgh, mgz, img, nii). Alternately, one can synthesize\n");
printf("white gaussian noise with --synth and --synth-frames.\n");
printf("\n");
printf("--subject subject (--s)\n");
printf("\n");
printf("Subject whose surface the input is defined on. Can use --s instead of\n");
printf("--subject.\n");
printf("\n");
printf("--surf surfname\n");
printf("\n");
printf("Compute AR1 on surface surfname. Default is white.\n");
printf("\n");
printf("--mask maskfile\n");
printf("\n");
printf("Compute AR1 only over voxels in the given mask. Format can be anything\n");
printf("accepted by mri_convert. See also --label.\n");
printf("\n");
printf("--mask-inv\n");
printf("\n");
printf("Invert mask, ie, compute AR1 only over voxels outside the given mask.\n");
printf("\n");
printf("--label label\n");
printf("\n");
printf("Use label as a mask. Can be inverted with --mask-inv.\n");
printf("\n");
printf("--cortex\n");
printf("\n");
printf("Use hemi.cortex.label as a mask. Can be inverted with --mask-inv.\n");
printf("\n");
printf("--hemi hemi (--h)\n");
printf("\n");
printf("Hemifield that the input is defined on. Legal values are lh and rh.\n");
printf("Can use --h instead of --hemi.\n");
printf("\n");
printf("--X x.mat\n");
printf("\n");
printf("Detrend data with the matrix in x.mat. Ie, y = (I-inv(X'*X)*X')*y, where\n");
printf("y is the input. x.mat must be a matlab4 matrix.\n");
printf("\n");
printf("--detrend order\n");
printf("\n");
printf("Detrend data with polynomial regressors upto order. If no output is specified,\n");
printf("then order=0 by default. If an output is specified, then no detrending is done.\n");
printf("\n");
printf("--sum sumfile\n");
printf("\n");
printf("Prints ascii summary to sumfile.\n");
printf("\n");
printf("--fwhm fwhm\n");
printf("\n");
printf("Smooth input by fwhm mm.\n");
printf("\n");
printf("--niters-only <nitersfile>\n");
printf("\n");
printf("Only report the number of iterations needed to achieve the FWHM given\n");
printf("by fwhm. If nitersfile is specified, the number of iterations is \n");
printf("written to the file.\n");
printf("\n");
printf("--o outfile\n");
printf("\n");
printf("Save (possibly synthesized and/or smoothed) data to outfile. Automatically\n");
printf("detects format. Format must be one accepted as by mri_convert. Note: do\n");
printf("not use analyze or nifit as these cannot store more than 32k in a dimension.\n");
printf("mri_surf2surf can store surface data in those formats. When an output\n");
printf("is specified, detrending is turned off. Normally the mean is removed.\n");
printf("\n");
printf("--synth\n");
printf("\n");
printf("Synthesize input with white gaussian noise. Ten frames are used by default,\n");
printf("but this can be changed with --synth-frames.\n");
printf("\n");
printf("--synth-frames nframes\n");
printf("\n");
printf("Synthesize input with white gaussian noise with the given number of frames.\n");
printf("Implies --synth.\n");
printf("\n");
printf("\n");
  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void) {
  printf("%s\n", vcid) ;
  exit(1) ;
}
/* --------------------------------------------- */
static void check_options(void) {
  if (subject == NULL) {
    printf("ERROR: need to specify --subject\n");
    exit(1);
  }
  if (hemi == NULL) {
    printf("ERROR: need to specify --hemi\n");
    exit(1);
  }
  if (inpath == NULL && !synth) {
    printf("ERROR: need to specify --in or --synth\n");
    exit(1);
  }
  if (maskpath && labelpath) {
    printf("ERROR: cannot specify both --label and --mask\n");
    exit(1);
  }
  if (infwhm == 0 && nitersonly) {
    printf("ERROR: must specify --fwhm with --niters-only\n");
    exit(1);
  }
  if(X != NULL && DetrendOrder > 0){
    printf("ERROR: cannot --X and --detrend\n");
    exit(1);
  }
  if(X == NULL && DetrendOrder < 0 && DoDetrend) DetrendOrder = 0;
  if(SmoothOnly && outpath == 0){
    printf("ERROR: must spec output with --smooth-only\n");
    exit(1);
  }
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if(SUBJECTS_DIR == NULL) {
    printf("ERROR: SUBJECTS_DIR not defined in environment\n");
    exit(1);
  }
  if(UseCortexLabel){
    if(labelpath != NULL){
      printf("ERROR: cannot spec --label and --cortex\n");
      exit(1);
    }
    sprintf(tmpstr,"%s/%s/label/%s.cortex.label",SUBJECTS_DIR,subject,hemi);
    labelpath = strcpyalloc(tmpstr);
  }
  return;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n",vcid);
  fprintf(fp,"%s\n",Progname);
  fprintf(fp,"FREESURFER_HOME %s\n",getenv("FREESURFER_HOME"));
  fprintf(fp,"SUBJECTS_DIR %s\n",getenv("SUBJECTS_DIR"));
  fprintf(fp,"cwd       %s\n",cwd);
  fprintf(fp,"cmdline   %s\n",cmdline);
  fprintf(fp,"timestamp %s\n",VERcurTimeStamp());
  fprintf(fp,"sysname   %s\n",uts.sysname);
  fprintf(fp,"hostname  %s\n",uts.nodename);
  fprintf(fp,"machine   %s\n",uts.machine);
  fprintf(fp,"user      %s\n",VERuser());
  fprintf(fp,"subject %s\n",subject);
  fprintf(fp,"hemi     %s\n",hemi);
  fprintf(fp,"surfname %s\n",surfname);
  if (inpath)   fprintf(fp,"input  %s\n",inpath);
  if (maskpath) fprintf(fp,"mask  %s\n",maskpath);
  fprintf(fp,"infwhm %lf\n",infwhm);
  fprintf(fp,"ingstd %lf\n",ingstd);
  fprintf(fp,"synth %d\n",synth);
  if (synth) {
    fprintf(fp,"synth-frames %d\n",nframes);
    fprintf(fp,"seed  %d\n",SynthSeed);
  }
  if (Xfile) fprintf(fp,"xfile %s\n",Xfile);
  if (sumfile) fprintf(fp,"sumfile  %s\n",sumfile);
  if (outpath) fprintf(fp,"out  %s\n",outpath);

  return;
}

/*-------------------------------------------------------------
  FixGroupAreaTest() - the purpose of this function is to see
  how the fwhm changes when accounting for differences in area
  of group subjects (such as fsaverage). The output file will
  have 3 columns: (1) niterations, (2) corresponding fwhm 
  when not fixing, (3) corresponding fwhm  when fixing.
  This should have no effect any more (4/9/10)
  -------------------------------------------------------------*/
int FixGroupAreaTest(MRIS *surf, char *outfile)
{
  double fwhm, fwhmfix;
  int niters;
  FILE *fp;

  fp = fopen(outfile,"w");
  if(fp == NULL){
    printf("ERROR: cannot open %s\n",outfile);
    return(1);
  }

  // fwhm is for unfixed
  for(fwhm = 1; fwhm < 50; fwhm ++){
    unsetenv("FIX_VERTEX_AREA");
    niters = MRISfwhm2niters(fwhm,surf);
    setenv("FIX_VERTEX_AREA","1",1);
    fwhmfix = MRISniters2fwhm(niters, surf);
    printf("%3d %g %g\n",niters,fwhm,fwhmfix);
    fprintf(fp,"%4d %6.2lf %6.2lf\n",niters,fwhm,fwhmfix);
    fflush(fp);
  }
  unsetenv("FIX_VERTEX_AREA");
  fclose(fp);

  return(0);
}

/*
  DHiters2fwhm(MRIS *surf, int vtxno, int niters, char *outfile) This
  is a function that replicates how Don Hagler computes FWHM from
  number of iterations, namely by smoothing a delta function on the
  surface, then counting the area of the vertices above half the max.
  Also fits to fwhm = beta*sqrt(k) model.
*/
double DHiters2fwhm(MRIS *surf, int vtxno, int niters, char *outfile)
{
  int k, nhits,c;
  MRI *mri;
  double XtX, Xty, vXty, b, bv;
  double f, fn, areasum, fwhmRet;
  double fwhm[1000], fwhmv[1000], fn2sum[1000], fwhmdng;
  FILE *fp;

  mri  = MRIalloc(surf->nvertices,1,1,MRI_FLOAT);
  MRIsetVoxVal(mri,vtxno,0,0,0, 100);
  XtX = 0;
  Xty = 0;
  vXty = 0;
  for(k = 0; k < niters; k++){
    MRISsmoothMRI(surf, mri, 1, NULL, mri);
    f = MRIgetVoxVal(mri,vtxno,0,0,0); // = max
    nhits = 0; // number of vertices over max/2
    areasum = 0.0; // area of vertices over max/2
    fn2sum[k] = 0;
    for(c=0; c < surf->nvertices; c++){
      fn = MRIgetVoxVal(mri,c,0,0,0); 
      fn2sum[k] += (fn*fn);
      if(fn > f/2.0){
	nhits++;
	areasum += surf->vertices[c].area;
      }
    }
    fn2sum[k] /= (100.0*100.0);
    fwhm[k]  = 2*sqrt(areasum/M_PI); // fwhm in mm
    fwhmv[k] = 2*sqrt(nhits/M_PI);   // fwhm in vertices
    if(k > 3){
      // Accumulate for fitting
      XtX += (k+1); 
      Xty += (fwhm[k]*sqrt((double)(k+1)));
      vXty += (fwhmv[k]*sqrt((double)(k+1)));
    }
  }
  fwhmRet = fwhm[k-1];

  // Fit
  b  = Xty/XtX;
  bv = vXty/XtX;
  printf("#DH %6d %7.4f %7.4f %lf %lf\n",
	 vtxno,b,bv,surf->total_area,surf->avg_vertex_dist);

  if(outfile != NULL){
    // Iteration   MeasFWHM FitFWHM  DNGfwhm MeasFWHMv FitFWHMv VRF
    fp = fopen(outfile,"w");
    fprintf(fp,"#DH %6d %7.4f %7.4f %lf %lf\n",
	   vtxno,b,bv,surf->total_area,surf->avg_vertex_dist);
    fflush(fp);
    for(k = 0; k < niters; k++){
      //nitersdng = MRISfwhm2niters(fwhm[k],surf);
      fwhmdng = MRISniters2fwhm(k+1,surf);
      fprintf(fp,"%3d  %7.3f %7.3f  %7.3f   %7.3f %7.3f %8.3f\n",
	      k+1,
	      fwhm[k],sqrt(k+1.0)*b,fwhmdng,
	      fwhmv[k],sqrt(k+1.0)*bv,
	      1.0/fn2sum[k]);
      fflush(fp);
    }
    fclose(fp);
  }

  return(fwhmRet);
}

#if 0
double DHfwhm(MRIS *surf, MRI *mri)
{
  double fwhm,dv,dv2sum,v0,vn;
  int vtxno, nbrvtxno, nNNbrs, nthnbr, ndv;

  dv2sum = 0.0;
  ndv = 0;
  for(vtxno = 0; vtxno < surf->nvertices; vtxno++){
    v0 = MRIgetVoxVal(mri,vtxno,0,0,0);
    nNNbrs = SphSurf->vertices[vtxno].vnum;

    for (nthnbr = 0; nthnbr < nNNbrs; nthnbr++)    {
      nbrvtxno = surf->vertices[vtxno].v[nthnbr];
      if (surf->vertices[nbrvtx].ripflag) continue;
      vn = MRIgetVoxVal(mri,nbrvtxno,0,0,0);
      dv2sum += ((v0-vn)*(v0-vn));
      ndv++;
    }
  }
  return(0.0);
}
#endif

int SetHop(int CenterVtx, MRI_SURFACE *Surf, int HopNo, int MaxHopNo)
{
  int nbr, nbr_vtx, nbr_hopno, nbr_has_been_center;

  if(HopNo >= MaxHopNo) return(0);

  Surf->vertices[CenterVtx].valbak = 1; // has been a center

  for (nbr=0; nbr < Surf->vertices[CenterVtx].vnum; nbr++) {
    nbr_vtx   = Surf->vertices[CenterVtx].v[nbr];
    nbr_hopno = Surf->vertices[nbr_vtx].undefval;
    if(nbr_hopno != 0 && nbr_hopno < HopNo + 1) continue;
    Surf->vertices[nbr_vtx].undefval = HopNo + 1;
  }

  for (nbr=0; nbr < Surf->vertices[CenterVtx].vnum; nbr++) {
    nbr_vtx = Surf->vertices[CenterVtx].v[nbr];
    nbr_has_been_center = Surf->vertices[nbr_vtx].valbak;
    if(nbr_has_been_center) continue;
    SetHop(nbr_vtx, Surf, HopNo+1, MaxHopNo);
  }

  return(0);
}


SURFHOPLIST *SetSurfHopListAlloc(MRI_SURFACE *Surf, int nHops)
{
  SURFHOPLIST *shl;
  int nthhop;
  shl = (SURFHOPLIST*) calloc(sizeof(SURFHOPLIST),1);
  shl->nhops = nHops;
  shl->hit = (char *)calloc(sizeof(char),Surf->nvertices);
  shl->nperhop = (int *)calloc(sizeof(int),nHops);
  shl->nperhop_alloced = (int *)calloc(sizeof(int),nHops);
  shl->vtxlist = (int **)calloc(sizeof(int *),nHops);
  for(nthhop=0; nthhop < nHops; nthhop++){
    shl->nperhop_alloced[nthhop] = 100;
    shl->vtxlist[nthhop] = (int *)calloc(sizeof(int),100);
  }
  return(shl);
}

int SurfHopListFree(SURFHOPLIST **shl0)
{
  SURFHOPLIST *shl = *shl0;  
  int nthhop;
  for(nthhop=0; nthhop < shl->nhops; nthhop++)
    free(shl->vtxlist[nthhop]);
  free(shl->vtxlist);
  free(shl->nperhop_alloced);
  free(shl->nperhop);
  free(shl->hit);
  free(shl);
  *shl0 = NULL;
  return(0);
}

/*!
  \fn SURFHOPLIST *SetSurfHopList(int CenterVtx, MRI_SURFACE *Surf, int nHops)
  \brief Fills in a SURFHOPLIST structure. This is a structure that indicates
  which vertices are a given number of links (hops) away. This can be used to
  set neighborhoods or compute the spatial autocorrelation function.
*/

SURFHOPLIST *SetSurfHopList(int CenterVtx, MRI_SURFACE *Surf, int nHops)
{
  SURFHOPLIST *shl;
  int nthhop, nhits, nthvtx, vtxno, nthnbr, nbr_vtxno;

  shl = SetSurfHopListAlloc(Surf, nHops);
  shl->cvtx = CenterVtx;

  // 0 hop is the center vertex
  shl->nperhop[0] = 1;
  shl->vtxlist[0][0] = CenterVtx;
  shl->hit[CenterVtx] = 1;

  // go through all hops
  for(nthhop = 1; nthhop < nHops; nthhop++){
    nhits = 0;
    // go through each vertex of the previous hop
    for(nthvtx = 0; nthvtx < shl->nperhop[nthhop-1]; nthvtx++){
      vtxno = shl->vtxlist[nthhop-1][nthvtx];
      // go through the neighbors of this vertex
      for(nthnbr=0; nthnbr < Surf->vertices[vtxno].vnum; nthnbr++) {
	nbr_vtxno  = Surf->vertices[vtxno].v[nthnbr];
	if(shl->hit[nbr_vtxno]) continue; // ignore if it has been hit already
	shl->hit[nbr_vtxno] = 1;
	if(nhits >= shl->nperhop_alloced[nthhop]){
	  // realloc if need to
	  shl->nperhop_alloced[nthhop] += 100;
	  shl->vtxlist[nthhop] = (int *)realloc(shl->vtxlist[nthhop],sizeof(int)*shl->nperhop_alloced[nthhop]);
	}
	// assign the vertex
	shl->vtxlist[nthhop][nhits] = nbr_vtxno;
	nhits ++;
      }
    }
    shl->nperhop[nthhop] = nhits;
  }

  // This assigns the hop number to the undefval, good for debugging
  for(nthhop = 0; nthhop < nHops; nthhop++){
    //printf("nper hop %d %6d\n",nthhop,shl->nperhop[nthhop]);
    for(nthvtx = 0; nthvtx < shl->nperhop[nthhop]; nthvtx++){
      vtxno = shl->vtxlist[nthhop][nthvtx];
      Surf->vertices[vtxno].undefval = nthhop;
    }
  }

  return(shl);
}

/*-----------------------------------------------------------------------
  MRISarN() - computes spatial autocorrelation function at each vertex
  by averaging the ARs within the neighborhood of a vertex. Note: does
  not try to take into account different distances between
  neighbors. N is the number of hops.  arN will have N frames where
  frame 0 is always 1, frame 1 is the average AR between the vertex
  and the vertices 1 hop away, etc.
  -----------------------------------------------------------------------*/
MRI *MRISarN(MRIS *surf, MRI *src, MRI *mask, MRI *arN, int N)
{
  int **crslut, nvox, vtx;

  nvox = src->width * src->height * src->depth;
  if (surf->nvertices != nvox){
    printf("ERROR: MRISarN: Surf/Src dimension mismatch.\n");
    return(NULL);
  }

  if(arN == NULL){
    arN = MRIcloneBySpace(src, MRI_FLOAT, N);
    if (arN == NULL){
      printf("ERROR: could not alloc\n");
      return(NULL);
    }
  }

  //Build LUT to map from col,row,slice to vertex
  crslut = MRIScrsLUT(surf,src);

  #ifdef _OPENMP
  #pragma omp parallel for 
  #endif
  for (vtx = 0; vtx < surf->nvertices; vtx++) {
    int nnbrs, frame, nbrvtx, nthnbr, c,r,s;
    int cnbr, rnbr,snbr, nnbrs_actual;
    double valvtx, valnbr, arsum, sumsqvtx, vtxvar, sumsqnbr, sumsqx, nbrvar;
    SURFHOPLIST *shl;
    int nthhop;

    if(surf->vertices[vtx].ripflag) continue;
    c = crslut[0][vtx];
    r = crslut[1][vtx];
    s = crslut[2][vtx];
    if(mask) if (MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;

    // Compute variance for this vertex
    sumsqvtx = 0;
    for(frame = 0; frame < src->nframes; frame ++) {
      valvtx = MRIFseq_vox(src,c,r,s,frame);
      sumsqvtx += (valvtx*valvtx);
    }
    if(sumsqvtx == 0) continue;  // exclude voxels with 0 variance
    vtxvar = sumsqvtx/src->nframes;

    // Zero hop is always 1
    MRIFseq_vox(arN,c,r,s,0) = 1;

    // Create structure to manage the multiple hops for this vertex
    shl = SetSurfHopList(vtx, surf, N);

    // loop through hops
    for(nthhop = 1; nthhop < N; nthhop++){
      nnbrs = shl->nperhop[nthhop];
      nnbrs_actual = 0;
      arsum = 0;
      // loop through the neighbors nthhop links away
      for(nthnbr = 0; nthnbr < nnbrs; nthnbr++){
	nbrvtx = shl->vtxlist[nthhop][nthnbr];
	if(surf->vertices[nbrvtx].ripflag) continue;
	cnbr = crslut[0][nbrvtx];
	rnbr = crslut[1][nbrvtx];
	snbr = crslut[2][nbrvtx];
	if(mask) if(MRIgetVoxVal(mask,cnbr,rnbr,snbr,0) < 0.5) continue;
	sumsqnbr = 0;
	sumsqx   = 0;
	for(frame = 0; frame < src->nframes; frame ++){
	  valvtx = MRIFseq_vox(src,c,r,s,frame);
	  valnbr = MRIFseq_vox(src,cnbr,rnbr,snbr,frame);
	  sumsqnbr += (valnbr*valnbr);
	  sumsqx   += (valvtx*valnbr);
	}
	if(sumsqnbr==0) continue;
	nbrvar = sumsqnbr/src->nframes;
	arsum += (sumsqx/src->nframes)/sqrt(vtxvar*nbrvar);
	nnbrs_actual ++;
      }/* end loop over hop neighborhood */
      
      if(nnbrs_actual != 0) MRIFseq_vox(arN,c,r,s,nthhop) = (arsum/nnbrs_actual);
    } /* end loop over hop */

    SurfHopListFree(&shl);

  } /* end loop over vertex */

  MRIScrsLUTFree(crslut);

  return(arN);
}

/*-----------------------------------------------------------------------
  MRISsmoothKernel() - 
  kernel = ACF^2
  -----------------------------------------------------------------------*/
MRI *MRISsmoothKernel(MRIS *surf, MRI *src, MRI *mask, MRI *mrikern, MATRIX *globkern, int SqrFlag, MRI *out)
{
  int vtx, **crslut, nvox;
  double *kern;
  int n,nhops;

  if(mrikern && globkern){
    printf("ERROR: MRISsmoothKernel(): cannot spec both mrikern and globkern\n");
    return(NULL);
  }

  nvox = src->width * src->height * src->depth;
  if (surf->nvertices != nvox){
    printf("ERROR: MRISsmoothKernel(): Surf/Src dimension mismatch.\n");
    return(NULL);
  }

  if(out == NULL){
    out = MRIcloneBySpace(src, MRI_FLOAT, src->nframes);
    if (out == NULL){
      printf("ERROR: MRISsmoothKernel(): could not alloc\n");
      return(NULL);
    }
  }
  else {
    if(out == src){
      printf("ERROR: MRISsmoothKernel(): cannot run in-place\n");
      return(NULL);
    }
    // Zero the output
  }

  if(mrikern)  nhops = mrikern->nframes;
  if(globkern) nhops = globkern->rows;
  printf("nhops = %d\n",nhops);
  kern = (double *) calloc(nhops,sizeof(double));
  if(globkern) {
    for(n = 0; n < nhops; n++) 
      kern[n] = globkern->rptr[n+1][1];
    if(SqrFlag){
      for(n = 0; n < nhops; n++) 
	kern[n] = kern[n]*kern[n];
    }
  }

  //Build LUT to map from col,row,slice to vertex
  crslut = MRIScrsLUT(surf,src);

  #ifdef _OPENMP
  #pragma omp parallel for 
  #endif
  for (vtx = 0; vtx < surf->nvertices; vtx++) {
    int nnbrs, frame, nbrvtx, nthnbr, c,r,s;
    int cnbr, rnbr,snbr, nnbrs_actual;
    double vtxval,*vkern,ksum,kvsum;
    SURFHOPLIST *shl;
    int nthhop;

    if(surf->vertices[vtx].ripflag) continue;
    c = crslut[0][vtx];
    r = crslut[1][vtx];
    s = crslut[2][vtx];
    if(mask) if (MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;

    if(mrikern){
      vkern = (double *) calloc(nhops,sizeof(double));
      for(nthhop = 0; nthhop < nhops; nthhop++) 
	vkern[nthhop] = MRIgetVoxVal(mrikern,vtx,0,0,nthhop);
      if(SqrFlag){
	for(nthhop = 0; nthhop < nhops; nthhop++) 
	  vkern[nthhop] = vkern[nthhop]*vkern[nthhop];
      }
    }
    else vkern = kern;

    // Create structure to manage the multiple hops for this vertex
    shl = SetSurfHopList(vtx, surf, nhops);

    for(frame = 0; frame < src->nframes; frame ++){
      // loop through hops and neighbors
      if(frame ==0) ksum = 0;
      kvsum = 0;
      nnbrs_actual = 0;
      for(nthhop = 0; nthhop < nhops; nthhop++){
	nnbrs = shl->nperhop[nthhop];
	// loop through the neighbors nthhop links away
	for(nthnbr = 0; nthnbr < nnbrs; nthnbr++){
	  nbrvtx = shl->vtxlist[nthhop][nthnbr];
	  if(surf->vertices[nbrvtx].ripflag) continue;
	  cnbr = crslut[0][nbrvtx];
	  rnbr = crslut[1][nbrvtx];
	  snbr = crslut[2][nbrvtx];
	  if(mask) if(MRIgetVoxVal(mask,cnbr,rnbr,snbr,0) < 0.5) continue;
	  vtxval = MRIFseq_vox(src,cnbr,rnbr,snbr,frame);
	  kvsum += (vtxval*kern[nthhop]);
	  if(frame ==0) ksum += kern[nthhop];
	  nnbrs_actual ++;
	} /* end loop over hop neighborhood */
      } /* end loop over hop */
      if(nnbrs_actual != 0) MRIFseq_vox(out,c,r,s,frame) = (kvsum/ksum);
      //printf("%5d %d %d %d %d%g\n",vtx,c,r,s,nnbrs_actual,ksum);
    } // end loop over frame
    SurfHopListFree(&shl);

    if(mrikern) free(vkern);
  } /* end loop over vertex */

  free(kern);
  MRIScrsLUTFree(crslut);

  return(out);
}

