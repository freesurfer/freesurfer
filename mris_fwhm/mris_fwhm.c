/**
 * @file  mris_fwhm.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2012/12/13 18:18:27 $
 *    $Revision: 1.37 $
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

Estimates the smoothness of a surface-based data set.

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

--sum sumfile

Prints ascii summary to sumfile.

--fwhm fwhm

Smooth by fwhm mm before estimating the fwhm. This is mainly good for
debuggging. But with --out can also be used to smooth data on the
surface (but might be better to use mri_surf2surf for this).

--niters-only <nitersfile>

Only report the number of iterations needed to achieve the FWHM given
by fwhm. If nitersfile is specified, the number of iterations is 
written to the file.

--out outfile

Save (possibly synthesized and/or smoothed) data to outfile. Automatically
detects format. Format must be one accepted as by mri_convert. Note: do
not use analyze or nifit as these cannot store more than 32k in a dimension.
mri_surf2surf can store surface data in those formats.

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

static char vcid[] = "$Id: mris_fwhm.c,v 1.37 2012/12/13 18:18:27 greve Exp $";
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

int FixGroupAreaTest(MRIS *surf, char *outfile);
char *GroupAreaTestFile = NULL;
char *nitersfile = NULL;
double DHiters2fwhm(MRIS *surf, int vtxno, int niters, char *outfile);
int DHvtxno=0, DHniters=0;
char *DHfile = NULL;
int UseCortexLabel = 0;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs, niters=0, Ntp, n, err;
  double fwhm = 0, ar1mn, ar1std, ar1max, avgvtxarea,ftmp, fwhmDH;
  double InterVertexDistAvg, InterVertexDistStdDev;
  MRI *ar1=NULL;
  FILE *fp;

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
    } else if (!strcasecmp(option, "--surf")) {
      if (nargc < 1) CMDargNErr(option,1);
      surfname = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--i")) {
      if (nargc < 1) CMDargNErr(option,1);
      inpath = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--mask")) {
      if (nargc < 1) CMDargNErr(option,1);
      maskpath = pargv[0];
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
    } else if (!strcasecmp(option, "--fwhm")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&infwhm);
      ingstd = infwhm/sqrt(log(256.0));
      nargsused = 1;
    } else if (!strcasecmp(option, "--synth-frames")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nframes);
      synth = 1;
      nargsused = 1;
    } else if (!strcasecmp(option, "--o")) {
      if (nargc < 1) CMDargNErr(option,1);
      outpath = pargv[0];
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
  printf("\n");
  printf("   --i input\n");
  printf("   --subject subject (--s)\n");
  printf("   --hemi hemi (--h)\n");
  printf("   --surf surf <white>\n");
  printf("   --label labelfile\n");
  printf("   --cortex : used hemi.cortex.label\n");
  printf("   --mask maskfile\n");
  printf("   --X x.mat : matlab4 detrending matrix\n");
  printf("   --detrend order : polynomial detrending (default 0)\n");
  printf("   --smooth-only : only smooth (implies --no-detrend)\n");
  printf("   --no-detrend : turn of poly detrending \n");
  printf("   --sqr : compute square of input before smoothing\n");
  printf("   --sum sumfile\n");
  printf("   --dat datfile (only contains fwhm)\n");
  printf("   --ar1dat ar1datfile (contains ar1mean ar1std)\n");
  printf("   --ar1 ar1vol : save spatial ar1 as an overlay\n");
  printf("   \n");
  printf("   --fwhm fwhm : apply before measuring\n");
  printf("   --niters-only <niters> : only report on niters for fwhm\n");
  printf("   --o output\n");
  printf("\n");
  printf("   --sd SUBJECTS_DIR \n");
  printf("   --synth \n");
  printf("   --synth-frames nframes : default is 10 \n");
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
printf("Estimates the smoothness of a surface-based data set.\n");
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
printf("--sum sumfile\n");
printf("\n");
printf("Prints ascii summary to sumfile.\n");
printf("\n");
printf("--fwhm fwhm\n");
printf("\n");
printf("Smooth by fwhm mm before estimating the fwhm. This is mainly good for\n");
printf("debuggging. But with --out can also be used to smooth data on the\n");
printf("surface (but might be better to use mri_surf2surf for this).\n");
printf("\n");
printf("--niters-only <nitersfile>\n");
printf("\n");
printf("Only report the number of iterations needed to achieve the FWHM given\n");
printf("by fwhm. If nitersfile is specified, the number of iterations is \n");
printf("written to the file.\n");
printf("\n");
printf("--out outfile\n");
printf("\n");
printf("Save (possibly synthesized and/or smoothed) data to outfile. Automatically\n");
printf("detects format. Format must be one accepted as by mri_convert. Note: do\n");
printf("not use analyze or nifit as these cannot store more than 32k in a dimension.\n");
printf("mri_surf2surf can store surface data in those formats.\n");
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
