/*
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


// mri_voldiff.c

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
#include "mri2.h"
#include "fio.h"
#include "version.h"
#include "label.h"
#include "annotation.h"
#include "fmriutils.h"
#include "cmdargs.h"
#include "fsglm.h"

// Exit codes. See dump_exit_codes.
#define DIMENSION_EC  2
#define PRECISION_EC  3
#define RESOLUTION_EC 4
#define VOX2RAS_EC    5
#define PIXEL_EC      6

static void dump_exit_codes(FILE *fp);

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);

int main(int argc, char *argv[]) ;
const char *Progname = NULL;

std::string vol1File, vol2File;
MRI *vol1=NULL, *vol2=NULL;
MATRIX *vox2ras1,*vox2ras2;

int debug = 0, checkoptsonly = 0;
char tmpstr[2000];

int AllowResolution = 0;
int AllowPrecision = 0;
int AllowVox2RAS = 0;

double vox2ras_thresh = 0;
double pixdiff_thresh = 0;

/*--------------------------------------------------*/
int main(int argc, char **argv) {
  int nargs,r,c;
  struct utsname uts;
  char *cmdline, cwd[2000];
  double maxdiff,d;
  int cmax,rmax,smax,fmax;

  nargs = handleVersionOption(argc, argv, "mri_voldiff");
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

  printf("\n");
  printf("%s\n",getVersion().c_str());
  printf("cwd %s\n",cwd);
  printf("cmdline %s\n",cmdline);
  printf("sysname  %s\n",uts.sysname);
  printf("hostname %s\n",uts.nodename);
  printf("machine  %s\n",uts.machine);

  dump_options(stdout);

  vol1 = MRIread(vol1File.c_str());
  if (vol1 == NULL) exit(1);
  vol2 = MRIread(vol2File.c_str());
  if (vol2 == NULL) exit(1);

  vox2ras1 = MRIxfmCRS2XYZ(vol1,0);
  vox2ras2 = MRIxfmCRS2XYZ(vol2,0);

  if (vol1->xsize   != vol2->xsize ||
      vol1->ysize   != vol2->ysize ||
      vol1->zsize   != vol2->zsize) {
    printf("volumes differ in resolution\n");
    if (!AllowResolution) exit(RESOLUTION_EC);
    printf("  but continuing\n");
  }

  if (vol1->type != vol2->type) {
    printf("volumes differ in precision\n");
    if (!AllowPrecision) exit(PRECISION_EC);
    printf("  but continuing\n");
  }

  for (c=1; c <= 4; c++) {
    for (r=1; r <= 4; r++) {
      d = fabs(vox2ras1->rptr[c][r] - vox2ras2->rptr[c][r]);
      if (fabs(d) > vox2ras_thresh) {
        printf("volumes differ in vox2ras %d %d %g %g\n",c,r,
               vox2ras1->rptr[c][r],vox2ras2->rptr[c][r]);
        if (!AllowVox2RAS) exit(VOX2RAS_EC);
        printf("  but continuing\n");
        c=4;
        r=4;
      }
    }
  }

  if (vol1->width   != vol2->width ||
      vol1->height  != vol2->height ||
      vol1->depth   != vol2->depth ||
      vol1->nframes != vol2->nframes) {
    printf("volumes differ in dimension (%d %d %d %d) (%d %d %d %d) \n",
           vol1->width,vol1->height,vol1->depth,vol1->nframes,
           vol2->width,vol2->height,vol2->depth,vol2->nframes);
    exit(DIMENSION_EC);
  }

  maxdiff = MRImaxAbsDiff(vol1,vol2,&cmax,&rmax,&smax,&fmax);
  printf("pixdiff %g at %d %d %d %d\n",maxdiff,cmax,rmax,smax,fmax);

  if (maxdiff > pixdiff_thresh) {
    printf("volumes differ in pixel data\n");
    exit(PIXEL_EC);
  }

  printf("volumes are consistent\n");
  printf("mri_voldiff done\n");
  return(0);
  exit(0);

}
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/

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
    else if (!strcasecmp(option, "--allow-res")) AllowResolution = 1;
    else if (!strcasecmp(option, "--allow-prec")) AllowPrecision = 1;
    else if (!strcasecmp(option, "--allow-vox2ras")) AllowVox2RAS = 1;

    else if (!strcmp(option, "--v1")) {
      if (nargc < 1) CMDargNErr(option,1);
      vol1File = fio_fullpath(pargv[0]);
      nargsused = 1;
    } else if (!strcmp(option, "--v2")) {
      if (nargc < 1) CMDargNErr(option,1);
      vol2File = fio_fullpath(pargv[0]);
      nargsused = 1;
    } else if (!strcmp(option, "--pix")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&pixdiff_thresh);
      nargsused = 1;
    } else if (!strcmp(option, "--vox2ras")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&vox2ras_thresh);
      nargsused = 1;
    } else {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if (CMDsingleDash(option))
        fprintf(stderr,"       Did you really mean -%s ?\n",option);
      exit(1);
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
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --v1 first input volume \n");
  printf("   --v2 second input volume \n");
  printf("\n");
  printf("   --vox2ras thresh \n");
  printf("   --pix thresh \n");
  printf("\n");
  printf("   --allow-prec\n");
  printf("   --allow-res\n");
  printf("   --allow-vox2ras\n");
  printf("\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  std::cout << getVersion() << std::endl;
  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf(
    "Determines whether two volumes are different. The difference\n"
    "can be in pixel data or in the dimension, precision, resolution\n"
    "or geometry. If there are no errors of differences in the volumes,\n"
    "then it exits with 0. If an error occurs, then it exits with 1.\n"
    "If the volumes are different, it will exit with one of the \n"
    "following codes:\n");

  printf("\n");
  dump_exit_codes(stdout);
  printf("\n");

  printf(
    "Some of the exit conditions can be allowed with --allow-prec,\n"
    "--allow-res, --allow-vox2ras.\n"
    "\n");


  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/* --------------------------------------------- */
static void check_options(void) {
  if (vol1File.size() == 0) {
    printf("ERROR: must specify a vol1 file\n");
    exit(1);
  }
  if (vol2File.size() == 0) {
    printf("ERROR: must specify a vol2 file\n");
    exit(1);
  }
  return;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"vol1    %s\n",vol1File.c_str());
  fprintf(fp,"vol2    %s\n",vol2File.c_str());
  fprintf(fp,"pix thresh  %g\n",pixdiff_thresh);
  fprintf(fp,"vox2ras thresh %g\n",vox2ras_thresh);
  return;
}

/* --------------------------------------------- */
static void dump_exit_codes(FILE *fp) {
  fprintf(fp,"dimensions inconsistent   %d\n",DIMENSION_EC);
  fprintf(fp,"precision  inconsistent   %d\n",PRECISION_EC);
  fprintf(fp,"resolution inconsistent   %d\n",RESOLUTION_EC);
  fprintf(fp,"vox2ras    inconsistent   %d\n",VOX2RAS_EC);
  fprintf(fp,"pixel      inconsistent   %d\n",PIXEL_EC);
}
