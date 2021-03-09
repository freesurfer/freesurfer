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


/*
   Name: mkxsubjreg.c
   Author: DouglasN. Greve
   Date: 8/24/03
   Purpose: Create a new registration matrix that will map a functional
   to the orig of another subject.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <sys/time.h>

#include "error.h"
#include "diag.h"
#include "proto.h"
#include "matrix.h"
#include "mri.h"
#include "registerio.h"
#include "fio.h"
#include "version.h"
#include "resample.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static void dump_options(FILE *fp);
static int  singledash(char *flag);

int main(int argc, char *argv[]) ;

const char *Progname = NULL;

char *srcregpath  = NULL;
char *targregpath = NULL;
const char *targsubj    = NULL;
const char *xfmrname    = "talairach.xfm";
char *subjects_dir = NULL;
char *fvolid = NULL;

MRI *SrcMRI, *TargMRI, *FuncMRI;

int fixtkreg = 1;
int debug = 0;
/*---------------------------------------------------------------*/
int main(int argc, char **argv) {
  char *srcsubj;
  float betplaneres, inplaneres, intensity;
  MATRIX *R, *Xsrc, *invXsrc, *Xtarg, *Rtarg;
  int float2int, err;

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if (argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();
  dump_options(stdout);

  /* Load the registration matrix, fix if necessary */
  err = regio_read_register(srcregpath, &srcsubj, &inplaneres,
                            &betplaneres, &intensity, &R,
                            &float2int);
  if (err) exit(1);
  printf("---- Input registration matrix --------\n");
  MatrixPrint(stdout,R);

  if (float2int == FLT2INT_TKREG && fixtkreg) {
    if (fvolid == NULL) {
      printf("ERROR: the input registration file requires that you "
             "supply an example functional volume with --fvol\n");
      exit(1);
    }
    FuncMRI = MRIreadHeader(fvolid,MRI_VOLUME_TYPE_UNKNOWN);
    if (FuncMRI==NULL) exit(1);
    printf("INFO: making tkreg matrix compatible with round\n");
    R = MRIfixTkReg(FuncMRI,R);
    printf("---- Fixed input registration matrix --------\n");
    MatrixPrint(stdout,R);
  }
  float2int = FLT2INT_ROUND;

  /* Load the source subject xfm */
  Xsrc = DevolveXFM(srcsubj,NULL,xfmrname);
  if (Xsrc == NULL) exit(1);
  invXsrc = MatrixInverse(Xsrc,NULL);

  /* Load the target subject xfm */
  Xtarg = DevolveXFM(targsubj,NULL,xfmrname);
  if (Xtarg == NULL) exit(1);

  /* Rtarg = R*inv(Xsrc)*Xtarg */
  Rtarg = MatrixMultiply(R,invXsrc,NULL);
  Rtarg = MatrixMultiply(Rtarg,Xtarg,Rtarg);

  printf("---- New registration matrix --------\n");
  MatrixPrint(stdout,Rtarg);

  err = regio_write_register(targregpath, targsubj, inplaneres,
                             betplaneres, intensity, Rtarg, float2int);

  if (err) {
    printf("ERROR: could not write to %s\n",targregpath);
    exit(1);
  }

  return(0);
  exit(0);
}
/* --------------------------------------------- */
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("Required Arguments\n");
  printf("   --srcreg   srcreg.dat\n");
  printf("   --targreg  targreg.dat\n");
  printf("\n");
  printf("Optional Arguments\n");
  printf("   --targsubj subjid : default is talairach\n");
  printf("   --xfm xfmrname    : xfm file name relative to transforms\n");
  printf("   --sd subjects_dir : default is env SUBJECTS_DIR\n");
  printf("   --fvol funcvol    : path to example functional volume\n");
  printf("   --help\n");
  printf("   --version\n");
  printf("\n");
}
/* ------------------------------------------------------------------ */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;
  int nargs;

  if (argc < 1) usage_exit();

  nargs = handleVersionOption(argc, argv, "mkxsubjreg");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

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

    else if (!strcmp(option, "--srcreg")) {
      if (nargc < 1) argnerr(option,1);
      srcregpath = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--targreg")) {
      if (nargc < 1) argnerr(option,1);
      targregpath = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--targsubj")) {
      if (nargc < 1) argnerr(option,1);
      targsubj = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--xfm")) {
      if (nargc < 1) argnerr(option,1);
      xfmrname = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--fvol")) {
      if (nargc < 1) argnerr(option,1);
      fvolid = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--sd")) {
      if (nargc < 1) argnerr(option,1);
      subjects_dir = pargv[0];
      nargsused = 1;
    } else {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if (singledash(option))
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
static void print_help(void) {
  print_usage() ;
  printf("Creates new registration matrix that maps from the functional \n");
  printf("volume of the source subject to the orig of the target subject\n");
  printf("through the talairach transform.\n");
  printf("\n");

  exit(1) ;
}
/* --------------------------------------------- */
static void check_options(void) {
  if (srcregpath == NULL) {
    printf("A source registration file must be supplied\n");
    exit(1);
  }
  if (targregpath == NULL) {
    printf("A target/output registration file must be supplied\n");
    exit(1);
  }
  if (targsubj == NULL) {
    targsubj = "talairach";
    printf("No target subject supplied, assuming %s\n",targsubj);
  }
  if (subjects_dir == NULL) {
    subjects_dir = getenv("SUBJECTS_DIR") ;
    if (subjects_dir==NULL) {
      printf("ERROR: SUBJECTS_DIR not defined\n");
      exit(1);
    }
  }

  return;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"srcregpath     %s\n",srcregpath);
  fprintf(fp,"targregpath    %s\n",targregpath);
  fprintf(fp,"targsubject    %s\n",targsubj);
  fprintf(fp,"xfm            %s\n",xfmrname);
  fprintf(fp,"subjects_dir   %s\n",subjects_dir);
  fprintf(fp,"Diag Level     %d\n",Gdiag_no);
  fprintf(fp,"%s\n", getVersion().c_str());
  return;
}


/*---------------------------------------------------------------*/
/*---------------------------------------------------------------*/
static int singledash(char *flag) {
  int len;
  len = strlen(flag);
  if (len < 2) return(0);
  if (flag[0] == '-' && flag[1] != '-') return(1);
  return(0);
}
/* --------------------------------------------- */
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/* --------------------------------------------- */
static void argnerr(char *option, int n) {
  if (n==1)
    fprintf(stderr,"ERROR: %s flag needs %d argument\n",option,n);
  else
    fprintf(stderr,"ERROR: %s flag needs %d arguments\n",option,n);
  exit(-1);
}
