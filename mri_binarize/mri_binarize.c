/**
 * @file  mri_binarize.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:04 $
 *    $Revision: 1.7 $
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


// $Id: mri_binarize.c,v 1.7 2006/12/29 02:09:04 nicks Exp $

/*
  BEGINHELP

Program to binarize a volume (or volume-encoded surface file). Can also
be used to merge with other binarizations.

--i invol

Input volume to be binarized.

--min min
--max max

Minimum and maximum thresholds. If the value at a voxel is >= min and
<= max, then its value in the output will be 1 (or --binval), otherwise
it is 0 (or --binvalnot) or the value of the merge volume at that voxel.
By default, min = -infinity and max = +infinity, but you must set one
of the thresholds.

--o outvol

Path to output volume.

--binval    binval
--binvalnot binvalnot

Value to use for those voxels that are in the threshold range
(--binval) or out of the range (--binvalnot). These must be integer
values. binvalnot only applies when a merge volume is not specified.

--frame frameno

Use give frame of the input. 0-based. Default is 0.

--merge mergevol

Merge binarization with the mergevol. If the voxel is within the threshold
range, then its value will be binval. If not, then it will inherit its
value from the value at that voxel in mergevol. mergevol must be the same
dimension as the input volume. Combining this with --binval allows you
to construct crude segmentations.

--mask maskvol
--mask-thresh thresh

Mask input with mask. The mask volume is itself binarized at thresh
(default is 0.5). If a voxel is not in the mask, then it will be assigned
binvalnot or the value from the merge volume.

  ENDHELP
*/

/*
  BEGINUSAGE

  ENDUSAGE
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
#include "mri2.h"
#include "fio.h"
#include "version.h"
#include "label.h"
#include "matrix.h"
#include "annotation.h"
#include "fmriutils.h"
#include "cmdargs.h"
#include "fsglm.h"
#include "pdf.h"
#include "fsgdf.h"
#include "timer.h"
#include "matfile.h"
#include "volcluster.h"
#include "surfcluster.h"


static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_binarize.c,v 1.7 2006/12/29 02:09:04 nicks Exp $";
char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

char *InVolFile=NULL;
char *OutVolFile=NULL;
char *MergeVolFile=NULL;
char *MaskVolFile=NULL;
double MinThresh, MaxThresh;
int MinThreshSet=0, MaxThreshSet=0;
int BinVal=1;
int BinValNot=0;
int frame=0;
int DoAbs=0;

MRI *InVol,*OutVol,*MergeVol,*MaskVol;
double MaskThresh = 0.5;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs, c, r, s, nhits, InMask;
  double val,outputval,maskval;

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
  dump_options(stdout);

  // Load the input volume
  InVol = MRIread(InVolFile);
  if (InVol==NULL) exit(1);
  if (frame >= InVol->nframes) {
    printf("ERROR: requested frame=%d >= nframes=%d\n",
           frame,InVol->nframes);
    exit(1);
  }

  if (DoAbs) {
    printf("Removing sign from input\n");
    MRIabs(InVol,InVol);
  }

  // Load the merge volume (if needed)
  if (MergeVolFile) {
    MergeVol = MRIread(MergeVolFile);
    if (MergeVol==NULL) exit(1);
    if (MergeVol->width != InVol->width) {
      printf("ERROR: dimension mismatch between input and merge volumes\n");
      exit(1);
    }
    if (MergeVol->height != InVol->height) {
      printf("ERROR: dimension mismatch between input and merge volumes\n");
      exit(1);
    }
    if (MergeVol->depth != InVol->depth) {
      printf("ERROR: dimension mismatch between input and merge volumes\n");
      exit(1);
    }
  }

  // Load the mask volume (if needed)
  if (MaskVolFile) {
    MaskVol = MRIread(MaskVolFile);
    if (MaskVol==NULL) exit(1);
    if (MaskVol->width != InVol->width) {
      printf("ERROR: dimension mismatch between input and mask volumes\n");
      exit(1);
    }
    if (MaskVol->height != InVol->height) {
      printf("ERROR: dimension mismatch between input and mask volumes\n");
      exit(1);
    }
    if (MaskVol->depth != InVol->depth) {
      printf("ERROR: dimension mismatch between input and mask volumes\n");
      exit(1);
    }
  }

  // Prepare the output volume
  OutVol = MRIalloc(InVol->width,InVol->height,InVol->depth,MRI_INT);
  if (OutVol == NULL) exit(1);
  MRIcopyHeader(InVol, OutVol);

  // Binarize
  InMask = 1;
  nhits = 0;
  for (c=0; c < InVol->width; c++) {
    for (r=0; r < InVol->height; r++) {
      for (s=0; s < InVol->depth; s++) {
        val = MRIgetVoxVal(InVol,c,r,s,frame);
        if (MaskVol) {
          maskval = MRIgetVoxVal(MaskVol,c,r,s,0);
          if (maskval > MaskThresh) InMask = 1;
          else                     InMask = 0;
        }
        if ((MinThreshSet && (val < MinThresh)) ||
            (MaxThreshSet && (val > MaxThresh)) ||
            !InMask ) {
          // Not in the Range
          if (MergeVol)
            outputval = MRIgetVoxVal(MergeVol,c,r,s,0);
          else
            outputval = BinValNot;
        } else {
          outputval = BinVal;
          nhits ++;
        }
        MRIsetVoxVal(OutVol,c,r,s,0,outputval);
      }
    }
  }
  printf("Found %d values in range\n",nhits);

  // Save output
  MRIwrite(OutVol,OutVolFile);

  printf("mri_binarize done\n");
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
    else if (!strcasecmp(option, "--abs")) DoAbs = 1;

    else if (!strcasecmp(option, "--i")) {
      if (nargc < 1) CMDargNErr(option,1);
      InVolFile = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--o")) {
      if (nargc < 1) CMDargNErr(option,1);
      OutVolFile = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--merge")) {
      if (nargc < 1) CMDargNErr(option,1);
      MergeVolFile = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--mask")) {
      if (nargc < 1) CMDargNErr(option,1);
      MaskVolFile = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--mask-thresh")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&MaskThresh);
      nargsused = 1;
    } else if (!strcasecmp(option, "--min")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&MinThresh);
      MinThreshSet = 1;
      nargsused = 1;
    } else if (!strcasecmp(option, "--max")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&MaxThresh);
      MaxThreshSet = 1;
      nargsused = 1;
    } else if (!strcasecmp(option, "--binval")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&BinVal);
      nargsused = 1;
    } else if (!strcasecmp(option, "--binvalnot")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&BinValNot);
      nargsused = 1;
    } else if (!strcasecmp(option, "--frame")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&frame);
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
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --i invol  : input volume \n");
  printf("   --min min  : min thresh (def is -inf)\n");
  printf("   --max max  : max thresh (def is +inf)\n");
  printf("   --o outvol : output volume \n");
  printf("   \n");
  printf("   --binval    val    : set vox within thresh to val (default is 1) \n");
  printf("   --binvalnot notval : set vox outside range to notval (default is 0) \n");
  printf("   --frame frameno    : use 0-based frame of input (default is 0) \n");
  printf("   --merge mergevol   : merge with mergevolume \n");
  printf("   --mask maskvol       : must be within mask \n");
  printf("   --mask-thresh thresh : set thresh for mask (def is 0.5) \n");
  printf("   --abs : take abs of invol first (ie, make unsigned)\n");
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
  printf("Program to binarize a volume (or volume-encoded surface file). Can also\n");
  printf("be used to merge with other binarizations.\n");
  printf("\n");
  printf("--i invol\n");
  printf("\n");
  printf("Input volume to be binarized. \n");
  printf("\n");
  printf("--min min\n");
  printf("--max max\n");
  printf("\n");
  printf("Minimum and maximum thresholds. If the value at a voxel is >= min and\n");
  printf("<= max, then its value in the output will be 1 (or --binval), otherwise\n");
  printf("it is 0 (or --binvalnot) or the value of the merge volume at that voxel.\n");
  printf("By default, min = -infinity and max = +infinity, but you must set one\n");
  printf("of the thresholds.\n");
  printf("\n");
  printf("--o outvol\n");
  printf("\n");
  printf("Path to output volume.\n");
  printf("\n");
  printf("--binval    binval\n");
  printf("--binvalnot binvalnot\n");
  printf("\n");
  printf("Value to use for those voxels that are in the threshold range\n");
  printf("(--binval) or out of the range (--binvalnot). These must be integer\n");
  printf("values. binvalnot only applies when a merge volume is not specified.\n");
  printf("\n");
  printf("--frame frameno\n");
  printf("\n");
  printf("Use give frame of the input. 0-based. Default is 0.\n");
  printf("\n");
  printf("--merge mergevol\n");
  printf("\n");
  printf("Merge binarization with the mergevol. If the voxel is within the threshold\n");
  printf("range, then its value will be binval. If not, then it will inherit its\n");
  printf("value from the value at that voxel in mergevol. mergevol must be the same\n");
  printf("dimension as the input volume. Combining this with --binval allows you\n");
  printf("to construct crude segmentations.\n");
  printf("\n");
  printf("--mask maskvol\n");
  printf("--mask-thresh thresh\n");
  printf("\n");
  printf("Mask input with mask. The mask volume is itself binarized at thresh \n");
  printf("(default is 0.5). If a voxel is not in the mask, then it will be assigned\n");
  printf("binvalnot or the value from the merge volume.\n");
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
  if (InVolFile == NULL) {
    printf("ERROR: must specify input volume\n");
    exit(1);
  }
  if (OutVolFile == NULL) {
    printf("ERROR: must specify output volume\n");
    exit(1);
  }
  if (MinThreshSet == 0 && MaxThreshSet == 0) {
    printf("ERROR: must specify minimum and/or maximum threshold\n");
    exit(1);
  }
  if (MaxThreshSet && MinThreshSet && MaxThresh < MinThresh) {
    printf("ERROR: max thresh = %g < min thresh = %g\n",
           MaxThresh,MinThresh);
    exit(1);
  }
  return;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n",vcid);
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  fprintf(fp,"\n");
  fprintf(fp,"input      %s\n",InVolFile);
  fprintf(fp,"frame      %d\n",frame);
  fprintf(fp,"output     %s\n",OutVolFile);
  if (MinThreshSet)
    fprintf(fp,"min        %g\n",MinThresh);
  else
    fprintf(fp,"min        -infinity\n");
  if (MaxThreshSet)
    fprintf(fp,"max        %g\n",MaxThresh);
  else
    fprintf(fp,"max        +infinity\n");
  fprintf(fp,"binval        %d\n",BinVal);
  fprintf(fp,"binvalnot     %d\n",BinValNot);
  if (MergeVolFile)
    fprintf(fp,"merge      %s\n",MergeVolFile);
  if (MaskVolFile) {
    fprintf(fp,"mask       %s\n",MaskVolFile);
    fprintf(fp,"maskthresh %lf\n",MaskThresh);
  }
  return;
}
