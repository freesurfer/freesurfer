/**
 * @brief Computes and merge differences in segmentation.
 *
 * This program computes and merges differences in segmentation volumes
 * (eg, aseg.auto.mgz and aseg.mgz) primarily for the purpose of
 * managing manual edits to aseg.mgz.
 */
/*
 * Original Author: greve
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
  BEGINHELP

  This program computes and merges differences in segmentation volumes
  (eg, aseg.auto.mgz and aseg.mgz) primarily for the purpose of
  managing manual edits to aseg.mgz. 

  When computing a difference, it compares seg1 to seg2 at each
  voxel. If they are the same, then the diff volume voxel gets a value
  of 256 (Voxel-Unchanged in FreeSurferColorLUT.txt). If they are
  different, then it takes the value of seg2. The diff seg can be
  loaded as a segmentation in tkmedit. If there is no difference,
  then the difference volume will not be written unless you use
  --diff-force.

  When merging a difference, a voxel in the merged seg will take the
  value of the input seg if the diff-in has a value of 256. If the
  diff-in value is something other than 256, then the merged value
  will be the diff-in value.

  EXAMPLES:

  # Determine manual edits 
  mri_seg_diff --seg1 aseg.auto.mgz --seg2 aseg.mgz --diff aseg.manedits.mgz

  # Merge manual edits 
  mri_seg_diff --seg aseg.auto.mgz --merged aseg.mgz --diff-in aseg.manedits.mgz

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

const char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

char *Seg1File=NULL;
char *Seg2File=NULL;
char *DiffFile=NULL;
char *InDiffFile=NULL;
char *MergedFile=NULL;
int ForceDiff = 0;

char *subject, *SUBJECTS_DIR;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int nargs, DiffFlag=0;
  MRI *seg1, *seg2, *diff;

  nargs = handleVersionOption(argc, argv, "mri_seg_diff");
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

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if (SUBJECTS_DIR == NULL) {
    printf("ERROR: SUBJECTS_DIR not defined in environment\n");
    exit(1);
  }

  seg1 = MRIread(Seg1File);
  if(seg1 == NULL) exit(1);

  // Compute diff of segs
  if(DiffFile != NULL){
    printf("Computing difference between segmentations\n");
    seg2 = MRIread(Seg2File);
    if(seg2 == NULL) exit(1);
    diff = MRIsegDiff(seg1,seg2,&DiffFlag);
    if(diff == NULL) exit(1);
    if(DiffFlag == 0) {
      printf("No difference found.\n");
      if(! ForceDiff) exit(0);
      printf(" ... but saving diff file anyway.\n");
    } else printf("A difference found, saving.\n");

    MRIwrite(diff,DiffFile);
    exit(0);
  }

  printf("Merging difference segmentation\n");
  diff = MRIread(InDiffFile);
  if(diff == NULL) exit(1);

  seg2 = MRIsegMergeDiff(seg1, diff);
  if(seg2 == NULL) exit(1);

  MRIwrite(seg2,MergedFile);

  exit(0);
  return(0);
}
/*-------------------------------------------------------*/
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
    else if (!strcasecmp(option, "--diff-force"))  ForceDiff = 1;

    else if (!strcasecmp(option, "--seg1") || !strcasecmp(option, "--seg")) {
      if (nargc < 1) CMDargNErr(option,1);
      Seg1File = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--seg2")) {
      if (nargc < 1) CMDargNErr(option,1);
      Seg2File = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--diff")) {
      if (nargc < 1) CMDargNErr(option,1);
      DiffFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--diff-in")) {
      if (nargc < 1) CMDargNErr(option,1);
      InDiffFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--merged")) {
      if (nargc < 1) CMDargNErr(option,1);
      MergedFile = pargv[0];
      nargsused = 1;
    } 
    else {
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
/*-------------------------------------------------------*/
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/*-------------------------------------------------------*/
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("Options for creating a diff\n");
  printf("   --seg1 seg1 : first segmentation (eg, unedited)\n");
  printf("   --seg2 seg2 : second segmentation (eg, edited)\n");
  printf("   --diff  outdifffile : output diff seg volume \n");
  printf("   --diff-force : force creation of a diff even if no diff\n");
  printf("\n");
  printf("Options for merging with a diff\n");
  printf("   --seg seg  : source seg (eg, unedited) \n");
  printf("   --diff-in indifffile : input diff seg volume \n");
  printf("   --merged mergefile   : unedited merged with diff \n");
  printf("\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  std::cout << getVersion() << std::endl;
  printf("\n");
}
/*-------------------------------------------------------*/
static void print_help(void) {
  print_usage() ;
printf("\n");
printf("  This program computes and merges differences in segmentation volumes\n");
printf("  (eg, aseg.auto.mgz and aseg.mgz) primarily for the purpose of\n");
printf("  managing manual edits to aseg.mgz. \n");
printf("\n");
printf("  When computing a difference, it compares seg1 to seg2 at each\n");
printf("  voxel. If they are the same, then the diff volume voxel gets a value\n");
printf("  of 256 (Voxel-Unchanged in FreeSurferColorLUT.txt). If they are\n");
printf("  different, then it takes the value of seg2. The diff seg can be\n");
printf("  loaded as a segmentation in tkmedit. If there is no difference,\n");
printf("  then the difference volume will not be written unless you use\n");
printf("  --diff-force.\n");
printf("\n");
printf("  When merging a difference, a voxel in the merged seg will take the\n");
printf("  value of the input seg if the diff-in has a value of 256. If the\n");
printf("  diff-in value is something other than 256, then the merged value\n");
printf("  will be the diff-in value.\n");
printf("\n");
printf("  EXAMPLES:\n");
printf("\n");
printf("  # Determine manual edits \n");
printf("  mri_seg_diff --seg1 aseg.auto.mgz --seg2 aseg.mgz --diff aseg.manedits.mgz\n");
printf("\n");
printf("  # Merge manual edits \n");
printf("  mri_seg_diff --seg aseg.auto.mgz --merged aseg.mgz --diff-in aseg.manedits.mgz\n");
printf("\n");
  exit(1) ;
}
/*-------------------------------------------------------*/
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/*-------------------------------------------------------*/
static void check_options(void) {
  if(Seg1File == NULL){
    printf("ERROR: need an input segmentation\n");
    exit(1);
  }
  if(DiffFile == NULL && InDiffFile == NULL){
    printf("ERROR: no diff file specified\n");
    exit(1);
  }
  if(DiffFile != NULL && InDiffFile != NULL){
    printf("ERROR: cannot specify --diff and --diff-in\n");
    exit(1);
  }
  if(DiffFile != NULL){
    if(Seg2File == NULL){
      printf("ERROR: --seg2 with --diff\n");
      exit(1);
    }
    if(MergedFile != NULL){
      printf("ERROR: cannot spec --merged with --diff\n");
      exit(1);
    }
  }
  if(InDiffFile != NULL){
    if(MergedFile == NULL){
      printf("ERROR: must spec --merged with --diff-in\n");
      exit(1);
    }
    if(Seg2File != NULL){
      printf("ERROR: cannot spec --seg2 with --diff-in\n");
      exit(1);
    }
  }

  return;
}
/*-------------------------------------------------------*/
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n", getVersion().c_str());
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  fprintf(fp,"Seg1     %s\n",Seg1File);
  fprintf(fp,"Seg2     %s\n",Seg2File);
  fprintf(fp,"Diff     %s\n",DiffFile);
  fprintf(fp,"InDiff   %s\n",InDiffFile);
  fprintf(fp,"Merged   %s\n",MergedFile);
  fprintf(fp,"ForceDiff %d\n",ForceDiff);

  return;
}
