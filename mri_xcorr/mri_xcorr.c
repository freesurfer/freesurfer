/**
 * @file  dummy.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2014/01/14 21:28:11 $
 *    $Revision: 1.1 $
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



/*!
\file dummy.c
\brief Example c file that can be used as a template.
\author Douglas Greve

*/


// $Id: mri_xcorr.c,v 1.1 2014/01/14 21:28:11 greve Exp $

/*
  BEGINHELP

  ENDHELP
*/

/*
  BEGINUSAGE

  ENDUSAGE
*/


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/utsname.h>

#include "macros.h"
#include "utils.h"
#include "fio.h"
#include "version.h"
#include "cmdargs.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "fmriutils.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_xcorr.c,v 1.1 2014/01/14 21:28:11 greve Exp $";
const char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

char *v1File=NULL;
char *v2File=NULL;
char *maskFile=NULL;
char *outFile=NULL;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs;
  MRI *v1=NULL, *v2=NULL, *mask=NULL, *xcorr=NULL;

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

  printf("Reading v1 %s\n",v1File);fflush(stdout);
  v1 = MRIread(v1File);
  if(v1 == NULL) exit(1);

  if(v2File){
    printf("Reading v2 %s\n",v2File);fflush(stdout);
    v2 = MRIread(v2File);
    if(v2 == NULL) exit(1);
  }

  if(maskFile){
    printf("Reading mask %s\n",maskFile); fflush(stdout);
    mask = MRIread(maskFile);
    if(mask == NULL) exit(1);
  }

  printf("Computing xcorr\n");fflush(stdout);
  xcorr = fMRIxcorr(v1, v2, mask, NULL);

  printf("Writing xcorr to %s\n",outFile);
  MRIwrite(xcorr,outFile);

  printf("mri_xcorr done\n");
  return(0);
  exit(0);
}
/*--------------------------------------------*/
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

    else if (!strcasecmp(option, "--v1")) {
      if (nargc < 1) CMDargNErr(option,1);
      v1File = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--v2")) {
      if (nargc < 1) CMDargNErr(option,1);
      v2File = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--mask")) {
      if (nargc < 1) CMDargNErr(option,1);
      maskFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--o")) {
      if (nargc < 1) CMDargNErr(option,1);
      outFile = pargv[0];
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
/*--------------------------------------------*/
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/*--------------------------------------------*/
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --v1 vol1file : volume1 \n");
  printf("   --v2 vol2file : volume2 (or leave dont spec for lrrev v1) \n");
  printf("   --mask maskfile : limit computations \n");
  printf("   --o outfile : output \n");
  printf("\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  printf("%s\n", vcid) ;
  printf("\n");
}
/*--------------------------------------------*/
static void print_help(void) {
  print_usage() ;
  printf("Computes voxelwise temporal correlation coefficient between two volumes. If v2 is not specified, then v2 is formed from v1 by left-right  reversing the columns.\n");
  exit(1) ;
}
/*--------------------------------------------*/
static void print_version(void) {
  printf("%s\n", vcid) ;
  exit(1) ;
}
/*--------------------------------------------*/
static void check_options(void) 
{
  if(v1File == NULL){
    printf("ERROR: need to spec v1 file\n");
    exit(1);
  }
  if(outFile == NULL){
    printf("ERROR: need to spec out file\n");
    exit(1);
  }
  if(v2File == NULL)
    printf("v2 not specified, reversing columns on v1 to make v2\n");
  return;
}
/*--------------------------------------------*/
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n",vcid);
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  return;
}
