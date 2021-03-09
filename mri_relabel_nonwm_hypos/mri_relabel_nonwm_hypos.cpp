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



/*!
\file dummy.c
\brief Example c file that can be used as a template.
\author Douglas Greve

*/



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
#include "mri2.h"
#include "colortab.h"
#include "fsenv.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

int PrintSegIds(int nsegs, int *segidlist, int *outsegidlist);
int DefaultSegIds(int *segidlist, int *outsegidlist);

const char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

char *InputSegFile=NULL,*OutputSegFile=NULL;
int nsegs;
int segidlist[1000],outsegidlist[1000];

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs,err;
  MRI *seg, *newseg;

  nargs = handleVersionOption(argc, argv, "mri_relabel_nonwm_hypos");
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

  seg = MRIread(InputSegFile);
  if(seg == NULL) exit(1);

  newseg = MRIrelabelNonWMHypos(seg, segidlist, nsegs, outsegidlist);
  if(newseg == NULL) exit(1);

  err = MRIwrite(newseg,OutputSegFile);
  if(err) exit(1);

  printf("mri_relabel_nonwm_hypos done\n");
  exit(0);
}
/* ---------------------------------------------------------------*/
/* ---------------------------------------------------------------*/
/* ---------------------------------------------------------------*/
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

    else if (!strcasecmp(option, "--i")) {
      if (nargc < 1) CMDargNErr(option,1);
      InputSegFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--o")) {
      if (nargc < 1) CMDargNErr(option,1);
      OutputSegFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--seg")) {
      if (nargc < 2) CMDargNErr(option,2);
      sscanf(pargv[0],"%d",&segidlist[nsegs]);
      sscanf(pargv[1],"%d",&outsegidlist[nsegs]);
      nsegs++;
      nargsused = 2;
    } 
    else if (!strcasecmp(option, "--seg-default")) 
      nsegs = DefaultSegIds(segidlist, outsegidlist);
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
  printf("   --i inputseg  : input segmentation with non-wm-hypos labeled as 80, 81, or 82\n");
  printf("   --o outputseg : segmentation with non-wm-hypos relabeled \n");
  printf("   --seg seg newseg : relabel hypos adjacent to seg as newseg (can use multiple --seg args)\n");
  printf("   --seg-default : use default relabeling scheme\n");
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
  printf("This program relabels non-WM hypointensities (80, 81, 82) based\n");
  printf("on proximity to a nearby label. The set of nearby labels is \n");
  printf("specified by the user using either --seg or --seg default.\n");
  printf("Eg, --seg 11 270 would replace non-wm hypos next to left caudate with \n");
  printf("the value 270. \nAlternatively, --seg-default will create the following list\n");
  printf("\n");
  nsegs = DefaultSegIds(segidlist, outsegidlist);
  PrintSegIds(nsegs, segidlist, outsegidlist);
  printf("\n");
  printf("The user will need to create a custom color to view/analyze the new segmentation\n");
  printf("The color table must be formated like $FREESURFER_HOME/FreeSurferColorLUT.txt, eg\n");
  printf("   270 left-caudate-nonwm-hypos   120 120 0 0\n");
  printf("This can be used with tkmedit, like tkmedit subject nu.mgz -seg newseg.mgz newcolortable\n");
  printf("or when running mri_segstats (passing the colortable with --ctab)\n");
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
  if(InputSegFile == NULL){
    printf("ERROR: must spec input segmentation\n");
    exit(1);
  }
  if(OutputSegFile == NULL){
    printf("ERROR: must spec output segmentation\n");
    exit(1);
  }
  if(nsegs == 0){
    printf("ERROR: must spec relabling scheme with --seg or --seg-default\n");
    exit(1);
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
  PrintSegIds(nsegs, segidlist, outsegidlist);
  return;
}

/*-------------------------------------------------------*/
int DefaultSegIds(int *segidlist, int *outsegidlist)
{
  int n;
  n = 0;
  segidlist[n] = 11; outsegidlist[n] = 270; n++;
  segidlist[n] = 12; outsegidlist[n] = 271; n++;
  segidlist[n] = 13; outsegidlist[n] = 272; n++;
  segidlist[n] = 17; outsegidlist[n] = 273; n++;
  segidlist[n] = 18; outsegidlist[n] = 274; n++;
  segidlist[n] = 26; outsegidlist[n] = 275; n++;
  segidlist[n] = 10; outsegidlist[n] = 276; n++;

  segidlist[n] = 50; outsegidlist[n] = 280; n++;
  segidlist[n] = 51; outsegidlist[n] = 281; n++;
  segidlist[n] = 52; outsegidlist[n] = 282; n++;
  segidlist[n] = 53; outsegidlist[n] = 283; n++;
  segidlist[n] = 54; outsegidlist[n] = 284; n++;
  segidlist[n] = 58; outsegidlist[n] = 285; n++;
  segidlist[n] = 49; outsegidlist[n] = 286; n++;
  
  return(n);
}

/*-------------------------------------------------------*/
int PrintSegIds(int nsegs, int *segidlist, int *outsegidlist)
{
  FSENV *fsenv;
  COLOR_TABLE *ct;
  char tmpstr[50000];
  int n,segid;

  fsenv = FSENVgetenv();
  sprintf(tmpstr,"%s/FreeSurferColorLUT.txt",fsenv->FREESURFER_HOME);
  ct = CTABreadASCII(tmpstr);

  for(n=0; n < nsegs; n++){
    segid = segidlist[n];
    printf("%4d %-25s will be relabeled as %4d\n",segid,ct->entries[segid]->name,outsegidlist[n]);
  }

  return(0);
}

