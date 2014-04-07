/**
 * @file  mri_gtmseg.c
 * @brief Creates a segmentation used with the geometric transfer matrix.
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Douglas N. Greve
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2014/04/07 19:50:09 $
 *    $Revision: 1.2 $
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
// $Id: mri_gtmseg.c,v 1.2 2014/04/07 19:50:09 greve Exp $

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
#include "timer.h"
#include "mri_identify.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#include "gtm.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_gtmseg.c,v 1.2 2014/04/07 19:50:09 greve Exp $";
char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

GTMSEG *gtmseg;
char *OutVolFile=NULL;
char *SUBJECTS_DIR;
int GTMdefaultSegReplacmentList(GTMSEG *gtmseg);
int nthreads;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs,err;
  char tmpstr[5000],*stem;
  struct timeb  mytimer, timer;

  nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd,2000);
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");

  gtmseg = (GTMSEG *) calloc(sizeof(GTMSEG),1);
  gtmseg->USF = 2;
  gtmseg->dmax = 5.0;
  gtmseg->KeepHypo = 0;
  gtmseg->KeepCC = 0;
  gtmseg->apasfile = "apas+head.mgz";
  gtmseg->ctxannotfile = "aparc.annot";
  gtmseg->SubSegWM = 0;
  gtmseg->ctxlhbase = 1000;
  gtmseg->ctxrhbase = 2000;
  if(gtmseg->SubSegWM){
    gtmseg->wmannotfile = "lobes.annot";
    gtmseg->wmlhbase =  3200;
    gtmseg->wmrhbase =  4200;
  }
  else  gtmseg->wmannotfile = NULL;

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  if (argc == 0) usage_exit();
  parse_commandline(argc, argv);
  check_options();
  if (checkoptsonly) return(0);

  GTMdefaultSegReplacmentList(gtmseg);
  dump_options(stdout);

  TimerStart(&timer);
  printf("Starting MRIgtmSeg()\n"); fflush(stdout);
  err = MRIgtmSeg(gtmseg);
  if(err) exit(1);

  sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,gtmseg->subject,OutVolFile);
  printf("Writing output file to %s\n",tmpstr);
  err = MRIwrite(gtmseg->seg,tmpstr);
  if(err) exit(1);

  stem = IDstemFromName(OutVolFile);
  sprintf(tmpstr,"%s/%s/mri/%s.lta",SUBJECTS_DIR,gtmseg->subject,stem);
  printf("Writing lta file to %s\n",tmpstr);
  err=LTAwrite(gtmseg->anat2seg,tmpstr);
  if(err) exit(1);

  printf("mri_gtmseg finished in %g minutes\n",TimerStop(&mytimer)/60000.0);
  printf("mri_gtmseg done\n");fflush(stdout);

  exit(0);  return 0;
}
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/

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

    else if (!strcasecmp(option, "--o")) {
      if (nargc < 1) CMDargNErr(option,1);
      OutVolFile = pargv[0];
      nargsused = 1;
    } 

    else if(!strcasecmp(option, "--s")) {
      if(nargc < 1) CMDargNErr(option,1);
      gtmseg->subject = pargv[0];
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--usf")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&gtmseg->USF);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--apas")) {
      if (nargc < 1) CMDargNErr(option,1);
      gtmseg->apasfile = pargv[0];
      nargsused = 1;
    } 

    else if(!strcasecmp(option, "--ctx-annot")) {
      if(nargc < 3) CMDargNErr(option,3);
      gtmseg->ctxannotfile = pargv[0];
      sscanf(pargv[1],"%d",&gtmseg->ctxlhbase);
      sscanf(pargv[2],"%d",&gtmseg->ctxrhbase);
      nargsused = 3;
    } 

    else if(!strcasecmp(option, "--subseg-wm"))    gtmseg->SubSegWM = 1;
    else if(!strcasecmp(option, "--no-subseg-wm")) gtmseg->SubSegWM = 0;
    else if(!strcasecmp(option, "--wm-annot")) {
      if(nargc < 3) CMDargNErr(option,3);
      gtmseg->wmannotfile = pargv[0];
      sscanf(pargv[1],"%d",&gtmseg->wmlhbase);
      sscanf(pargv[2],"%d",&gtmseg->wmrhbase);
      gtmseg->SubSegWM = 1;
      nargsused = 3;
    } 
    else if(!strcasecmp(option, "--dmax")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%f",&gtmseg->dmax);
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--keep-hypo"))    gtmseg->KeepHypo = 1;
    else if(!strcasecmp(option, "--no-keep-hypo")) gtmseg->KeepHypo = 0;
    else if(!strcasecmp(option, "--keep-cc"))    gtmseg->KeepCC = 1;
    else if(!strcasecmp(option, "--no-keep-cc")) gtmseg->KeepCC = 0;

    else if (!strcmp(option, "--sd") || !strcmp(option, "-SDIR")) {
      if(nargc < 1) CMDargNErr(option,1);
      setenv("SUBJECTS_DIR",pargv[0],1);
      SUBJECTS_DIR = getenv("SUBJECTS_DIR");
      nargsused = 1;
    } 

    else if (!strcasecmp(option, "--threads")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nthreads);
      #ifdef _OPENMP
      omp_set_num_threads(nthreads);
      #endif
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--max-threads")){
      nthreads = 1;
      #ifdef _OPENMP
      nthreads = omp_get_max_threads();
      omp_set_num_threads(nthreads);
      #endif
    } 
    else if (!strcasecmp(option, "--max-threads-1")){
      nthreads = 1;
      #ifdef _OPENMP
      nthreads = omp_get_max_threads()-1;
      if(nthreads < 0) nthreads = 1;
      omp_set_num_threads(nthreads);
      #endif
    } 

    else {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if(CMDsingleDash(option))
        fprintf(stderr,"       Did you really mean -%s ?\n",option);
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/*-----------------------------------------------------------------------------------*/
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/*-----------------------------------------------------------------------------------*/
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --o outvol  : output volume (output will be subject/mri/outvol\n");
  printf("   --s subject : source subject \n");
  printf("   --usf USF : upsampling factor (default %d)\n",gtmseg->USF);
  printf("   --apas apasfile : defines extra-cerebral segmentations (%s)\n",gtmseg->apasfile);
  printf("   --ctx-annot annot lhbase rhbase : use annot to segment cortex (%s,%d,%d)\n",
	 gtmseg->ctxannotfile,gtmseg->ctxlhbase,gtmseg->ctxrhbase);
  printf("   --subseg-wm : turn on segmenting of WM into smaller parts (off by default)\n");
  printf("       sets wmannot to lobes.annot and lhbase,rhbase to 3200,4200\n");
  printf("   --wm-annot annot lhbase rhbase : use annot to subseg wm\n");
  printf("   --dmax dmax : distance from ctx for wmseg to be considered 'unsegmented' (%f)\n",gtmseg->dmax);
  printf("   --keep-hypo : do not convert WM hypointensities to a white matter label \n");
  printf("   --keep-cc : do not convert corpus callosum to a white matter label \n");
  printf("\n");
  #ifdef _OPENMP
  printf("   --threads N : use N threads (with Open MP)\n");
  printf("   --threads-max : use the maximum allowable number of threads for this computer\n");
  printf("   --threads-max-1 : use one less than the maximum allowable number of threads for this computer\n");
  #endif
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  printf("%s\n", vcid) ;
  printf("\n");
}
/*-----------------------------------------------------------------------------------*/
static void print_help(void) {
  print_usage() ;
  printf("WARNING: this program is not yet tested!\n");
  exit(1) ;
}
/*-----------------------------------------------------------------------------------*/
static void print_version(void) {
  printf("%s\n", vcid) ;
  exit(1) ;
}
/*-----------------------------------------------------------------------------------*/
static void check_options(void) {
  return;
}
/*-----------------------------------------------------------------------------------*/
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n",vcid);
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  GTMSEGprint(gtmseg, stdout);
  return;
}

/*-----------------------------------------------------------------------------------*/
int MakeGTMSeg(char *subject, int USF, int SubSegWM, int LabelCCAsWM, int LabelHypoAsWM, float dmax, char *outsegfile)
{
  GTMSEG *gtmseg;
  int nlist, *srclist, *targlist, err;
  char *SUBJECTS_DIR,tmpstr[5000];

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");

  gtmseg = (GTMSEG *) calloc(sizeof(GTMSEG),1);
  gtmseg->subject = strcpyalloc(subject);
  gtmseg->apasfile = "apas+head.mgz";
  gtmseg->ctxannotfile = "aparc.annot";
  gtmseg->ctxlhbase = 1000;
  gtmseg->ctxrhbase = 2000;
  gtmseg->SubSegWM = SubSegWM;
  gtmseg->KeepHypo = !LabelHypoAsWM;
  gtmseg->KeepCC = !LabelCCAsWM;
  gtmseg->dmax = dmax;
  gtmseg->USF = USF;
  if(gtmseg->SubSegWM){
    gtmseg->wmannotfile = "lobes.annot";
    gtmseg->wmlhbase =  3200;
    gtmseg->wmrhbase =  4200;
  }
  else  gtmseg->wmannotfile = NULL;

  srclist  = &(gtmseg->srclist[0]);
  targlist = &(gtmseg->targlist[0]);

  nlist = 0;
  srclist[nlist] = 1033; targlist[nlist] = 1030; nlist++; // temppole=stg
  srclist[nlist] = 2033; targlist[nlist] = 2030; nlist++; // temppole=stg
  srclist[nlist] = 1034; targlist[nlist] = 1030; nlist++; // transtemp=stg
  srclist[nlist] = 2034; targlist[nlist] = 1030; nlist++; // transtemp=stg
  srclist[nlist] = 1001; targlist[nlist] = 1015; nlist++; // bankssts=mtg
  srclist[nlist] = 2001; targlist[nlist] = 2015; nlist++; // bankssts=mtg
  srclist[nlist] = 1032; targlist[nlist] = 1027; nlist++; // frontpole=rmf
  srclist[nlist] = 2032; targlist[nlist] = 2027; nlist++; // frontpole=rmf
  //srclist[nlist] = 1016; targlist[nlist] = 1006; nlist++; // parahip=entorhinal ?
  //srclist[nlist] = 2016; targlist[nlist] = 2006; nlist++; // parahip=entorhinal ?

  // There should not be any cortex unknown after MRIannot2CorticalSeg()
  srclist[nlist] = 1000; targlist[nlist] =    0; nlist++; // cortex unknown
  srclist[nlist] = 2000; targlist[nlist] =    0; nlist++; // cortex unknown

  // Should I replace subcorts before hires seg?
  srclist[nlist] =   85; targlist[nlist] =    0; nlist++; // optic chiasm
  srclist[nlist] =    4; targlist[nlist] =   24; nlist++; // LLatVent
  srclist[nlist] =    5; targlist[nlist] =   24; nlist++; // LInfLatVent
  srclist[nlist] =   14; targlist[nlist] =   24; nlist++; // 3rd
  srclist[nlist] =   15; targlist[nlist] =   24; nlist++; // 4th
  srclist[nlist] =   72; targlist[nlist] =   24; nlist++; // 5th
  srclist[nlist] =   31; targlist[nlist] =   24; nlist++; // LChoroidP ?
  srclist[nlist] =   43; targlist[nlist] =   24; nlist++; // RLatVent
  srclist[nlist] =   44; targlist[nlist] =   24; nlist++; // RInfLatVent
  srclist[nlist] =   63; targlist[nlist] =   24; nlist++; // RChoroidP ?
  srclist[nlist] =   30; targlist[nlist] =   24; nlist++; // LVessel ?
  srclist[nlist] =   62; targlist[nlist] =   24; nlist++; // RVessel ?
  srclist[nlist] =   80; targlist[nlist] =   24; nlist++; // non-WM-hypo ?

  /* Repace CC segments with one CC if not unsegmenting CC */
  if(gtmseg->KeepCC){
    printf(" Relabeling CC as a single segmentation\n");
    srclist[nlist] =  251; targlist[nlist] =  192; nlist++; 
    srclist[nlist] =  252; targlist[nlist] =  192; nlist++; 
    srclist[nlist] =  253; targlist[nlist] =  192; nlist++; 
    srclist[nlist] =  254; targlist[nlist] =  192; nlist++; 
    srclist[nlist] =  255; targlist[nlist] =  192; nlist++; 
  }

  gtmseg->nlist = nlist;

  GTMSEGprint(gtmseg, stdout);
  MRIgtmSeg(gtmseg);
  sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,gtmseg->subject,outsegfile);
  printf("Writing output file to %s\n",tmpstr);
  err = MRIwrite(gtmseg->seg,tmpstr);
  if(err) return(1);

  char *stem = IDstemFromName(outsegfile);
  sprintf(tmpstr,"%s/%s/mri/%s.lta",SUBJECTS_DIR,gtmseg->subject,stem);
  printf("Writing lta file to %s\n",tmpstr);
  err=LTAwrite(gtmseg->anat2seg,tmpstr);
  if(err) return(1);

  printf("MakeGTMSeg() done\n");fflush(stdout);

  return(0);
}


/*-----------------------------------------------------------------------------*/
int GTMdefaultSegReplacmentList(GTMSEG *gtmseg)
{

  int nlist, *srclist, *targlist;

  srclist  = &(gtmseg->srclist[0]);
  targlist = &(gtmseg->targlist[0]);

  nlist = 0;
  srclist[nlist] = 1033; targlist[nlist] = 1030; nlist++; // temppole=stg
  srclist[nlist] = 2033; targlist[nlist] = 2030; nlist++; // temppole=stg
  srclist[nlist] = 1034; targlist[nlist] = 1030; nlist++; // transtemp=stg
  srclist[nlist] = 2034; targlist[nlist] = 1030; nlist++; // transtemp=stg
  srclist[nlist] = 1001; targlist[nlist] = 1015; nlist++; // bankssts=mtg
  srclist[nlist] = 2001; targlist[nlist] = 2015; nlist++; // bankssts=mtg
  srclist[nlist] = 1032; targlist[nlist] = 1027; nlist++; // frontpole=rmf
  srclist[nlist] = 2032; targlist[nlist] = 2027; nlist++; // frontpole=rmf
  //srclist[nlist] = 1016; targlist[nlist] = 1006; nlist++; // parahip=entorhinal ?
  //srclist[nlist] = 2016; targlist[nlist] = 2006; nlist++; // parahip=entorhinal ?

  // There should not be any cortex unknown after MRIannot2CorticalSeg()
  srclist[nlist] = 1000; targlist[nlist] =    0; nlist++; // cortex unknown
  srclist[nlist] = 2000; targlist[nlist] =    0; nlist++; // cortex unknown

  // Should I replace subcorts before hires seg?
  srclist[nlist] =   85; targlist[nlist] =    0; nlist++; // optic chiasm
  srclist[nlist] =    4; targlist[nlist] =   24; nlist++; // LLatVent
  srclist[nlist] =    5; targlist[nlist] =   24; nlist++; // LInfLatVent
  srclist[nlist] =   14; targlist[nlist] =   24; nlist++; // 3rd
  srclist[nlist] =   15; targlist[nlist] =   24; nlist++; // 4th
  srclist[nlist] =   72; targlist[nlist] =   24; nlist++; // 5th
  srclist[nlist] =   31; targlist[nlist] =   24; nlist++; // LChoroidP ?
  srclist[nlist] =   43; targlist[nlist] =   24; nlist++; // RLatVent
  srclist[nlist] =   44; targlist[nlist] =   24; nlist++; // RInfLatVent
  srclist[nlist] =   63; targlist[nlist] =   24; nlist++; // RChoroidP ?
  srclist[nlist] =   30; targlist[nlist] =   24; nlist++; // LVessel ?
  srclist[nlist] =   62; targlist[nlist] =   24; nlist++; // RVessel ?
  srclist[nlist] =   80; targlist[nlist] =   24; nlist++; // non-WM-hypo ?

  /* Repace CC segments with one CC if not unsegmenting CC */
  if(gtmseg->KeepCC){
    srclist[nlist] =  251; targlist[nlist] =  192; nlist++; 
    srclist[nlist] =  252; targlist[nlist] =  192; nlist++; 
    srclist[nlist] =  253; targlist[nlist] =  192; nlist++; 
    srclist[nlist] =  254; targlist[nlist] =  192; nlist++; 
    srclist[nlist] =  255; targlist[nlist] =  192; nlist++; 
  }

  gtmseg->nlist = nlist;
  return(0);
}
