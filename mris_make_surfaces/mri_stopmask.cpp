/**
 * @brief Creates a "stop mask" used to stop mri_place_surface when searching for the max gradient,
 * eg, in areas that may be dark such as lesions, VR spaces, ventricles, and we don't want the
 * surface to wander down into that regsion.
 */
/*
 * Original Author: Douglas N Greve 
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
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/utsname.h>
#include <errno.h>
#include "timer.h"
#include "utils.h"
#include "fio.h"
#include "mri2.h"
#include "error.h"
#include "diag.h"
#include "cmdargs.h"
#include "version.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);

struct utsname uts;
char *cmdline, cwd[2000];
int debug = 0, checkoptsonly = 0;

int main(int argc, char *argv[]) ;

const char *Progname = "mris_place_surfaces";
SCMstopMask sm;
char *stopmaskpath=NULL;
char *SUBJECTS_DIR=NULL;
char *subject=NULL;

/*--------------------------------------------------*/
int main(int argc, char **argv) 
{
  int nargs, i, msec;
  char *cmdline2, cwd[2000];

  nargs = handleVersionOption(argc, argv, "mris_place_surface");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd,2000);
  cmdline2 = argv2cmdline(argc,argv);

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  Gdiag |= DIAG_SHOW ;

  if(argc == 0) usage_exit();
  parse_commandline(argc, argv);
  check_options();
  if(checkoptsonly) return(0);

  MRI *stopmask = sm.getmask();
  MRIwrite(stopmask,stopmaskpath);

  printf("#VMPC# mri_stopmask VmPeak  %d\n",GetVmPeak());
  printf("mri_stopmask done\n");

  return(0);

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
    else if(!strcasecmp(option, "--version")) print_version() ;
    else if(!strcasecmp(option, "--debug"))   debug = 1;
    else if(!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if(!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;
    else if(!strcasecmp(option, "--aseg")){
      if(nargc < 1) CMDargNErr(option,1);
      sm.aseg  = MRIread(pargv[0]);
      if(!sm.aseg) exit(1);
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--wm")){
      if(nargc < 1) CMDargNErr(option,1);
      sm.wm  = MRIread(pargv[0]);
      if(!sm.wm) exit(1);
      sm.DoWM255=1;
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--no-wm")) sm.DoWM255=0;
    else if(!strcasecmp(option, "--bfs")){
      if(nargc < 1) CMDargNErr(option,1);
      sm.bfs  = MRIread(pargv[0]);
      if(!sm.bfs) exit(1);
      sm.DoBFS255=1;
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--no-bfs")) sm.DoBFS255=0;
    else if(!strcasecmp(option, "--filled")){
      if(nargc < 2) CMDargNErr(option,2);
      sm.filledauto  = MRIread(pargv[0]);
      if(!sm.filledauto) exit(1);
      sm.filled = MRIread(pargv[1]);
      if(!sm.filled) exit(1);
      sm.DoFilled=1;
      nargsused = 2;
    } 
    else if(!strcasecmp(option, "--no-filled")) sm.DoFilled=0;
    else if(!strcasecmp(option, "--lv"))  sm.DoLV=1;
    else if(!strcasecmp(option, "--no-lv"))  sm.DoLV=0;
    else if(!strcasecmp(option, "--wmsa")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&sm.WMSAErodeMM);
      sm.DoWMSA=1;
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--no-wmsa")) sm.DoWMSA=0;
    else if(!strcasecmp(option, "--o")){
      if(nargc < 1) CMDargNErr(option,1);
      stopmaskpath = pargv[0];
      nargsused = 1;
    } 
    else if(!strcmp(option, "--sd")){
      if(nargc < 1) CMDargNErr(option,1);
      printf("using %s as SUBJECTS_DIR...\n", pargv[0]) ;
      setenv("SUBJECTS_DIR",pargv[0],1);
      nargsused = 1;
    }
    else if(!strcmp(option, "--s")){
      if(nargc < 1) CMDargNErr(option,1);
      subject = pargv[0];
      nargsused = 1;
    }
    else {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if(CMDsingleDash(option)) fprintf(stderr," Did you really mean -%s ?\n",option);
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* --------------------------------------------- */
static void check_options(void) 
{
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if(stopmaskpath == NULL) {
    printf("ERROR: must spec an output file with --o\n");
    exit(1);
  }
  if(sm.DoLV && ! sm.aseg){
    printf("ERROR: --lv requires aseg\n");
    exit(1);
  }
  if(sm.DoWMSA && ! sm.aseg){
    printf("ERROR: --wmsa requires aseg\n");
    exit(1);
  }
  if(!sm.DoBFS255 && !sm.DoFilled && !sm.DoLV && !sm.DoWM255&& !sm.DoWMSA){
    printf("ERROR: nothing to do\n");
    exit(1);
  }

  return;
}

#include "mri_stopmask.help.xml.h"
static void print_usage(void){
  outputHelpXml(mri_stopmask_help_xml,mri_stopmask_help_xml_len);
}

static void
print_help(void)
{
  print_usage() ;
  exit(1) ;
}
/* ------------------------------------------------------ */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n", getVersion().c_str());
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  return;
}
