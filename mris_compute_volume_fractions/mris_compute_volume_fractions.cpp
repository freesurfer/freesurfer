/**
 * @brief Computes accurate volume fractions remaining within a surface
 *
 * This program computes an accurate estimate of the fraction of the volume 
 * remaining wihtin a surface. 
 */
/*
 * Original Author: Ender Konukoglu
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
#include <sys/utsname.h>

#include "macros.h"
#include "utils.h"
#include "fio.h"
#include "version.h"
#include "cmdargs.h"
#include "error.h"
#include "diag.h"
#include "mrisurf.h"
#include "mris_compVolFrac.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);

const char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

char *VolFile=NULL;
char *SurfFile = NULL;
char *OutFile = NULL;
double Accuracy = -1000.0;
int main(int argc, char *argv[]) {
  
  printf("working!\n");
  int nargs;

  nargs = handleVersionOption(argc, argv, "mris_compute_volume_fractions");

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
  volFraction v;v.frac = 0.0; v.err = 0.0;
  MRI* mri = MRIread(VolFile);
  MRI_SURFACE* mris = MRISread(SurfFile);    
  printf("running the computation...\n");
  MRI* mri_fractions = MRIcomputeVolumeFractionFromSurface(mris, Accuracy, mri, NULL);
  printf("computation is finished...\n");
  MRIwrite(mri_fractions, OutFile);

  return 0;
}

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

    else if (!strcasecmp(option, "--vol")) {
      if (nargc < 1) CMDargNErr(option,1);
      VolFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--surf")){
      if (nargc < 1) CMDargNErr(option,1);
      SurfFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--acc")){
      if (nargc < 1) CMDargNErr(option,1);
      Accuracy = atof(pargv[0]); 
      nargsused = 1; 
    }
    else if (!strcasecmp(option, "--out")){
      if (nargc < 1) CMDargNErr(option,1);
      OutFile = pargv[0];
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
/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void usage_exit(void)
\brief Prints usage and exits
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void print_usage(void)
\brief Prints usage and returns (does not exit)
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --vol volume_file : volume \n");
  printf("   --surf surface_file: surface\n");
  printf("   --acc accuracy: required accuracy\n");
  printf("   --out out_file: output volume file for the fractions\n");
  printf("\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  std::cout << getVersion() << std::endl;
  printf("\n");
}
/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void print_help(void)
\brief Prints help and exits
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void print_help(void) {
  print_usage() ;
  printf("WARNING: this program is not yet tested!\n");
  exit(1) ;
}
/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void print_version(void)
\brief Prints version and exits
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void check_options(void)
\brief Checks command-line options
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void check_options(void) {
  if(VolFile == NULL || SurfFile == NULL || Accuracy < 0 || OutFile == NULL)
    {
      print_usage(); 
      exit(1);
    }
    
  return;
}

/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void dump_options(FILE *fp)
\brief Prints command-line options to the given file pointer
\param FILE *fp - file pointer
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n", getVersion().c_str());
  fprintf(fp,"Working Directory: %s\n",cwd);
  fprintf(fp,"cmdline: %s\n",cmdline);
  /*
  fprintf(fp,"sysname:  %s\n",uts.sysname);
  fprintf(fp,"hostname: %s\n",uts.nodename);
  fprintf(fp,"machine:  %s\n",uts.machine);
  fprintf(fp,"user:     %s\n",VERuser());
  */
  return;
}

