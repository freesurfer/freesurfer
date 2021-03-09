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
\file mri_z2p.c
\brief Converts a z to a -log10(p)
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
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <unistd.h>

#include "macros.h"
#include "utils.h"
#include "mrisutils.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mri2.h"
#include "fio.h"
#include "version.h"
#include "cmdargs.h"
#include "randomfields.h"
#include "mri_identify.h"

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

char *ZVolFile=NULL;
char *PVolFile=NULL;
char *Log10PVolFile=NULL;
char *MaskVolFile=NULL;
int TwoSided = 1;

MRI *z,*p,*sig,*mask=NULL;
char *featdir = NULL;
const char *fmt;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs, nthzstat;
  char tmpstr[1000];

  nargs = handleVersionOption(argc, argv, "mri_z2p");
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
  if (checkoptsonly) exit(0);
  dump_options(stdout);

  //---------------------------------------------------------
  if(featdir == NULL){
    z = MRIread(ZVolFile);
    if (z==NULL) exit(1);

    if (MaskVolFile) {
      printf("Reading mask %s\n",MaskVolFile);
      mask = MRIread(MaskVolFile);
      if (mask==NULL) exit(1);
    }

    p = RFz2p(z, mask, TwoSided, NULL);
    if(PVolFile) MRIwrite(p,PVolFile);
    if(Log10PVolFile == NULL) exit(0);
    sig  = MRIlog10(p,mask,sig,1);
    MRIwrite(sig,Log10PVolFile);
    printf("mri_z2p: done\n");
    exit(0);
  }

  //---------------------------------------------------------
  // feat dir -------------------
  printf("Processing data in feat dir %s\n",featdir);
  if(fmt == NULL){
    if (strcmp(getenv("FSLOUTPUTTYPE"),"NIFTI")==0)    fmt = "nii";
    if (strcmp(getenv("FSLOUTPUTTYPE"),"NIFTI_GZ")==0) fmt = "nii.gz";
    if (strcmp(getenv("FSLOUTPUTTYPE"),"ANALYZE")==0)  fmt = "img";
    if(fmt == NULL){
      printf("ERROR: cannot determine FSLOUTPUTTYPE \n");
      exit(1);
    }
  }
  printf("Using %s as FSL/FEAT format extension\n",fmt);

  sprintf(tmpstr,"%s/mask",featdir);
  printf("Mask from %s\n",tmpstr);
  MaskVolFile = IDnameFromStem(tmpstr);
  if(MaskVolFile == NULL) exit(1);
  mask = MRIread(MaskVolFile);
  if (mask==NULL) exit(1);

  printf("Processing z stats\n");
  nthzstat = 1;
  while(1){
    sprintf(tmpstr,"%s/stats/zstat%d",featdir,nthzstat);
    ZVolFile = IDnameFromStem(tmpstr);
    if(ZVolFile == NULL) break;

    printf("%2d Reading z from %s\n",nthzstat,ZVolFile);
    z = MRIread(ZVolFile);
    if (z==NULL) exit(1);

    TwoSided = 1; // Use 2-sided for t
    p   = RFz2p(z, mask, TwoSided, NULL);
    sig = MRIlog10(p,mask,sig,1);
    sprintf(tmpstr,"%s/stats/zsig%d.%s",featdir,nthzstat,fmt);
    printf("Writing sig to %s\n",tmpstr);
    MRIwrite(sig,tmpstr);

    MRIfree(&z);
    MRIfree(&p);
    MRIfree(&sig);
    free(ZVolFile);
    nthzstat ++;
  }

  printf("Processing zf stats\n");
  nthzstat = 1;
  while(1){
    sprintf(tmpstr,"%s/stats/zfstat%d",featdir,nthzstat);
    ZVolFile = IDnameFromStem(tmpstr);
    if(ZVolFile == NULL) break;

    printf("%2d Reading z from %s\n",nthzstat,ZVolFile);
    z = MRIread(ZVolFile);
    if (z==NULL) exit(1);

    TwoSided = 0; // Use 1-sided for F
    p   = RFz2p(z, mask, TwoSided, NULL);
    sig = MRIlog10(p,mask,sig,1);
    sprintf(tmpstr,"%s/stats/zfsig%d.%s",featdir,nthzstat,fmt);
    printf("Writing sig to %s\n",tmpstr);
    MRIwrite(sig,tmpstr);

    MRIfree(&z);
    MRIfree(&p);
    MRIfree(&sig);
    free(ZVolFile);
    nthzstat ++;
  }

  printf("mri_z2p: done\n");
  exit(0);
}
/* ------ Doxygen markup starts on the line below ---- */
/*!
\fn int parse_commandline(int argc, char **argv)
\brief Parses the command-line arguments
\param argc - number of command line arguments
\param argv - pointer to a character pointer
*/
/* ------ Doxygen markup ends on the line above ---- */
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
    else if (!strcasecmp(option, "--two-sided"))   TwoSided = 1;
    else if (!strcasecmp(option, "--unsigned"))    TwoSided = 1;
    else if (!strcasecmp(option, "--one-sided"))   TwoSided = 0;
    else if (!strcasecmp(option, "--signed"))      TwoSided = 0;

    else if (!strcasecmp(option, "--z")) {
      if (nargc < 1) CMDargNErr(option,1);
      ZVolFile = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--p")) {
      if (nargc < 1) CMDargNErr(option,1);
      PVolFile = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--log10p")) {
      if (nargc < 1) CMDargNErr(option,1);
      Log10PVolFile = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--mask")) {
      if (nargc < 1) CMDargNErr(option,1);
      MaskVolFile = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--feat")) {
      if (nargc < 1) CMDargNErr(option,1);
      featdir = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--featfmt")) {
      if (nargc < 1) CMDargNErr(option,1);
      fmt = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--nii")) fmt = "nii";
    else if (!strcasecmp(option, "--nii.gz")) fmt = "nii.gz";
    else if (!strcasecmp(option, "--mgh")) fmt = "mgh";
    else if (!strcasecmp(option, "--mgz")) fmt = "mgz";
    else if (!strcasecmp(option, "--img")) fmt = "img";
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
  printf("   --z zvolfile : z volume \n");
  printf("   --p pvolfile : p volume \n");
  printf("   --log10p sigvolfile : sig volume \n");
  printf("\n");
  printf("   --mask maskfile : mask volume \n");
  printf("   --two-sided : (default) assume a two-sided, unsigned test (keeps sign of input)\n");
  printf("   --one-sided : assume a one-sided, signed test\n");
  printf("   --signed : two-sided/signed pvalue (p = 2*(1-p))\n");
  printf("\n");
  printf("   --feat featdir : convert all zstats and zfstats to sigs\n");
  printf("   --featfmt extension : use given format (eg, nii, nii.gz, mgh, etc)\n");
  printf("   --nii : use nii output format\n");
  printf("   --nii.gz : use nii.gz output format\n");
  printf("   --mgh : use mgh output format\n");
  printf("   --mgz : use mgz output format\n");
  printf("   --img : use img output format (analyze)\n");
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
  printf("Have not gotten to this yet:)\n");
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
  if(featdir == NULL){
    if (ZVolFile == NULL) {
      printf("ERROR: need zvol\n");
      exit(1);
    }
    if (PVolFile == NULL && Log10PVolFile == NULL) {
      printf("ERROR: need output vol (either --p or --log10p\n");
      exit(1);
    }
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
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  fprintf(fp,"twosided   %d\n",TwoSided);

  return;
}
