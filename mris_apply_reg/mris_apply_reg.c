/**
 * @file  mris_apply_reg.c
 * @brief Applies multiple surface registrations
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Douglas N. Greve
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2012/10/04 18:01:53 $
 *    $Revision: 1.4 $
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
\file mris_apply_reg.c
\brief Example c file that can be used as a template.
\author Douglas Greve

*/


// $Id: mris_apply_reg.c,v 1.4 2012/10/04 18:01:53 greve Exp $

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
#include "mrisurf.h"
#include "resample.h"
#include "pdf.h"
#include "icosahedron.h"


static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
void usage_message(FILE *stream);
void usage(FILE *stream);
int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mris_apply_reg.c,v 1.4 2012/10/04 18:01:53 greve Exp $";
char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

char *SrcValFile=NULL;
char *TrgValFile=NULL;
char *SurfRegFile[100];
int ReverseMapFlag = 1;
int DoJac = 0;
int UseHash = 1;
int nsurfs = 0;
int DoSynthRand = 0;
int DoSynthOnes = 0;
int SynthSeed = -1;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs,n;
  MRIS *SurfReg[100];
  MRI *SrcVal, *TrgVal;
  char *base;

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

  printf("Loading %s\n",SrcValFile);
  SrcVal = MRIread(SrcValFile);
  if(SrcVal==NULL) exit(1);

  if(DoSynthRand) {
    if (SynthSeed < 0) SynthSeed = PDFtodSeed();
    printf("INFO: synthesizing, seed = %d\n",SynthSeed);
    srand48(SynthSeed);
    MRIrandn(SrcVal->width, SrcVal->height, SrcVal->depth,
             SrcVal->nframes,0, 1, SrcVal);
  }
  if(DoSynthOnes != 0) {
    printf("INFO: filling input with all 1s\n");
    MRIconst(SrcVal->width, SrcVal->height, SrcVal->depth,
             SrcVal->nframes, 1, SrcVal);
  }

  for(n=0; n<nsurfs;n++){
    printf("%d Loading %s\n",n+1,SurfRegFile[n]);
    base = fio_basename(SurfRegFile[n],".tri");
    if(strcmp(base,"ic7")==0){
      // Have to do it this way to rescale. Need to find a better more robust way.
      printf("   reading as ico 7, rescaling radius to 100\n");
      SurfReg[n] = ReadIcoByOrder(7, 100);
    }
    else
      SurfReg[n] = MRISread(SurfRegFile[n]);
    free(base);
    if(SurfReg[n]==NULL) exit(1);
  }

  TrgVal = MRISapplyReg(SrcVal, SurfReg, nsurfs,ReverseMapFlag,DoJac,UseHash);
  if(TrgVal == NULL) exit(1);

  printf("Writing %s\n",TrgValFile);
  MRIwrite(TrgVal,TrgValFile);
  
  printf("mris_apply_reg done\n");
  return 0;
}
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
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
    else if (!strcasecmp(option, "--norev")) ReverseMapFlag = 0;
    else if (!strcasecmp(option, "--no-rev")) ReverseMapFlag = 0;
    else if (!strcasecmp(option, "--nnf")) ReverseMapFlag = 0;
    else if (!strcasecmp(option, "--nnfr")) ReverseMapFlag = 1;
    else if (!strcasecmp(option, "--no-hash")) UseHash = 0;
    else if (!strcasecmp(option, "--jac")) DoJac = 1;
    else if (!strcasecmp(option, "--no-jac")) DoJac = 0;
    else if (!strcasecmp(option, "--randn")) DoSynthRand = 1;
    else if (!strcasecmp(option, "--ones")) DoSynthOnes = 1;

    else if (!strcasecmp(option, "--src") || !strcasecmp(option, "--sval")) {
      if (nargc < 1) CMDargNErr(option,1);
      SrcValFile = pargv[0];
      if(!fio_FileExistsReadable(SrcValFile)){
	printf("ERROR: %s does not exist or is not readable by you\n",SrcValFile);
	exit(1);
      }
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--trg") || !strcasecmp(option, "--tval") 
	     || !strcasecmp(option, "--o")) {
      if (nargc < 1) CMDargNErr(option,1);
      TrgValFile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--streg")) {
      if (nargc < 2) CMDargNErr(option,2);
      SurfRegFile[nsurfs] = pargv[0];
      nsurfs++;
      SurfRegFile[nsurfs] = pargv[1];
      nsurfs++;
      nargsused = 2;
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
/*--------------------------------------------------------------*/
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/*--------------------------------------------------------------*/
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --src srcvalfile : source values\n");
  printf("   --trg trgvalfile : output target values (--o)\n");
  printf("   --streg srcreg1 trgreg1 : source and target reg files\n");
  printf("     --streg srcreg2 trgreg2 : source and target reg files ...\n");
  printf("\n");
  printf("   --jac : use jacobian correction\n");
  printf("   --no-rev : do not do reverse mapping\n");
  printf("   --randn : replace input with WGN\n");
  printf("   --ones  : replace input with ones\n");
  printf("\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  printf("%s\n", vcid) ;
  printf("\n");
}
/*--------------------------------------------------------------*/
static void print_help(void) {
  print_usage() ;
  printf("WARNING: this program is not yet tested!\n");
  usage(stdout);
  exit(1) ;
}
/*--------------------------------------------------------------*/
static void print_version(void) {
  printf("%s\n", vcid) ;
  exit(1) ;
}
/*--------------------------------------------------------------*/
static void check_options(void) {
  int n;
  if(SrcValFile == NULL){
    printf("ERROR: need to specify source value file\n");
    exit(1);
  }
  if(TrgValFile == NULL){
    printf("ERROR: need to specify target value file\n");
    exit(1);
  }
  if(nsurfs == 0){
    printf("ERROR: must specify at least one source:target registration pair\n");
    exit(1);
  }
  for(n=0; n<nsurfs;n++){
    if(!fio_FileExistsReadable(SurfRegFile[n])){
      printf("ERROR: %s does not exist or is not readable by you\n",SurfRegFile[n]);
      exit(1);
    }
  }
  return;
}
/*--------------------------------------------------------------*/
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n",vcid);
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  fprintf(fp,"srcvalfile  %s\n",SrcValFile);
  fprintf(fp,"trgvalfile  %s\n",TrgValFile);
  fprintf(fp,"nsurfs  %d\n",nsurfs);
  fprintf(fp,"jac  %d\n",DoJac);
  fprintf(fp,"revmap  %d\n",ReverseMapFlag);
  return;
}

#include "mris_apply_reg.help.xml.h"
void usage(FILE *stream)
{
  outputHelpXml(mris_apply_reg_help_xml,mris_apply_reg_help_xml_len);
} /* end usage() */

/* EOF */
