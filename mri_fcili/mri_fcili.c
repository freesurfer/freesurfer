/**
 * @file  mri_fcili.c
 * @brief Computes intrinsic laterality index (iLI) based on supplied waveforms.
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Douglas N. Greve
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2013/11/22 19:41:44 $
 *    $Revision: 1.3 $
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


// $Id: mri_fcili.c,v 1.3 2013/11/22 19:41:44 greve Exp $

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
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#include "macros.h"
#include "utils.h"
#include "fio.h"
#include "version.h"
#include "cmdargs.h"
#include "error.h"
#include "diag.h"
#include "mri.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;
MRI *MRIfcIntrinsicLI(MRI *lh, MRI *rh, double DenThresh);

static char vcid[] = "$Id: mri_fcili.c,v 1.3 2013/11/22 19:41:44 greve Exp $";
const char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;
char tmpstr[4000];

char *outdir=NULL;
char *lhfile=NULL, *rhfile=NULL;
char *outfmt = "nii.gz";
double DenThresh = 0.2;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs,err;
  MRI *lh, *rh;
  MRI *lhn, *rhn, *iLI, *vLL, *vLR, *vRL, *vRR;
  int c,f, nrois, roi1, roi2, nframes;
  double v,ss, v1,v2,den,num, LL, LR, RL, RR;

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

  lh = MRIread(lhfile);
  if(lh == NULL) exit(1);
  rh = MRIread(rhfile);
  if(rh == NULL) exit(1);

  if(lh->nframes != rh->nframes){
    printf("ERROR: mri_fcili: frame mismatch\n");
    exit(1);
  }
  nframes = lh->nframes;

  if(lh->width != rh->width){
    printf("ERROR: mri_fcili: roi mismatch\n");
    exit(1);
  }
  nrois = lh->width;

  printf("nrois = %d\n",nrois);
  printf("nframes = %d\n",nframes);
  fflush(stdout);

  err = mkdir(outdir,0777);
  if (err != 0 && errno != EEXIST) {
    printf("ERROR: creating directory %s\n",outdir);
    perror(NULL);
    exit(1);
  }

  printf("Allocating\n");fflush(stdout);
  iLI = MRIalloc(nrois,nrois,1,MRI_FLOAT);
  if(iLI == NULL){
    printf("ERROR: mri_fcili: could not alloc %d\n",nrois);
    exit(1);
  }
  vLL = MRIalloc(nrois,nrois,1,MRI_FLOAT);
  if(vLL == NULL){
    printf("ERROR: mri_fcili: could not alloc vLL %d\n",nrois);
    exit(1);
  }
  vLR = MRIalloc(nrois,nrois,1,MRI_FLOAT);
  if(vLR == NULL){
    printf("ERROR: mri_fcili: could not alloc vLR %d\n",nrois);
    exit(1);
  }
  vRL = MRIalloc(nrois,nrois,1,MRI_FLOAT);
  if(vRL == NULL){
    printf("ERROR: mri_fcili: could not alloc vRL %d\n",nrois);
    exit(1);
  }
  vRR = MRIalloc(nrois,nrois,1,MRI_FLOAT);
  if(vRR == NULL){
    printf("ERROR: mri_fcili: could not alloc vRR %d\n",nrois);
    exit(1);
  }

  // Temporal normalization
  printf("Normalizing\n");fflush(stdout);
  lhn = MRIcopy(lh,NULL);
  rhn = MRIcopy(rh,NULL);
  for(c=0; c < nrois; c++){
    ss = 0;
    for(f=0; f < nframes; f++){
      v = MRIgetVoxVal(lh,c,0,0,f);
      ss += (v*v);
    }
    ss = sqrt(ss);
    for(f=0; f < nframes; f++){
      v = MRIgetVoxVal(lh,c,0,0,f);
      MRIsetVoxVal(lhn,c,0,0,f, v/ss);
    }
    ss = 0;
    for(f=0; f < nframes; f++){
      v = MRIgetVoxVal(rh,c,0,0,f);
      ss += (v*v);
    }
    ss = sqrt(ss);
    for(f=0; f < nframes; f++){
      v = MRIgetVoxVal(rh,c,0,0,f);
      MRIsetVoxVal(rhn,c,0,0,f, v/ss);
    }
  }
  //MRIwrite(lhn,"lhn.nii");
  //MRIwrite(rhn,"rhn.nii");

  printf("Computing iLI\n");fflush(stdout);
  for(roi1 = 0; roi1 < nrois; roi1++){
    for(roi2 = 0; roi2 < nrois; roi2++){

      LL = 0;
      for(f=0; f < nframes; f++){
	v1 = MRIgetVoxVal(lhn,roi1,0,0,f);
	v2 = MRIgetVoxVal(lhn,roi2,0,0,f);
	LL += (v1*v2);
      }
      LR = 0;
      for(f=0; f < nframes; f++){
	v1 = MRIgetVoxVal(lhn,roi1,0,0,f);
	v2 = MRIgetVoxVal(rhn,roi2,0,0,f);
	LR += (v1*v2);
      }
      RL = 0; // not the same as LR (matrices are transposes)
      for(f=0; f < nframes; f++){
	v1 = MRIgetVoxVal(rhn,roi1,0,0,f);
	v2 = MRIgetVoxVal(lhn,roi2,0,0,f);
	RL += (v1*v2);
      }
      RR = 0;
      for(f=0; f < nframes; f++){
	v1 = MRIgetVoxVal(rhn,roi1,0,0,f);
	v2 = MRIgetVoxVal(rhn,roi2,0,0,f);
	RR += (v1*v2);
      }

      num = ((LL-RL)-(RR-LR));
      den = (fabs(LL)+fabs(LR)+fabs(RR)+fabs(RL));
      if(den>DenThresh) v = num/den;
      else              v = 0.0;
      MRIsetVoxVal(iLI,roi1,roi2,0,0, v);
      MRIsetVoxVal(vLL,roi1,roi2,0,0, LL);
      MRIsetVoxVal(vLR,roi1,roi2,0,0, LR);
      MRIsetVoxVal(vRL,roi1,roi2,0,0, RL);
      MRIsetVoxVal(vRR,roi1,roi2,0,0, RR);

    } // roi2
  } // roi1

  MRIfree(&lhn);
  MRIfree(&rhn);

  printf("Saving ...\n");fflush(stdout);
  sprintf(tmpstr,"%s/iLI.%s",outdir,outfmt);
  printf("writing iLI to %s\n",tmpstr);
  err = MRIwrite(iLI,tmpstr);
  if(err) exit(1);

  sprintf(tmpstr,"%s/LL.%s",outdir,outfmt);
  printf("writing LL to %s\n",tmpstr);
  err = MRIwrite(vLL,tmpstr);
  if(err) exit(1);

  sprintf(tmpstr,"%s/LR.%s",outdir,outfmt);
  printf("writing LR to %s\n",tmpstr);
  err = MRIwrite(vLR,tmpstr);
  if(err) exit(1);

  sprintf(tmpstr,"%s/RL.%s",outdir,outfmt);
  printf("writing RL to %s\n",tmpstr);
  err = MRIwrite(vRL,tmpstr);
  if(err) exit(1);

  sprintf(tmpstr,"%s/RR.%s",outdir,outfmt);
  printf("writing RR to %s\n",tmpstr);
  err = MRIwrite(vRR,tmpstr);
  if(err) exit(1);

  printf("mri_fcili done\n");
  exit(0);
}
/* ------------------------------------------------*/
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
      outdir = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--fmt")) {
      if (nargc < 1) CMDargNErr(option,1);
      outfmt = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--den-thresh")) {
      if (nargc < 1) CMDargNErr(option,1);
      scanf(pargv[0],"%ld",&DenThresh);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--lh")) {
      if (nargc < 1) CMDargNErr(option,1);
      lhfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--rh")) {
      if (nargc < 1) CMDargNErr(option,1);
      rhfile = pargv[0];
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
/*-----------------------------------------------------------*/
static void check_options(void) {
  if(outdir == NULL){
    printf("ERROR: must spec output dir\n");
    exit(1);
  }
  if(lhfile == NULL){
    printf("ERROR: must spec lh time course file\n");
    exit(1);
  }
  if(rhfile == NULL){
    printf("ERROR: must spec rh time course file\n");
    exit(1);
  }
  return;
}
/* ------------------------------------------------*/
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/*-----------------------------------------------------------*/
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --o outputdir \n");
  printf("   --lh lhtimecourses (nii, mgh, etc) \n");
  printf("   --rh rhtimecourses (nii, mgh, etc) \n");
  printf("   --den-thresh thesh : denominator threshold (%g) \n",DenThresh);
  printf("   --fmt extension : output format (%s) \n",outfmt);
  printf("\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  printf("%s\n", vcid) ;
  printf("\n");
}
/*-----------------------------------------------------------*/
static void print_help(void) {
  print_usage() ;
  printf("This program computes the intrinsic laterality index (iLI)\n");
  printf("based on left hemisphere and right hemisphere registered\n");
  printf("waveforms. Based on Liu, et al, PNAS, 2009, vol 106\n");
  printf("Evidence from intrinsic activity that asymmetry of the\n");
  printf("human brain is controlled by multiple factors.\n");
  printf("\n");
  exit(1) ;
}
/*-----------------------------------------------------------*/
static void print_version(void) {
  printf("%s\n", vcid) ;
  exit(1) ;
}
/*-----------------------------------------------------------*/
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n",vcid);
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  fprintf(fp,"lhfile   %s\n",lhfile);
  fprintf(fp,"rhfile   %s\n",rhfile);
  fprintf(fp,"outdir   %s\n",outdir);
  fprintf(fp,"outfmt   %s\n",outfmt);
  fprintf(fp,"DenThresh  %g\n",DenThresh);

  return;
}

/*-----------------------------------------------------------*/

MRI *MRIfcIntrinsicLI(MRI *lh, MRI *rh, double DenThresh)
{
  MRI *lhn, *rhn, *iLI;
  int c,f, nrois, roi1, roi2, nframes;
  double v,ss, v1,v2,den,num, LL, LR, RL, RR;

  if(lh->nframes != rh->nframes){
    printf("ERROR: MRIfcIntrinsicLI(): frame mismatch\n");
    return(NULL);
  }
  nframes = lh->nframes;

  if(lh->width != rh->width){
    printf("ERROR: MRIfcIntrinsicLI(): roi mismatch\n");
    return(NULL);
  }
  nrois = lh->width;

  iLI = MRIalloc(nrois,nrois,1,MRI_FLOAT);
  if(iLI == NULL){
    printf("ERROR: MRIfcIntrinsicLI(): could not alloc %d\n",nrois);
    return(NULL);
  }

  lhn = MRIcopy(lh,NULL);
  rhn = MRIcopy(rh,NULL);

  for(c=0; c < nrois; c++){
    ss = 0;
    for(f=0; f < nframes; f++){
      v = MRIgetVoxVal(lh,c,0,0,f);
      ss += (v*v);
    }
    ss = sqrt(ss);
    for(f=0; f < nframes; f++){
      v = MRIgetVoxVal(lh,c,0,0,f);
      MRIsetVoxVal(lhn,c,0,0,f, v/ss);
    }
    ss = 0;
    for(f=0; f < nframes; f++){
      v = MRIgetVoxVal(rh,c,0,0,f);
      ss += (v*v);
    }
    ss = sqrt(ss);
    for(f=0; f < nframes; f++){
      v = MRIgetVoxVal(rh,c,0,0,f);
      MRIsetVoxVal(rhn,c,0,0,f, v/ss);
    }
  }
  //MRIwrite(lhn,"lhn.nii");
  //MRIwrite(rhn,"rhn.nii");


  for(roi1 = 0; roi1 < nrois; roi1++){
    for(roi2 = 0; roi2 < nrois; roi2++){

      LL = 0;
      for(f=0; f < nframes; f++){
	v1 = MRIgetVoxVal(lhn,roi1,0,0,f);
	v2 = MRIgetVoxVal(lhn,roi2,0,0,f);
	LL += (v1*v2);
      }
      LR = 0;
      for(f=0; f < nframes; f++){
	v1 = MRIgetVoxVal(lhn,roi1,0,0,f);
	v2 = MRIgetVoxVal(rhn,roi2,0,0,f);
	LR += (v1*v2);
      }
      RL = 0; // not the same as LR (matrices are transposes)
      for(f=0; f < nframes; f++){
	v1 = MRIgetVoxVal(rhn,roi1,0,0,f);
	v2 = MRIgetVoxVal(lhn,roi2,0,0,f);
	RL += (v1*v2);
      }
      RR = 0;
      for(f=0; f < nframes; f++){
	v1 = MRIgetVoxVal(rhn,roi1,0,0,f);
	v2 = MRIgetVoxVal(rhn,roi2,0,0,f);
	RR += (v1*v2);
      }

      num = ((LL-RL)-(RR-LR));
      den = (fabs(LL)+fabs(LR)+fabs(RR)+fabs(RL));
      if(den>DenThresh) v = num/den;
      else              v = 0.0;
      MRIsetVoxVal(iLI,roi1,roi2,0,0, v);

    } // roi2
  } // roi1

  MRIfree(&lhn);
  MRIfree(&rhn);
  return(iLI);
}

