/**
 * @file  mri_fieldsign
 * @brief Computes retinotopic field sign
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Douglas N. Greve
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2007/09/28 22:08:17 $
 *    $Revision: 1.3 $
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



/*!
\file mri_fieldsign.c
\brief Computes retinotopy field sign
\author Douglas Greve

*/


// $Id: mri_fieldsign.c,v 1.3 2007/09/28 22:08:17 greve Exp $

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
double round(double x);
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <unistd.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "utils.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "mri.h"
#include "mri2.h"
#include "fio.h"
#include "version.h"
#include "cmdargs.h"
#include "fsenv.h"
#include "retinotopy.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
MRI *SFA2MRI(MRI *eccen, MRI *polar, int SFATrue);

int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_fieldsign.c,v 1.3 2007/09/28 22:08:17 greve Exp $";
char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

char *FieldSignFile = NULL;

int DoSFA = 0;
char *EccenSFAFile=NULL,*PolarSFAFile=NULL;

int DoComplex = 0;
char *EccenRealFile=NULL, *EccenImagFile=NULL;
char *PolarRealFile=NULL, *PolarImagFile=NULL;

char *subject, *hemi, *SUBJECTS_DIR;
char *PatchFile = NULL;
double fwhm = -1;
int nsmooth = -1;
char tmpstr[2000];
int ReverseSign = 0;
int SFATrue = 0;
int UseSphere = 0;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs, err, reshapefactor;
  MRIS *surf;
  MRI *eccensfa, *polarsfa, *mri, *mritmp, *mritmp2;
  MRI *eccenreal,*eccenimag,*polarreal,*polarimag;

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

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if (SUBJECTS_DIR == NULL) {
    printf("ERROR: SUBJECTS_DIR not defined in environment\n");
    exit(1);
  }

  // Load the surface
  sprintf(tmpstr,"%s/%s/surf/%s.sphere",SUBJECTS_DIR,subject,hemi);
  printf("Reading %s\n",tmpstr);
  surf = MRISread(tmpstr);
  if(!surf) exit(1);

  // Load the patch
  if(PatchFile){
    sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,subject,hemi,PatchFile);
    printf("Reading %s\n",tmpstr);
    err = MRISreadPatchNoRemove(surf, tmpstr) ;
    if(err) exit(1);
  } else {
    printf("Using spherical coordinates\n");
    MRISsphericalCoords(surf);
  }

  if(DoSFA){
    eccensfa = MRIread(EccenSFAFile);
    if(eccensfa == NULL) exit(1);
    polarsfa = MRIread(PolarSFAFile);
    if(polarsfa == NULL) exit(1);
    mri = SFA2MRI(eccensfa, polarsfa, SFATrue);
    MRIfree(&eccensfa);
    MRIfree(&polarsfa);
  }

  if(DoComplex){
    eccenreal = MRIread(EccenRealFile);
    if(eccenreal == NULL) exit(1);
    eccenimag = MRIread(EccenImagFile);
    if(eccenimag == NULL) exit(1);
    polarreal = MRIread(PolarRealFile);
    if(polarreal == NULL) exit(1);
    polarimag = MRIread(PolarImagFile);
    if(polarimag == NULL) exit(1);

    mritmp  = MRIconcatenateFrames(eccenreal, eccenimag, NULL);
    mritmp2 = MRIconcatenateFrames(polarreal, polarimag, NULL);
    mri = MRIconcatenateFrames(mritmp,mritmp2,NULL);

    MRIfree(&eccenreal);
    MRIfree(&eccenimag);
    MRIfree(&polarreal);
    MRIfree(&polarimag);
    MRIfree(&mritmp);
    MRIfree(&mritmp2);
  }

  if (mri->height != 1 || mri->depth != 1) {
    reshapefactor = mri->height * mri->depth;
    printf("Reshaping %d\n",reshapefactor);
    mritmp = mri_reshape(mri, reshapefactor*mri->width,
                           1, 1, mri->nframes);
    MRIfree(&mri);
    mri = mritmp;
    reshapefactor = 0; /* reset for output */
  }

  printf("Ripping Zeros\n");
  MRISripZeros(surf,mri);

  if(fwhm > 0) {
    nsmooth = MRISfwhm2nitersSubj(fwhm,subject,hemi,"white");
    if(nsmooth == -1) exit(1);
    printf("Approximating gaussian smoothing of target with fwhm = %lf,\n"
           "with %d iterations of nearest-neighbor smoothing\n",
           fwhm,nsmooth);
  }

  if(nsmooth > 0){
    printf("Smoothing %d steps\n",nsmooth);
    mritmp = MRISsmoothMRI(surf, mri, nsmooth, NULL, NULL);
    MRIfree(&mri);
    mri = mritmp;
  }

  MRIScopyMRI(surf, mri, 0, "val");    // eccen real
  MRIScopyMRI(surf, mri, 1, "val2");   // eccen imag
  MRIScopyMRI(surf, mri, 2, "valbak"); // polar real
  MRIScopyMRI(surf, mri, 3, "val2bak");// polar imag

  RETcompute_angles(surf);
  RETcompute_fieldsign(surf);

  mritmp = MRIcopyMRIS(NULL, surf, 0, "fieldsign");
  if(ReverseSign){
    printf("Reversing sign\n");
    MRIscalarMul(mritmp, mritmp, -1.0);
  }

  MRIwrite(mritmp,FieldSignFile);

  return 0;
}
/*---------------------------------------------------------*/
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;

  if (argc < 1) usage_exit();

  nargc = argc;
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
    else if (!strcasecmp(option, "--occip")) PatchFile = "occip.patch.flat";
    else if (!strcasecmp(option, "--rev")) ReverseSign = 1;
    else if (!strcasecmp(option, "--sfa-true")) SFATrue = 1;
    else if (!strcasecmp(option, "--sphere")) UseSphere = 1;

    else if (!strcasecmp(option, "--eccen-sfa")) {
      if (nargc < 1) CMDargNErr(option,1);
      EccenSFAFile = pargv[0];
      if(!fio_FileExistsReadable(EccenSFAFile)){
	printf("ERROR: cannot find %s\n",EccenSFAFile);
	exit(1);
      }
      DoSFA = 1;
      nargsused = 1;
    } else if (!strcasecmp(option, "--polar-sfa")) {
      if (nargc < 1) CMDargNErr(option,1);
      PolarSFAFile = pargv[0];
      if(!fio_FileExistsReadable(PolarSFAFile)){
	printf("ERROR: cannot find %s\n",PolarSFAFile);
	exit(1);
      }
      DoSFA = 1;
      nargsused = 1;
    } else if (!strcasecmp(option, "--eccen")) {
      if (nargc < 2) CMDargNErr(option,2);
      EccenRealFile = pargv[0];
      EccenImagFile = pargv[1];
      if(!fio_FileExistsReadable(EccenRealFile)){
	printf("ERROR: cannot find %s\n",EccenRealFile);
	exit(1);
      }
      DoComplex = 1;
      nargsused = 1;
    } else if (!strcasecmp(option, "--polar")) {
      if (nargc < 2) CMDargNErr(option,2);
      PolarRealFile = pargv[0];
      PolarImagFile = pargv[1];
      if(!fio_FileExistsReadable(PolarRealFile)){
	printf("ERROR: cannot find %s\n",PolarRealFile);
	exit(1);
      }
      DoComplex = 1;
      nargsused = 1;
    } else if (!strcasecmp(option, "--s")) {
      if (nargc < 1) CMDargNErr(option,1);
      subject = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--hemi")) {
      if (nargc < 1) CMDargNErr(option,1);
      hemi = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--patch")) {
      if (nargc < 1) CMDargNErr(option,1);
      PatchFile = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--fwhm")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&fwhm);
      nargsused = 1;
    } else if (!strcasecmp(option, "--nsmooth")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nsmooth);
      nargsused = 1;
    } else if (!strcmp(option, "--sd")) {
      if (nargc < 1) CMDargNErr(option,1);
      FSENVsetSUBJECTS_DIR(pargv[0]);
      nargsused = 1;
    } else if (!strcmp(option, "--fs")) {
      if (nargc < 1) CMDargNErr(option,1);
      FieldSignFile = pargv[0];
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
/*---------------------------------------------------------*/
static void usage_exit(void) {
print_usage() ;
  exit(1) ;
}
/*------------------------------------------------*/
static void print_usage(void) {
  printf("%s \n",Progname) ;
  printf("   --fs fieldsignfile : output\n");
  printf("\n");
  printf("   --eccen-sfa sfafile : eccen selfreqavg file \n");
  printf("   --polar-sfa sfafile : polar selfreqavg file \n");
  printf("   --sfa-true          : use true real and imag\n");
  printf("\n");
  printf("   --s subject \n");
  printf("   --hemi hemi \n");
  printf("   --patch patchfile : without hemi \n");
  printf("   --occip : patchfile = occip.patch.flat\n");
  printf("   --sphere : use spherical surface instead of patch\n");
  printf("\n");
  printf("   --fwhm fwhm_mm\n");
  printf("   --nsmooth nsmoothsteps\n");
  printf("\n");
  printf("   --rev : reverse sign\n");
  printf("\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  printf("%s\n", vcid) ;
  printf("\n");
}
/*--------------------------------------------------*/
static void print_help(void) {
  print_usage() ;
  printf("WARNING: this program is not yet tested!\n");
  exit(1) ;
}
/*--------------------------------------------------*/
static void print_version(void) {
  printf("%s\n", vcid) ;
  exit(1) ;
}
/*--------------------------------------------------*/
static void check_options(void) 
{
  if(!DoSFA){
    printf("Need SFA files\n");
    exit(1);
  }
  if(FieldSignFile == NULL){
    printf("Need output field sign file\n");
    exit(1);
  }
  if(subject == NULL){
    printf("Need subject\n");
    exit(1);
  }
  if(hemi == NULL){
    printf("Need hemi\n");
    exit(1);
  }
  if(fwhm > 0 && nsmooth > 0){
    printf("Cannot --fwhm and --nsmooth\n");
    exit(1);
  }
  if(PatchFile == NULL && ! UseSphere) {
    printf("ERROR: must spec --patch or --sphere\n");
    exit(1);
  }

  return;
}
/*--------------------------------------------------*/
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n",vcid);
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  if(DoSFA){
    fprintf(fp,"eccen-sfa  %s\n",EccenSFAFile);
    fprintf(fp,"polar-sfa  %s\n",PolarSFAFile);
    fprintf(fp,"sfa-true   %d\n",SFATrue);
  }
  if(DoComplex){
    fprintf(fp,"eccen real  %s\n",EccenRealFile);
    fprintf(fp,"eccen imag  %s\n",EccenImagFile);
    fprintf(fp,"polar real  %s\n",PolarRealFile);
    fprintf(fp,"polar imag  %s\n",PolarImagFile);
  }
  fprintf(fp,"patch     %s\n",PatchFile);
  fprintf(fp,"subject     %s\n",subject);
  fprintf(fp,"hemi        %s\n",hemi);
  fprintf(fp,"fwhm        %lf\n",fwhm);
  fprintf(fp,"nsmooth     %d\n",nsmooth);
  fprintf(fp,"fieldsign   %s\n",FieldSignFile);
  return;
}
/*--------------------------------------------------------------------
  SFA2MRI(MRI *eccen, MRI *polar) - pack two SFAs int a single MRI. An
  SFA is the output of selfreqavg. Must be sampled to the surface. 
  ----------------------------------------------------------------*/
MRI *SFA2MRI(MRI *eccen, MRI *polar, int SFATrue)
{
  MRI *mri;
  int c,r,s;
  double v;

  mri = MRIallocSequence(eccen->width, eccen->height, eccen->depth, 
			 MRI_FLOAT, 4);
  MRIcopyHeader(eccen,mri);

  for(c=0; c < eccen->width; c++){
    for(r=0; r < eccen->height; r++){
      for(s=0; s < eccen->depth; s++){

	if(SFATrue){
	  // Use pure real and imag
	  v = MRIgetVoxVal(eccen,c,r,s,7); // eccen-real
	  MRIsetVoxVal(mri,c,r,s,0, v);
	  v = MRIgetVoxVal(eccen,c,r,s,8); // eccen-imag
	  MRIsetVoxVal(mri,c,r,s,1, v);
	  v = MRIgetVoxVal(polar,c,r,s,7); // polar-real
	  MRIsetVoxVal(mri,c,r,s,2, v);
	  v = MRIgetVoxVal(polar,c,r,s,8); // polar-imag
	  MRIsetVoxVal(mri,c,r,s,3, v);
	} else {
	  // Use real and imag weighted by log10(p)
	  // This corresponds to Marty's original code
	  v = MRIgetVoxVal(eccen,c,r,s,1); // eccen-real
	  MRIsetVoxVal(mri,c,r,s,0, v);
	  v = MRIgetVoxVal(eccen,c,r,s,2); // eccen-imag
	  MRIsetVoxVal(mri,c,r,s,1, v);
	  v = MRIgetVoxVal(polar,c,r,s,1); // polar-real
	  MRIsetVoxVal(mri,c,r,s,2, v);
	  v = MRIgetVoxVal(polar,c,r,s,2); // polar-imag
	  MRIsetVoxVal(mri,c,r,s,3, v);
	}

      }
    }
  }

  return(mri);
}
