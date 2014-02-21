/**
 * @file  mri_compute_volume_fraction.c
 * @brief compute the % of gm, wm and CSF in each voxel in a volume
 *
 * Program to compute partial volume fractions for every voxel in e.g. an EPI image using
 * the aseg and the surfaces. Uses a high-resolution internal representation to compute the
 * cortical fractions.
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2014/02/21 18:47:36 $
 *    $Revision: 1.15 $
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
#include "colortab.h"
#include "transform.h"
#include "mri2.h"
#include "mrisutils.h"
#include "cma.h"
#include "timer.h"
#include "fmriutils.h"


int main(int argc, char *argv[]) ;
static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);

static char vcid[] = "$Id: mri_compute_volume_fractions.c,v 1.15 2014/02/21 18:47:36 greve Exp $";
char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

char *TempVolFile=NULL;
char *subject, *subjectoverride=NULL, *SUBJECTS_DIR;
char *regfile=NULL;
char *outstem=NULL,*stackfile=NULL,*gmfile=NULL;
int USF=2;
char *wsurf = "white";
char *psurf = "pial";
char *asegfile="aseg.mgz";
char *csfmaskfile = NULL;
COLOR_TABLE *ct=NULL;
char tmpstr[5000];
char *fmt = "mgz";
LTA *aseg2vol=NULL;
int regtype;
int nOptUnknown = 0;
int UseAseg = 1;
int FillCSF=1;
int nDil=3;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs,err,tt,nTT;
  MRI *aseg,**pvf,*mritmp,*stack,*gm,*csfmask;
  MRIS *lhw, *lhp, *rhw, *rhp;
  struct timeb start ;

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

  TimerStart(&start) ;
  nTT = ct->ctabTissueType->nentries-1; // -1 to exclude unknown

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if (SUBJECTS_DIR == NULL) {
    printf("ERROR: SUBJECTS_DIR not defined in environment\n");
    exit(1);
  }
  printf("Reading in aseg and surfs from %s/%s\n",SUBJECTS_DIR,subject);

  if(UseAseg){
    sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,subject,asegfile);
    printf("Loading %s\n",tmpstr);
    aseg = MRIread(tmpstr);
    if(aseg==NULL) exit(1);
    if(FillCSF){
      printf("Filling empty voxels with extracerebral CSF, nDil=%d\n",nDil);
      if(aseg->type == MRI_UCHAR){
	mritmp = MRIchangeType(aseg,MRI_INT,0,0,0);
	MRIfree(&aseg);
	aseg = mritmp;
      }
      aseg = MRIaddExtraCerebralCSF(aseg, nDil, aseg);
    }
  } else aseg=NULL;

  sprintf(tmpstr,"%s/%s/surf/lh.%s",SUBJECTS_DIR,subject,wsurf);
  lhw = MRISread(tmpstr);
  if(lhw==NULL) exit(1);

  sprintf(tmpstr,"%s/%s/surf/lh.%s",SUBJECTS_DIR,subject,psurf);
  lhp = MRISread(tmpstr);
  if(lhp==NULL) exit(1);

  sprintf(tmpstr,"%s/%s/surf/rh.%s",SUBJECTS_DIR,subject,wsurf);
  rhw = MRISread(tmpstr);
  if(rhw==NULL) exit(1);

  sprintf(tmpstr,"%s/%s/surf/rh.%s",SUBJECTS_DIR,subject,psurf);
  rhp = MRISread(tmpstr);
  if(rhp==NULL) exit(1);

  if(regtype == REGISTER_DAT){
    mritmp = MRIreadHeader(TempVolFile,MRI_VOLUME_TYPE_UNKNOWN);
    if(mritmp==NULL) exit(1);
    getVolGeom(mritmp, &aseg2vol->xforms[0].src);
    getVolGeom(aseg, &aseg2vol->xforms[0].dst);
    MRIfree(&mritmp);
    if(debug) LTAprint(stdout,aseg2vol);
  }

  printf("  t = %g\n",TimerStop(&start)/1000.0) ;
  printf("Computing PVF (USF=%d)\n",USF);
  pvf = MRIpartialVolumeFractionAS(aseg2vol,aseg,lhw,lhp,rhw,rhp,USF,ct);
  if(pvf == NULL) exit(1);

  printf("  t = %g\n",TimerStop(&start)/1000.0) ;

  if(csfmaskfile){
    printf("Masking CSF with %s\n",csfmaskfile);
    sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,subject,csfmaskfile);
    csfmask = MRIread(tmpstr);
    if(csfmask==NULL) exit(1);
    MRIbinarize(csfmask, csfmask, 0.5, 0, 1);
    MRImask(pvf[3], csfmask, pvf[3], 0, 0);
  }
  printf("  t = %g\n",TimerStop(&start)/1000.0) ;

  printf("Writing results nTT=%d\n",nTT);

  for(tt=0; tt < nTT; tt++){
    sprintf(tmpstr,"%s.%s.%s",outstem,ct->ctabTissueType->entries[tt+1]->name,fmt);
    err = MRIwrite(pvf[tt],tmpstr);
    if(err) exit(1);
  }

  if(stackfile){
    stack = MRIallocSequence(pvf[0]->width,pvf[0]->height,pvf[0]->depth,MRI_FLOAT,nTT);
    MRIcopyHeader(pvf[0],stack);
    MRIcopyPulseParameters(pvf[0],stack);
    for(tt=0; tt < nTT; tt++) fMRIinsertFrame(pvf[tt], 0, stack, tt);
    err=MRIwrite(stack,stackfile);
    if(err) exit(1);
  }

  if(gmfile){
    gm = MRIadd(pvf[0],pvf[1],NULL);
    MRIcopyPulseParameters(pvf[0],gm);
    err=MRIwrite(gm,gmfile);
    if(err) exit(1);
  }

  printf("#@# CVF-Run-Time-Sec %g\n",TimerStop(&start)/1000.0) ;
  printf("mri_compute_volume_fraction done\n");
  exit(0);
}
/*---------------------------------------------*/
/*---------------------------------------------*/
/*---------------------------------------------*/
// Note: some of the command-line parsing is needed for backwards compatibility
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
    else if (!strcasecmp(option, "--no-aseg")) UseAseg = 0;
    else if (!strcasecmp(option, "--fill-csf")) FillCSF = 1;
    else if (!strcasecmp(option, "--no-fill-csf")) FillCSF = 0;
    else if (!strcasecmp(option, "--ttype+head")) 
      ct = TissueTypeSchema(NULL,"default-jan-2014+head");

    else if (!strcasecmp(option, "--reg")) {
      if(nargc < 1) CMDargNErr(option,1);
      regfile = pargv[0];
      nargsused = 1;
      aseg2vol = LTAread(regfile);
      if(aseg2vol == NULL) exit(1);
      regtype = TransformFileNameType(regfile);
      if(regtype == REGISTER_DAT){
	if(nargc < 2 || CMDisFlag(pargv[1])){
	  printf("ERROR: registration file %s is not an LTA so --reg requires a target template volume\n",regfile);
	  exit(1);
	}
	TempVolFile = pargv[1];
	nargsused++;
      }
    } 
    else if (!strcasecmp(option, "--o")) {
      if (nargc < 1) CMDargNErr(option,1);
      outstem = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--stack")) {
      if (nargc < 1) CMDargNErr(option,1);
      stackfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--gm")) {
      if (nargc < 1) CMDargNErr(option,1);
      gmfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--csf-mask")) {
      if (nargc < 1) CMDargNErr(option,1);
      csfmaskfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--seg")) {
      if (nargc < 1) CMDargNErr(option,1);
      asegfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--wsurf")) {
      if (nargc < 1) CMDargNErr(option,1);
      wsurf = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--psurf")) {
      if (nargc < 1) CMDargNErr(option,1);
      psurf = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "-s") || !strcasecmp(option, "--s")) {
      if (nargc < 1) CMDargNErr(option,1);
      subjectoverride = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "-nii.gz")) fmt = "nii.gz";
    else if (!strcasecmp(option, "-nii")) fmt = "nii";
    else if (!strcasecmp(option, "-mgh")) fmt = "mgh";
    else if (!strcasecmp(option, "-mgz")) fmt = "mgz";
    else if (!strcasecmp(option, "--nii.gz")) fmt = "nii.gz";
    else if (!strcasecmp(option, "--nii")) fmt = "nii";
    else if (!strcasecmp(option, "--mgh")) fmt = "mgh";
    else if (!strcasecmp(option, "--mgz")) fmt = "mgz";
    else if (!strcmp(option, "--sd") || !strcmp(option, "-SDIR")) {
      if(nargc < 1) CMDargNErr(option,1);
      setenv("SUBJECTS_DIR",pargv[0],1);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--usf")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&USF);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--ndil")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nDil);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "-r") || !strcasecmp(option, "--r")) {
      double resolution;
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&resolution);
      USF = round(1/resolution);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--diag-no")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&Gdiag_no);
      nargsused = 1;
    } 
    else {
      // Make backwards compatible for flagless inputs
      if(nOptUnknown == 0){
	regfile = option;
	nargsused = 0;
	aseg2vol = LTAread(regfile);
	if(aseg2vol == NULL) exit(1);
	regtype = TransformFileNameType(regfile);
	printf("Setting registration file to %s\n",regfile);
      }
      else if(nOptUnknown == 1){
	TempVolFile = option;
	nargsused = 0;
	printf("Setting template volume to %s\n",TempVolFile);
      }
      else if(nOptUnknown == 2){
	outstem = option;
	nargsused = 0;
	printf("Setting outstem to %s\n",outstem);
      }
      else {
	fprintf(stderr,"ERROR: Option %s unknown\n",option);
	if (CMDsingleDash(option))
	  fprintf(stderr,"       Did you really mean -%s ?\n",option);
	exit(-1);
      }
      nOptUnknown++;
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/*---------------------------------------------*/

static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/*---------------------------------------------*/
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --reg regfile : can be LTA or reg.dat\n");
  printf("   --o   outstem : output will be oustem.{cortex,subcort_gm,wm,csf}.mgz\n");
  printf("\n");
  printf("   --usf USF : upsample factor (default %d)\n",USF);
  printf("   --r res   : resolution, sets USF = round(1/res)\n");
  printf("   --seg  segfile : use segfile instead of %s\n",asegfile);
  printf("   --wsurf wsurf : white surface (default is %s)\n",wsurf);
  printf("   --psurf psurf : pial surface (default is %s)\n",psurf);
  printf("   --no-aseg : do not include aseg (good for testing)\n");
  printf("   --stack stackfile : put ctx,subcortgm,wm,csf into a single multi-frame file\n");
  printf("   --gm gmfile : put ctx+subcortgm into a single-frame file\n");
  printf("   --no-fill-csf : do not attempt to fill voxels surrounding seg with the extracerebral CSF segmetation \n");
  printf("     Note: when the fill is done, there is no attempt to actually segment xCSF voxels.\n");
  printf("     The passed segmentation is dilated and the new voxels become xCSF\n");
  printf("     Note: if the passed seg already has the CSF_ExtraCerebral seg, nothing will be done\n");
  printf("   --dil N : for xCSF fill, dilate by N (default is %d); use -1 to fill the entire volume \n",nDil);
  printf("\n");
  printf("   --mgz    : use mgz format (default)\n");
  printf("   --mgh    : use mgh format\n");
  printf("   --nii    : use nii format\n");
  printf("   --nii.gz : use nii.gz format\n");
  printf("   --ttype+head : use default+head instead of default tissue type info for seg\n");
  printf("\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  printf("%s\n", vcid) ;
  printf("\n");
}
/*---------------------------------------------*/
static void print_help(void) {
  print_usage() ;
  printf("Computes partial volume fractions for cortex, subcortical GM, WM and CSF\n");
  exit(1) ;
}
/*---------------------------------------------*/
static void print_version(void) {
  printf("%s\n", vcid) ;
  exit(1) ;
}
/*---------------------------------------------*/
static void check_options(void) {

  if(regfile == NULL){
    printf("ERROR: no registration file specified\n");
    exit(1);
  }
  if(outstem == NULL){
    printf("ERROR: no output stem specified\n");
    exit(1);
  }
  if(subjectoverride == NULL) subject = aseg2vol->subject;
  else                        subject = subjectoverride;
  if(regtype == REGISTER_DAT && TempVolFile==NULL){
    printf("ERROR: template volume needed with register.dat file\n");
    exit(1);
  }
  if(ct==NULL) ct = TissueTypeSchema(NULL,"default-jan-2014");

  return;
}
/*---------------------------------------------*/
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n",vcid);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  fprintf(fp,"setenv SUBJECTS_DIR %s\n",SUBJECTS_DIR);
  fprintf(fp,"cd %s\n",cwd);
  fprintf(fp,"%s\n",cmdline);
  fprintf(fp,"outstem %s\n",outstem);
  fprintf(fp,"regfile %s\n",regfile);
  fprintf(fp,"regtype %d\n",regtype);
  fprintf(fp,"segfile %s\n",asegfile);
  fprintf(fp,"wsurf %s\n",wsurf);
  fprintf(fp,"psurf %s\n",psurf);
  if(TempVolFile) fprintf(fp,"TempVolFile %s\n",regfile);
  fprintf(fp,"USF %d\n",USF);
  if(subjectoverride) fprintf(fp,"subjectoverride %s\n",subjectoverride );
  return;
}

