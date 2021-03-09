/**
 * @brief compute the % of gm, wm and CSF in each voxel in a volume
 *
 * Program to compute partial volume fractions for every voxel in e.g. an EPI image using
 * the aseg and the surfaces. Uses a high-resolution internal representation to compute the
 * cortical fractions.
 */
/*
 * Original Author: Bruce Fischl
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
#include "colortab.h"
#include "transform.h"
#include "mri2.h"
#include "mrisutils.h"
#include "cma.h"
#include "timer.h"
#include "fmriutils.h"
#include "registerio.h"


int main(int argc, char *argv[]) ;
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

char *TempVolFile=NULL;
char *subject, *subjectoverride=NULL, *SUBJECTS_DIR;
char *regfile=NULL;
char *outstem=NULL,*stackfile=NULL,*gmfile=NULL;
int USF=2;
const char *wsurf = "white";
const char *psurf = "pial";
const char *asegfile="aseg.mgz";
char *outsegfile=NULL,*ttsegfile=NULL,*ttsegctabfile=NULL;
char *csfmaskfile = NULL;
COLOR_TABLE *ct=NULL;
char tmpstr[5000];
const char *fmt = "mgz";
LTA *aseg2vol=NULL;
int regtype;
int nOptUnknown = 0;
int UseAseg = 1;
int FillCSF=1;
int nDil=3;
int RegHeader=0;
double resmm=0;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs,err,tt,nTT;
  MRI *aseg,*pvf,*mritmp,*csfmask;
  MRIS *lhw, *lhp, *rhw, *rhp;
  Timer start ;

  nargs = handleVersionOption(argc, argv, "mri_compute_volume_fractions");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd,2000);
  vg_isEqual_Threshold = 10e-4;

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  if (argc == 0) usage_exit();

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if (SUBJECTS_DIR == NULL) {
    printf("ERROR: SUBJECTS_DIR not defined in environment\n");
    exit(1);
  }

  parse_commandline(argc, argv);
  check_options();
  if (checkoptsonly) return(0);
  dump_options(stdout);

  start.reset() ;
  nTT = ct->ctabTissueType->nentries-1; // -1 to exclude unknown

  printf("Reading in aseg and surfs from %s/%s\n",SUBJECTS_DIR,subject);

  if(ttsegctabfile)
    CTABwriteFileASCII(ct->ctabTissueType, ttsegctabfile);

  if(UseAseg){
    sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,subject,asegfile);
    printf("Loading %s\n",tmpstr);
    aseg = MRIread(tmpstr);
    if(aseg==NULL) exit(1);
    if(aseg->type == MRI_UCHAR){
      printf("Changing type of seg from UCHAR to INT\n");
      mritmp = MRIchangeType(aseg, MRI_INT, 0, 1, 1);
      MRIfree(&aseg);
      aseg = mritmp;
    }
    if(FillCSF){
      printf("Filling empty voxels with extracerebral CSF (if not there already), nDil=%d\n",nDil);
      if(aseg->type == MRI_UCHAR){
	mritmp = MRIchangeType(aseg,MRI_INT,0,0,0);
	MRIfree(&aseg);
	aseg = mritmp;
      }
      aseg = MRIaddExtraCerebralCSF(aseg, nDil, aseg);
    }
    if(outsegfile){
      err = MRIwrite(aseg,outsegfile);
      if(err) exit(1);
    }
    if(ttsegfile){    
      MRI *ttseg;
      ttseg = MRIseg2TissueType(aseg, ct, NULL);
      if(ttseg == NULL) exit(1);
      err = MRIwrite(ttseg,ttsegfile);
      if(err) exit(1);
      MRIfree(&ttseg);
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

  if(RegHeader){
    MRI *orig;
    MRI *temp;
    sprintf(tmpstr,"%s/%s/mri/orig.mgz",SUBJECTS_DIR,subject);
    orig = MRIreadHeader(tmpstr,MRI_VOLUME_TYPE_UNKNOWN);
    if(orig==0) exit(1);
    temp = MRIreadHeader(TempVolFile,MRI_VOLUME_TYPE_UNKNOWN);
    if(temp==0) exit(1);
    aseg2vol = TransformRegDat2LTA(orig, temp, NULL);
    printf("Computing registration from header\n");
    MatrixPrint(stdout,aseg2vol->xforms[0].m_L);
    MRIfree(&orig);
    MRIfree(&temp);
  }

  printf("  t = %g\n",start.seconds()) ;
  if(resmm == 0){
    if(aseg) resmm = aseg->xsize/(USF);
    else     resmm = 1.0/(USF);
  }
  printf("Computing PVF (USF=%d, resmm=%lf)\n",USF,resmm);
  pvf = MRIpartialVolumeFractionAS(aseg2vol,aseg,lhw,lhp,rhw,rhp,USF,resmm,ct,NULL);
  if(pvf == NULL) exit(1);

  printf("  t = %g\n",start.seconds()) ;

  if(csfmaskfile){
    printf("Masking CSF with %s\n",csfmaskfile);
    sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,subject,csfmaskfile);
    csfmask = MRIread(tmpstr);
    if(csfmask==NULL) exit(1);
    MRIbinarize(csfmask, csfmask, 0.5, 0, 1);
    MRImask(pvf, csfmask, pvf, 0, 0);
  }
  printf("  t = %g\n",start.seconds()) ;

  printf("Writing results nTT=%d\n",nTT);

  if(outstem){
    for(tt=0; tt < nTT; tt++){
      sprintf(tmpstr,"%s.%s.%s",outstem,ct->ctabTissueType->entries[tt+1]->name,fmt);
      mritmp = fMRIframe(pvf,tt,NULL);
      err = MRIwrite(mritmp,tmpstr);
      if(err) exit(1);
      MRIfree(&mritmp);
    }
  }

  if(stackfile){
    err=MRIwrite(pvf,stackfile);
    if(err) exit(1);
  }

  if(gmfile){
    MRI *ctx, *subctxgm, *gm;
    ctx = fMRIframe(pvf,0,NULL);
    subctxgm = fMRIframe(pvf,1,NULL);
    gm = MRIadd(ctx,subctxgm,NULL);
    MRIcopyPulseParameters(pvf,gm);
    err=MRIwrite(gm,gmfile);
    if(err) exit(1);
    MRIfree(&ctx);
    MRIfree(&subctxgm);
    MRIfree(&gm);
  }

  printf("#@# CVF-Run-Time-Sec %g\n",start.seconds()) ;
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
    else if (!strcasecmp(option, "--regheader")) {
      if(nargc < 2) CMDargNErr(option,2);
      subject = pargv[0];
      TempVolFile = pargv[1];
      RegHeader = 1;
      nargsused = 2;
    } 
    else if(!strcasecmp(option, "--o") || !strcasecmp(option, "--pvf")) {
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
    else if (!strcasecmp(option, "--out-seg")) {
      if (nargc < 1) CMDargNErr(option,1);
      outsegfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--ttseg")) {
      if (nargc < 1) CMDargNErr(option,1);
      ttsegfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--ttseg-ctab")) {
      if (nargc < 1) CMDargNErr(option,1);
      ttsegctabfile = pargv[0];
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
      SUBJECTS_DIR = getenv("SUBJECTS_DIR");
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--usf")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&USF);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--resmm")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&resmm);
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
      USF = nint(1/resolution);
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--diag-no")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&Gdiag_no);
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--vg-thresh")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&vg_isEqual_Threshold);
      nargsused = 1;
    }
    else {
      // Make backwards compatible for flagless inputs
      if(nOptUnknown == 0){
	regfile = option;
	nargsused = 0;
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
  printf("   --o   outstem : output will be oustem.{cortex,subcort_gm,wm,csf}.mgz\n");
  printf("   --reg regfile : can be LTA or reg.dat (if reg.dat, then need template volume)\n");
  printf("   --regheader subject\n");
  printf("\n");
  printf("   --usf USF : set anatomical upsample factor (default %d)\n",USF);
  printf("   --r res   : resolution, sets USF = round(1/res)\n");
  printf("   --resmm resmm : set functional upsampling resolution (default is aseg->xsize/(USF))\n");
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
  printf("   --out-seg outseg : save seg (after adding xcsf voxels)\n");
  printf("   --ttseg ttseg : save tissue type segmentation (probably not that useful)\n");
  printf("   --ttseg-ctab ctab : save tissue type segmentation ctab (probably not that useful)\n");
  printf("\n");
  printf("   --mgz    : use mgz format (default)\n");
  printf("   --mgh    : use mgh format\n");
  printf("   --nii    : use nii format\n");
  printf("   --nii.gz : use nii.gz format\n");
  printf("   --ttype+head : use default+head instead of default tissue type info for seg\n");
  printf("\n");
  printf("   --vg-thresh thrshold : threshold for  'ERROR: LTAconcat(): LTAs 0 and 1 do not match'\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  std::cout << getVersion() << std::endl;
  printf("\n");
}
/*---------------------------------------------*/
static void print_help(void) {
  print_usage() ;
  printf("Computes partial volume fractions for cortex, subcortical GM, WM and CSF\n");
  printf("\n");
  printf("Example\n");
  printf("   mri_compute_volume_fractions --nii.gz --reg reg.lta --usf 3 --o pvf.\n");
  printf("This will create files called pvf.{cortex,subcort_gm,wm,csf}.nii.gz\n");
  printf("The value at each voxel will be between 0 and 1 and represent the fraction of \n");
  printf("the given tissue type\n");
  exit(1) ;
}
/*---------------------------------------------*/
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/*---------------------------------------------*/
static void check_options(void) {

  if(outstem == NULL && stackfile == NULL){
    printf("ERROR: no output specified\n");
    exit(1);
  }
  if(regfile == NULL && !RegHeader){
    printf("ERROR: no registration file specified\n");
    exit(1);
  }
  if(RegHeader){
    if(regfile){
      printf("ERROR: cannot spec both --regheader and registration file\n");
      exit(1);
    }
  }
  if((regtype == REGISTER_DAT || RegHeader) && TempVolFile==NULL){
    printf("ERROR: template volume needed with register.dat file or --reg-header\n");
    exit(1);
  }
  if (!RegHeader)
  {
    if(regtype == REGISTER_DAT){
      char *subjecttmp;
      subjecttmp = regio_read_subject(regfile);
      if(subjecttmp==NULL) exit(1);
      sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,subjecttmp,asegfile);
      aseg2vol = ltaReadRegisterDat(regfile, TempVolFile, tmpstr);
      if(aseg2vol == NULL){
	printf("ERROR: reading %s\n",regfile);
	printf("  maybe something wrong with %s\n",tmpstr);
	exit(1);
      }
      free(subjecttmp);
    }
    else aseg2vol = LTAread(regfile);
    if(debug) LTAprint(stdout,aseg2vol);
  }

  if(!RegHeader){
    if(subjectoverride == NULL) subject = aseg2vol->subject;
    else                        subject = subjectoverride;
  }
  if(ct==NULL) ct = TissueTypeSchema(NULL,"default-jan-2014");

  return;
}
/*---------------------------------------------*/
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n", getVersion().c_str());
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

