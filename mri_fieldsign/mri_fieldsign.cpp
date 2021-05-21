/**
 * @brief Computes retinotopic field sign
 *
 */
/*
 * Original Author: Douglas N. Greve
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
\file mri_fieldsign.c
\brief Computes retinotopy field sign
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
#include "fsglm.h"
#include "surfcluster.h"
#include "timer.h"


static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
MRI *SFA2MRI(MRI *eccen, MRI *polar, int SFATrue);

int main(int argc, char *argv[]) ;

const char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

char *FieldSignFile = NULL;
char *CMFFile = NULL;
char *NNbrFile = NULL;
char *RVarFile = NULL;

int DoSFA = 0;
char *EccenSFAFile=NULL,*PolarSFAFile=NULL;

int DoComplex = 0;
char *EccenRealFile=NULL, *EccenImagFile=NULL;
char *PolarRealFile=NULL, *PolarImagFile=NULL;
char *EccenOut=NULL, *PolarOut=NULL;

char *subject, *hemi, *SUBJECTS_DIR;
const char *PatchFile = NULL;
double fwhm = -1;
int nsmooth = -1;
char tmpstr[2000];
int ReverseSign = 0;
int SFATrue = 0;
int UseSphere = 0;
int usenew = 1;
double EccenRotAngle = 0;
double PolarRotAngle = 0;

int RETcompute_fieldsign2(MRIS *mris);

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs, err, reshapefactor, r, c, s;
  double v;
  MRIS *surf;
  MRI *eccensfa, *polarsfa, *mri, *mritmp;
  MRI *eccenreal,*eccenimag,*polarreal,*polarimag;

  nargs = handleVersionOption(argc, argv, "mri_fieldsign");
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

  sprintf(tmpstr,"%s/%s/label/%s.aparc.annot",SUBJECTS_DIR,subject,hemi);
  printf("Reading %s\n",tmpstr);
  err = MRISreadAnnotation(surf, tmpstr);
  if(err) exit(1);

  if(DoSFA){
    printf("SFA\n");
    eccensfa = MRIread(EccenSFAFile);
    if(eccensfa == NULL) exit(1);
    polarsfa = MRIread(PolarSFAFile);
    if(polarsfa == NULL) exit(1);
    mri = SFA2MRI(eccensfa, polarsfa, SFATrue);
    //MRIwrite(mri,"epc.mgh");
    MRIfree(&eccensfa);
    MRIfree(&polarsfa);
  }

  if(DoComplex){
    printf("Complex\n");
    eccenreal = MRIread(EccenRealFile);
    if(eccenreal == NULL) exit(1);
    eccenimag = MRIread(EccenImagFile);
    if(eccenimag == NULL) exit(1);
    polarreal = MRIread(PolarRealFile);
    if(polarreal == NULL) exit(1);
    polarimag = MRIread(PolarImagFile);
    if(polarimag == NULL) exit(1);

    mri = MRIallocSequence(eccenreal->width, eccenreal->height, 
			   eccenreal->depth, MRI_FLOAT, 4);
    for(c=0; c < mri->width; c++){
      for(r=0; r < mri->height; r++){
	for(s=0; s < mri->depth; s++){
	  v = MRIgetVoxVal(eccenreal,c,r,s,0);
	  MRIsetVoxVal(mri,c,r,s,0, v);
	  v = MRIgetVoxVal(eccenimag,c,r,s,0);
	  MRIsetVoxVal(mri,c,r,s,1, v);
	  v = MRIgetVoxVal(polarreal,c,r,s,0);
	  MRIsetVoxVal(mri,c,r,s,2, v);
	  v = MRIgetVoxVal(polarimag,c,r,s,0);
	  MRIsetVoxVal(mri,c,r,s,3, v);
	}
      }
    }
    MRIfree(&eccenreal);
    MRIfree(&eccenimag);
    MRIfree(&polarreal);
    MRIfree(&polarimag);
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
  err = MRISripZeros(surf,mri);
  if(err) exit(1);

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

  // Note: angle is also in 10th frame (f0=9) of SFA
  printf("Rot (rad): %g %g\n",EccenRotAngle,PolarRotAngle);
  RETcompute_angles(surf,EccenRotAngle,PolarRotAngle);
  if(EccenOut){
    mritmp = MRIcopyMRIS(NULL, surf, 0, "val");
    MRIwrite(mritmp,EccenOut);
    MRIfree(&mritmp);
  }
  if(PolarOut){
    mritmp = MRIcopyMRIS(NULL, surf, 0, "valbak");
    MRIwrite(mritmp,PolarOut);
    MRIfree(&mritmp);
  }

  if(usenew) RETcompute_fieldsign2(surf);
  else       RETcompute_fieldsign(surf);

  mritmp = MRIcopyMRIS(NULL, surf, 0, "fieldsign");
  if(ReverseSign){
    printf("Reversing sign\n");
    RETreverseSign(mritmp);
  }
  MRIwrite(mritmp,FieldSignFile);
  MRIfree(&mritmp);

  if(CMFFile){
    mritmp = MRIcopyMRIS(NULL, surf, 0, "stat");
    MRIwrite(mritmp,CMFFile);
    MRIfree(&mritmp);
  }

  if(NNbrFile){
    mritmp = MRIcopyMRIS(NULL, surf, 0, "K");
    MRIwrite(mritmp,NNbrFile);
    MRIfree(&mritmp);
  }
  if(RVarFile){
    mritmp = MRIcopyMRIS(NULL, surf, 0, "H");
    MRIwrite(mritmp,RVarFile);
    MRIfree(&mritmp);
  }

  printf("mri_fieldsign done\n");
  exit(0);
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
    else if (!strcasecmp(option, "--new")) usenew = 1;
    else if (!strcasecmp(option, "--old")) usenew = 0;

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
    } else if (!strcasecmp(option, "--sfa")) {
      if(nargc < 1) CMDargNErr(option,1);
      sprintf(tmpstr,"%s/eccen/h.nii",pargv[0]);
      EccenSFAFile = strcpyalloc(tmpstr);
      sprintf(tmpstr,"%s/polar/h.nii",pargv[0]);
      PolarSFAFile = strcpyalloc(tmpstr);
      if(!fio_FileExistsReadable(EccenSFAFile)){
	printf("ERROR: cannot find %s\n",EccenSFAFile);
	exit(1);
      }
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
      DoSFA = 0;
      nargsused = 2;
    } else if (!strcasecmp(option, "--polar")) {
      if (nargc < 2) CMDargNErr(option,2);
      PolarRealFile = pargv[0];
      PolarImagFile = pargv[1];
      if(!fio_FileExistsReadable(PolarRealFile)){
	printf("ERROR: cannot find %s\n",PolarRealFile);
	exit(1);
      }
      DoComplex = 1;
      DoSFA = 0;
      nargsused = 2;
    } else if (!strcasecmp(option, "--eccen-out")) {
      if(nargc < 1) CMDargNErr(option,1);
      EccenOut = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--polar-out")) {
      if(nargc < 1) CMDargNErr(option,1);
      PolarOut = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--eccen-rot")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&EccenRotAngle); //degrees
      EccenRotAngle *= M_PI/180.0;
      nargsused = 1;
    } else if (!strcasecmp(option, "--polar-rot")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&PolarRotAngle); //degrees
      PolarRotAngle *= M_PI/180.0;
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
    } else if (!strcmp(option, "--cmf")) {
      if (nargc < 1) CMDargNErr(option,1);
      CMFFile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--nnbr")) {
      if (nargc < 1) CMDargNErr(option,1);
      NNbrFile = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--rvar")) {
      if (nargc < 1) CMDargNErr(option,1);
      RVarFile = pargv[0];
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
  printf("   --fs fieldsignfile : output file\n");
  printf("   --eccen real imag\n");
  printf("   --polar real imag\n");
  printf("\n");
  printf("   --s subject \n");
  printf("   --hemi hemi \n");
  printf("   --patch patchfile : without hemi \n");
  printf("   --occip : patchfile = occip.patch.flat\n");
  printf("   --sphere : use spherical surface instead of patch\n");
  printf("   --fwhm fwhm_mm\n");
  printf("   --nsmooth nsmoothsteps\n");
  printf("   --rev : reverse sign\n");
  printf("   --old : use old FS estimation code (default is new)\n");
  printf("\n");
  printf("   --eccen-rot rotangle : rotate eccen by rotangle (degrees)\n");
  printf("   --polar-rot rotangle : rotate polar by rotangle (degrees)\n");
  printf("   --eccen-out eccenangle : output\n");
  printf("   --polar-out polarangle : output\n");
  printf("\n");
  printf("   --eccen-sfa sfafile : eccen selfreqavg file \n");
  printf("   --polar-sfa sfafile : polar selfreqavg file \n");
  printf("   --sfa sfadir :  \n");
  printf("   --sfa-true          : use true real and imag (only affects when smoothing)\n");
  printf("\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  std::cout << getVersion() << std::endl;
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
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/*--------------------------------------------------*/
static void check_options(void) 
{
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
  fprintf(fp,"%s\n", getVersion().c_str());
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
  fprintf(fp,"usenew      %d\n",usenew);
  return;
}
/*--------------------------------------------------------------------
  SFA2MRI(MRI *eccen, MRI *polar) - pack two SFAs int a single MRI. An
  SFA is the output of selfreqavg. Must be sampled to the surface. Note:
  the actual phase is in the 10th frame (f0=9) of the SFA.
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
	  v = MRIgetVoxVal(eccen,c,r,s,7); // eccen-real (8th)
	  MRIsetVoxVal(mri,c,r,s,0, v);
	  v = MRIgetVoxVal(eccen,c,r,s,8); // eccen-imag (9th)
	  MRIsetVoxVal(mri,c,r,s,1, v);
	  v = MRIgetVoxVal(polar,c,r,s,7); // polar-real (8th)
	  MRIsetVoxVal(mri,c,r,s,2, v);
	  v = MRIgetVoxVal(polar,c,r,s,8); // polar-imag (9th)
	  MRIsetVoxVal(mri,c,r,s,3, v);
	} else {
	  // Use real and imag weighted by log10(p)
	  // This corresponds to Marty's original code
	  v = MRIgetVoxVal(eccen,c,r,s,2); // eccen-real (3rd)
	  MRIsetVoxVal(mri,c,r,s,0, v);
	  v = MRIgetVoxVal(eccen,c,r,s,1); // eccen-image (2nd)
	  MRIsetVoxVal(mri,c,r,s,1, v);
	  v = MRIgetVoxVal(polar,c,r,s,2); // polar-real (3rd)
	  MRIsetVoxVal(mri,c,r,s,2, v);
	  v = MRIgetVoxVal(polar,c,r,s,1); // polar-imag (2nd)
	  MRIsetVoxVal(mri,c,r,s,3, v);
	}

      }
    }
  }

  return(mri);
}

int MRISextendedHopNeighbors(MRIS *surf,int TargVtxNo, int CurVtxNo,
			     int nhops, int *XNbrVtxNo, int *nXNbrs)
{
  static int nthhop = 0;
  VERTEX *vcur, *vnbr;
  int nNNbrs, n, NbrVtxNo;

  // Get the current vertex
  vcur = &surf->vertices[CurVtxNo];
  // Return if this vertex has been hit
  if( (int)vcur->val2bak == TargVtxNo ) return(0);
  // Return if this vertex is ripped
  if(vcur->ripflag) return(0);

  // Init the number of hops
  if(CurVtxNo == TargVtxNo){
    *nXNbrs = 0;
    nthhop = 0;
    XNbrVtxNo[*nXNbrs] = CurVtxNo;
    (*nXNbrs)++;
    vcur->val2bak = TargVtxNo; // record a hit
  }

  // Increment the number of hops
  nthhop++;

  // Stopping criteria
  if(nthhop > nhops) return(0);

  // Add nearest neighbors of current vertex
  nNNbrs = surf->vertices_topology[CurVtxNo].vnum;
  for (n = 0; n < nNNbrs; n++) {
    NbrVtxNo = surf->vertices_topology[CurVtxNo].v[n];
    vnbr = &surf->vertices[NbrVtxNo];
    if(vnbr->ripflag) continue;
    if(vnbr->val2bak == TargVtxNo) continue;
    // record a hit
    XNbrVtxNo[*nXNbrs] = NbrVtxNo;
    (*nXNbrs)++;
    vnbr->val2bak = TargVtxNo; 
  }

  // Now, loop over the current nearest neighbors
  nNNbrs = surf->vertices_topology[CurVtxNo].vnum;
  for (n = 0; n < nNNbrs; n++) {
    NbrVtxNo = surf->vertices_topology[CurVtxNo].v[n];
    MRISextendedHopNeighbors(surf, TargVtxNo, NbrVtxNo, nhops,
			     XNbrVtxNo, nXNbrs);
  }

  return(0);
}
/*------------------------------------------------------------*/
int RETcompute_fieldsign2(MRIS *mris)
{
  int k,n, knbr, shape;
  VERTEX *v, *vnbr;
  double dthresh = 1.5; //mm
  int *vtxlist, nlist=0;
  double *uctx, *vctx, d, deccen, dpolar, det=0;
  MATRIX *X, *y, *J;
  GLMMAT *glm;
  int msecTime, msecTot, nhit;
  Timer mytimer;
  int annot, annotindex, ok;
  FILE *fp;

  MRISremoveTriangleLinks(mris) ;
  printf("surfer: compute_fieldsign2()\n");

  J = MatrixAlloc(2,2,MATRIX_REAL);
  vtxlist = (int *) calloc(mris->nvertices,sizeof(int));
  uctx = (double *) calloc(mris->nvertices,sizeof(double));
  vctx = (double *) calloc(mris->nvertices,sizeof(double));

  if(mris->patch){
    // Assume flat map if a patch
    shape = 0;
    for(k=0;k<mris->nvertices;k++) {
      v = &(mris->vertices[k]);
      uctx[k] = v->x;
      vctx[k] = v->y;
    }
  }
  else{
    shape = SPHERICAL_COORDS;
    v = &(mris->vertices[0]);
    for(k=0;k<mris->nvertices;k++) {
      v = &(mris->vertices[k]);
      uctx[k] = atan2(v->y,v->x);
      d = sqrt(v->x * v->x + v->y * v->y);
      vctx[k] = atan2(d,v->z);
    }
  }

  printf("dthresh = %g\n",dthresh);
  printf("shape = %d\n",shape);
  printf("nvertices = %d\n",mris->nvertices);
  msecTot = 0;
  mytimer.reset() ;
  nhit = 0;
  fp = NULL;
  fp = fopen("tmp.dat","w");
  for (k=0; k < mris->nvertices; k++) {
    //printf("%d ------------------------\n",k);
    v = &(mris->vertices[k]);
    if(v->ripflag)      continue;
    if(v->val == 0)     continue;
    if(v->valbak == 0)  continue;
    if(v->val2 == 0)    continue;
    if(v->val2bak == 0) continue;

    if(0){
      // This can speed things up when testing
      annot = v->annotation;
      CTABfindAnnotation(mris->ct, annot, &annotindex);
      ok = 0;
      if(annotindex == 44) ok = 1;
      if(annotindex == 18) ok = 1;
      if(annotindex ==  5) ok = 1;
      if(annotindex == 42) ok = 1;
      if(annotindex == 53) ok = 1;
      if(annotindex == 66) ok = 1;
      if(annotindex == 14) ok = 1;
      if(annotindex == 16) ok = 1;
      if(!ok) continue;
    }

    nlist = sclustGrowByDist(mris, k, dthresh, shape, -1, vtxlist);
    X = MatrixAlloc(nlist,3,MATRIX_REAL);
    y = MatrixAlloc(nlist,2,MATRIX_REAL);
    //printf("nlist = %d\n",nlist);
    for(n = 0; n < nlist; n++){
      knbr = vtxlist[n];
      vnbr = &(mris->vertices[knbr]);
      deccen = RETcircsubtract(v->val,   vnbr->val);
      dpolar = RETcircsubtract(v->valbak,vnbr->valbak);
      y->rptr[n+1][1] = deccen;
      y->rptr[n+1][2] = dpolar;
      X->rptr[n+1][1] = uctx[k] - uctx[knbr];
      X->rptr[n+1][2] = vctx[k] - vctx[knbr];
      X->rptr[n+1][3] = 1;
      if(0){
	fprintf(fp,"%6d %6d  %g %g  %g %g %g %g\n",n,k,
		deccen,dpolar,
		uctx[k],uctx[knbr],
		vctx[k],uctx[knbr]);
      }
    }
    v->K = nlist; // size of neighborhood

    // Solve GLM
    glm = GLMalloc();
    glm->X = X;
    glm->y = y;
    GLMxMatrices(glm);
    GLMfit(glm);

    if(! glm->ill_cond_flag){
      v->H = glm->rvar; // residual variance

      // Fill Jacobian matrix
      J->rptr[1][1] = glm->beta->rptr[1][1]; // dR/dU
      J->rptr[1][2] = glm->beta->rptr[1][2]; // dR/dV
      J->rptr[2][1] = glm->beta->rptr[2][1]; // dTheta/dU
      J->rptr[2][2] = glm->beta->rptr[2][2]; // dTheta/dV
      
      if(0){
	printf("J -----------\n");
	MatrixPrint(stdout,J);
	printf("J -----------\n");
      }
      
      det = MatrixDeterminant(J);
      v->stat = 1000*fabs(det); // 1000 to make range a little easier
      if(det < 0) v->fieldsign = -1.0;
      else        v->fieldsign = +1.0;
      v->fsmask = sqrt(v->val2*v->val2bak);  /* geom mean of r,th power */
    }

    GLMfree(&glm); // Also frees X and y

    if(nhit%1000 == 0 || nhit == 0) {
      msecTime = mytimer.milliseconds() ;
      msecTot += msecTime;
      printf("%5d %5d %3d %6.2f  %g\n",k,nhit,nlist,msecTot/1000.0,det);
      mytimer.reset() ;
    }
    nhit ++;

    //if(nhit > 10000) break;

  } // loop over vertices
  msecTime = mytimer.milliseconds() ;
  msecTot += msecTime;
  printf("done: %5d %g\n",nhit,msecTot/1000.0);

  free(vtxlist);
  free(uctx);
  free(vctx);
  MatrixFree(&J);

  return(0);
}

