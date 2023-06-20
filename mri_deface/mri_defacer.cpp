/**
 * @brief Part of a defacing algorithm (see the defacer script)
 *
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
#include "mri2.h"
#include "timer.h"
#include "fmriutils.h"
#include "cma.h"
#include "mrimorph.h"
#include "resample.h"
#include "numerics.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "cma.h"
#include "mri_identify.h"
#include "mris_sphshapepvf.h"
#include "mideface.h"

#ifdef _OPENMP
#include "romp_support.h"
#endif

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

/*!
  MRI *MRISpaintSphere(MRIS *surf, LABEL *label, MRI *img)
  This function is not used in this binary, this is just an ok place
  to park it until I get around to improving it.  The way it works is
  that you pass it a sphere with a label on the sphere. It then
  samples the img in the area of the label, scaling to fit.  It uses
  the phi/theta coords as an approx to a uniform grid.  It uses the
  tkreg xyz to compute phi/theta which makes it dependent upon the
  slicing of the volume used to create the sphere. Eg, the label may
  be at the pole, which would be bad. This was used to create the
  MIDEFACE label used as a watermark.
 */
MRI *MRISpaintSphere(MRIS *surf, LABEL *label, MRI *img)
{
  BasicSpherePVF bsph;

  double tmin=10e10, tmax=-10e10, pmin=10e10, pmax=-10e10;
  std::array<double,3> rtp;  
  for(int n=0; n < label->n_points; n++){
    int vtxno = label->lv[n].vno;
    VERTEX *v = &(surf->vertices[vtxno]);
    std::array<double,3> xyz = {v->x,v->z,v->y};
    rtp = bsph.XYZ2RTP(xyz);
    if(tmin > rtp[1]) tmin = rtp[1];
    if(tmax < rtp[1]) tmax = rtp[1];
    if(pmin > rtp[2]) pmin = rtp[2];
    if(pmax < rtp[2]) pmax = rtp[2];
  }
  printf("range %g %g  %g %g\n",tmin,tmax,pmin,pmax);
  MATRIX *vox2ras = MatrixIdentity(3,NULL);
  vox2ras->rptr[1][1] = (pmax-pmin)/img->width;
  vox2ras->rptr[2][2] = (tmax-tmin)/img->height;
  vox2ras->rptr[1][3] = pmin;
  vox2ras->rptr[2][3] = tmin;
  MatrixPrint(stdout,vox2ras);
  MATRIX *ras2vox = MatrixInverse(vox2ras,NULL);

  MATRIX *ras = MatrixAlloc(3,1,MATRIX_REAL);
  MATRIX *vox = NULL;
  ras->rptr[3][1] = 1;
  MRI *overlay = MRIallocSequence(surf->nvertices,1,1,MRI_FLOAT,img->nframes);
  int vtxnodb = -1;
  for(int n=0; n < label->n_points; n++){
    int vtxno = label->lv[n].vno;
    VERTEX *v = &(surf->vertices[vtxno]);
    std::array<double,3> xyz = {v->x,v->z,v->y};
    rtp = bsph.XYZ2RTP(xyz);
    ras->rptr[1][1] = rtp[2];
    ras->rptr[2][1] = rtp[1];
    vox = MatrixMultiply(ras2vox,ras,vox);
    int c = nint(vox->rptr[1][1]);
    int r = nint(vox->rptr[2][1]);
    if(vtxno == vtxnodb){
      printf("%d (%g,%g,%g) %g %g  %d %d\n",vtxno,v->x,v->y,v->z,rtp[2],rtp[1],c,r);
      fflush(stdout);
    }
    if(c < 0 || c >= img->width || r < 0 || r >= img->height){
      MRIsetVoxVal(overlay,vtxno,0,0,0,0);
      continue;
    }
    for(int f=0; f < img->nframes; f++){
      double val = MRIgetVoxVal(img,c,r,0,f);
      if(vtxno == vtxnodb){
	printf("%d %g %g  %d %d val=%g\n",vtxno,rtp[2],rtp[1],c,r,val);
	fflush(stdout);
      }
      MRIsetVoxVal(overlay,vtxno,0,0,f,val);
    }
  }

  return(overlay);
}

char *involpath=NULL, *headmaskpath=NULL, *tempsurfpath=NULL;
char *regpath=NULL,*xmaskpath=NULL;
char *templabelpathlist[100];
int ntemplabelpathlist=0;
char *outvolpath=NULL, *facesegpath=NULL, *minsurfpath=NULL, *maxsurfpath=NULL;
char *distdatpath=NULL, *distboundspath=NULL, *distoverlaypath=NULL, *statspath=NULL;
char *watermarkpath = NULL;
double dwatermark = 1;
MiDeface defacer;
double DistInMinList[200],DistInMaxList[200],DistInMin=2,DistInMax=20;
double DistInList[200];
char *outtempsurfpath=NULL;
int embed = 1;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int nargs,err=0;

  srand48(53);
  nargs = handleVersionOption(argc, argv, "mri_defacer");
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

  defacer.invol = MRIread(involpath);
  if(defacer.invol==NULL) exit(1);
  defacer.tempsurf = MRISread(tempsurfpath);
  if(defacer.tempsurf==NULL) exit(1);
  defacer.tempsurf->hemisphere = 3;
  defacer.headmask = MRIread(headmaskpath);
  if(defacer.headmask==NULL) exit(1);
  err = MRIdimMismatch(defacer.invol, defacer.headmask, 0);
  if(err){
    printf("ERROR: dimension mismatch with headmask\n");
    exit(1);
  }
  if(xmaskpath){
    defacer.xmask = MRIread(xmaskpath);
    if(defacer.xmask==NULL) exit(1);
    err = MRIdimMismatch(defacer.invol, defacer.xmask, 0);
    if(err){
      printf("ERROR: dimension mismatch with xmask\n");
      exit(1);
    }
  }

  if(watermarkpath){
    // Raise the template surface in the watermark
    LABEL *watermark = LabelRead("",watermarkpath);
    if(watermark == NULL) exit(1);
    printf("Applying watermark d=%g\n",dwatermark);
    defacer.watermark(watermark,dwatermark);
    LabelFree(&watermark);
  }

  if(regpath){
    // Generally not used
    printf("Applying reg %s to the template surface\n",regpath);
    LTA *lta = LTAread(regpath);
    if(lta==NULL) exit(1);
    int err = MRISltaMultiply(defacer.tempsurf, lta);
    if(err) exit(1);
    MRIScomputeMetricProperties(defacer.tempsurf);
    MRISsmoothSurfaceNormals(defacer.tempsurf,10);
  }

  // Create the output volume
  defacer.outvol = MRIallocSequence(defacer.invol->width, defacer.invol->height, 
    defacer.invol->depth, defacer.invol->type, defacer.invol->nframes);
  MRIcopyHeader(defacer.invol, defacer.outvol);
  MRIcopyPulseParameters(defacer.invol, defacer.outvol);
  // Allocate the min and max surfaces (only used for display)
  defacer.minsurf = MRISclone(defacer.tempsurf);
  defacer.maxsurf = MRISclone(defacer.tempsurf);

  defacer.SetDeltaDist();
  defacer.PrintParams(stdout);
  printf("\n");

  FILE *fpDLP = NULL;
  if(distboundspath) fpDLP = fopen(distboundspath,"w");

  for(int n=0; n < ntemplabelpathlist; n++){
    // Set the distance bounds for each label separately
    // Labels should be mutally exclusive
    if(defacer.templabel) LabelFree(&defacer.templabel);
    printf("===============================================\n");
    printf("Label %d %s %g %g\n",n,templabelpathlist[n],DistInMinList[n],DistInMaxList[n]);fflush(stdout);
    defacer.templabel = LabelRead("",templabelpathlist[n]);
    if(defacer.DoRipple) {
      // Apply ripple to the template surface in this label
      printf("Applying ripple %g %g\n",defacer.rippleamp,defacer.rippleperiod);
      defacer.ripple(defacer.templabel);
    }
    defacer.DistInMin = DistInMinList[n];
    defacer.DistInMax = DistInMaxList[n];
    if(defacer.templabel==NULL) exit(1);
    defacer.DistanceBounds();
    DistInList[n] = defacer.DistIn;
    printf("\n");
    if(distboundspath){
      // Write the bounds into a file
      fprintf(fpDLP,"%2d %5d %6.4f %6.4f %6.4f %6.4f\n",n+1,defacer.templabel->n_points,
	      defacer.DistInRaw,defacer.DistOutRaw,defacer.DistIn,defacer.DistOut);
      fflush(fpDLP);
    }
    if(distdatpath){
      // This saves the distance for each vertex in a text file, probably not useful
      char tmpstr[2000];
      sprintf(tmpstr,"%s.label%02d.dat",distdatpath,n+1);
      FILE *fp = fopen(tmpstr,"w");
      for(int nthp=0; nthp < defacer.templabel->n_points; nthp++){
	int vtxno = defacer.templabel->lv[nthp].vno;
	fprintf(fp,"%4d %6d %g %g\n",nthp,vtxno,defacer.DistInList[nthp],defacer.DistOutList[nthp]);
      }
      fclose(fp);
    }
    // Now segment the face mask for this label using the distance bounds found above
    defacer.SegFace();
    fflush(stdout);
  }
  if(distboundspath) fclose(fpDLP);
  printf("===============================================\n\n");

  // Compute stats in each compartment
  defacer.FaceIntensityStats();

  // Deface by filling in the face mask segmentatoin
  defacer.Deface();

  // Save some handy stats
  defacer.PrintStats(stdout);
  if(statspath){
    FILE *fp = fopen(statspath,"w");
    defacer.PrintStats(fp);
    for(int n=0; n < ntemplabelpathlist; n++)
      fprintf(fp,"DistIn.%d %g\n",n+1,DistInList[n]);
    fclose(fp);
  }

  // Write out a surface overlay of the distance at each vertex
  if(distoverlaypath){
    MRI *distoverlay = MRIcopyMRIS(NULL, defacer.tempsurf, 1, "valbak");
    distoverlay = MRIcopyMRIS(distoverlay, defacer.tempsurf, 0, "val");
    err = MRIwrite(distoverlay,distoverlaypath);
    if(err) exit(1);
  }

  if(embed){
    printf("Embedding code %s into output\n",defacer.embedded_code);
    defacer.EmbedCode(defacer.outvol,defacer.outvol);
  }

  err = MRIwrite(defacer.outvol,outvolpath);
  if(err) exit(1);

  if(facesegpath){
    err = MRIwrite(defacer.faceseg,facesegpath);
    if(err) exit(1);
  }
  if(minsurfpath){
    err = MRISwrite(defacer.minsurf,minsurfpath);
    if(err) exit(1);
  }
  if(maxsurfpath){
    err = MRISwrite(defacer.maxsurf,maxsurfpath);
    if(err) exit(1);
  }
  if(outtempsurfpath){
    err = MRISwrite(defacer.tempsurf,outtempsurfpath);
    if(err) exit(1);
  }

  printf("mri_defacer done\n");
  return(0);
  exit(0);
} // end of main
/*--------------------------------------------------------------------*/
/*---------------------------------------------------------------*/
/*---------------------------------------------------------------*/
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;

  if (argc < 1) usage_exit();

  nargc   = argc;
  pargv = argv;
  while (nargc > 0) {

    option = pargv[0];
    if(debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if(!strcasecmp(option, "--help"))  print_help() ;
    else if(!strcasecmp(option, "--version")) print_version() ;
    else if(!strcasecmp(option, "--debug"))   debug = 1;
    else if(!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if(!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;
    else if(!strcasecmp(option, "--no-ripple")) defacer.DoRipple = 0;
    else if(!strcasecmp(option, "--check-code")){
      // --check-code volume <outfile>
      if(nargc < 1) CMDargNErr(option,1);
      MRI *vol = MRIread(pargv[0]);
      if(vol==NULL) exit(1);
      int check = defacer.EmbedCodeCheck(vol);
      printf("%d\n",check);
      if(CMDnthIsArg(nargc, pargv, 1)){
	FILE *fp = fopen(pargv[1],"w");
	if(fp==NULL){
	  printf("ERROR: opening %s\n",pargv[1]);
	  exit(1);
	}
	fprintf(fp,"%d\n",check);
	fclose(fp);
      }
      exit(0);
    }
    else if(!strcasecmp(option, "--ripple")){
      if(nargc < 2) CMDargNErr(option,2);
      defacer.DoRipple = 1;
      sscanf(pargv[0],"%lf",&defacer.rippleamp);
      sscanf(pargv[1],"%lf",&defacer.rippleperiod);
      nargsused = 2;
    }
    else if(!strcasecmp(option, "--i")){
      if(nargc < 1) CMDargNErr(option,1);
      involpath = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--hm")){
      if(nargc < 1) CMDargNErr(option,1);
      headmaskpath = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--ts")){
      if(nargc < 1) CMDargNErr(option,1);
      tempsurfpath = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--ots")){
      if(nargc < 1) CMDargNErr(option,1);
      outtempsurfpath = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--ripple-center")){
      if(nargc < 3) CMDargNErr(option,3);
      for(int k=0; k<3;k++) sscanf(pargv[k],"%lf",&defacer.ripplecenter[k]);
      nargsused = 3;
    }
    else if(!strcasecmp(option, "--apply-ripple")){
      // stand-alone: surf axis amp period lable out
      // Change the center using --rippple-center before --apply-ripple
      // For atlas surf ripple, use center = (0,-50,-73),axis=1,rippleamp=2,rippleperiod=30
      if(nargc < 6) CMDargNErr(option,6);
      defacer.DoRipple = 1;
      defacer.tempsurf = MRISread(pargv[0]);
      sscanf(pargv[1],"%d",&defacer.rippleaxis);//1 or 2
      sscanf(pargv[2],"%lf",&defacer.rippleamp);
      sscanf(pargv[3],"%lf",&defacer.rippleperiod);
      LABEL *label = LabelRead("",pargv[4]);
      defacer.ripple(label);
      MRISwrite(defacer.tempsurf,pargv[5]);
      exit(0);
    }
    else if(!strcasecmp(option, "--l")){
      if(nargc < 1) CMDargNErr(option,1);
      templabelpathlist[ntemplabelpathlist] = pargv[0];
      DistInMinList[ntemplabelpathlist] = DistInMin;
      DistInMaxList[ntemplabelpathlist] = DistInMax;
      nargsused = 1;
      if(CMDnthIsArg(nargc, pargv, 1)) {
        sscanf(pargv[1],"%lf",&DistInMinList[ntemplabelpathlist]);
        nargsused ++;
	if(CMDnthIsArg(nargc, pargv, 2)) {
	  sscanf(pargv[2],"%lf",&DistInMaxList[ntemplabelpathlist]);
	  nargsused ++;
	} 
      } 
      ntemplabelpathlist++;
    }
    else if(!strcasecmp(option, "--w")){
      if(nargc < 2) CMDargNErr(option,2);
      watermarkpath = pargv[0];
      sscanf(pargv[1],"%lf",&dwatermark);
      nargsused = 2;
    }
    else if(!strcasecmp(option, "--reg")){
      if(nargc < 1) CMDargNErr(option,1);
      regpath = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--min")){
      if(nargc < 1) CMDargNErr(option,1);
      minsurfpath = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--max")){
      if(nargc < 1) CMDargNErr(option,1);
      maxsurfpath = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--fill-const")){
      if(nargc < 2) CMDargNErr(option,2);
      sscanf(pargv[0],"%lf",&defacer.FillConstIn);
      sscanf(pargv[1],"%lf",&defacer.FillConstOut);
      defacer.FillType=2;
      nargsused = 2;
    }
    else if(!strcasecmp(option, "--dist-in-frac")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&defacer.DistInFrac);
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--dist-in-min")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&DistInMin);
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--dist-in-max")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&DistInMax);
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--dist-out-frac")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&defacer.DistOutFrac);
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--dist-out-min")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&defacer.DistOutMin);
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--dist-out-max")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&defacer.DistOutMax);
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--distbounds")){
      if(nargc < 1) CMDargNErr(option,1);
      distboundspath = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--distdat")){
      if(nargc < 1) CMDargNErr(option,1);
      distdatpath = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--distoverlay")){
      if(nargc < 1) CMDargNErr(option,1);
      distoverlaypath = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--stats")){
      if(nargc < 1) CMDargNErr(option,1);
      statspath = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--o")){
      if(nargc < 1) CMDargNErr(option,1);
      outvolpath = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--xmask")){
      if(nargc < 1) CMDargNErr(option,1);
      xmaskpath = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--m")){
      if(nargc < 1) CMDargNErr(option,1);
      facesegpath = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--apply")){
      // input facemask reg out
      if(nargc < 4) CMDargNErr(option,4);
      defacer.invol   = MRIread(pargv[0]);
      defacer.faceseg = MRIread(pargv[1]);
      int resampleneeded = 0;
      MRIvol2VolLTA v2v;
      if(strcmp(pargv[2],"regheader") != 0){
	v2v.lta = LTAread(pargv[2]);
	if(v2v.lta==NULL) exit(1);
	resampleneeded = 1;
      } 
      else { // regheader
	v2v.targ = defacer.invol;
	VOL_GEOM vginvol;
	getVolGeom(defacer.invol, &vginvol);
	VOL_GEOM vgfaceseg;
	getVolGeom(defacer.faceseg, &vgfaceseg);
	vg_isEqual_Threshold = 10e-4;
	if(!vg_isEqual(&vginvol, &vgfaceseg)) resampleneeded = 1;
      }
      if(resampleneeded){
	v2v.InterpCode = SAMPLE_NEAREST;
	v2v.mov = defacer.faceseg;
	MRI *mritmp = v2v.vol2vol(NULL);
	if(mritmp==NULL) exit(1);
	MRIfree(&defacer.faceseg);
	defacer.faceseg = mritmp;
      }
      defacer.outvol = MRIallocSequence(defacer.invol->width, defacer.invol->height, 
					defacer.invol->depth, defacer.invol->type, 
					defacer.invol->nframes);
      MRIcopyHeader(defacer.invol, defacer.outvol);
      MRIcopyPulseParameters(defacer.invol, defacer.outvol);
      defacer.FaceIntensityStats();
      defacer.Deface();
      if(embed){
	printf("Embedding code %s into output\n",defacer.embedded_code);
	defacer.EmbedCode(defacer.outvol,defacer.outvol);
      }
      int err = MRIwrite(defacer.outvol,pargv[3]);
      exit(err);
      nargsused = 4;
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
/*---------------------------------------------------------------*/
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/*---------------------------------------------------------------*/
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --i  inputvol \n");
  printf("   --hm headmask \n");
  printf("   --ts tempsurf \n");
  printf("   --l templabel1 <--l templabel2> \n");
  printf("   --w watermark d \n");
  printf("   --o  defacedvol \n");
  printf("   --m  facemask \n");
  printf("   --fill-const ConstIn ConstOut\n");
  printf("   --xmask xmask : exclude anything in this mask from defacing\n");
  printf("\n");
  printf("   --reg tempreg.lta : apply to surface\n");
  printf("   --min minsurfpath : output 'minimum surface'\n");
  printf("   --max maxsurfpath : output 'maximum surface'\n");
  printf("   --distbounds distboundspath : text file with info about distance bounds for each label\n");
  printf("   --distoverlay dist.overlay.mgz : overlay of distance for each vertex\n");
  printf("   --distdat distdatpath : text file with distances for each vertex\n");
  printf("   --stats statspath : has info about nxmask and means and modes\n");
  printf("   --ots outputtempsurf : after any watermark and/or ripple\n");
  printf("\n");
  printf("   --apply vol facemask reg output : apply to another volume (use regheader if no reg needed)\n");
  printf("   --ripple-center R A S\n");
  printf("   --apply-ripple insurf axis amp period label outsurf\n");
  printf("\n");
  printf("   --gdiag diagno : set diagnostic level\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("   --check-code vol <outfile> : determine whether mideface code is in volume\n");
  printf("\n");
  std::cout << getVersion() << std::endl;
  printf("\n");
}
/*---------------------------------------------------------------*/
static void print_help(void) {
  print_usage() ;
  exit(1) ;
}
/*---------------------------------------------------------------*/
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/*---------------------------------------------------------------*/
static void check_options(void) 
{
  if(involpath == NULL){
    printf("ERROR: need input --i\n");
    exit(1);
  }
  if(outvolpath == NULL){
    printf("ERROR: need output --o\n");
    exit(1);
  }
  if(headmaskpath == NULL){
    printf("ERROR: need headmask --hm\n");
    exit(1);
  }
  if(ntemplabelpathlist == 0){
    printf("ERROR: need at least one template label --l\n");
    exit(1);
  }
  if(tempsurfpath == NULL){
    printf("ERROR: need template surface --tl\n");
    exit(1);
  }

  return;
}
/*---------------------------------------------------------------*/
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n", getVersion().c_str());
  fprintf(fp,"cd %s\n",cwd);
  fprintf(fp,"%s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  return;
}

