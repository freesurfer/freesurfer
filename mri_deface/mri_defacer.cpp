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

class FsDefacer {
public:
  MRI *invol=NULL, *headmask=NULL, *xmask=NULL;
  MRIS *tempsurf=NULL;
  LABEL *templabel=NULL;
  MRI *tempmask=NULL;
  MRI *outvol=NULL, *faceseg=NULL;
  MRIS *minsurf=NULL,*maxsurf=NULL;
  int nface1vox=0,nface2vox=0;
  double gmean1=0, gstddev1=0,gmean2=0, gstddev2=0;
  double DistIn=0,DistOut=0,DeltaDist=-1; 
  double DistInMax=100,DistOutMax=100; 
  int cPad=5, rPad=5, sPad=5;
  int StatsComputed=0;
  int SegFace(void);
  int FaceIntensityStats(void);
  int SetDeltaDist(void);
  int DistanceBounds(void);
  int Deface(void);
  int VoxOutOfBounds(int c, int r, int s);
  int YeaYea(void);
};


int FsDefacer::YeaYea(void)
{
  printf("FsDefacer::YeaYea() DeltaDist=%g\n",DeltaDist);

  //outvol = MRIallocSequence(invol->width, invol->height, invol->depth, invol->type, invol->nframes);
  outvol = MRIcopy(invol,NULL);
  MRIcopyHeader(invol, outvol);
  MRIcopyPulseParameters(invol, outvol);
  MRIS_SurfRAS2VoxelMap* sras2v_map = MRIS_makeRAS2VoxelMap(headmask, tempsurf);
  int nthp,vtxno;
  for(nthp=0; nthp < templabel->n_points; nthp++){
    vtxno = templabel->lv[nthp].vno;
    VERTEX *v = &(tempsurf->vertices[vtxno]);

    double d,dmax;
    double x=0,y=0,z=0, c,r,s;

    dmax = drand48()*3+1;
    d = -dmax;
    while(d < dmax){
      // project xyz inward to this distance
      x = v->x - v->nx * d;
      y = v->y - v->ny * d;
      z = v->z - v->nz * d;
      MRIS_useRAS2VoxelMap(sras2v_map, headmask, x, y, z, &c, &r, &s);
      int OutOfBounds = FsDefacer::VoxOutOfBounds(c,r,s);
      if(OutOfBounds) break;
      int ic=nint(c);
      int ir=nint(r);
      int is=nint(s);
      if(xmask){
	int m = MRIgetVoxVal(xmask,ic,ir,is,0);
	if(m > 0.5) continue; // hit the exclusion mask
      }
      double val = (drand48()-0.5)*gstddev1 + gmean1;
      MRIsetVoxVal(outvol,ic,ir,is,0,val);
      d += DeltaDist;
    }
  }// vertex

  return(0);
}

int FsDefacer::Deface(void)
{
  printf("FsDefacer::Deface()\n");
  outvol = MRIallocSequence(invol->width, invol->height, invol->depth, invol->type, invol->nframes);
  MRIcopyHeader(invol, outvol);
  MRIcopyPulseParameters(invol, outvol);
  int c;
  for(c=0; c < invol->width; c++){
    int r,s,f;
    for(r=0; r < invol->height; r++){
      for(s=0; s < invol->depth; s++){
	int m = MRIgetVoxVal(faceseg,c,r,s,0);
	for(f=0; f < invol->nframes; f++){
	  double v=0;
	  if(m == 0) v = MRIgetVoxVal(invol,c,r,s,0);// copy input
	  else{
	    if(m==1 || m==2) v = (drand48()-0.5)*gstddev1 + gmean1;
	    if(m==3 || m==4) v = (drand48()-0.5)*gstddev2 + gmean2;
	  }
	  MRIsetVoxVal(outvol,c,r,s,0,v);
	}
      }
    }
  }
  return(0);
}

int FsDefacer::FaceIntensityStats(void)
{
  if(faceseg == NULL) FsDefacer::SegFace();
  printf("FsDefacer::FaceIntensityStats() StatsComputed=%d\n",StatsComputed);
  if(StatsComputed) return(0);

  int c;
  double sum1=0, sumsq1=0, sum2=0, sumsq2=0;
  for(c=0; c < invol->width; c++){
    int r,s;
    for(r=0; r < invol->height; r++){
      for(s=0; s < invol->depth; s++){
	int m = MRIgetVoxVal(faceseg,c,r,s,0);
	if(m == 1){ // inside template surface and in the head mask
	  double v = MRIgetVoxVal(invol,c,r,s,0);
	  sum1 += v;
	  sumsq1 += (v * v);
	  nface1vox++;
	}
	if(m == 4){// outside template surface and not in the head mask
	  double v = MRIgetVoxVal(invol,c,r,s,0);
	  sum2 += v;
	  sumsq2 += (v * v);
	  nface2vox++;
	}
      }
    }
  }
  gmean1 = sum1 / nface1vox;
  gstddev1 = sqrt(sumsq1 / nface1vox - (gmean1) * (gmean1));
  gmean2 = sum2 / nface2vox;
  gstddev2 = sqrt(sumsq2 / nface2vox - (gmean2) * (gmean2));
  printf("nface1vox %d  gmean %g  gstddev %g\n",nface1vox,gmean1,gstddev1);
  printf("nface2vox %d  gmean %g  gstddev %g\n",nface2vox,gmean2,gstddev2);
  StatsComputed=1;
  return(0);
}

int FsDefacer::VoxOutOfBounds(int c, int r, int s)
{
  if(c < cPad) return(1);
  if(c >= invol->width - cPad) return(1);
  if(r < rPad) return(1);
  if(r >= invol->height - rPad) return(1);
  if(s < sPad) return(1);
  if(s >= invol->depth - sPad) return(1);
  return(0);
}

int FsDefacer::SetDeltaDist(void)
{
  DeltaDist = MIN(MIN(headmask->xsize/2.0,headmask->ysize/2.0),headmask->zsize/2.0);
  printf("DeltaDist = %g\n",DeltaDist);
  return(0);
}
int FsDefacer::DistanceBounds(void)
{
  printf("FsDefacer::DistanceBounds() DeltaDist=%g\n",DeltaDist);

  MRIS_SurfRAS2VoxelMap* sras2v_map = MRIS_makeRAS2VoxelMap(headmask, tempsurf);

  minsurf = MRISclone(tempsurf);
  maxsurf = MRISclone(tempsurf);
  DistIn = 0; 
  DistOut = 0;
  int nthp,vtxno;
  for(nthp=0; nthp < templabel->n_points; nthp++){
    vtxno = templabel->lv[nthp].vno;
    VERTEX *v = &(tempsurf->vertices[vtxno]);

    double d;
    double x=0,y=0,z=0, c,r,s;

    // Project inward until reached max or hit headmask or hit xmask
    d = 0;
    while(d < DistInMax){
      // project xyz inward to this distance
      x = v->x - v->nx * d;
      y = v->y - v->ny * d;
      z = v->z - v->nz * d;
      MRIS_useRAS2VoxelMap(sras2v_map, headmask, x, y, z, &c, &r, &s);
      int OutOfBounds = FsDefacer::VoxOutOfBounds(c,r,s);
      if(OutOfBounds) break;
      int ic=nint(c);
      int ir=nint(r);
      int is=nint(s);
      int m = MRIgetVoxVal(headmask,ic,ir,is,0);
      if(m > 0.5) break; // hit the headmask
      if(xmask){
	int m = MRIgetVoxVal(xmask,ic,ir,is,0);
	if(m > 0.5) break; // hit the exclusion mask
      }
      d += DeltaDist;
    }
    if(DistIn < d) DistIn=d;

    if(0){
      minsurf->vertices[vtxno].x = x;
      minsurf->vertices[vtxno].y = y;
      minsurf->vertices[vtxno].z = z;
    }

    // Project outward until reached max or go outside of headmask or hit xmask
    // Can it go out and then back in? Eg, nostril?
    d=0;
    while(d < DistOutMax){
      // project xyz outward to this distance
      x = v->x + v->nx * d;
      y = v->y + v->ny * d;
      z = v->z + v->nz * d;
      MRIS_useRAS2VoxelMap(sras2v_map, headmask, x, y, z, &c, &r, &s);
      int OutOfBounds = FsDefacer::VoxOutOfBounds(c,r,s);
      if(OutOfBounds) break;
      int ic=nint(c);
      int ir=nint(r);
      int is=nint(s);
      int m = MRIgetVoxVal(headmask,ic,ir,is,0);
      if(m < 0.5) break; // now out of headmask
      if(xmask){
	int m = MRIgetVoxVal(xmask,ic,ir,is,0);
	if(m > 0.5) break; // hit the exclusion mask
      }
      d += DeltaDist;
    }
    if(DistOut < d) DistOut=d;

    if(0){
      maxsurf->vertices[vtxno].x = x;
      maxsurf->vertices[vtxno].y = y;
      maxsurf->vertices[vtxno].z = z;
    }

  }// vertex

  printf("DistIn = %g, DistOut = %g\n",DistIn,DistOut);

#if 1
  for(nthp=0; nthp < templabel->n_points; nthp++){
    vtxno = templabel->lv[nthp].vno;
    VERTEX *v0 = &(tempsurf->vertices[vtxno]);
    VERTEX *v;
    // project min xyz inward to this distance
    v = &(minsurf->vertices[vtxno]);
    v->x = v->x - v0->nx * DistIn;
    v->y = v->y - v0->ny * DistIn;
    v->z = v->z - v0->nz * DistIn;
    // project max xyz outward to this distance
    v = &(maxsurf->vertices[vtxno]);
    v->x = v->x + v0->nx * DistOut;
    v->y = v->y + v0->ny * DistOut;
    v->z = v->z + v0->nz * DistOut;
  }
#endif

  MRIScomputeMetricProperties(minsurf);
  MRISstoreMetricProperties(minsurf);
  MRIScomputeNormals(minsurf);
  MRIScomputeMetricProperties(maxsurf);
  MRISstoreMetricProperties(maxsurf);
  MRIScomputeNormals(maxsurf);

  return(0);
}

int FsDefacer::SegFace(void)
{
  printf("FsDefacer::SegFace()\n");

  MRIS_SurfRAS2VoxelMap* sras2v_map = MRIS_makeRAS2VoxelMap(invol, tempsurf);
  faceseg = MRIallocSequence(invol->width, invol->height, invol->depth, MRI_INT, 1);
  MRIcopyHeader(invol, faceseg);
  MRIcopyPulseParameters(invol, faceseg);

  int nthp;
  for(nthp=0; nthp < templabel->n_points; nthp++){
    int vtxno = templabel->lv[nthp].vno;
    VERTEX *v = &(tempsurf->vertices[vtxno]);
    double d;
    for(d = -DistIn; d <= DistOut; d+= DeltaDist){ 
      double x,y,z, c,r,s;
      // project xyz to this distance
      x = v->x + v->nx * d;
      y = v->y + v->ny * d;
      z = v->z + v->nz * d;
      MRIS_useRAS2VoxelMap(sras2v_map, invol, x, y, z, &c, &r, &s);
      int OutOfBounds = MRIindexNotInVolume(invol,c,r,s);
      if(OutOfBounds) continue;
      int ic=nint(c);
      int ir=nint(r);
      int is=nint(s);
      if(xmask){
	int m = MRIgetVoxVal(xmask,ic,ir,is,0);
	if(m>0.5)continue;
      }
      int hm = MRIgetVoxVal(headmask,ic,ir,is,0);
      if(d<=0) {
	if(hm) MRIsetVoxVal(faceseg,ic,ir,is,0,1);// inside surf AND in mask
	else   MRIsetVoxVal(faceseg,ic,ir,is,0,2);// inside surf AND out mask
      }
      else{
	if(hm) MRIsetVoxVal(faceseg,ic,ir,is,0,3); // outside surf AND in mask
	else   MRIsetVoxVal(faceseg,ic,ir,is,0,4); // outside surf AND out mask
      }
    }
  }
  return(0);
}

char *involpath=NULL, *headmaskpath=NULL, *tempsurfpath=NULL, *templabelpath=NULL;
char *regpath=NULL;
char *outvolpath=NULL, *facesegpath=NULL;
/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int nargs,err=0;
  FsDefacer defacer;

  nargs = handleVersionOption(argc, argv, "mri_gtmpvc");
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
  defacer.templabel = LabelRead("",templabelpath);
  if(defacer.templabel==NULL) exit(1);
  defacer.headmask = MRIread(headmaskpath);
  if(defacer.headmask==NULL) exit(1);

  if(regpath){
    printf("Applying reg %s to the template surface\n",regpath);
    LTA *lta = LTAread(regpath);
    if(lta==NULL) exit(1);
    int err = MRISltaMultiply(defacer.tempsurf, lta);
    if(err) exit(1);
    MRIScomputeMetricProperties(defacer.tempsurf);
  }

  defacer.SetDeltaDist();
  defacer.gmean1 = 100;
  defacer.gstddev1 = 10;
  defacer.YeaYea();
  if(0){
  defacer.DistanceBounds();
  defacer.SegFace();
  defacer.FaceIntensityStats();
  defacer.Deface();
  }
  err = MRIwrite(defacer.outvol,outvolpath);
  if(err) exit(1);

  if(facesegpath){
    //err = MRIwrite(defacer.faceseg,facesegpath);
    //if(err) exit(1);
  }
  //err = MRISwrite(defacer.minsurf,"rh.min.surf");
  //err = MRISwrite(defacer.maxsurf,"rh.max.surf");

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
    else if(!strcasecmp(option, "--tl")){
      if(nargc < 1) CMDargNErr(option,1);
      templabelpath = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--reg")){
      if(nargc < 1) CMDargNErr(option,1);
      regpath = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--o")){
      if(nargc < 1) CMDargNErr(option,1);
      outvolpath = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--m")){
      if(nargc < 1) CMDargNErr(option,1);
      facesegpath = pargv[0];
      nargsused = 1;
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
  printf("   --i   inputvol \n");
  printf("   --hm headmask \n");
  printf("   --ts tempsurf \n");
  printf("   --tl templabel \n");
  printf("   --reg tempreg.lta \n");
  printf("   --o defacedvol \n");
  printf("   --m facemask \n");
  printf("\n");
  printf("   --gdiag diagno : set diagnostic level\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
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
    printf("ERROR: need template label --hm\n");
    exit(1);
  }
  if(templabelpath == NULL){
    printf("ERROR: need template label --tl\n");
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


