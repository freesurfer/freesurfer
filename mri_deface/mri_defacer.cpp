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
  MRIS *headsurf=NULL;
  MRI *facesurfmask=NULL;
  MRI *outvol=NULL, *faceseg=NULL;
  MRIS *facesurffit=NULL,*minsurf=NULL,*maxsurf=NULL;
  MRI *facesurffitnorm=NULL;
  int nface1vox=0, nfacevtxs=0;
  double gmean=0, gstddev=0;
  double DistIn=0,DistOut=0,DeltaDist=-1; 
  double DistInMax=100,DistOutMax=100; 
  int cPad=10, rPad=10, sPad=10;
  int StatsComputed=0;
  int SegFace(void);
  int SegFace2(void);
  int FaceIntensityStats(void);
  int FitFace(void);
  int DistanceBounds(void);
  int Deface(void);
  int VoxOutOfBounds(int c, int r, int s);
};

int FsDefacer::FitFace(void)
{
  printf("Fitting face\n");
  int n;
  nfacevtxs = 0;
  for(n=0; n< headsurf->nvertices;n++)
    if(MRIgetVoxVal(facesurfmask,n,0,0,0)>0.5) nfacevtxs++;

  GLMMAT *glm = GLMalloc();
  glm->ncontrasts = 0;

  glm->X = MatrixAlloc(nfacevtxs,6,MATRIX_REAL);
  glm->y = MatrixAlloc(nfacevtxs,1,MATRIX_REAL);
  int k = 0;
  for(n=0; n< headsurf->nvertices;n++){
    if(MRIgetVoxVal(facesurfmask,n,0,0,0)<0.5) continue;
    VERTEX *v = &(headsurf->vertices[n]);
    glm->y->rptr[k+1][1] = v->z;
    glm->X->rptr[k+1][1] = 1;
    glm->X->rptr[k+1][2] = v->x;
    glm->X->rptr[k+1][3] = v->y;
    glm->X->rptr[k+1][4] = v->x * v->x;
    glm->X->rptr[k+1][5] = v->y * v->y;
    glm->X->rptr[k+1][6] = v->x * v->y;
    k = k + 1;
  }

  GLManalyze(glm);

  //MatrixWriteTxt("X.mtx", glm->X);
  //MatrixWriteTxt("y.mtx", glm->y);
  //MatrixWriteTxt("yhat.mtx", glm->yhat);

  if(facesurffit) MRISfree(&facesurffit);
  facesurffit = MRISclone(headsurf);

  k = 0;
  for(n=0; n< headsurf->nvertices;n++){
    if(MRIgetVoxVal(facesurfmask,n,0,0,0)<0.5) continue;
    VERTEX *v = &(facesurffit->vertices[n]);
    v->z = glm->yhat->rptr[k+1][1] ;
    k = k + 1;
  }

  MRIScomputeMetricProperties(facesurffit);
  MRISstoreMetricProperties(facesurffit);
  MRIScomputeNormals(facesurffit);

  // While here, compute the normal to the fit surface
  // Don't rely on the vertex calcs
  facesurffitnorm = MRIallocSequence(headsurf->nvertices, 1, 1, MRI_FLOAT, 3);
  MATRIX *norm = MatrixAlloc(3,1,MATRIX_REAL);
  norm->rptr[3][1] = -1;
  // F = X*beta - z;
  for(n=0; n< headsurf->nvertices;n++){
    if(MRIgetVoxVal(facesurfmask,n,0,0,0)<0.5) continue;
    VERTEX *v = &(headsurf->vertices[n]);
    // dF/dx = b2 + 2*b4*x + b6*y
    norm->rptr[1][1] = glm->beta->rptr[2][1] + 2*glm->beta->rptr[4][1]*v->x + glm->beta->rptr[6][1]*v->y;
    // dF/dy = b3 + 2*b5*y + b6*x
    norm->rptr[2][1] = glm->beta->rptr[3][1] + 2*glm->beta->rptr[5][1]*v->y + glm->beta->rptr[6][1]*v->x;
    // dF/dz = -1
    VectorNormalize(norm,norm);
    for(int k=0; k<3; k++) MRIsetVoxVal(facesurffitnorm,n,0,0,k, norm->rptr[k+1][1]);
  }

  GLMfree(&glm);

  return(0);
}


int FsDefacer::Deface(void)
{
  FsDefacer::FitFace();
  FsDefacer::DistanceBounds();
  FsDefacer::SegFace();
  FsDefacer::FaceIntensityStats();

  printf("FsDefacer::Deface()\n");
  outvol = MRIallocSequence(invol->width, invol->height, invol->depth, MRI_INT, 1);
  MRIcopyHeader(invol, outvol);
  MRIcopyPulseParameters(invol, outvol);
  int c;
  for(c=0; c < invol->width; c++){
    int r,s;
    for(r=0; r < invol->height; r++){
      for(s=0; s < invol->depth; s++){
	int m = MRIgetVoxVal(faceseg,c,r,s,0);
	double v;
	if(m == 0) v = MRIgetVoxVal(invol,c,r,s,0);
	else       v = (drand48()-0.5)*gstddev + gmean;
	MRIsetVoxVal(outvol,c,r,s,0,v);
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
  double sum=0, sumsq=0;
  for(c=0; c < invol->width; c++){
    int r,s;
    for(r=0; r < invol->height; r++){
      for(s=0; s < invol->depth; s++){
	int m = MRIgetVoxVal(faceseg,c,r,s,0);
	if(m != 1) continue;
	double v = MRIgetVoxVal(invol,c,r,s,0);
	sum += v;
	sumsq += (v * v);
	nface1vox++;
      }
    }
  }
  gmean = sum / nface1vox;
  gstddev = sqrt(sumsq / nface1vox - (gmean) * (gmean));
  printf("nface1vox %d  gmean %g  gstddev %g\n",nface1vox,gmean,gstddev);
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

int FsDefacer::DistanceBounds(void)
{
  printf("FsDefacer::DistanceBounds()\n");

  MRIS_SurfRAS2VoxelMap* sras2v_map = MRIS_makeRAS2VoxelMap(headmask, headsurf);

  DeltaDist = MIN(MIN(headmask->xsize/2.0,headmask->ysize/2.0),headmask->zsize/2.0);
  printf("DeltaDist = %g\n",DeltaDist);

  minsurf = MRISclone(facesurffit);
  maxsurf = MRISclone(facesurffit);
  DistIn = 0; 
  DistOut = 0;
  int nthvtx;
  for(nthvtx = 0; nthvtx < headsurf->nvertices; nthvtx++){
    if(MRIgetVoxVal(facesurfmask,nthvtx,0,0,0) < 0.5) continue;
    VERTEX *v;
    v = &(facesurffit->vertices[nthvtx]);

    double d;
    double x=0,y=0,z=0, c,r,s;

    // Project inward until reached max or hit headmask or hit xmask
    d = 0;
    while(d < DistInMax){
      // project xyz inward to this distance
      double nx, ny, nz;
      nx = MRIgetVoxVal(facesurffitnorm,nthvtx,0,0,0);
      ny = MRIgetVoxVal(facesurffitnorm,nthvtx,0,0,1);
      nz = MRIgetVoxVal(facesurffitnorm,nthvtx,0,0,2);
      x = v->x - nx * d;
      y = v->y - ny * d;
      z = v->z - nz * d;
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
      minsurf->vertices[nthvtx].x = x;
      minsurf->vertices[nthvtx].y = y;
      minsurf->vertices[nthvtx].z = z;
    }

    // Project outward until reached max or go outside of headmask or hit xmask
    // Can it go out and then back in? Eg, nostril?
    d=0;
    while(d < DistOutMax){
      // project xyz outward to this distance
      double nx, ny, nz;
      nx = MRIgetVoxVal(facesurffitnorm,nthvtx,0,0,0);
      ny = MRIgetVoxVal(facesurffitnorm,nthvtx,0,0,1);
      nz = MRIgetVoxVal(facesurffitnorm,nthvtx,0,0,2);
      x = v->x + nx * d;
      y = v->y + ny * d;
      z = v->z + nz * d;
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
      maxsurf->vertices[nthvtx].x = x;
      maxsurf->vertices[nthvtx].y = y;
      maxsurf->vertices[nthvtx].z = z;
    }

  }// vertex

  printf("DistIn = %g, DistOut = %g\n",DistIn,DistOut);

#if 1
  for(nthvtx = 0; nthvtx < headsurf->nvertices; nthvtx++){
    if(MRIgetVoxVal(facesurfmask,nthvtx,0,0,0) < 0.5) continue;
    VERTEX *v;
    // project min xyz inward to this distance
    v = &(minsurf->vertices[nthvtx]);
    double nx,ny,nz;
    nx = MRIgetVoxVal(facesurffitnorm,nthvtx,0,0,0);
    ny = MRIgetVoxVal(facesurffitnorm,nthvtx,0,0,1);
    nz = MRIgetVoxVal(facesurffitnorm,nthvtx,0,0,2);
    v->x = v->x - nx * DistIn;
    v->y = v->y - ny * DistIn;
    v->z = v->z - nz * DistIn;
    // project max xyz outward to this distance
    v = &(maxsurf->vertices[nthvtx]);
    nx = MRIgetVoxVal(facesurffitnorm,nthvtx,0,0,0);
    ny = MRIgetVoxVal(facesurffitnorm,nthvtx,0,0,1);
    nz = MRIgetVoxVal(facesurffitnorm,nthvtx,0,0,2);
    v->x = v->x + nx * DistOut;
    v->y = v->y + ny * DistOut;
    v->z = v->z + nz * DistOut;
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

  MRIS_SurfRAS2VoxelMap* sras2v_map = MRIS_makeRAS2VoxelMap(invol, headsurf);
  faceseg = MRIallocSequence(invol->width, invol->height, invol->depth, MRI_INT, 1);
  MRIcopyHeader(invol, faceseg);
  MRIcopyPulseParameters(invol, faceseg);

  int nthvtx;
  for(nthvtx = 0; nthvtx < headsurf->nvertices; nthvtx++){
    if(MRIgetVoxVal(facesurfmask,nthvtx,0,0,0) < 0.5) continue;
    VERTEX *v = &(facesurffit->vertices[nthvtx]);
    double d;
    for(d = -DistIn; d <= DistOut; d+= DeltaDist){ 
      double x,y,z, c,r,s;
      // project xyz to this distance
      double nx, ny, nz;
      nx = MRIgetVoxVal(facesurffitnorm,nthvtx,0,0,0);
      ny = MRIgetVoxVal(facesurffitnorm,nthvtx,0,0,1);
      nz = MRIgetVoxVal(facesurffitnorm,nthvtx,0,0,2);
      x = v->x + nx * d;
      y = v->y + ny * d;
      z = v->z + nz * d;
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
      if(d<=0) MRIsetVoxVal(faceseg,ic,ir,is,0,1);
      else     MRIsetVoxVal(faceseg,ic,ir,is,0,2);
    }
  }
  return(0);
}

int FsDefacer::SegFace2(void)
{
  printf("FsDefacer::SegFace2()\n");
  DeltaDist = MIN(MIN(headmask->xsize/2.0,headmask->ysize/2.0),headmask->zsize/2.0);
  printf("DeltaDist = %g\n",DeltaDist);

  MRIS_SurfRAS2VoxelMap* sras2v_map = MRIS_makeRAS2VoxelMap(invol, headsurf);
  faceseg = MRIallocSequence(invol->width, invol->height, invol->depth, MRI_INT, 1);
  MRIcopyHeader(invol, faceseg);
  MRIcopyPulseParameters(invol, faceseg);

  int nthvtx;
  for(nthvtx = 0; nthvtx < headsurf->nvertices; nthvtx++){
    if(MRIgetVoxVal(facesurfmask,nthvtx,0,0,0) < 0.5) continue;
    double x0,y0,z0,dx,dy,dz,dtot;
    VERTEX *v0 = &(headsurf->vertices[nthvtx]);
    VERTEX *v  = &(facesurffit->vertices[nthvtx]);
    x0 = v0->x;
    y0 = v0->y;
    z0 = v0->z;
    dx = v->x - v0->x;
    dy = v->y - v0->y;
    dz = v->z - v0->z;
    dtot = sqrt(dx*dx + dy*dy + dz*dz);
    double deltax, deltay, deltaz;
    deltax = x0 + DeltaDist*dx/dtot;
    deltay = y0 + DeltaDist*dy/dtot;
    deltaz = z0 + DeltaDist*dz/dtot;

    // Determine whether going inside or outside the mean surface
    double dot = deltax * v0->nx + deltay * v0->ny + deltaz * v0->nz;
    int seg;
    if(dot <= 0) seg = 1; // inside the mean surface
    if(dot >  0) seg = 2; // outside the mean surface

    double d;
    while(d<dtot){
      double x,y,z, c,r,s;
      // project xyz to this distance
      x = x0 + DeltaDist*dx/dtot;
      y = y0 + DeltaDist*dy/dtot;
      z = z0 + DeltaDist*dz/dtot;
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
      MRIsetVoxVal(faceseg,ic,ir,is,0, seg);
    }
  }
  return(0);
}

char *involpath=NULL, *headmaskpath=NULL, *headsurfpath=NULL, *facesurfmaskpath=NULL;
char *outvolpath=NULL, *facesegpath=NULL, *facefitsurfpath=NULL;
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
  defacer.headmask = MRIread(headmaskpath);
  if(defacer.headmask==NULL) exit(1);
  defacer.headsurf = MRISread(headsurfpath);
  if(defacer.headsurf==NULL) exit(1);
  defacer.facesurfmask = MRIread(facesurfmaskpath);
  if(defacer.facesurfmask==NULL) exit(1);

  defacer.Deface();

  err = MRIwrite(defacer.outvol,outvolpath);
  if(err) exit(1);

  if(facesegpath){
    err = MRIwrite(defacer.faceseg,facesegpath);
    if(err) exit(1);
  }

  if(facefitsurfpath){
    err = MRISwrite(defacer.facesurffit,facefitsurfpath);
    if(err) exit(1);
  }
  err = MRISwrite(defacer.minsurf,"min.surf");
  err = MRISwrite(defacer.maxsurf,"max.surf");

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
    else if(!strcasecmp(option, "--hs")){
      if(nargc < 1) CMDargNErr(option,1);
      headsurfpath = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--fsm")){
      if(nargc < 1) CMDargNErr(option,1);
      facesurfmaskpath = pargv[0];
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
    else if(!strcasecmp(option, "--facefit")){
      if(nargc < 1) CMDargNErr(option,1);
      facefitsurfpath = pargv[0];
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
  printf("   --hs headsurf \n");
  printf("   --fsm facesurfmask \n");
  printf("   --o defacedvol \n");
  printf("   --m facemask \n");
  printf("   --facefit surf : surf with face fit to poly");
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
    printf("ERROR: need head mask --hm\n");
    exit(1);
  }
  if(headsurfpath == NULL){
    printf("ERROR: need head surface --hs\n");
    exit(1);
  }
  if(facesurfmaskpath == NULL){
    printf("ERROR: need face surface mask --fsm\n");
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


