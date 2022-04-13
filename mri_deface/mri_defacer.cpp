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
#include "mrishash.h"

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
  MRI *outvol=NULL, *faceseg=NULL;
  MRIS *minsurf=NULL,*maxsurf=NULL;
  int nface1vox=0,nface2vox=0;
  double gmean1=0, gstddev1=0,gmean2=0, gstddev2=0;
  double DistIn=0,DistOut=0,DeltaDist=-1,dL=-1;
  float *DistInList=NULL,*DistOutList=NULL;
  double DistInMax=100,DistOutMax=100; 
  double DistInFrac=0.9, DistOutFrac=1.0;
  int cPad=5, rPad=5, sPad=5;
  int SegFace(void);
  int FaceIntensityStats(void);
  int SetDeltaDist(void);
  int DistanceBounds(void);
  int Deface(void);
  int VoxOutOfBounds(int c, int r, int s);
  MRI *Surf2VolProjFill(MRI *vol, double Dist, double FillVal);
};

/*!
  \fn MRI *FsDefacer::Surf2VolProjFill(MRI *vol, double Dist, double FillVal)
  \brief This function fills in voxels within Dist of the template
  surface.  It goes through all the vertices in the template label,
  then fills the voxels intersecting with each neighboring face. This
  is done at each projection depth from 0 to Dist step DeltaDist
  (where Dist can be negative). The triangle filling is done by
  sampling the triangle in barycentric space with sampling distance
  dL, which must be sufficiently small to prevent missing of
  voxels. Simply filling along a projected normal does not work very
  well because it leaves a lot of voxels unsampled (worse for
  projecting out than in). Finding voxels inside a projected surface
  does not well because the projection often creates folds and
  intersections rendering the normal inaccurate.
 */
MRI *FsDefacer::Surf2VolProjFill(MRI *vol, double Dist, double FillVal)
{
  if(vol == NULL) {
    vol = MRIcopy(invol,NULL);
    MRIcopyHeader(invol, vol);
    MRIcopyPulseParameters(invol, vol);
    MRIsetValues(vol,0);
  }

  MATRIX *tkras = MatrixAlloc(4,1,MATRIX_REAL);
  tkras->rptr[4][1] = 1;
  MATRIX *crs = NULL;
  MATRIX *Vox2TkRAS = MRIxfmCRS2XYZtkreg(vol);
  MATRIX *TkRAS2Vox = MatrixInverse(Vox2TkRAS,NULL);
  MatrixFree(&Vox2TkRAS);

  // Get min and max absolute depth
  double dmin, dmax;
  dmin = MIN(0,Dist);
  dmax = MAX(0,Dist);

  double p[3][3];
  int nthp,vtxno;
  for(nthp=0; nthp < templabel->n_points; nthp++){
    vtxno = templabel->lv[nthp].vno;
    VERTEX_TOPOLOGY *vt = &(tempsurf->vertices_topology[vtxno]);
    for(int nthface=0; nthface < vt->num; nthface++){
      int faceno = vt->f[nthface];
      FACE *face = &(tempsurf->faces[faceno]);
      for(double d=dmin; d<=dmax; d += DeltaDist){
	// Load the corner xyz's at this depth into the p[][] array
	for(int n=0; n < 3; n++){
	  int vnof = face->v[n];
	  VERTEX *vf = &(tempsurf->vertices[vnof]);
	  tkras->rptr[1][1] = vf->x + d*vf->nx;
	  tkras->rptr[2][1] = vf->y + d*vf->ny;
	  tkras->rptr[3][1] = vf->z + d*vf->nz;
	  crs = MatrixMultiplyD(TkRAS2Vox,tkras,crs);
	  for(int k=0; k<3; k++) p[n][k] = crs->rptr[k+1][1];
	}
	// Fill this triangle
	MRIfillTriangle(vol, p[0],p[1],p[2], dL, FillVal);
      } // dist
    } // face
  } //nthpoint

  MatrixFree(&tkras);
  MatrixFree(&crs);
  MatrixFree(&TkRAS2Vox);
  return(vol);
}

int FsDefacer::Deface(void)
{
  printf("FsDefacer::Deface()\n");
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
  printf("FsDefacer::FaceIntensityStats()\n");

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
  // DeltaDist and dL control how finely the surf2vol sampling is done.
  // It is probably ok if both are MinVoxSize/2
  DeltaDist = MIN(MIN(headmask->xsize,headmask->ysize),headmask->zsize)/3.0;
  dL = DeltaDist/3.0;
  printf("DeltaDist = %g, dL = %g\n",DeltaDist,dL);
  return(0);
}

/*!
  \fn int FsDefacer::DistanceBounds(void)
  \brief Computes the min and max distance from template face needed
  to cover individual's face. At each face point, the headmask is
  sampled over a range of depths. When projecting inward, the distance
  for a vertex is how far it has to project before it hits the head
  mask.  When projecting outward, the distance for a vertex is how far
  it has to project before it leaves the head mask. The final (fixed)
  distance will encompas DistInFrac or DistOutFrac fraction of the
  vertices.  There is little harm in setting DistOutFrac to 1. Setting
  of DistInFrac is trickier. The closer it is to 1, the more tissue
  will be removed.  The closer DistInFrac is to 0, the more of the
  subject's face will be unmasked. The risk of a break is still very
  low as the unmasked part will likely be behind some mask, and any
  problems are usually around the chin.
*/
int FsDefacer::DistanceBounds(void)
{
  printf("FsDefacer::DistanceBounds() DeltaDist=%g\n",DeltaDist);

  MRIS_SurfRAS2VoxelMap* sras2v_map = MRIS_makeRAS2VoxelMap(headmask, tempsurf);

  DistInList  = (float *) calloc(templabel->n_points,sizeof(float));
  DistOutList = (float *) calloc(templabel->n_points,sizeof(float));

  DistIn = 0; 
  DistOut = 0;
  int nthp,vtxno;
  double d;
  double x=0,y=0,z=0, c,r,s;
  for(nthp=0; nthp < templabel->n_points; nthp++){
    vtxno = templabel->lv[nthp].vno;
    VERTEX *v = &(tempsurf->vertices[vtxno]);

    // Project inward until reached max or hit headmask
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
      d += DeltaDist;
    }
    DistInList[nthp] = d;
    v->val = d;

    // Project outward until reached max or go outside of headmask
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
      d += DeltaDist;
    }
    DistOutList[nthp] = d;
    v->valbak = d;

  }// vertex
  
  // Sort the inward distances to get the right value
  float *DistInListSorted  = (float *) calloc(templabel->n_points,sizeof(float));
  memcpy(DistInListSorted,DistInList,templabel->n_points*sizeof(float));
  qsort(DistInListSorted,  templabel->n_points, sizeof(float), compare_floats);
  int k = round((templabel->n_points-1) * DistInFrac);
  DistIn  = DistInListSorted[k];
  // Sort the outward distances to get the right value
  float *DistOutListSorted = (float *) calloc(templabel->n_points,sizeof(float));
  memcpy(DistOutListSorted,DistOutList,templabel->n_points*sizeof(float));
  qsort(DistOutListSorted, templabel->n_points, sizeof(float), compare_floats);
  k = round((templabel->n_points-1)* DistOutFrac);
  DistOut = DistOutListSorted[k];

  printf("DistIn = %g, DistOut = %g\n",DistIn,DistOut);
  free(DistInListSorted);
  free(DistOutListSorted);

  // Compute min and max surfaces 
  // Used only for display/qa/debug
  // Don't have to do it here
  for(nthp=0; nthp < templabel->n_points; nthp++){
    vtxno = templabel->lv[nthp].vno;
    VERTEX *v0 = &(tempsurf->vertices[vtxno]);
    if(v0->val    < DistIn)  v0->val = DistIn;
    if(v0->valbak < DistOut) v0->valbak = DistOut;
    VERTEX *v;
    // project min xyz inward to this distance
    v = &(minsurf->vertices[vtxno]);
    v->x = v0->x - v0->nx * v0->val;
    v->y = v0->y - v0->ny * v0->val;
    v->z = v0->z - v0->nz * v0->val;
    // project max xyz outward to this distance
    v = &(maxsurf->vertices[vtxno]);
    v->x = v0->x + v0->nx * v0->valbak;
    v->y = v0->y + v0->ny * v0->valbak;
    v->z = v0->z + v0->nz * v0->valbak;
  }

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

  // First, fill the inner and outer face mask areas
  faceseg = Surf2VolProjFill(faceseg, -DistIn,  1); // inside  template face
  faceseg = Surf2VolProjFill(faceseg, +DistOut, 4); // outside template face

  // Now change depending upon xmask and whether the head is there or not
  int c;
  for(c=0; c < invol->width; c++){
    int r,s,m,xm,v;
    for(r=0; r < invol->height; r++){
      for(s=0; s < invol->depth; s++){
	m = MRIgetVoxVal(faceseg,c,r,s,0);
	if(m < 0.5) continue; // not in the face mask, skip
	if(xmask){
	  xm = MRIgetVoxVal(xmask,c,r,s,0);
	  if(xm>0.5) {
	    // this vox is in the exclusion mask so zero it
	    MRIsetVoxVal(faceseg,c,r,s,0, 0); 
	    continue;
	  }
	}
	v = MRIgetVoxVal(headmask,c,r,s,0);
	if(m == 1){ // inside template face
	  if(!v) MRIsetVoxVal(faceseg,c,r,s,0, 2); // outside the head mask
	}
	if(m == 4){// outside template face
	  if(v) MRIsetVoxVal(faceseg,c,r,s,0, 3); // inside the head mask
	}
      }
    }
  }

  return(0);
}

char *involpath=NULL, *headmaskpath=NULL, *tempsurfpath=NULL, *templabelpath=NULL;
char *regpath=NULL,*xmaskpath=NULL;
char *outvolpath=NULL, *facesegpath=NULL, *minsurfpath=NULL, *maxsurfpath=NULL;
char *distdatpath=NULL;
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
  defacer.tempsurf->hemisphere = 3;
  defacer.headmask = MRIread(headmaskpath);
  if(defacer.headmask==NULL) exit(1);
  if(xmaskpath){
    defacer.xmask = MRIread(xmaskpath);
    if(defacer.xmask==NULL) exit(1);
  }

  if(regpath){
    printf("Applying reg %s to the template surface\n",regpath);
    LTA *lta = LTAread(regpath);
    if(lta==NULL) exit(1);
    int err = MRISltaMultiply(defacer.tempsurf, lta);
    if(err) exit(1);
    MRIScomputeMetricProperties(defacer.tempsurf);
    MRISsmoothSurfaceNormals(defacer.tempsurf,10);
  }

  defacer.outvol = MRIallocSequence(defacer.invol->width, defacer.invol->height, defacer.invol->depth, defacer.invol->type, defacer.invol->nframes);
  MRIcopyHeader(defacer.invol, defacer.outvol);
  MRIcopyPulseParameters(defacer.invol, defacer.outvol);
  defacer.minsurf = MRISclone(defacer.tempsurf);
  defacer.maxsurf = MRISclone(defacer.tempsurf);

  defacer.SetDeltaDist();

  defacer.templabel = LabelRead("",templabelpath);
  if(defacer.templabel==NULL) exit(1);

  defacer.DistanceBounds();
  defacer.SegFace();
  defacer.FaceIntensityStats();
  defacer.Deface();


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
  if(distdatpath){
    FILE *fp = fopen(distdatpath,"w");
    for(int nthp=0; nthp < defacer.templabel->n_points; nthp++){
      int vtxno = defacer.templabel->lv[nthp].vno;
      fprintf(fp,"%4d %6d %g %g\n",nthp,vtxno,defacer.DistInList[nthp],defacer.DistOutList[nthp]);
    }
    fclose(fp);
  }

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
    else if(!strcasecmp(option, "--distdat")){
      if(nargc < 1) CMDargNErr(option,1);
      distdatpath = pargv[0];
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
  printf("   --tl templabel \n");
  printf("   --o  defacedvol \n");
  printf("   --m  facemask \n");
  printf("\n");
  printf("   --reg tempreg.lta \n");
  printf("   --min minsurfpath \n");
  printf("   --max maxsurfpath \n");
  printf("   --distdat distdatpath \n");
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


