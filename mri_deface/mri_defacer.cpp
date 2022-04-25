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

class FsDefacer {
public:
  MRI *invol=NULL, *headmask=NULL, *xmask=NULL;
  MRIS *tempsurf=NULL;
  LABEL *templabel=NULL;
  MRI *outvol=NULL, *faceseg=NULL;
  MRIS *minsurf=NULL,*maxsurf=NULL;
  int nface1vox=0,nface2vox=0;
  double gmean1=0, gstddev1=0,gmean2=0, gstddev2=0;
  double min1=0, max1=0, min2=0, max2=0, mode2=0;
  double DistIn=0,DistOut=0,DeltaDist=-1,dL=-1;
  double DistInRaw=0,DistOutRaw=0;
  float *DistInList=NULL,*DistOutList=NULL;
  double DistInMin=2,DistInMax=100; // these will be set in the caller
  double DistOutMin=2,DistOutMax=100; 
  double DistInFrac=0.9, DistOutFrac=1.0;
  int FillType=1;
  double FillConstIn=0,FillConstOut=0;
  int cPad=5, rPad=5, sPad=5;
  int nxmask = 0;
  int DoRipple=0,rippleaxis=1;
  double *ripplecenter=NULL,rippleperiod=20,rippleamp=1;
  int SegFace(void);
  int FaceIntensityStats(void);
  int SetDeltaDist(void);
  int DistanceBounds(void);
  int Deface(void);
  int VoxOutOfBounds(int c, int r, int s);
  MRI *Surf2VolProjFill(MRI *vol, double Dist, double FillVal);
  int PrintParams(FILE *fp);
  int PrintStats(FILE *fp);
  int ripple(LABEL *label);
  int watermark(LABEL *watermark, double dwatermark);
};

int FsDefacer::PrintParams(FILE *fp)
{
  fprintf(fp,"DistInMin   %g\n",DistInMin);
  fprintf(fp,"DistInMax   %g\n",DistInMax);
  fprintf(fp,"DistInFrac  %g\n",DistInFrac);
  fprintf(fp,"DistOutMin  %g\n",DistOutMin);
  fprintf(fp,"DistOutMax  %g\n",DistOutMax);
  fprintf(fp,"DistOutFrac %g\n",DistOutFrac);
  fprintf(fp,"DeltaDist   %g\n",DeltaDist);
  fprintf(fp,"dL          %g\n",dL);
  fprintf(fp,"Pad %d %d %d\n",cPad,rPad,sPad);
  fprintf(fp,"FillType    %d\n",FillType);
  fprintf(fp,"FillConstIn  %g\n",FillConstIn);
  fprintf(fp,"FillConstOut %g\n",FillConstOut);
  fprintf(fp,"DoRipple     %d\n",DoRipple);
  if(DoRipple){
    fprintf(fp,"RippleAmp    %g\n",rippleamp);
    fprintf(fp,"RipplePeriod %g\n",rippleperiod);
  }
  fflush(fp);
  return(0);
}
int FsDefacer::PrintStats(FILE *fp)
{
  fprintf(fp,"nface1vox  %d\n",nface1vox);
  fprintf(fp,"gmean1     %g\n",gmean1);
  fprintf(fp,"gstddev1   %g\n",gstddev1);
  fprintf(fp,"min1       %g\n",min1);
  fprintf(fp,"max1       %g\n",max1);
  fprintf(fp,"nface2vox  %d\n",nface2vox);
  fprintf(fp,"gmean2     %g\n",gmean2);
  fprintf(fp,"gstddev2   %g\n",gstddev2);
  fprintf(fp,"mode2      %g\n",mode2);
  fprintf(fp,"min2       %g\n",min2);
  fprintf(fp,"max2       %g\n",max2);
  fprintf(fp,"gmeanratio %g %g\n",gmean1/gmean2,gmean1/mode2);
  if(xmask) fprintf(fp,"nxmask     %d\n",nxmask);
  else      fprintf(fp,"nxmask     %d\n",-1);
  fflush(fp);
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
	    nxmask ++;
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
	    if(FillType==1){
	      if(m==1 || m==2) v = (drand48()-0.5)*gstddev1 + gmean1;
	      if(m==3 || m==4) v = (drand48()-0.5)*gstddev2 + mode2;
	    }
	    if(FillType==2){
	      if(m==1 || m==2) v = FillConstIn;
	      if(m==3 || m==4) v = FillConstOut;
	    }
	  }
	  v = fabs(v);
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
  min1=0;  max1=0;
  min2=0;  max2=0;
  nface1vox=0;
  nface2vox=0;
  for(c=0; c < invol->width; c++){
    int r,s;
    for(r=0; r < invol->height; r++){
      for(s=0; s < invol->depth; s++){
	int m = MRIgetVoxVal(faceseg,c,r,s,0);
	if(m == 1){ 
	  // inside template surface and in the head mask (ie, in tissue)
	  double v = MRIgetVoxVal(invol,c,r,s,0);
	  sum1 += v;
	  sumsq1 += (v * v);
	  nface1vox++;
	  if(min1 > v) min1 = v;
	  if(max1 < v) max1 = v;
	}
	if(m == 4){
	  // outside template surface and not in the head mask (ie, in background)
	  double v = MRIgetVoxVal(invol,c,r,s,0);
	  sum2 += v;
	  sumsq2 += (v * v);
	  nface2vox++;
	  if(min2 > v) min2 = v;
	  if(max2 < v) max2 = v;
	}
      }
    }
  }
  gmean1 = sum1 / nface1vox;
  gstddev1 = sqrt(sumsq1 / nface1vox - (gmean1) * (gmean1));
  gmean2 = sum2 / nface2vox;
  gstddev2 = sqrt(sumsq2 / nface2vox - (gmean2) * (gmean2));

  // Compute the mode of the outside voxels because there can 
  // be a lot of crud out there.
  HISTO *h2 = HISTOalloc((int)ceil(max2)+1) ;
  HISTOinit(h2, h2->nbins, 0, ceil(max2)) ;
  float bin_size = (h2->max - h2->min)/((float)h2->nbins - 1);
  for(c=0; c < invol->width; c++){
    int r,s;
    for(r=0; r < invol->height; r++){
      for(s=0; s < invol->depth; s++){
	int m = MRIgetVoxVal(faceseg,c,r,s,0);
	if(m == 4){
	  // outside template surface and not in the head mask (in background)
	  double v = MRIgetVoxVal(invol,c,r,s,0);
	  int bn = nint((v - h2->min) / bin_size);
	  h2->counts[bn]++;
	  // Note: should use HISTOaddSample() but it appears to be broken when binsize=1
	}
      }
    }
  }
  int max_count=0, b, bmax=0;
  for (b=0; b < h2->nbins; b++) {
    int center_val = h2->counts[b];
    if(center_val > max_count) {
      max_count = center_val;
      bmax = b;
    }
  }
  printf("Mode2 %d %g  %g\n",bmax,h2->bins[bmax],h2->counts[bmax]);
  mode2 = h2->bins[bmax];
  HISTOfree(&h2);
  printf("nface1vox %d  gmean %g  gstddev %g  min %g max %g\n",nface1vox,gmean1,gstddev1,min1,max1);
  printf("nface2vox %d  gmean %g  gstddev %g  min %g max %g mode %g\n",nface2vox,gmean2,gstddev2,min2,max2,mode2);
  printf("gmeanratio %g %g\n",gmean1/gmean2,gmean1/mode2);

  fflush(stdout);
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

  // This will generate a memory leak with multiple labels
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
  free(DistInListSorted);
  free(DistOutListSorted);

  DistInRaw  = DistIn;
  DistOutRaw = DistOut;
  printf("Raw DistIn = %g, DistOut = %g\n",DistIn,DistOut);
  if(DistIn < DistInMin){
    printf("DistIn=%g < DistInMin=%g, resetting\n",DistIn,DistInMin);
    DistIn = DistInMin;
  }
  if(DistOut < DistOutMin){
    printf("DistOut=%g < DistOutMin=%g, resetting\n",DistOut,DistOutMin);
    DistOut = DistOutMin;
  }
  printf("PostMinCheck DistIn = %g, DistOut = %g\n",DistIn,DistOut);

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
int FsDefacer::ripple(LABEL *label)
{
  BasicSpherePVF sph;
  MRIS *surf = tempsurf;
  int nmax = surf->nvertices;
  if(label) nmax = label->n_points;
  for(int n=0; n < nmax; n++){
    int vtxno;
    if(label) vtxno = label->lv[n].vno;
    else      vtxno = n;
    VERTEX *v = &(surf->vertices[vtxno]);
    double dx = v->x;
    double dy = v->y;
    double dz = v->z;
    if(ripplecenter){
      dx -= ripplecenter[0];
      dy -= ripplecenter[1];
      dz -= ripplecenter[2];
    }
    std::array<double,3> xyz = {dx,dy,dz};
    std::array<double,3> rtp = sph.XYZ2RTP(xyz);
    rtp[0] = rtp[0] + rippleamp*sin(2*M_PI*rtp[0]*rtp[rippleaxis]/rippleperiod);
    xyz = sph.RTP2XYZ(rtp);
    v->x = xyz[0];
    v->y = xyz[1];
    v->z = xyz[2];
    if(ripplecenter){
      v->x += ripplecenter[0];
      v->y += ripplecenter[1];
      v->z += ripplecenter[2];
    }
  }
  return(0);
}

int FsDefacer::watermark(LABEL *watermark, double dwatermark)
{
  printf("Applying watermark %d\n",watermark->n_points);
  for(int n=0; n < watermark->n_points; n++){
    int vtxno = watermark->lv[n].vno;
    VERTEX *v = &(tempsurf->vertices[vtxno]);
    v->x += dwatermark*v->nx;
    v->y += dwatermark*v->ny;
    v->z += dwatermark*v->nz;
  }
  return(0);
}

char *involpath=NULL, *headmaskpath=NULL, *tempsurfpath=NULL;
char *regpath=NULL,*xmaskpath=NULL;
char *templabelpathlist[100];
int ntemplabelpathlist=0;
char *outvolpath=NULL, *facesegpath=NULL, *minsurfpath=NULL, *maxsurfpath=NULL;
char *distdatpath=NULL, *distboundspath=NULL, *distoverlaypath=NULL, *statspath=NULL;
char *watermarkpath = NULL;
double dwatermark = 1;
FsDefacer defacer;
double DistInMinList[200],DistInMaxList[200],DistInMin=2,DistInMax=20;
double DistInList[200];

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
  printf("\n");
  printf("   --reg tempreg.lta : apply to surface\n");
  printf("   --min minsurfpath : output 'minimum surface'\n");
  printf("   --max maxsurfpath : output 'maximum surface'\n");
  printf("   --distbounds distboundspath : text file with info about distance bounds for each label\n");
  printf("   --distoverlay dist.overlay.mgz : overlay of distance for each vertex\n");
  printf("   --distdat distdatpath : text file with distances for each vertex\n");
  printf("   --stats statspath : has info about nxmask and means and modes\n");
  printf("\n");
  printf("   --apply vol facemask reg output : apply to another volume (use regheader if no reg needed)\n");
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


