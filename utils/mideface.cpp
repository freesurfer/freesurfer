/*
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

// mideface.cpp - utilities for minimally invasive defacing tool
#undef X

int MiDeface::PrintParams(FILE *fp)
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
int MiDeface::PrintStats(FILE *fp)
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

int MiDeface::SegFace(void)
{
  printf("MiDeface::SegFace()\n");

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
  \fn MRI *MiDeface::Surf2VolProjFill(MRI *vol, double Dist, double FillVal)
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
MRI *MiDeface::Surf2VolProjFill(MRI *vol, double Dist, double FillVal)
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

int MiDeface::Deface(void)
{
  printf("MiDeface::Deface()\n");
  int c;
  for(c=0; c < invol->width; c++){
    int r,s,f;
    for(r=0; r < invol->height; r++){
      for(s=0; s < invol->depth; s++){
	int m = MRIgetVoxVal(faceseg,c,r,s,0);
	for(f=0; f < invol->nframes; f++){
	  double v=0;
	  if(m == 0) v = MRIgetVoxVal(invol,c,r,s,f);// copy input
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
	  MRIsetVoxVal(outvol,c,r,s,f,v);
	}
      }
    }
  }
  return(0);
}

int MiDeface::FaceIntensityStats(void)
{
  if(faceseg == NULL) MiDeface::SegFace();
  printf("MiDeface::FaceIntensityStats()\n");

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

int MiDeface::VoxOutOfBounds(int c, int r, int s)
{
  if(c < cPad) return(1);
  if(c >= invol->width - cPad) return(1);
  if(r < rPad) return(1);
  if(r >= invol->height - rPad) return(1);
  if(s < sPad) return(1);
  if(s >= invol->depth - sPad) return(1);
  return(0);
}

int MiDeface::SetDeltaDist(void)
{
  // DeltaDist and dL control how finely the surf2vol sampling is done.
  // It is probably ok if both are MinVoxSize/2
  DeltaDist = MIN(MIN(headmask->xsize,headmask->ysize),headmask->zsize)/3.0;
  dL = DeltaDist/3.0;
  printf("DeltaDist = %g, dL = %g\n",DeltaDist,dL);
  return(0);
}

/*!
  \fn int MiDeface::DistanceBounds(void)
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
int MiDeface::DistanceBounds(void)
{
  printf("MiDeface::DistanceBounds() DeltaDist=%g\n",DeltaDist);

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
      int OutOfBounds = MiDeface::VoxOutOfBounds(c,r,s);
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
      int OutOfBounds = MiDeface::VoxOutOfBounds(c,r,s);
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
int MiDeface::ripple(LABEL *label)
{
  printf("ripple (%g,%g,%g) amp=%g period=%g\n",
	 ripplecenter[0],ripplecenter[1],ripplecenter[2],
	 rippleamp,rippleperiod);

  BasicSpherePVF sph;
  MRIS *surf = tempsurf;
  MATRIX *M = VGtkreg2RAS(&(surf->vg),NULL);
  MATRIX *invM = MatrixInverse(M,NULL);
  MATRIX *tkras = MatrixAlloc(4,1,MATRIX_REAL);
  tkras->rptr[4][1] = 1;
  MATRIX *ras = MatrixAlloc(4,1,MATRIX_REAL);
  ras->rptr[4][1] = 1;

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
    tkras->rptr[1][1] = dx;
    tkras->rptr[2][1] = dy;
    tkras->rptr[3][1] = dz;
    ras = MatrixMultiplyD(M,tkras,ras);
    ras->rptr[1][1] -= ripplecenter[0];
    ras->rptr[2][1] -= ripplecenter[1];
    ras->rptr[3][1] -= ripplecenter[2];
    std::array<double,3> xyz = {ras->rptr[1][1],ras->rptr[2][1],ras->rptr[3][1]};
    std::array<double,3> rtp = sph.XYZ2RTP(xyz);
    rtp[0] = rtp[0] + rippleamp*sin(2*M_PI*rtp[0]*rtp[rippleaxis]/rippleperiod);
    xyz = sph.RTP2XYZ(rtp);
    ras->rptr[1][1] = xyz[0];
    ras->rptr[2][1] = xyz[1];
    ras->rptr[3][1] = xyz[2];
    ras->rptr[1][1] += ripplecenter[0];
    ras->rptr[2][1] += ripplecenter[1];
    ras->rptr[3][1] += ripplecenter[2];
    tkras = MatrixMultiplyD(invM,ras,tkras);
    v->x = tkras->rptr[1][1];
    v->y = tkras->rptr[2][1];
    v->z = tkras->rptr[3][1];
  }
  MatrixFree(&M);
  MatrixFree(&invM);
  MatrixFree(&tkras);
  MatrixFree(&ras);
  return(0);
}

int MiDeface::watermark(LABEL *watermark, double dwatermark)
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

/*!
  \fn MRI *MiDeface::EmbedCode(MRI *input, MRI *output)
  \brief Embeds the mideface code into the first few voxels of a volume
 */
MRI *MiDeface::EmbedCode(MRI *input, MRI *output)
{
  if(input == NULL) {
    output = MRIcopy(input,NULL);
    MRIcopyHeader(input, output);
    MRIcopyPulseParameters(input, output);
  }
  else{
    int err = MRIdimMismatch(input, output, 0);
    if(err){
      printf("ERROR: MiDeface::EmbedCode() dimension mismatch\n");
      return(NULL);
    }
  }
  int codelen = strlen(embedded_code) + strlen(version) + 1;
  if(output->width < codelen){
    printf("ERROR: MiDeface::EmbedCode() dimension not long enough\n");
    return(NULL);
  }
  for(int f=0; f < output->nframes; f++){
    int k=0;
    for(int n=0; n < strlen(embedded_code); n++)
      MRIsetVoxVal(output,k++,0,0,f, embedded_code[n]);
    MRIsetVoxVal(output,k++,0,0,f, '-');
    for(int n=0; n < strlen(version); n++)
      MRIsetVoxVal(output,k++,0,0,f, version[n]);
  }
  return(output);
}

/*!
  \fn int MiDeface::EmbedCodeCheck(const MRI *vol)
  \brief Checks to see if volume has the mideface code in it. The code 
  can start in any corner to handle the situation where the volume
  has been reoriented.
 */
int MiDeface::EmbedCodeCheck(const MRI *vol)
{
  int codelen = strlen(embedded_code) + strlen(version) + 1;
  if(embeddebug) printf("MiDeface::EmbedCodeCheck() %d %d\n",codelen,vol->width);
  if(vol->width < codelen) return(0);

  for(int direction = 0; direction < 2; direction++){
    for(int axis = 0; axis < 3; axis++){
      int ok=1;
      if(embeddebug) printf("%d %d  ",direction,axis);
      for(int n=0; n < strlen(embedded_code); n++){
	int val = 0;
	if(axis==0) {
	  if(direction==0) val = MRIgetVoxVal(vol,n,0,0,0);
	  else             val = MRIgetVoxVal(vol,vol->width-1-n,0,0,0);
	}
	if(axis==1) {
	  if(direction==0) val = MRIgetVoxVal(vol,0,n,0,0);
	  else             val = MRIgetVoxVal(vol,0,vol->height-1-n,0,0);
	}
	if(axis==2) {
	  if(direction==0) val = MRIgetVoxVal(vol,0,0,n,0);
	  else             val = MRIgetVoxVal(vol,0,0,vol->depth-1-n,0);
	}
	if(embeddebug) printf("%c",(char)val);
	if(val != embedded_code[n]) {
	  if(embeddebug) printf("\n");
	  ok = 0;
	  break;
	}
      }
      if(ok){
	if(embeddebug) {printf("\n");fflush(stdout);}
	return(1);
      }
    }
  }
  if(embeddebug) {printf("Match not found\n");fflush(stdout);}
  return(0);
}
