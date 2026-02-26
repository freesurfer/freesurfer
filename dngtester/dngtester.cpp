/**
 * @brief dougs super special test code
 *
 */
/*
 * Original Author: Doug Greve
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
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <float.h>
#include "mrisurf.h"
#include "mrisutils.h"
#include "geodesics.h"
#include "timer.h"
#include "utils.h"
#include "annotation.h"
#include "error.h"
#include "dmatrix.h"
#include "surfgrad.h"
#include "diag.h"
#include "DICOMRead.h"
#include "region.h"
#include "surfcluster.h"
#include "volcluster.h"
#include "mris_sphshapepvf.h"
#include "gcamorph.h"
#include "resample.h"
#include "mrisurf_metricProperties.h"
#include "surfcluster.h"
#include "fsglm.h"
#include "cma.h"
#include "mrishash.h"
#include "randomfields.h"
#include "fmriutils.h"

#include "romp_support.h"
#undef private
#include "tfce.h"

#if 0
int GCAMgeomMatch(VOL_GEOM *vg, GCAM *gcam)
{
  if(vg_isEqual(&gcam->atlas, vg)) return(1);
  if(vg_isEqual(&gcam->image, vg)) return(2);
  if(Gdiag_no){
    printf("GCAMgeomMatch(): not match\n");
    printf("Input ====================\n");
    LTAdumpVolGeom(stdout,vg);
    printf("Atlas ====================\n");
    LTAdumpVolGeom(stdout,&gcam->atlas);
    printf("Image ====================\n");
    LTAdumpVolGeom(stdout,&gcam->image);
    fflush(stdout);
  }
  return(0);
}
#endif

#if 0
MRI *Image2Sphere(MRIS *sph, MRI *imgsrc)
{
  printf("%g %g %g\n",MRIgetVoxVal(imgsrc,10,20,0,0),MRIgetVoxVal(imgsrc,10,20,0,1),MRIgetVoxVal(imgsrc,10,20,0,2));
  MRI *img;
  int nframes = imgsrc->nframes;
  int rgb = 0;
  if(imgsrc->type == MRI_RGB){
    nframes = 3;
    rgb = 1;
    img = MRIallocSequence(imgsrc->width,imgsrc->height,imgsrc->depth,MRI_FLOAT,3);
    for(int c=0; c<img->width; c++){
      for(int r=0; r<img->height; r++){
        for(int s=0; s<img->depth; s++){
          for(int f=0; f<3; f++) {
            double v = MRIgetVoxVal(imgsrc,c,r,s,f);
	    MRIsetVoxVal(img,c,r,s,f,v);
          }
        }
      }
    }
  }
  else img = MRIcopy(imgsrc,NULL);
  //MRI *img = imgsrc;

  printf("w=%d h=%d nframes = %d, type=%d\n",img->width,img->height,img->nframes,img->type);
  img->xsize = (2*M_PI)/(img->width-1);
  img->ysize = (M_PI)/(img->height-1);
  img->x_r = 1;
  img->x_a = 0;
  img->x_s = 0;
  img->y_r = 0;
  img->y_a = 1;
  img->y_s = 0;
  img->z_r = 0;
  img->z_a = 0;
  img->z_s = 1;

  MRIp0ToCRAS(img, 0,0,0);
  MATRIX *ras2vox = img->get_RAS2Vox();

  MATRIX *phitheta = MatrixAlloc(4,1,MATRIX_REAL);
  phitheta->rptr[3][1] = 0;
  phitheta->rptr[4][1] = 1;
  MATRIX *cr = MatrixAlloc(4,1,MATRIX_REAL);

  MRI *map = MRIallocSequence(sph->nvertices,1,1,MRI_FLOAT,nframes);
  for(int vno=0; vno < sph->nvertices; vno++){
    VERTEX *v = &(sph->vertices[vno]);
    double radius = sqrt(v->x*v->x + v->y*v->y + v->z*v->z);
    double theta  = acos(v->z/radius);
    double phi    = atan2(v->y,v->x);
    if(phi < 0) phi += 2*M_PI;
    phitheta->rptr[1][1] = phi;
    phitheta->rptr[2][1] = theta;
    cr = MatrixMultiply(ras2vox,phitheta,cr);
    int c = nint(cr->rptr[1][1]);
    int r = nint(cr->rptr[2][1]);
    if(vno == 43478) printf("vno=%d x=%g;y=%g;z=%g; rad=%g t=%g p=%g c=%d r=%d\n",vno,v->x,v->y,v->z,radius,theta,phi,c,r);
    if(c<0 || c > img->width-1)  continue;
    if(r<0 || r > img->height-1) continue;
    double vsum=0;
    for(int f=0; f < img->nframes; f++) {
      double v = MRIgetVoxVal(img,c,r,0,f);
      if(! rgb) MRIsetVoxVal(map,vno,0,0,f,v);
      vsum += v;
      if(vno == 43478) printf("f=%d  v=%g\n",f,v);
    }//frame
    if(rgb) MRIsetVoxVal(map,vno,0,0,0,vsum/3);
  }// vno
  return(map);
}
#endif
//int mrisComputeThicknessMinimizationTerm(MRI_SURFACE *mris, double l_thick_min, INTEGRATION_PARMS *parms);
int mrisComputeTargetLocationTerm(MRI_SURFACE *mris, double l_location, INTEGRATION_PARMS *parms);

#if 0
int SmoothSurf(MRIS *surf, int nsubiters)
{
  MRIScomputeMetricProperties(surf);

  INTEGRATION_PARMS parms;

  parms.fill_interior = 0 ;
  parms.projection = NO_PROJECTION ;
  parms.tol = .1; //1e-4 ;
  parms.dt = 5 ; // 0.5
  parms.base_dt = parms.dt ;

  parms.l_hinge = 1;
  parms.l_location = 1;
  parms.l_curv = 0.0 ;
  parms.l_intensity = 0.0 ;
  parms.l_tspring = 0.0 ;
  parms.l_nspring = 0.0 ;
  parms.l_spring = 0.0 ;
  parms.l_surf_repulse = 0.0 ;
  parms.l_spring_nzr = 0.0 ;
  parms.l_spring_nzr_len = 0.0 ;
  parms.l_tsmooth = 0;

  parms.niterations = 0 ;
  parms.write_iterations = 0 /*WRITE_ITERATIONS */;
  parms.integration_type = INTEGRATE_MOMENTUM ;
  parms.momentum = 8.0 /*0.8*/ ;
  parms.dt_increase = 1.0 /* DT_INCREASE */;
  parms.dt_decrease = 0.50 /* DT_DECREASE*/ ;
  parms.error_ratio = 50.0 /*ERROR_RATIO */;
  if(parms.momentum < 0.0) parms.momentum = 0.0 ;
  parms.niterations = 30;

  MRISsaveVertexPositions(surf,TARGET_VERTICES);

  if(surf->edges == NULL){
    printf("First pass, creating edges\n");
    MRISedges(surf);
  }

  MRI *involPS = MRIallocFromVolGeom(&(surf->vg), MRI_UCHAR, 1, 1);
  for(int subiter = 0; subiter < nsubiters; subiter++){
    printf("Positioning surface %d ==========================================\n",subiter);
    fflush(stdout);
    double hingecost = MRISedgeAngleCost(surf,parms.l_hinge, 1);
    double loccost   = 0; //mrisComputeRmsDistanceError(surf);
    printf("#@# %d %10.6lf %10.6lf\n",subiter,1000*hingecost,1000*loccost);

    MRISclearGradient(surf);
    mrisComputeTargetLocationTerm(surf, parms.l_location, &parms);
    MRISedgeMetric(surf, 2); //2=do hinge only

    //MRISpositionSurface(surf, involPS, involPS, &parms);
    for(int vno=0; vno < surf->nvertices; vno++){
      VERTEX *v = &(surf->vertices[vno]);
      v->x -= (v->dx*parms.dt);
      v->y -= (v->dy*parms.dt);
      v->z -= (v->dz*parms.dt);

    }
    //char tmpstr[1000];
    //sprintf(tmpstr,"lh.surf.i%03d",subiter);
    //MRISwrite(surf,tmpstr);
  }

  MRIfree(&involPS);
  return(0);
}
#endif

#if 0
MRI *MRIsetEdges(MRI *in, double SetVal, MRI *out)
{
  out = MRIcopy(in,out);
  if(!out) return(NULL);
  MRIcopyPulseParameters(in, out);
  if(in->ct) out->ct = CTABdeepCopy(in->ct);

  int c, r, s, f, k;
  for(k=0; k < 2; k++){
    if(k==0) c = 0;
    if(k==1) c = out->width-1;
    for(r=0; r < out->height; r++){
      for(s=0; s < out->depth; s++){
	for(f=0; f < out->nframes; f++){
	  MRIsetVoxVal(out,c,r,s,f,SetVal);
	}
      }
    }
    if(k==0) r = 0;
    if(k==1) r = out->height-1;
    for(c=0; c < out->width; c++){
      for(s=0; s < out->depth; s++){
	for(f=0; f < out->nframes; f++){
	  MRIsetVoxVal(out,c,r,s,f,SetVal);
	}
      }
    }
    if(k==0) s = 0;
    if(k==1) s = out->depth-1;
    for(c=0; c < out->width; c++){
      for(r=0; r < out->height; r++){
	for(f=0; f < out->nframes; f++){
	  MRIsetVoxVal(out,c,r,s,f,SetVal);
	}
      }
    }
  }

  return(out);
}
#endif

int sub2ind(MRI *mri, int c, int r, int s)
{
  int index = c + r*mri->width + s*mri->width*mri->height;
  return(index);
}

int CountFaces(MRI *mri)
{
  int nfaces = 0;
  for(int c=1; c < mri->width-1; c++){
    for(int r=1; r < mri->height-1; r++){
      for(int s=1; s < mri->depth-1; s++){
	double val = MRIgetVoxVal(mri,c,r,s,0);
	if(val < 0.5) continue;
        for(int dc = -1; dc < 2; dc++) {
          for(int dr = -1; dr < 2; dr++) {
            for(int ds = -1; ds < 2; ds++) {
	      int topo = fabs(dc)+fabs(dr)+fabs(ds);
	      if(topo != 1) continue;
	      double val2 = MRIgetVoxVal(mri,c+dc,r+dr,s+ds,0);
	      if(val2 > 0.5) continue;
	      nfaces++;
	      int eindices[4][4];
	      for(int e=0; e<2; e++){
		int de=0;
		if(e==0) de = +1;
		if(e==1) de = -1;
		if(dc != 0){
		  eindices[e][0]  = sub2ind(mri,c,r,s);
		  eindices[e][1]  = sub2ind(mri,c,r+de,s);
		  eindices[e][2]  = sub2ind(mri,c+dc,r+dr,s+ds);
		  eindices[e][3]  = sub2ind(mri,c+dc,r+dr+de,s+ds);
		  eindices[e+2][0]  = sub2ind(mri,c,r,s);
		  eindices[e+2][1]  = sub2ind(mri,c,r,s+de);
		  eindices[e+2][2]  = sub2ind(mri,c+dc,r+dr,s+ds);
		  eindices[e+2][3]  = sub2ind(mri,c+dc,r+dr,s+ds+de);
		}
		if(dr != 0){
		  eindices[e][0]  = sub2ind(mri,c,r,s);
		  eindices[e][1]  = sub2ind(mri,c+de,r,s);
		  eindices[e][2]  = sub2ind(mri,c+dc,r+dr,s+ds);
		  eindices[e][3]  = sub2ind(mri,c+dc+de,r+dr,s+ds);
		  eindices[e+2][0]  = sub2ind(mri,c,r,s);
		  eindices[e+2][1]  = sub2ind(mri,c,r,s+de);
		  eindices[e+2][2]  = sub2ind(mri,c+dc,r+dr,s+ds);
		  eindices[e+2][3]  = sub2ind(mri,c+dc,r+dr,s+ds+de);
		}
		if(ds != 0){
		  eindices[e][0]  = sub2ind(mri,c,r,s);
		  eindices[e][1]  = sub2ind(mri,c,r+de,s);
		  eindices[e][2]  = sub2ind(mri,c+dc,r+dr,s+ds);
		  eindices[e][3]  = sub2ind(mri,c+dc,r+dr+de,s+ds);
		  eindices[e+2][0]  = sub2ind(mri,c,r,s);
		  eindices[e+2][1]  = sub2ind(mri,c+de,r,s);
		  eindices[e+2][2]  = sub2ind(mri,c+dc,r+dr,s+ds);
		  eindices[e+2][3]  = sub2ind(mri,c+dc+de,r+dr,s+ds);
		}
	      } // e
	      for(int e=0; e < 4; e++){
		qsort(eindices[e], 4, sizeof(int), compare_ints);
		for(int k=0; k<4; k++) printf("%3d ",eindices[e][k]);
		printf("\n");
	      }
	      fflush(stdout);

	    } // ds
	  }
	}

      }
    }
  }
  //printf("nfaces %d\n",nfaces/2);
  return(nfaces);
}

MRI *MRIgetEditMask(MRI *bm, MRI *t1, MRI *mask, double thresh, int nhitsthresh)
{
  MRI *emask = MRIclone(mask,NULL);
  if(mask->ct) emask->ct = CTABdeepCopy(mask->ct);  

  int nhitslh = 0,nhitsrh = 0;
#ifdef HAVE_OPENMP
#pragma omp parallel for reduction(+ : nhitslh, nhitsrh)
#endif
  for(int c=0; c < mask->width; c++){
    for(int r=0; r < mask->height; r++){
      for(int s=0; s < mask->depth; s++){
	int m = MRIgetVoxVal(mask,c,r,s,0);
	if(m < 0.5) continue;
	double bmv = MRIgetVoxVal(bm,c,r,s,0);
	if(bmv == 1){
	  MRIsetVoxVal(emask,c,r,s,0,m);
	  if(m==1)  nhitslh ++;
	  if(m==2)  nhitsrh ++;
	  continue;
	}
	double t1v = MRIgetVoxVal(t1,c,r,s,0);
	if(bmv == 0 && t1v > thresh){
	  MRIsetVoxVal(emask,c,r,s,0,m);
	  if(m==1)  nhitslh ++;
	  if(m==2)  nhitsrh ++;
	  continue;
	}

      }
    }
  }
  printf("nhits  %d  %d\n",nhitslh,nhitsrh);
  if(nhitslh < nhitsthresh || nhitsrh < nhitsthresh ){
    printf("nhits below threshold of %d, exiting\n",nhitsthresh);
    exit(1010);
  }

  return(emask);
}

// takes mri = nvertices-by-height-by-3. Assigns point xyz to the
// [0,1,2] at each height point for the given vertex. Could set up so
// that each height point is an xyz from a projection.
fsPointSet GetPointSet(MRI *mri, int vno)
{
  fsPointSet  ps;
  for(int i=0; i < mri->height; i++){
    fsPointSet::Point p = fsPointSet::Point();
    p.index = i;
    p.value = 0;
    p.x = MRIgetVoxVal(mri,vno,i,0,0);
    p.y = MRIgetVoxVal(mri,vno,i,0,1);
    p.z = MRIgetVoxVal(mri,vno,i,0,2);
    ps.add(p);
  }
  return(ps);
}

// maps the annot into a volume with proj fraction
// vol, surf(s), segids, min, max, delta, reg
MRI *MRISmapAnnot(MRI *vol, MRIS *surf, LTA *reg, MRI *surfseg, std::vector<int> segids, 
		  double dinward, double doutward, double delta, MRI *inseg)
{
  int err=0;

  MRI *outseg=NULL;
  if(inseg) {
    outseg = inseg;
    vol = inseg;
  }
  if(reg){
    MRISsaveVertexPositions(surf, TMP_VERTICES);
    err = MRISltaMultiply(surf, reg);
    if(err) exit(1);
  }

  if(outseg == NULL){
    outseg = MRIallocSequence(vol->width, vol->height, vol->depth, MRI_INT, 1);
    MRIcopyHeader(vol,outseg);
    MRIcopyPulseParameters(vol,outseg);
  }
  err = MRIdimMismatch(vol, outseg,0);
  if(err) {
    printf("ERROR: MRImapAnnot(): dimension mismatch\n");
    return (NULL);
  }

  MATRIX *M = outseg->get_TkregRAS2Vox();
  MATRIX *xyz = MatrixAlloc(4,1,MATRIX_REAL);
  MATRIX *vox=NULL;
  xyz->rptr[4][1] = 1;
  for(int vno=0; vno < surf->nvertices; vno++){
    int segid = MRIgetVoxVal(surfseg,vno,0,0,0);

    VERTEX *v = &(surf->vertices[vno]);
    for(double dist = dinward; dist < doutward; dist += delta){
      xyz->rptr[1][1] = v->x + dist*v->nx;
      xyz->rptr[2][1] = v->y + dist*v->ny;
      xyz->rptr[3][1] = v->z + dist*v->nz;
      vox = MatrixMultiply(M,xyz,vox);
      int c = nint(vox->rptr[1][1]);
      int r = nint(vox->rptr[2][1]);
      int s = nint(vox->rptr[3][1]);
      if(c < 0 || c >= outseg->width)  continue;
      if(r < 0 || r >= outseg->height) continue;
      if(s < 0 || s >= outseg->depth)  continue;
      MRIsetVoxVal(outseg,c,r,s,0,segid+1000);
    }
  }

  if(reg) MRISrestoreVertexPositions(surf, TMP_VERTICES);

  return(outseg);
}

int MRISfaceCentroid(MRIS *surf, int fno, double centroid[3])
{
  FACE *f = &(surf->faces[fno]);
  for(int k=0; k<3; k++) centroid[k] = 0.0;
  for(int nthv = 0; nthv < 3; nthv++){
    int vno = f->v[nthv];
    VERTEX *const v = &(surf->vertices[vno]);
    centroid[0] += v->x;
    centroid[1] += v->y;
    centroid[2] += v->z;
  }
  for(int k=0; k<3; k++) centroid[k] /= 3.0;
  return(0);
}


int CheckFaces(MRIS *surf, int fno1, int fno2, double dthresh, double dotthresh)
{
  //printf("%d %d\n",fno1,fno2);
  double centroid1[3], centroid2[3];
  MRISfaceCentroid(surf,fno1,centroid1);
  MRISfaceCentroid(surf,fno2,centroid2);
  double rmsdist = 0;
  for(int k=0; k<3; k++) {
    double d = centroid1[k]-centroid2[k];
    rmsdist += (d*d);
  }
  rmsdist = sqrt(rmsdist);
  if(rmsdist > dthresh) return(0);

  FACE *f1 = &(surf->faces[fno1]);
  FACE *f2 = &(surf->faces[fno2]);
  double dot = 0;
  for(int k=1;k<=3;k++) dot += (f1->norm->rptr[k][1]*f2->norm->rptr[k][1]);
  if(dot > dotthresh) return(0);

  return(1);

}

int MarkFaces(MRIS *surf, MRI *mask, double dthresh, double dotthresh)
{
  int *facex = (int*)calloc(sizeof(int),surf->nfaces);
  for(int fno=0; fno < surf->nfaces; fno++){
    FACE *f = &surf->faces[fno];
    for(int nthv = 0; nthv < 3; nthv++){
      int vno = f->v[nthv];
      if(MRIgetVoxVal(mask,vno,0,0,0)>0){
	facex[fno] = 1;
	break;
      }
    }
  }
  MRISclearMarks(surf);
  int nhits=0;
  std::vector<std::vector<int>> facepairs;
#ifdef HAVE_OPENMP
#pragma omp parallel for reduction(+ : nhits)
#endif
  for(int fno1=0; fno1 < surf->nfaces-1; fno1++){
    if(facex[fno1]) continue;
    //if(fno1%1000==0)  printf("fno1 %d  %d\n",fno1,nhits); fflush(stdout);
    for(int fno2=fno1+1; fno2 < surf->nfaces; fno2++){
      if(facex[fno2]) continue;
      int m = CheckFaces(surf, fno1, fno2, dthresh, dotthresh);
      if(!m) continue;
      //printf("fno1 %d  fno2 %d  m=%d\n",fno1,fno2,m); fflush(stdout);      
      nhits ++;
      FACE *f1 = &surf->faces[fno1];
      for(int nthv = 0; nthv < 3; nthv++){
	int vno = f1->v[nthv];
	VERTEX *v = &surf->vertices[vno];
	v->marked = 1;
      }
      FACE *f2 = &surf->faces[fno2];
      for(int nthv = 0; nthv < 3; nthv++){
	int vno = f2->v[nthv];
	VERTEX *v = &surf->vertices[vno];
	v->marked = 1;
      }
      std::vector<int> fp;
      fp.push_back(fno1);
      fp.push_back(fno2);
      facepairs.push_back(fp);
    }
  }
  printf("npairs %d\n",(int) facepairs.size());
  return(0);
}

int DNGMRISmarkEdge(MRIS *surf, MRI *mask, int metricid, double thresh, int FillHoles)
{
  MRISfaceMetric(surf,0);
  if(surf->edges == NULL) MRISedges(surf);
  MRIScomputeMetricProperties(surf);
  MRISedgeMetric(surf,0);

  MRISclearMarks(surf);
  int nhits = 0;
  for(int edgeno = 0; edgeno < surf->nedges; edgeno++){
    MRI_EDGE *e = &(surf->edges[edgeno]);
    // Skip this edge if any vertex is ripped or masked out
    int skip = 0;
    int nvmax = 4;
    if(metricid == 0) nvmax = 2; // length does not need full hinge
    for(int nthv=0; nthv < nvmax; nthv++){
      int vno = e->vtxno[nthv];
      VERTEX  * const v = &(surf->vertices[vno]);
      if(v->ripflag) skip = 1;
      if(mask && MRIgetVoxVal(mask,vno,0,0,0) < 0.5) skip = 1;
    }
    if(skip) continue;
    int mark = 0, nmark=4;
    switch(metricid){
    case 0: if(e->len       < thresh) {mark=1; nmark=2;} break;
    case 1: if(fabs(e->dot) > thresh) mark=1; break;
    case 2: if(e->angle     > thresh) mark=1; break;
    default:
      printf("ERROR: MRISmaxEdgeStatToOverlay() metricid %d unrecognized\n",metricid);
      return(-1);
    }
    if(!mark) continue;
    for(int nthv=0; nthv < nmark; nthv++){ 
      int vno = e->vtxno[nthv];
      surf->vertices[vno].marked = 1;
    }
    nhits++;
  }

  if(FillHoles){
    // Note: this will not respect the mask
    int nfill = MRISfillHoles(surf, (char*) "marked", NULL, 0.5);
    nhits += nfill;
  }

  return(nhits);
}

int MRISremoveHinge(MRI_SURFACE *mris, double AngleThresh, int FillHoles, int ndil, int nsmoothiters, MRI *mask)
{
  int no_progress_max = 5;
  int kmax = 100;

  printf("MRISremoveHinge(): %g %d %d %d %d %d\n",AngleThresh,FillHoles,ndil,nsmoothiters,no_progress_max,kmax);

  MRISclearMarks(mris);
  int num = MRISmarkEdge(mris, mask, 2, AngleThresh, FillHoles);//2=edge angle metric id

  printf(" Found %d vertices with high angle hinges\n",num);
  if(num == 0) return (NO_ERROR);
  MRISdilateMarked(mris, ndil);
  MRISsaveVertexPositions(mris, TMP2_VERTICES);
  int k = 0, n=0, no_progress = 0, old_num = mris->nvertices;
  while (num > 0) {
    printf("k=%d num=%d =================================\n",k,num);fflush(stdout);
    if(num >= old_num){  // couldn't remove any
      no_progress++;
      printf("step %d with no progress (num=%d, old_num=%d)\n", no_progress, num, old_num); fflush(stdout);
      if(no_progress > no_progress_max) break;
    }
    else no_progress = 0;
    if(k > kmax) break; // don't let it go forever

    printf("  %03d: %d high angle hinges \n", n, num); fflush(stdout);
    MRISnotMarked(mris);  // turn off->on and on->off so soap bubble is correct (marked are fixed)
    MRISsoapBubbleVertexPositions(mris, nsmoothiters);

    MRISclearMarks(mris);
    old_num = num;
    num = MRISmarkEdge(mris, mask, 2, AngleThresh, FillHoles);//2=edge angle metric id
    k++;
  }

  printf("  terminating search with %d high angle hing vertices remaining\n", num);
  return(NO_ERROR);
}

// Look for two vertices very close together
int CheckSurf(MRIS *surf, double thresh)
{
  int nhits = 0;
#ifdef HAVE_OPENMP
#pragma omp parallel for reduction(+ : nhits)
#endif
  for(int vno1=0; vno1 < surf->nvertices-1; vno1++){
    VERTEX *v1 = &surf->vertices[vno1];
    for(int vno2=vno1+1; vno2 < surf->nvertices; vno2++){
      VERTEX *v2 = &surf->vertices[vno2];
      double dx = (v1->x - v2->x);
      double dy = (v1->y - v2->y);
      double dz = (v1->z - v2->z);
      double d2 = dx*dx + dy*dy + dz*dz;
      if(d2 > thresh) continue;
      nhits++;
      printf("%5d %6d %6d  %20.19lf\n",nhits,vno1,vno2,d2);
    }
  }
  printf("nhits = %d\n",nhits);
  return(nhits);
}

double NormalDispersion(MRIS *surf, MRI *seg, int segid)
{

  double nx=0, ny=0, nz=0;
  int nhits=0;
  for(int vtxno=0; vtxno < surf->nvertices; vtxno++){
    if(MRIgetVoxVal(seg,vtxno,0,0,0) != segid) continue;
    nhits++;
    VERTEX *v = &(surf->vertices[vtxno]);
    nx += v->nx;
    ny += v->ny;
    nz += v->nz;
  }
  if(nhits < 3) return(0);
  double mag = sqrt(nx*nx + ny*ny + nz*nz);
  nx /= mag;
  ny /= mag;
  nz /= mag;

  double *dotlist = (double*) calloc(nhits,sizeof(double));
  double dotsum=0;
  nhits=0;
  for(int vtxno=0; vtxno < surf->nvertices; vtxno++){
    if(MRIgetVoxVal(seg,vtxno,0,0,0) != segid) continue;
    VERTEX *v = &(surf->vertices[vtxno]);
    double dot = nx*v->nx + ny*v->ny + nz*v->nz;
    dotlist[nhits] = dot;
    nhits++;
    dotsum += dot;
  }
  double dotmean = dotsum/nhits;

  double s2 = 0;
  for(int n=0; n < nhits; n++){
    double ddiff = (dotlist[n]-dotmean);
    s2 += (ddiff*ddiff);
  }
  double stddev = sqrt(s2/(nhits-1));
  printf("segid %d  %g  %5d\n",segid,stddev,nhits);
  free(dotlist);
  return(stddev);
}

MRI *DoIt(MRIS *surf, int nhops, double k1thresh){
  //SURFHOPLIST **shlarray=NULL;
  MRISsetNeighborhoodSizeAndDist(surf,2); // neighborhood size 2 
  MRIScomputeSecondFundamentalFormDiscrete(surf, 0);

  MRI *dotstd = MRIalloc(surf->nvertices,1,1,MRI_FLOAT);

  for(int vno = 0; vno < surf->nvertices; vno++) {
    SURFHOPLIST *shl = SetSurfHopList(vno, surf, nhops);
    std::vector<int> vtxlist;
    for(int hop=0; hop < nhops; hop++){
      int nper = shl->nperhop[hop];
      for(int n=0; n < nper; n++) vtxlist.push_back(shl->vtxlist[hop][n]);
    }
    SurfHopListFree(&shl);

    double nx=0, ny=0, nz=0;
    for(int n=0; n < vtxlist.size(); n++){
      VERTEX *v = &(surf->vertices[vtxlist[n]]);
      nx += v->nx;
      ny += v->ny;
      nz += v->nz;
    }
    double mag = sqrt(nx*nx + ny*ny + nz*nz);
    nx /= mag;
    ny /= mag;
    nz /= mag;
    double dotsum=0;
    double dotsum2=0;
    for(int n=0; n < vtxlist.size(); n++){
      VERTEX *v = &(surf->vertices[vtxlist[n]]);
      double dot = nx*v->nx + ny*v->ny + nz*v->nz;
      dotsum  += dot;
      dotsum2 += (dot*dot);
    }
    double vdotstd = sum2stddev(dotsum, dotsum2, vtxlist.size());
    MRIsetVoxVal(dotstd,vno,0,0,0,vdotstd);
  }

  return(dotstd);
}

int RemoveLeaks(MRIS *surf, MRI *leakseg, int nmax, int nsmoothiters)
{
  printf("nmax = %d nsm = %d  nv = %d\n",nmax,nsmoothiters,surf->nvertices);
  int nhits=0;
  MRISclearMarks(surf);
  for(int vtxno=0; vtxno < surf->nvertices; vtxno++){
    int segid = MRIgetVoxVal(leakseg,vtxno,0,0,0);
    if(segid == 0 || segid > nmax) surf->vertices[vtxno].marked = 0;
    else surf->vertices[vtxno].marked = 1;
    if(surf->vertices[vtxno].marked) nhits++;
    //printf("%5d %3d %d %6d\n",vtxno,(int)MRIgetVoxVal(leakseg,vtxno,0,0,0),surf->vertices[vtxno].marked,nhits);
  }
  // try to remove holes
  MRISdilateMarked(surf, 1);
  MRISerodeMarked(surf, 1);
  MRISnotMarked(surf);  // turn off->on and on->off so soap bubble is correct (unmarked are fixed)
  printf("soap %d\n",nhits);
  // Does not smooth vertices that are marked, smooth unmarked
  MRISsoapBubbleVertexPositions(surf, nsmoothiters);
  MRISnotMarked(surf);  
  MRI *mri = MRIcopyMRIS(NULL,surf,0,"marked");
  MRIwrite(mri,"dng.marked.mgz");
  //MRISclearMarks(surf);

  return(0);
}

int MRISTfindAdjacentFace(MRIS *surf, int vno1, int vno2, int vno3)
{
  for(int fno=0; fno < surf->nfaces; fno++){
    int nhits = 0;
    for(int n = 0; n < VERTICES_PER_FACE; n++) {
      int cvno = surf->faces[fno].v[n]; // vtxno at this corner
      if(cvno == vno1) nhits++;
      if(cvno == vno2) nhits++;
      if(cvno == vno3) nhits--;
    }
    if(nhits == 2) return(fno);
  }
  return(-1);
}
int MRISTcountVertexNeighbors(MRIS *surf)
{
  // First, set the number to 0
  for(int vtxno=0; vtxno < surf->nvertices; vtxno++)
    surf->vertices_topology[vtxno].vnum=0;
  // Go through each corner of each face and increment the number of
  // neighbors at the corner vertex.
  for(int fno=0; fno < surf->nfaces; fno++){
    for(int n = 0; n < VERTICES_PER_FACE; n++) {
      int cvno = surf->faces[fno].v[n]; // vtxno at this corner
      VERTEX_TOPOLOGY *vtop = &(surf->vertices_topology[cvno]);
      vtop->vnum++;
    }
  }
  // Update (redundant) num and check vnum
  int err = 0;
  for(int vtxno=0; vtxno < surf->nvertices; vtxno++){
    surf->vertices_topology[vtxno].num = surf->vertices_topology[vtxno].vnum;
    if(surf->vertices_topology[vtxno].vnum < VERTICES_PER_FACE){
      // might be too aggressive, eg, with defects
      printf("WARNING: MRISTcountVertexNeighbors(): vertex %d has %d neighbors (expecting >=3)\n",
	     vtxno,surf->vertices_topology[vtxno].vnum);
      err++;
    }
  }

  return(err);
}
int MRISTbuildVertexTopology(MRIS *surf, int vno, int fno0)
{
  // Assumes MRISTcountVertexNeighbors() has been run so that vnum is valid
  VERTEX_TOPOLOGY *vtop = &(surf->vertices_topology[vno]);
  if(vtop->v) free(vtop->v);
  if(vtop->f) free(vtop->f);
  vtop->v = (int*)calloc(sizeof(int),vtop->vnum);
  vtop->f = (int*)calloc(sizeof(int),vtop->vnum);
  int fno = fno0;
  for(int n=0; n < vtop->vnum; n++){
    FACE *f = &(surf->faces[fno]);
    int c;
    for(c = 0; c < VERTICES_PER_FACE; c++) {
      int cvno = f->v[c];
      if(cvno == vno) break;
    }
    // Get the next corner index
    c++; if(c==VERTICES_PER_FACE) c=0; // wrap-around
    int vno2 = f->v[c];
    vtop->v[n] = vno2;
    vtop->f[n] = fno; // face neighbor

    int vno3;
    // Get the 3rd vertex in this face
    c++;  if(c==VERTICES_PER_FACE) c=0; // wrap-around
    vno3 = f->v[c];
    fno = MRISTfindAdjacentFace(surf, vno, vno3, vno2);
  }
  return(0);
}
int MRISTbuildTopology0(MRIS *surf)
{
  printf("   buildTopo\n");fflush(stdout);
  printf("     Counting vertices\n");fflush(stdout);
  int err = MRISTcountVertexNeighbors(surf);
  if(err) return(1);

  printf("     Generating face list\n");fflush(stdout);
  std::vector<int> facelist(surf->nvertices);
  for(int faceno=0; faceno < surf->nfaces; faceno++){
    FACE *f = &(surf->faces[faceno]);
    for(int c = 0; c < VERTICES_PER_FACE; c++) {
      int cvno = f->v[c];
      facelist[cvno] = faceno;
    }
  }

  printf("     Generating vertex topo\n");fflush(stdout);
  for(int vno=0; vno < surf->nvertices; vno++)
    MRISTbuildVertexTopology(surf, vno, facelist[vno]);

  printf("     Done building\n");fflush(stdout);
  return(0);
}

int FacesSame(MRIS *surf, int faceno1, int faceno2){
  FACE *f1 = &(surf->faces[faceno1]);
  FACE *f2 = &(surf->faces[faceno2]);
  int nmatch = 0;
  for(int c=0; c < VERTICES_PER_FACE; c++){ 
    for(int k=0; k < VERTICES_PER_FACE; k++){
      if(f1->v[c] == f2->v[k]) nmatch++;
    }
  }
  if(nmatch == 3) return(1);
  return(0);
}

int CheckFaceRep(MRIS *surf, int vno){
  std::vector<int> facelist;
  for(int fno=0; fno < surf->nfaces; fno++){
    for(int c = 0; c < VERTICES_PER_FACE; c++) {
      if(surf->faces[fno].v[c] == vno) facelist.push_back(fno);
    }
  }
  printf("CheckFaceRep(): vno=%d nfaces = %d\n",vno,(int)facelist.size());
  for(int nthf = 0; nthf < facelist.size()-1; nthf++){
    int nfno = facelist[nthf];
    FACE *nf = &(surf->faces[nfno]);
    for(int mthf = nthf+1; mthf < facelist.size(); mthf++){
      int mfno = facelist[mthf];
      FACE *mf = &(surf->faces[mfno]);
      int nmatch=0;
      for(int c=0; c < VERTICES_PER_FACE; c++){ 
	for(int k=0; k < VERTICES_PER_FACE; k++){
	  if(nf->v[c] == mf->v[k]) nmatch++;
	}
      }
      if(nmatch == 3) {
	printf(" CheckFaceRep(): vno = %d  f1 = %d f2 = %d\n",vno,nfno,mfno);
	return(1);
      }
    }
  }
  return(0);
}
int MRISconsistency(MRIS *surf)
{
  // Add check for repeat vtxs and/or faces
  int inconsistent = 0;

  for(int vtxno=0; vtxno < surf->nvertices; vtxno++){
    VERTEX_TOPOLOGY *vtop = &(surf->vertices_topology[vtxno]);
    // number of vertex neighbs must equal number of face neighbors
    if(vtop->vnum != vtop->num) {
      printf("vtx %d vnum = %d num = %d\n",vtxno, vtop->vnum,vtop->num);
      MRISprintVertexInfo(stdout, surf, vtxno);
      continue;
    }
    // For this vertex, visit all the neighbors
    for(int n=0; n < vtop->vnum; n++) {
      int bvtxno = vtop->v[n];
      if(bvtxno < 0 || bvtxno >= surf->nvertices){
	inconsistent++;
	printf("vtx %d nbr %d bvtx %d out of range\n",vtxno,n,bvtxno);
	MRISprintVertexInfo(stdout, surf, vtxno);
	continue;
      }
      // Check neighbors of neighbor vertex to make sure vtxno is one of them
      VERTEX_TOPOLOGY *btop = &(surf->vertices_topology[bvtxno]);
      int ok = 0;
      for(int m=0; m < btop->vnum; m++){
	if(btop->v[m] == vtxno){
	  ok = 1;
	  break;
	}
      }
      if(!ok){
	inconsistent++;
	printf("vtx %d nbr %d bvtxno  %d not a neighbor\n",vtxno,n,bvtxno);
	MRISprintVertexInfo(stdout, surf, vtxno);
      }
      // Check corners of the face to make sure that vtxno is one of them
      int bfaceno = vtop->f[n];
      if(bfaceno < 0 || bfaceno >= surf->nfaces){
	inconsistent++;
	printf("vtx %d nbr %d bfaceno %d out of range\n",vtxno,n,bfaceno);
	MRISprintVertexInfo(stdout, surf, vtxno);
	continue;
      }
      FACE *f = &surf->faces[bfaceno];
      for(int m=0; m < 3; m++){
	if(f->v[m] == vtxno){
	  ok = 1;
	  break;
	}
      }
      if(!ok){
	inconsistent++;
	printf("vtx %d nbr %d bfaceno %d not a neighbor\n",vtxno,n,bfaceno);
      }
    } // neighbor 
  }
  printf("  found %d inconsistencies\n",inconsistent);
  return(inconsistent);
}

int CheckFaces(MRIS *mris)
{
  int nerr = 0;
  for(int fno = 0; fno < mris->nfaces; fno++) {
    FACE *f = &mris->faces[fno];
    for(int n = 0; n < VERTICES_PER_FACE; n++) {
      if(f->v[n] >= mris->nvertices || f->v[n] < 0){
        printf("  f[%d]->v[%d] = %d >= %d - out of range!\n", fno, n, f->v[n],mris->nvertices);
	nerr++;
      }
    }
    // Make sure all vertices are unique
    for(int c1 = 0; c1 < 2; c1++){
      for(int c2 = c1+1; c2 < 3; c2++){
	if(f->v[c1] == f->v[c2] ){
	  printf("  face %d has non-unique vertices %d %d %d\n",fno,f->v[0],f->v[1],f->v[2]);
	  nerr++;
	}
      }
    }
  }
  printf(" CheckFaces nerr = %d\n",nerr);
  return(nerr);
}

int MRISTbuildTopology(MRIS *surf)
{
  printf("    buildTopo --------------\n");fflush(stdout);

  printf("     Checking Faces\n");fflush(stdout);
  CheckFaces(surf);

  printf("     Counting vertex neighbors\n");fflush(stdout);
  std::vector<int> nnbrs(surf->nvertices);
  for(int fno=0; fno < surf->nfaces; fno++){
    for(int c = 0; c < VERTICES_PER_FACE; c++) {
      int cvno = surf->faces[fno].v[c]; // vtxno at this corner
      nnbrs[cvno]++;
    }
  }
  int same = FacesSame(surf, 52044, 67404);
  if(same) printf("  same\n");

  printf("     Allocating vertex topos\n");fflush(stdout);
  for(int vno=0; vno < surf->nvertices; vno++){
    VERTEX_TOPOLOGY *vtop = &(surf->vertices_topology[vno]);
    if(vtop->v) free(vtop->v);
    if(vtop->f) free(vtop->f);
    vtop->v = (int*)calloc(sizeof(int),nnbrs[vno]);
    vtop->f = (int*)calloc(sizeof(int),nnbrs[vno]);
    vtop->vnum = 0;
    vtop->num  = 0;
    for(int n=0; n < nnbrs[vno]; n++) vtop->v[n] = -1; 
  }

  printf("     Generating vertex topo\n");fflush(stdout);
  std::vector<std::vector<int>> vtxnbrlist(surf->nvertices);
  for(int faceno=0; faceno < surf->nfaces; faceno++){
    FACE *f = &(surf->faces[faceno]);
    for(int c = 0; c < VERTICES_PER_FACE; c++) {
      int cvno = f->v[c];
      VERTEX_TOPOLOGY *cvtop = &(surf->vertices_topology[cvno]);
      // Add this face to the list of faces for this vertex
      cvtop->f[cvtop->num] = faceno;
      cvtop->num++;

      if(cvno == 72836){
	printf(" TTT faceno=%d corners: %d %d %d c=%d cvno=%d vnum=%d nbrs: ",faceno,f->v[0],f->v[1],f->v[2],c,cvno,cvtop->vnum);
	for(int k=0; k < cvtop->vnum; k++) printf("%d ",cvtop->v[k]);
	printf("\n");
      }

      // For this corner vertex, are the other two corners already neighbors?
      int ctry = c;
      for(int k=0; k < 2; k++){
	ctry++;
	if(ctry == VERTICES_PER_FACE) ctry = 0;
	int ctryvno = f->v[ctry];
	vtxnbrlist[cvno].push_back(ctryvno);
	int hit = 0;
	for(int n=0; n < cvtop->vnum; n++){
	  if(cvtop->v[n] == ctryvno){
	    hit = 1;
	    break;
	  }
	}
	if(!hit){
	  cvtop->v[cvtop->vnum] = ctryvno;
	  cvtop->vnum++;
	}
      }

    } // face corner
  } // face

  printf("     Checking vertex topo\n");fflush(stdout);
  int err = 0;
  for(int vno=0; vno < surf->nvertices; vno++){
    VERTEX_TOPOLOGY *vtop = &(surf->vertices_topology[vno]);
    if(vtop->vnum != vtop->num){
      printf("       vtx %5d  nv=%d nf=%d nnbrs=%d==================\n",vno,vtop->vnum,vtop->num,nnbrs[vno]);
      MRISprintVertexInfo(stdout, surf, vno);
      CheckFaceRep(surf, vno);
      MRISconsistency(surf);
      for(int k=0; k < vtxnbrlist[vno].size(); k++) printf("  k=%d  %d\n",k,vtxnbrlist[vno][k]);
      printf("     Checking Faces B\n");
      CheckFaces(surf);
      printf("     Face neighbors\n");
      for(int n=0; n < vtop->num; n++){
	int fno = vtop->f[n];
	FACE *f = &(surf->faces[fno]);
	printf("  %d %5d ",n,fno);
	for(int c=0; c < 3; c++) printf("%d ",f->v[c]);
	printf("\n");
      }
      fflush(stdout);
      err++;
    }
  }

  printf("     Done building topo nerr=%d\n",err);fflush(stdout);
  if(err) {
    MRISwrite(surf,"lh.remesh");
    exit(1);
  }
  return(err);
}


int ShareTwo(MRIS *surf,int kvtxno, int rvtxno)
{
  VERTEX_TOPOLOGY *ktop = &(surf->vertices_topology[kvtxno]);
  VERTEX_TOPOLOGY *rtop = &(surf->vertices_topology[rvtxno]);

  // Go through each face neighbor around k
  for(int nknf = 0; nknf < ktop->vnum; nknf++){
    int knfno = ktop->f[nknf];
    FACE *knf = &(surf->faces[knfno]);
    // Go through each face neighbor around 4
    for(int nrnf = 0; nrnf < rtop->vnum; nrnf++){
      int rnfno = rtop->f[nrnf];
      FACE *rnf = &(surf->faces[rnfno]);
      // For these two triangles, count how many vertices they
      // share that are not kvtxno or rvtxno
      int nhits=0;
      for(int ck = 0; ck < VERTICES_PER_FACE; ck++){
	int ckvno = knf->v[ck];
	if(ckvno == kvtxno) continue;
	for(int cr = 0; cr < VERTICES_PER_FACE; cr++){
	  int crvno = rnf->v[cr];
	  if(crvno == rvtxno) continue;
	  if(ckvno == crvno) nhits++;
	}
      }
      if(nhits == 2) return(1);
    }
  }
  return(0);
}


int MRIScollapseEdge(MRIS *surf, int kvtxno, int rvtxno, int rfaceno[2])
{

  int sh2 = ShareTwo(surf,kvtxno,rvtxno);
  if(sh2) return(2);

  VERTEX *kvtx = &(surf->vertices[kvtxno]);
  VERTEX *rvtx = &(surf->vertices[rvtxno]);
  VERTEX_TOPOLOGY *rtop = &(surf->vertices_topology[rvtxno]);

  int debug = 0;

  int ok = 0;
  int nthrface=0;
  for(int nthnbr = 0; nthnbr < rtop->vnum; nthnbr++){
    int bvtxno = rtop->v[nthnbr];
    if(bvtxno == kvtxno) ok = 1; // assure that r and k are neighbors
    int bfaceno = rtop->f[nthnbr];
    FACE *bface = &(surf->faces[bfaceno]);
    int removethisface = 0;
    if(debug) printf("nthnbr %d ----------------------\n",nthnbr);
    for(int n=0; n < VERTICES_PER_FACE; n++){
      if(debug) printf("  rvtxno=%d kvtxno=%d nthnbr=%d c=%d cvtxno=%d bfaceno=%d ",
		       rvtxno,kvtxno,nthnbr,n,bface->v[n],bfaceno);
      if(bface->v[n] == kvtxno){
	// This is a face that has both kvtx and rvtx, which
	// means that it will be removed in the collapse
	rfaceno[nthrface] = bfaceno;
	nthrface++;
	if(debug) printf("   removing face %d\n",nthrface);
	removethisface = 1;
	break;
      }
      else if(debug) printf("\n");
    }
    if(removethisface) {
      if(debug) printf("\n");
      int same = FacesSame(surf, 52044, 67404);
      if(same) printf("  Collapse1 %d bfaceno=%d same k=%d r=%d  rfs %d %d %d\n",nthnbr,bfaceno,kvtxno,rvtxno,nthrface,rfaceno[0],rfaceno[1]);
      continue;
    }
    if(debug) printf("  not removing face\n");
    // If it gets here, then this is a face that is not shared with kvtx, so
    // change rvtxno to kvtxno to remove rvtxno
    for(int n=0; n < VERTICES_PER_FACE; n++){
      if(debug) printf("  rvtxno=%d kvtxno=%d nthnbr=%d c=%d cvtxno=%d bfaceno=%d ",
		       rvtxno,kvtxno,nthnbr,n,bface->v[n],bfaceno);
      if(bface->v[n] == rvtxno){
	bface->v[n] = kvtxno;
	if(debug) printf(" setting bvtxno to kvtxno\n");
	break;
      }
      else {
	if(debug) printf(" NOT setting bvtxno to kvtxno\n");
      }
    }
    int same = FacesSame(surf, 52044, 67404);
    if(same) printf("  Collapse1 %d bfaceno=%d same k=%d r=%d  rfs %d %d %d\n",nthnbr,bfaceno,kvtxno,rvtxno,nthrface,rfaceno[0],rfaceno[1]);
  }
  if(debug) fflush(stdout);
  if(!ok){
    printf("ERROR: MRIScollapseEdge(): kvtxno=%d and rvtxno=%d are not neighbors\n",kvtxno,rvtxno);
    MRISprintVertexInfo(stdout, surf, kvtxno);
    MRISprintVertexInfo(stdout, surf, rvtxno);
    return(1);
  }
  if(nthrface != 2){
    printf("ERROR: MRIScollapseEdge(): kvtxno=%d and rvtxno=%d nthrface = %d\n",kvtxno,rvtxno,nthrface);
    return(1);
  }

  // Set the kept vertex xyz to the average of the two
  kvtx->x = (kvtx->x + rvtx->x)/2;
  kvtx->y = (kvtx->y + rvtx->y)/2;
  kvtx->z = (kvtx->z + rvtx->z)/2;

  return(0);
}

int MRISremoveVertex(MRIS *surf, int rvtxno)
{
  // last vertex in the list
  int qvtxno = surf->nvertices-1;

  // if removing the last vertex, then just dec the number of vertices
  if(rvtxno == qvtxno){
    surf->nvertices--;
    return(0);
  }

  // Vertex to be removed
  VERTEX *rvtx = &(surf->vertices[rvtxno]);
  // Replace the vertex in the list with the last vertex
  VERTEX *qvtx = &(surf->vertices[qvtxno]);
  rvtx->x = qvtx->x;
  rvtx->y = qvtx->y;
  rvtx->z = qvtx->z;

  // Go through all the triangles and replace instances of qvtxno with
  // rvtxno.  This is not so efficient, but can't rely on the vertex
  // neighbors being up-to-date.
  int nhits = 0;
  for(int n=0; n < surf->nfaces; n++){
    FACE *f = &(surf->faces[n]);
    for(int c=0; c < VERTICES_PER_FACE; c++){
      if(f->v[c] == qvtxno){
	f->v[c] = rvtxno;
	nhits ++;
	break;
      }
    }
  }
  surf->nvertices--;
  return(0);
}

int MRISremoveFace(MRIS *surf, int rfaceno)
{

  // last face in the list
  int qfaceno = surf->nfaces-1;

  // if removing the last face, then just dec the number of faces
  if(qfaceno == rfaceno){
    surf->nfaces--;
    return(0);
  }

  FACE *rface = &(surf->faces[rfaceno]);
  FACE *qface = &(surf->faces[qfaceno]);

  // copy face info from last face to rface
  //memcpy(rface,qface,sizeof(FACE));
  // rface and qface now have the same info
  for(int c=0; c < VERTICES_PER_FACE; c++) rface->v[c]=qface->v[c];

  // Don't have to update the vertex topologies because that will
  // be done by MRISTbuildTopology()
  surf->nfaces--;
  return(0);
}

int remesh(MRIS *surf, int nmax)
{
  MRISfaceMetric(surf,0);
  if(surf->edges == NULL) MRISedges(surf);
  MRIScomputeMetricProperties(surf);
  MRISedgeMetric(surf,0);

  printf(" remesh nv=%d nf=%d ne=%d euler=%d\n",
	 surf->nvertices,surf->nfaces,surf->nedges,surf->nvertices+surf->nfaces-surf->nedges);

  // Make a list of hinge angles at each edge
  std::vector< std::pair<float,int> > hingang;
  for(int edgeno = 0; edgeno < surf->nedges; edgeno++){
    MRI_EDGE *e = &(surf->edges[edgeno]);
    //pair<float,int> a = [e->angle,edgno];
    hingang.push_back(std::make_pair(e->angle,edgeno));
  }
  // Sort lowest to highest
  sort(hingang.begin(),hingang.end());

  std::vector<int> edgenolist, rvtxlist;
  for(int nthedge = 0; nthedge < nmax; nthedge++){
    int edgeno = hingang[nthedge].second;
    MRI_EDGE *e = &(surf->edges[edgeno]);
    int skip = 0;
    for(int n=0; n < 2; n++){
      // Check whether this vertex is already in the remove list
      int vtxno = e->vtxno[n];
      if(std::find(rvtxlist.begin(), rvtxlist.end(), vtxno) != rvtxlist.end()){
	skip = 1;
	break;
      }
      // Check whether a neighbor of this vertex is already in the remove list
      VERTEX_TOPOLOGY *btop = &(surf->vertices_topology[vtxno]);
      for(int nthnbr = 0; nthnbr < btop->vnum; nthnbr++){
	int bvtxno = btop->v[nthnbr];
	if(std::find(rvtxlist.begin(), rvtxlist.end(), bvtxno) != rvtxlist.end()){
	  skip = 1;
	  break;
	}
      }
    }
    if(skip) continue;
    rvtxlist.push_back(e->vtxno[1]);
    edgenolist.push_back(edgeno);
  }

  // Now go through all the edges
  printf(" Collapsing %d edges\n",(int)edgenolist.size()); fflush(stdout);
  rvtxlist.clear();
  std::vector<int> rfacelist;
  for(int nthedge = 0; nthedge < edgenolist.size(); nthedge++){
    int edgeno = edgenolist[nthedge];
    MRI_EDGE *e = &(surf->edges[edgeno]);
    int kvtxno = e->vtxno[0];
    int rvtxno = e->vtxno[1];
    int rfaceno[2];
    rfaceno[0] = -1;
    rfaceno[1] = -1;
    int err = MRIScollapseEdge(surf, kvtxno, rvtxno, rfaceno);
    if(err == 2) continue;
    rvtxlist.push_back(rvtxno);
    rfacelist.push_back(rfaceno[0]);
    rfacelist.push_back(rfaceno[1]);
    int same = FacesSame(surf, 52044, 67404);
    if(same) {
      printf("  nthedge = %d same k=%d r=%d  rfs %d %d\n",nthedge,kvtxno,rvtxno,rfaceno[0],rfaceno[1]);
      exit(1);
    }
  }
  CheckFaces(surf);

  // Free edges so it will be recreated later on with new topology
  free(surf->edges); 
  surf->edges = NULL;

  printf(" Found %d vertices and %d faces to remove\n",(int)rvtxlist.size(),(int)rfacelist.size());

  // have to reverse sort so that largest vtxno/facno is first
  std::sort(rvtxlist.begin(), rvtxlist.end(), std::greater<int>());
  std::sort(rfacelist.begin(), rfacelist.end(), std::greater<int>());

  printf(" Removing vertices\n");
  fflush(stdout);
  for(int nthvtx=0; nthvtx < rvtxlist.size(); nthvtx++){
    MRISremoveVertex(surf,rvtxlist[nthvtx]);
    int same = FacesSame(surf, 52044, 67404);
    if(same) printf("  nthvtx = %d same\n",nthvtx);
  }

  printf(" Removing faces\n"); fflush(stdout);
  for(int nthface=0; nthface < rfacelist.size(); nthface++){
    MRISremoveFace(surf,rfacelist[nthface]);
    int same = FacesSame(surf, 52044, 67404);
    if(same) printf("  nthface = %d same\n",nthface);
  }

  printf(" Building topology ============================\n");fflush(stdout);
  MRISTbuildTopology(surf);
  MRIScomputeMetricProperties(surf);

  printf(" Checking consistency ============================\n");fflush(stdout);
  MRISconsistency(surf);

  // have to rebuild edges
  //printf(" remesh post nv=%d nf=%d ne=%d euler=%d\n",
  //surf->nvertices,surf->nfaces,surf->nedges,surf->nvertices+surf->nfaces-surf->nedges);

  return((int)rvtxlist.size());
}

/*
  double Curv2(MRIS *surf, int vno)

  For about 99% of vertices, this function will generate the same
  result (within .001) as the MRIScomputeSecondFundamentalForm (2FF)
  functions. For those few remaining, the problem voxels appear to be
  in the fundus of a pial sulcus, or places where the surface is
  wonky, or where the matrix is ill-conditioned. The 2FF code has some
  error traps to prevent extreme values. The 2FF does not have a
  "keepself" option so make sure that is turned off when making
  comparisons, although probably a good idea to keep it on in
  general. Another diff between this and the 2FF is that this uses
  MatrixMultiplyD() which uses double precision in the
  multiplications. Before calling this function, make sure to run
  MRISsetNeighborhoodSizeAndDist(surf,N) and
  mrisComputeTangentPlanes(surf). Make sure to #include fsglm.h as
  well as other typical includes. This function should be thread safe.
 */
double Curv2(MRIS *surf, int vno)
{
  VERTEX *vtx = &(surf->vertices[vno]);
  VERTEX_TOPOLOGY *vtop = &(surf->vertices_topology[vno]);

  // Get the total number of vertices in the extended neighborhood.
  // This is set by a call to MRISsetNeighborhoodSizeAndDist(surf,2);
  // where 2, means 2 hops (I think)
  int nnbrs = vtop->vtotal;
  int keepself = 1;
  if(keepself) nnbrs++;

  int debug = 1;
  //if(vno == 138268) debug = 1;
  //if(vno == 156791) debug = 1;

  // Create a matrix of the delta xyz of the neighbors. 
  // Should self be included?
  MATRIX *xyz0 = MatrixAlloc(nnbrs,3,MATRIX_REAL);
  MATRIX *xyz = MatrixAlloc(nnbrs,3,MATRIX_REAL);
  if(keepself){
    xyz->rptr[1][1] = 0;
    xyz->rptr[1][2] = 0;
    xyz->rptr[1][3] = 0;
    xyz0->rptr[1][1] = vtx->x;
    xyz0->rptr[1][2] = vtx->y;
    xyz0->rptr[1][3] = vtx->z;
  }
  for(int i = keepself; i < nnbrs; i++) {
    int vnbno = vtop->v[i];
    VERTEX *vnb = &surf->vertices[vnbno];
    xyz->rptr[i+1][1] = vnb->x - vtx->x; 
    xyz->rptr[i+1][2] = vnb->y - vtx->y; 
    xyz->rptr[i+1][3] = vnb->z - vtx->z; 
    xyz0->rptr[i+1][1] = vnb->x;
    xyz0->rptr[i+1][2] = vnb->y;
    xyz0->rptr[i+1][3] = vnb->z;
  }

  // Create a matrix to rotate the xyz into the norm/tangent coord
  // system vtx->n is the normal vector vtx-e{1,2} are tangential
  // vectors; see mrisComputeTangentPlanes(surf) The tangent vectors
  // only need to be orthog to each other and the normal. Note: in
  // BF's code, the tangent vectors are replaced with the principal
  // directions.
  MATRIX *G = MatrixAlloc(3,3,MATRIX_REAL);
  G->rptr[1][1] = vtx->e1x;  G->rptr[2][1] = vtx->e1y;  G->rptr[3][1] = vtx->e1z;
  G->rptr[1][2] = vtx->e2x;  G->rptr[2][2] = vtx->e2y;  G->rptr[3][2] = vtx->e2z;
  G->rptr[1][3] = vtx->nx;   G->rptr[2][3] = vtx->ny;   G->rptr[3][3] = vtx->nz;

  // Rotate
  MatrixMultiplyD(xyz,G,xyz);

  // Create and alloc a GLM
  GLMMAT *glm = GLMalloc();
  GLMallocX(glm,nnbrs,3);
  GLMallocY(glm);
  // Load the GLM X and y matrices
  for(int i = 0; i < nnbrs; i++) {
    double x = xyz->rptr[i+1][1];
    double y = xyz->rptr[i+1][2];
    double z = xyz->rptr[i+1][3];
    // z = a*x*x + 2*b*x*y + c*y*y
    // same order as BF
    glm->X->rptr[i+1][1] =   x*x;
    glm->X->rptr[i+1][2] = 2*x*y;
    glm->X->rptr[i+1][3] =   y*y;
    glm->y->rptr[i+1][1] =     z;
  }

  // solve beta = inv(X'*X)*X'*y
  GLMxMatrices(glm);
  GLMfit(glm);
  if(glm->ill_cond_flag) {
    printf("vno = %d ill-cond %d\n",vno,nnbrs);
    MatrixFree(&xyz);
    MatrixFree(&G);
    GLMfree(&glm);
    vtx->k1 = -1;
    vtx->k2 = -1;
    return(-1);
  }

  //MatrixPrint(stdout,glm->X);
  //MatrixPrint(stdout,glm->beta);

  // The Hessian matrix is the second deriv, ie,
  // H = [d2z/dx2 d2z/dxdy; d2z/dxdy d2z/dy2] 
  // Where z = a*x*x + 2*b*x*y + c*y*y
  // d2z/dx2 = 2a
  // d2z/dy2 = 2c
  // d2z/dxdy = 2b
  // So H = [2a 2b; 2b 2c];
  double a = 2*glm->beta->rptr[1][1];
  double b = 2*glm->beta->rptr[3][1];
  double c = 2*glm->beta->rptr[2][1];

  // k1 and k2 are the eigen values of the Hessian. For a 2x2 matrix,
  // the eigvals can be computed without decomposition
  double d = sqrt(4*c*c + (a-b)*(a-b));
  vtx->k1 = (a+b+d)/2;
  vtx->k2 = (a+b-d)/2;
  if(fabs(vtx->k1) < fabs(vtx->k2)){
    vtx->k1 = (a+b-d)/2;
    vtx->k2 = (a+b+d)/2;
  }

  double meancurv = (vtx->k1 + vtx->k2)/2; // = (a+b)/2
  //printf("B vno %d %g %g %g\n",vno,vtx->k1,vtx->k2,meancurv);

  if(debug){
    MatrixWriteTxt("G.txt", G);
    MatrixWriteTxt("xyz.txt", xyz0);
    MatrixWriteTxt("X.txt", glm->X);
    MatrixWriteTxt("y.txt", glm->y);
    MatrixWriteTxt("beta.txt", glm->beta);
    printf("  k1=%g k2=%g  km=%g\n",vtx->k1,vtx->k2,meancurv);
  }

  MatrixFree(&xyz);
  MatrixFree(&G);
  GLMfree(&glm);

  return(meancurv);
}

int MRIScomputeSecondFundamentalFormAtVertex2(MRIS *mris, int vno)
{
  int i, n, nbad = 0;
  VERTEX *vertex, *vnb;
  MATRIX *m_U, *m_Ut, *m_tmp1, *m_tmp2, *m_inverse;
  VECTOR *v_z;
  static MATRIX *m_Q, *m_eigen;
  static VECTOR *v_c = NULL, *v_n, *v_e1, *v_e2, *v_yi;
  float k1, k2, evalues[3], a11, a12, a21, a22, cond_no, rsq, k, kmin, kmax;
  double ui, vi;

  VERTEX_TOPOLOGY *vtop = &(mris->vertices_topology[vno]);
  int vnum = vtop->vtotal;

  if (v_c == NULL) {
    v_c = VectorAlloc(3, MATRIX_REAL);
    v_n = VectorAlloc(3, MATRIX_REAL);
    v_e1 = VectorAlloc(3, MATRIX_REAL);
    v_e2 = VectorAlloc(3, MATRIX_REAL);
    v_yi = VectorAlloc(3, MATRIX_REAL);
    m_Q = MatrixAlloc(2, 2, MATRIX_REAL); /* the quadratic form */
    m_eigen = MatrixAlloc(2, 2, MATRIX_REAL);
  }

  vertex = &mris->vertices[vno];
  if (vertex->ripflag) {
    return (ERROR_BADPARM);
  }

  if (vno == 142915) {
    DiagBreak();
  }
  VECTOR_LOAD(v_n, vertex->nx, vertex->ny, vertex->nz);
  VECTOR_LOAD(v_e1, vertex->e1x, vertex->e1y, vertex->e1z);
  VECTOR_LOAD(v_e2, vertex->e2x, vertex->e2y, vertex->e2z);
  // v_n is the normal vector
  // e1 is the 1st principle direction
  // e2 is the 2nd principle direction
  // these three vectors form a coordinate system. 
  // e1/e2 are the tangent plane
  //MatrixPrint(stdout,v_e1);
  //MatrixPrint(stdout,v_e2);
  //MatrixPrint(stdout,v_n);

  if (vnum <= 0) {
    return (ERROR_BADPARM);
  }

  m_U = MatrixAlloc(vnum, 3, MATRIX_REAL);
  v_z = VectorAlloc(vnum, MATRIX_REAL);

  if (vno == Gdiag_no) {
    DiagBreak();
  }

  /* fit a quadratic form to the surface at this vertex */
  // where v_z = A*ui^2 + B*2*ui*vi + C*vi^2
  kmin = 10000.0f;
  kmax = -kmin;
  n = 0;
  //FILE *fp = fopen("xyz.dat","w");
  //FILE *fp2 = fopen("uvz.dat","w");
  for (i = 0; i < vnum; i++) {
    int vnbno = vtop->v[i];
    vnb = &mris->vertices[vnbno];
    //fprintf(fp,"%7.4f %7.4f %7.4f \n",vnb->x-vertex->x,vnb->y-vertex->y,vnb->z-vertex->z);
    if (vnb->ripflag) {
      continue;
    }
    /* Calculate the projection of this vertex onto the local tangent plane */
    // v_yi = the difference from the cener vertex to the neighbor
    VECTOR_LOAD(v_yi, vnb->x - vertex->x, vnb->y - vertex->y, vnb->z - vertex->z);
    ui = V3_DOT(v_yi, v_e1); // distance along the e1 axis
    vi = V3_DOT(v_yi, v_e2); // distance along the e2 axis
    *MATRIX_RELT(m_U, n + 1, 1) = ui * ui;
    *MATRIX_RELT(m_U, n + 1, 2) = 2 * ui * vi;
    *MATRIX_RELT(m_U, n + 1, 3) = vi * vi;
    VECTOR_ELT(v_z, n + 1) = V3_DOT(v_n, v_yi); /* height above TpS */
    //fprintf(fp2,"%7.4f %7.4f %7.4f \n",ui,vi,V3_DOT(v_n, v_yi));

    // distance^2 from center vertex to nbr vertex with in the tangent plane
    rsq = ui * ui + vi * vi; 
    if (!FZERO(rsq)) {
      // This appears to be some curvature measure to use when ill-conditioned
      k = VECTOR_ELT(v_z, n + 1) / rsq;
      if (k > kmax) {
        kmax = k; // The max will be k1 when ill-cond
      }
      if (k < kmin) {
        kmin = k; // The min will be k2 when ill-cond
      }
    }
    n++;
  }
  //fclose(fp);
  //fclose(fp2);
  // v_z = m_U*v_c, v_c = inv(m_U'*m_U)*m_U'*v_z
  //MatrixPrint(stdout,m_U);
  //MatrixWriteTxt("mu.dat",m_U);

  m_Ut = MatrixTranspose(m_U, NULL);        /* Ut */
  m_tmp2 = MatrixMultiply(m_Ut, m_U, NULL); /* Ut U */
  cond_no = MatrixConditionNumber(m_tmp2);
  m_inverse = MatrixSVDInverse(m_tmp2, NULL); /* (Ut U)^-1 */
  if (!m_inverse) /* singular matrix - must be planar?? */
  {
    nbad++;
    evalues[0] = evalues[1] = 0.0;
  }
  else {
    m_tmp1 = MatrixMultiply(m_Ut, v_z, NULL); /* Ut z */
    MatrixMultiply(m_inverse, m_tmp1, v_c);   /* (Ut U)^-1 Ut z */

    /* Build 2x2 sym Hessian matrix for the quadratic formula above. 
       The Hessian is the second derivative matrix, so
       m_Q(1,1) = d2z/du2 = 2*A (A = v_c(1))
       First princ curv k1 will be the first eigen val of Q
       Second princ curv k2 will be the second eigen val of Q
       ev1 = (v_c(1)+v_c(2) + d)/2
       ev2 = (v_c(1)+v_c(2) - d)/2
       d = sqrt(4*v_c(3)^2 + (vc_(1)-vc_(2).^2))
       Mean Curv = (ev1+ev2)/2 = (v_c(1)+v_c(2))/2 (v_c(3) not there!)
     */
    *MATRIX_RELT(m_Q, 1, 1) = 2 * VECTOR_ELT(v_c, 1);
    *MATRIX_RELT(m_Q, 1, 2) = *MATRIX_RELT(m_Q, 2, 1) = 2 * VECTOR_ELT(v_c, 2);
    *MATRIX_RELT(m_Q, 2, 2) = 2 * VECTOR_ELT(v_c, 3);

    //printf("beta0\n");
    //MatrixPrint(stdout,v_c);
    //printf("Q\n");
    //MatrixPrint(stdout,m_Q);

    /* the columns of m_eigen will be the eigenvectors of m_Q */
    if (MatrixEigenSystem(m_Q, evalues, m_eigen) == NULL) {
      nbad++;
      MatrixSVDEigenValues(m_Q, evalues);
      vertex->k1 = k1 = evalues[0];
      vertex->k2 = k2 = evalues[1];
      vertex->K = k1 * k2;
      vertex->H = (k1 + k2) / 2;
      MatrixFree(&m_Ut);
      MatrixFree(&m_tmp2);
      MatrixFree(&m_U);
      VectorFree(&v_z);
      MatrixFree(&m_tmp1);
      MatrixFree(&m_inverse);
      return (ERROR_BADPARM);
    }

    MatrixFree(&m_tmp1);
    MatrixFree(&m_inverse);
  } // end not illcond
  k1 = evalues[0];
  k2 = evalues[1];
  vertex->k1 = k1;
  vertex->k2 = k2;
  //printf("A vno %d  k1 %g k2 %g  %g\n",vno,k1,k2,(k1+k2)/2);

  vertex->K = k1 * k2;
  vertex->H = (k1 + k2) / 2;
  if (vno == Gdiag_no && (Gdiag & DIAG_SHOW))
    fprintf(stdout, "v %d: k1=%2.3f, k2=%2.3f, K=%2.3f, H=%2.3f\n", vno, vertex->k1, vertex->k2, vertex->K, vertex->H);
  if (vertex->K < mris->Kmin) {
    mris->Kmin = vertex->K;
  }
  if (vertex->H < mris->Hmin) {
    mris->Hmin = vertex->H;
  }
  if (vertex->K > mris->Kmax) {
    mris->Kmax = vertex->K;
  }
  if (vertex->H > mris->Hmax) {
    mris->Hmax = vertex->H;
  }
  mris->Ktotal += (double)k1 * (double)k2 * (double)vertex->area;

  /* now update the basis vectors to be the principal directions */
  a11 = *MATRIX_RELT(m_eigen, 1, 1);
  a12 = *MATRIX_RELT(m_eigen, 1, 2);
  a21 = *MATRIX_RELT(m_eigen, 2, 1);
  a22 = *MATRIX_RELT(m_eigen, 2, 2);
  if (V3_LEN(v_e1) < 0.5) {
    DiagBreak();
  }
  vertex->e1x = V3_X(v_e1) * a11 + V3_X(v_e2) * a21;
  vertex->e1y = V3_Y(v_e1) * a11 + V3_Y(v_e2) * a21;
  vertex->e1z = V3_Z(v_e1) * a11 + V3_Z(v_e2) * a21;
  vertex->e2x = V3_X(v_e1) * a12 + V3_X(v_e2) * a22;
  vertex->e2y = V3_Y(v_e1) * a12 + V3_Y(v_e2) * a22;
  vertex->e2z = V3_Z(v_e1) * a12 + V3_Z(v_e2) * a22;
  if (SQR(vertex->e1x) + SQR(vertex->e1y) + SQR(vertex->e1z) < 0.5) {
    DiagBreak();
  }

  MatrixFree(&m_Ut);
  MatrixFree(&m_tmp2);
  MatrixFree(&m_U);
  VectorFree(&v_z);

  if (Gdiag & DIAG_SHOW && (nbad > 0)) {
    fprintf(stdout, "%d ill-conditioned points\n", nbad);
  }
  return (NO_ERROR);
}

MRI *SimAtrophy(MRIS *surf, int vno, int nhops, int navgs)
{

  // Soap Bubble will ignore all marked=1, so mark everything
  // and then unmark the places we want to smooth
  MRISsetMarks(surf,1);

  VERTEX *vtx = &(surf->vertices[vno]);
  vtx->marked = 0;
  
  MRI *mri = MRIalloc(surf->nvertices,1,1,MRI_FLOAT);
  SURFHOPLIST *shl = SetSurfHopList(vno, surf, nhops);
  for(int hop=0; hop < nhops; hop++){
    int nper = shl->nperhop[hop];
    for(int n=0; n < nper; n++) {
      int nvno = shl->vtxlist[hop][n];
      VERTEX *vtx = &(surf->vertices[nvno]);
      vtx->marked = 0;
      vtx->x -= 0.5*vtx->nx;
      vtx->y -= 0.5*vtx->ny;
      vtx->z -= 0.5*vtx->nz;
      MRIsetVoxVal(mri,nvno,0,0,0,1);
    }
  }
  //MRISsoapBubbleVertexPositions(surf, navgs);

  return(mri);
}

MRI *SimAtrophyLabel(MRIS *surf, LABEL *label, double dist)
{

  MRI *mri = MRIalloc(surf->nvertices,1,1,MRI_FLOAT);
  for(int n = 0; n < label->n_points; n++) {
    int vno = label->lv[n].vno;
    VERTEX *vtx = &(surf->vertices[vno]);
    vtx->x -= dist*vtx->nx;
    vtx->y -= dist*vtx->ny;
    vtx->z -= dist*vtx->nz;
    MRIsetVoxVal(mri,vno,0,0,0,1);
  }

  return(mri);
}

MRI *MRIaseg2acqueduct(MRI *aseg)
{
  MATRIX *vox2ras = aseg->get_Vox2RAS(0);
  MATRIX *crs = MatrixAlloc(4,1,MATRIX_REAL);
  crs->rptr[4][1] = 1;
  MATRIX *ras = NULL;

  // LIA
  //col = R->L
  //row = S->I
  //slice = P->A

  // Below uses vox2ras but the rest assumes LIA. It would be good
  // to make it independent of orientation someday.
  // Find the most superior voxel (row) of the 4th vent
  int crssup[3]={0,0,0},nhits=0;
  double smax = -10e10;
  int sant = 0; // most anterior slice
  int rsup = 10e5; // most superior row
  for(int c=0; c < aseg->width; c++){
    for(int r=0; r < aseg->height; r++){
      for(int s=0; s < aseg->depth; s++){
	int segid = MRIgetVoxVal(aseg,c,r,s,0);
	if(segid != 15) continue; //15 = 4th vent
	nhits++;
	if(r < rsup) rsup = r; // most superior row
	if(s > sant) sant = s; // most anterior slice
	crs->rptr[1][1] = c;
	crs->rptr[2][1] = r;
	crs->rptr[3][1] = s;
	ras = MatrixMultiplyD(vox2ras,crs,ras);
	if(smax < ras->rptr[3][1]) {
	  smax = ras->rptr[3][1];
	  crssup[0] = c;
	  crssup[1] = r;
	  crssup[2] = s;
	}
      }
    }
  }
  int spost = round(sant - 6/aseg->zsize);//6mm in paper, but not in matlab
  printf("nhits = %d, crs = %d %d %d  rsup=%d sant=%d spost=%d %g\n",
	 nhits,crssup[0],crssup[1],crssup[2],rsup,sant,spost,smax);

  MRI *aq = MRIalloc(aseg->width,aseg->width,aseg->depth,MRI_INT);
  MRIcopyHeader(aseg,aq);
  MRIcopyPulseParameters(aseg,aq);
  MRIcopy(aseg,aq);

  // Starting at the most superior row. Examine each slice to see
  // where the acqueduct is not totally surrounded by brainstem.
  int r = crssup[1];
  int done = 0;
  while(!done){
    // Get single slice binarizerd to non-brainstem. The 4thvent/aqd
    // will show up as an island
    MRI *slice = MRIalloc(aseg->width,1,aseg->depth,MRI_INT);
    int nhits = 0;
    for(int c=0; c < aseg->width; c++){
      for(int s=0; s < aseg->depth; s++){
	int segid = MRIgetVoxVal(aseg,c,r,s,0);
	int val=0;
	if(segid != 16) {//not brainstem
	  val = 1;
	  nhits++;
	}
	MRIsetVoxVal(slice,c,0,s,0,val);
      }
    }
    if(nhits==0){ // No BS voxels, probably does not happen
      MRIfree(&slice);
      break;
    }
    // Clusterize. 
    int nClusters;
    VOLCLUSTER **ClusterList = clustGetClusters(slice, 0, 0.5, 0, 0, 0, NULL, &nClusters, NULL);
    if(nClusters > 1){
      // Look through all the clusters and relable the 4th vent voxels as acqueduct
      for(int cno = 1; cno < nClusters; cno++){
	VOLCLUSTER *vc = ClusterList[cno];
	for(int n = 0; n < vc->nmembers; n++){
	  int segid = MRIgetVoxVal(aseg,vc->col[n],r,vc->slc[n],0);
	  if(segid != 15) continue;
	  MRIsetVoxVal(aq, vc->col[n], r, vc->slc[n], 0, 540); // acqueduct
	}
      }
    }
    else done = 1;
    MRIfree(&slice);
    printf("r=%d nhits=%d nc=%d\n",r,nhits,nClusters);
    clustFreeClusterList(&ClusterList, nClusters);
    r++;
  }
  MRIwrite(aq,"aq.mgz");
  int rinf = r;
  printf("rsup = %d rinf = %d\n",rsup,rinf);

  // DR search space is c=0:width,r=rsup:rinf,s=spost:sant
  // MR search space is c=0:width,r=rinf+1:limit,s=spost:sant

  MatrixFree(&vox2ras);
  return(aq);
}

int SetNextMax(MRI *seg, int segid0, int cmin, int cmax, int rmin, int rmax, int smin, int smax, MRI *stat)
{
  int Ncmax=-1,Nrmax=-1,Nsmax=-1;
  double statmax=-10e10;
  for(int c=cmin; c <= cmax; c++){
    for(int r=rmin; r <= rmax; r++){
      for(int s=smin; s < smax; s++){
	int segid = MRIgetVoxVal(seg,c,r,s,0);
	if(segid != segid0) continue;
	for(int dc = -1; dc < 2; dc ++){
	  for(int dr = -1; dr < 2; dr ++){
	    for(int ds = -1; ds < 2; ds ++){
	      segid = MRIgetVoxVal(seg,c+dc,r+dr,s+ds,0);
	      if(segid == segid0) continue;
	      if(segid != 16) continue; // not in brainstem
	      if(segid == 540) continue; // not in acqueduct
	      double v = MRIgetVoxVal(stat,c+dc,r+dr,s+ds,0);
	      if(statmax < v){
		statmax = v;
		Ncmax = c+dc;
		Nrmax = r+dr;
		Nsmax = s+ds;
	      }
	    }
	  }
	}
      }
    }
  }
  MRIsetVoxVal(seg,Ncmax,Nrmax,Nsmax,0,segid0);
  return(0);
}
MRI *MRIsegRaphe(MRI *asegaq, MRI *dasb)
{

  double drmax=0, mrmax=0;
  int drseed[3]={0,0,0}, mrseed[3]={0,0,0};
  for(int c=0; c < asegaq->width; c++){
    for(int r=0; r < asegaq->height; r++){
      for(int s=0; s < asegaq->depth; s++){
	int segid = MRIgetVoxVal(asegaq,c,r,s,0);
	if(segid != 15 && segid != 540) continue;
	for(int ss=s; ss < 256; ss++){
	  int segid2 = MRIgetVoxVal(asegaq,c,r,ss,0);
	  if(segid2 != 16) continue; // must be in brainstem
	  double v = MRIgetVoxVal(dasb,c,r,ss,0);
	  if(segid == 540 && drmax < v) {
	    drmax = v;
	    drseed[0] = c;
	    drseed[1] = r;
	    drseed[2] = ss;
	  }
	  if(segid == 15 && mrmax < v) {
	    mrmax = v;
	    mrseed[0] = c;
	    mrseed[1] = r;
	    mrseed[2] = ss;
	  }
	} // ss
      } // s
    } // r
  } // c
  MRIsetVoxVal(asegaq,drseed[0],drseed[1],drseed[2],0,118);
  MRIsetVoxVal(asegaq,mrseed[0],mrseed[1],mrseed[2],0,119);
  printf("DR %d %d %d  %g\n",drseed[0],drseed[1],drseed[2],drmax);
  printf("MR %d %d %d  %g\n",mrseed[0],mrseed[1],mrseed[2],mrmax);
  MRIwrite(asegaq,"raphe1.mgz");

#if 0
  for(int n = 0; n<115; n++) {
    printf("DR %d\n",n);
    SetNextMax(asegaq, 118, dasb);
  }
  for(int n = 0; n<64; n++){
    SetNextMax(asegaq, 119, dasb);
    printf("MR %d\n",n);
  }
  MRIwrite(asegaq,"raphe2.mgz");
#endif
  return(asegaq);
}

class RapheSeg 
{
public:
  MRI *aseg0=NULL, *dasb0=NULL;
  int AqInfSlice=1000, AqSupSlice=-1;
  int AqPostRow=1000, AqAntRow=-1;
  int BsLeftCol=1000, BsRightCol=-1;
  MRI *dr0=NULL, *mr0=NULL;
  MRI *doit(void);
  int RapheDice(MRI *aseg);
};
MRI *RapheSeg::doit(void)
{
  double drtargvol = 150.3; // vs 115 mm3
  double mrtargvol =  99.6; // vs  64 mm3
  double voxvol = aseg0->xsize*aseg0->ysize*aseg0->zsize;
  double DRAPmm = 6.0; // VB used 6
  double MRAPmm = 11.0; // VB used 6

  // Make a copy so that they can be reoriented if needed
  MRI *aseg = MRIcopy(aseg0,NULL);
  MRI *dasb = MRIcopy(dasb0,NULL);

  // Reorient to RAS
  char ostr[5];
  ostr[4] = '\0';
  MRIdircosToOrientationString(aseg,ostr);
  printf("Aseg Orientation %s\n",ostr);
  if(strcmp(ostr,"RAS")){
    printf("Reorientating Aseg to RAS\n");
    MRI *rastemplate = MRIcopy(aseg,NULL);
    MRIorientationStringToDircos(rastemplate, "RAS");
    MRI *tmp = MRIresample(aseg, rastemplate, SAMPLE_NEAREST);
    MRIfree(&aseg);
    aseg = tmp;
    MRIfree(&rastemplate);
  } 
  MRIdircosToOrientationString(dasb,ostr);
  printf("DASB Orientation %s\n",ostr);
  if(strcmp(ostr,"RAS")){
    printf("Reorientating DASB to RAS\n");
    MRI *rastemplate = MRIcopy(dasb,NULL);
    MRIorientationStringToDircos(rastemplate, "RAS");
    MRI *tmp = MRIresample(dasb, rastemplate, SAMPLE_NEAREST);
    MRIfree(&dasb);
    dasb = tmp;
    MRIfree(&rastemplate);
  } 
  int DRAPvox = ceil(DRAPmm/aseg->ysize);
  int MRAPvox = ceil(MRAPmm/aseg->ysize);

  int err = MRIdimMismatch(aseg, dasb, 0);
  if(err) {
    printf("ERROR: MRIrapheSeg(): dimension mismatch\n");
    return (NULL);
  }

  // Go through the sagittal plane to find the largest row 
  std::vector<int> supslice;
  for(int r=0; r < aseg->height; r++){
    int aqfound=0;
    for(int s=0; s < aseg->depth;  s++){
      // At this row/slice, search the cols for aqueduct
      for(int c=0; c < aseg->width; c++){
	int segid = MRIgetVoxVal(aseg,c,r,s,0);
	if(segid == 540){
	  aqfound=1;
	  break;
	}
      }
      if(aqfound){
	supslice.push_back(s-1);
	aqfound=1;
	printf("row %d %d\n",r,s-1);
	break;
      }
    } // slice
    if(!aqfound) supslice.push_back(0);
  } // row
  int naq=0, nbs=0;
  double dasbaqsum = 0, dasbbssum = 0;
  for(int c=0; c < aseg->width; c++){
    for(int r=0; r < aseg->height; r++){
      for(int s=0; s < aseg->depth;  s++){
	int segid = MRIgetVoxVal(aseg,c,r,s,0);
	if(segid == 540){
	  naq++;
	  dasbaqsum += MRIgetVoxVal(dasb,c,r,s,0);
	  if(AqInfSlice > s) AqInfSlice = s;
	  if(AqSupSlice < s) AqSupSlice = s;
	  if(AqPostRow > r)  AqPostRow = r;
	  if(AqAntRow  < r)  AqAntRow = r;
	}
	if(segid == 16){
	  nbs++;
	  dasbbssum += MRIgetVoxVal(dasb,c,r,s,0);
	  if(BsLeftCol  > c) BsLeftCol=c;
	  if(BsRightCol < c) BsRightCol=c;
	}
      }
    }
  }
  // give it another index
  AqInfSlice--;
  AqPostRow--;

  double dasbaqmn = dasbaqsum/naq;
  double dasbbsmn = dasbbssum/naq;
  printf("nbs = %d dasbbs = %g  naq = %d dasbaq = %g\n",nbs,dasbbsmn,naq,dasbaqmn);
  printf("Bounding Box: c = (%d,%d), r=(%d,%d), s=(%d,%d) DRAPVox=%d MRAPVox=%d\n",
	 BsLeftCol,BsRightCol,AqPostRow,AqAntRow,AqInfSlice,AqSupSlice,DRAPvox,MRAPvox);

  // DR ==============================================================
  int drseedc=0, drseedr=0, drseeds=0;
  double drmax = 0;
  for(int r=AqAntRow-DRAPvox; r <= AqAntRow; r++){
    for(int s=AqInfSlice; s <= supslice[r]; s++){
      for(int c=BsLeftCol; c <= BsRightCol; c++){
	int segid = MRIgetVoxVal(aseg,c,r,s,0);
	double val = MRIgetVoxVal(dasb,c,r,s,0);
	//printf("%d %d %d   %d  %g\n",c,r,s,segid,val);
	if(segid != 16) continue;
	if(val > drmax){
	  drmax = val;
	  drseedc = c;
	  drseedr = r;
	  drseeds = s;
	}
      }
    }
  }
  printf("DR: max=%g, seed = (%d,%d,%d)\n",drmax,drseedc,drseedr,drseeds);
  MRIsetVoxVal(aseg,drseedc,drseedr,drseeds,0,118);
  int drvolvoxtarg = round(drtargvol/voxvol);
  printf("voxvol = %g, drvolvoxtarg = %d\n",voxvol,drvolvoxtarg);
  std::vector<std::vector<int>> drcrs;
  std::vector<int>crs = {drseedc,drseedr,drseeds};
  drcrs.push_back(crs);
  while(1){
    double dasbmax = -1;
    int cnbr=0, rnbr=0, snbr=0;
    for(int n=0; n < drcrs.size(); n++){
      for(int cd = -1; cd < 2; cd++){
	for(int rd = -1; rd < 2; rd++){
	  for(int sd = -1; sd < 2; sd++){
	    if(abs(cd)+abs(rd)+abs(sd) != 1) continue; // face nbrs only
	    int c = drcrs[n][0]+cd;
	    int r = drcrs[n][1]+rd;
	    int s = drcrs[n][2]+sd;
	    int segid = MRIgetVoxVal(aseg,c,r,s,0);
	    if(segid != 16) continue;
	    if(s > supslice[r]) continue;
	    if(s < AqInfSlice ) continue;
	    if(r < AqAntRow-DRAPvox || r > AqAntRow) continue;
	    double val = MRIgetVoxVal(dasb,c,r,s,0);
	    if(val > dasbmax){
	      dasbmax = val;
	      cnbr = c;
	      rnbr = r;
	      snbr = s;
	    }
	  }//sd
	}//rd
      }//cd
    } // loop over voxels in current DR
    if(dasbmax < 0){
      printf("ERROR: could not fill DR up to %d (%d)\n",drvolvoxtarg,(int)drcrs.size());
      break;
    }
    crs.clear();
    crs = {cnbr,rnbr,snbr};
    drcrs.push_back(crs);
    MRIsetVoxVal(aseg,cnbr,rnbr,snbr,0,118);
    if(drcrs.size() == drvolvoxtarg) break;
  }
  printf("#DR# Found %d voxels\n",(int)drcrs.size());

  // MR ==============================================================
  int mrseedc=0, mrseedr=0, mrseeds=0;
  double mrmax = 0;
  for(int r=AqAntRow-MRAPvox; r <= AqAntRow; r++){
    for(int s=0; s <= AqInfSlice; s++){
      for(int c=BsLeftCol; c <= BsRightCol; c++){
	int segid = MRIgetVoxVal(aseg,c,r,s,0);
	double val = MRIgetVoxVal(dasb,c,r,s,0);
	//printf("%d %d %d   %d  %g\n",c,r,s,segid,val);
	if(segid != 16) continue;
	if(val > mrmax){
	  mrmax = val;
	  mrseedc = c;
	  mrseedr = r;
	  mrseeds = s;
	}
      }
    }
  }
  printf("MR: max=%g, seed = (%d,%d,%d)\n",mrmax,mrseedc,mrseedr,mrseeds);
  MRIsetVoxVal(aseg,mrseedc,mrseedr,mrseeds,0,119);
  int mrvolvoxtarg = round(mrtargvol/voxvol);
  printf("voxvol = %g, mrvolvoxtarg = %d\n",voxvol,mrvolvoxtarg);
  std::vector<std::vector<int>> mrcrs;
  crs.clear();
  crs = {mrseedc,mrseedr,mrseeds};
  mrcrs.push_back(crs);
  while(1){
    double dasbmax = -1;
    int cnbr=0, rnbr=0, snbr=0;
    for(int n=0; n < mrcrs.size(); n++){
      for(int cd = -1; cd < 2; cd++){
	for(int rd = -1; rd < 2; rd++){
	  for(int sd = -1; sd < 2; sd++){
	    if(abs(cd)+abs(rd)+abs(sd) != 1) continue; // face nbrs only
	    int c = mrcrs[n][0]+cd;
	    int r = mrcrs[n][1]+rd;
	    int s = mrcrs[n][2]+sd;
	    int segid = MRIgetVoxVal(aseg,c,r,s,0);
	    if(segid != 16) continue;
	    if(s > supslice[r]) continue;
	    if(r < AqAntRow-MRAPvox || r > AqAntRow) continue;
	    double val = MRIgetVoxVal(dasb,c,r,s,0);
	    if(val > dasbmax){
	      dasbmax = val;
	      cnbr = c;
	      rnbr = r;
	      snbr = s;
	    }
	  }//sd
	}//rd
      }//cd
    } // loop over voxels in current MR
    if(dasbmax < 0){
      printf("ERROR: could not fill MR up to %d (%d)\n",mrvolvoxtarg,(int)mrcrs.size());
      break;
    }
    crs.clear();
    crs = {cnbr,rnbr,snbr};
    mrcrs.push_back(crs);
    MRIsetVoxVal(aseg,cnbr,rnbr,snbr,0,119);
    if(mrcrs.size() == mrvolvoxtarg) break;
  }
  printf("#MR# Found %d voxels\n",(int)mrcrs.size());

  if(strcmp(ostr,"RAS")){
    printf("Reorientating back to %s\n",ostr);
    MRI *tmp = MRIresample(aseg, aseg0, SAMPLE_NEAREST);
    MRIfree(&aseg);
    aseg = tmp;
  }

  // Update the ctab
  if(!aseg->ct) aseg->ct = CTABreadDefault();
  CTE *cte = aseg->ct->entries[118];
  if(cte == NULL) {
    cte = (CTE *)malloc(sizeof(CTE));
    aseg->ct->entries[118] = cte;
  }
  strcpy(cte->name,"Dorsal-Raphe");
  cte->ri=255; cte->gi=20; cte->bi=147; cte->ai=0;
  cte->rf=cte->ri/255.0; cte->gf=cte->gi/255.0; cte->bf=cte->gi/255.0;

  cte = aseg->ct->entries[119];
  if(cte == NULL) {
    cte = (CTE *)malloc(sizeof(CTE));
    aseg->ct->entries[118] = cte;
  }
  strcpy(cte->name,"Median-Raphe");
  cte->ri=205; cte->gi=179; cte->bi=139; cte->ai=0;
  cte->rf=cte->ri/255.0; cte->gf=cte->gi/255.0; cte->bf=cte->gi/255.0;

  cte = aseg->ct->entries[540];
  if(cte == NULL) {
    cte = (CTE *)malloc(sizeof(CTE));
    aseg->ct->entries[540] = cte;
  }
  strcpy(cte->name,"Cerebral-Aqueduct");
  cte->ri=170; cte->gi=85; cte->bi=255; cte->ai=0;
  cte->rf=cte->ri/255.0; cte->gf=cte->gi/255.0; cte->bf=cte->gi/255.0;

  if(dr0 || mr0) RapheDice(aseg);

  return(aseg);
}
int RapheSeg::RapheDice(MRI *aseg)
{
  // Dice ===============================================
  char ostr0[5];
  ostr0[4] = '\0';
  MRIdircosToOrientationString(aseg,ostr0);
  printf("Dice: Aseg Orientation %s\n",ostr0);
  if(dr0 || mr0){
    char ostr[5];
    ostr[4] = '\0';
    if(dr0){
      MRIdircosToOrientationString(dr0,ostr);
      printf("DR0 Orientation %s\n",ostr);
      if(strcmp(ostr,ostr0)){
	printf("Reorientating DR0 to RAS\n");
	MRI *rastemplate = MRIcopy(dr0,NULL);
	MRIorientationStringToDircos(rastemplate, ostr0);
	MRI *tmp = MRIresample(dr0, rastemplate, SAMPLE_NEAREST);
	MRIfree(&dr0);
	dr0 = tmp;
	MRIfree(&rastemplate);
      } 
    }
    if(mr0){
      MRIdircosToOrientationString(mr0,ostr);
      printf("MR0 Orientation %s\n",ostr);
      if(strcmp(ostr,ostr0)){
	printf("Reorientating MR0 to RAS\n");
	MRI *rastemplate = MRIcopy(mr0,NULL);
	MRIorientationStringToDircos(rastemplate, ostr0);
	MRI *tmp = MRIresample(mr0, rastemplate, SAMPLE_NEAREST);
	MRIfree(&mr0);
	mr0 = tmp;
	MRIfree(&rastemplate);
      } 
    }
    int ndr0=0, nmr0=0, ndr=0, nmr=0, ndrboth=0, nmrboth=0;
    for(int c=0; c < aseg->width; c++){
      for(int r=0; r < aseg->height; r++){
	for(int s=0; s < aseg->depth;  s++){
	  int segid = MRIgetVoxVal(aseg,c,r,s,0);
	  if(segid == 118) ndr++;
	  if(dr0){
	    double val = MRIgetVoxVal(dr0,c,r,s,0);
	    if(val > 0.01){
	      ndr0++;
	      if(segid == 118) ndrboth++;
	    }
	  }
	  if(segid == 119) nmr++;
	  if(mr0){
	    double val = MRIgetVoxVal(mr0,c,r,s,0);
	    if(val > 0.01){
	      nmr0++;
	      if(segid == 119) nmrboth++;
	    }
	  }
	}
      }
    }
    if(dr0){
      double dice = (double)ndrboth/((ndr0+ndr)/2.0);
      printf("DRdice: %3d %3d %3d %6.4lf\n",ndr0,ndr,ndrboth,dice);
    }
    if(mr0){
      double dice = (double)nmrboth/((nmr0+nmr)/2);
      printf("MRdice: %3d %3d %3d %6.4lf\n",nmr0,nmr,nmrboth,dice);
    }
  }
  return(0);
}
class MRISflat2mri {
public:
  MRIS *surf=NULL;
  MRIS *surfreg=NULL;
  //LABEL *label=NULL;
  double lhfsasph1xyz[3] = {11.40, 15.42, -98.14}; //vno 136722 
  double lhfsasph2xyz[3] = {-4.52, 30.08, -95.26}; //vno 52903
  double lhfsasph3xyz[3] = {2.92, 36.70, -92.98}; //vno 136674
  int fsavno1=-1, fsavno2=-1;
  int nvox[2] = {-1,-1};
  double voxsize[2] = {0,0};
  double dthresh = 2;
  MRI *flat2mri(MRI *ov);
};

MRI *MRISflat2mri::flat2mri(MRI *ov)
{
  // Get a list of the xy coords of the vertices in the patch
  std::vector<double> vx,vy;
  std::vector<int> vnolist;
  for(int vno=0; vno < surf->nvertices; vno++){
    VERTEX *v = &(surf->vertices[vno]);    
    if(v->ripflag) continue;
    vx.push_back(v->x);
    vy.push_back(v->y);
    vnolist.push_back(vno);
    // get xmin/max ymin/max too
  }
  printf("Found %d vertices in patch\n",(int)vx.size());

  // Three points to define the coordinate system
  float dmin;
  int vno1 = MRISfindClosestVertex(surfreg, lhfsasph1xyz[0], lhfsasph1xyz[1], lhfsasph1xyz[2], &dmin, CURRENT_VERTICES);
  int vno2 = MRISfindClosestVertex(surfreg, lhfsasph2xyz[0], lhfsasph2xyz[1], lhfsasph2xyz[2], &dmin, CURRENT_VERTICES);
  int vno3 = MRISfindClosestVertex(surfreg, lhfsasph3xyz[0], lhfsasph3xyz[1], lhfsasph3xyz[2], &dmin, CURRENT_VERTICES);

  MRI *mri = MRIallocSequence(nvox[0],nvox[1],1,ov->type,ov->nframes);
  MRIcopyHeader(ov,mri);
  MRIcopyPulseParameters(ov,mri);
  if(ov->ct) mri->ct = CTABdeepCopy(ov->ct);  
  mri->valid = 1;
  mri->xsize = voxsize[0];
  mri->ysize = voxsize[1];
  mri->zsize = 1; // flat in this dim

  // vno1 and vno2 are two vertices that dictate the direction of the columns/width
  VERTEX *v1 = &(surf->vertices[vno1]);
  VERTEX *v2 = &(surf->vertices[vno2]);
  VERTEX *v3 = &(surf->vertices[vno3]);
  double d12 = sqrt( pow(v1->x-v2->x,2.0) + pow(v1->y-v2->y,2.0) + pow(v1->z-v2->z,2.0));

  // The center is defined midway between these two vertices
  mri->c_r  = (v1->x+v2->x)/2;
  mri->c_a  = (v1->y+v2->y)/2;
  mri->c_s  = 0;

  // Set the geometry to be axial. When loading the slice into FV, it
  // will appear that you are looking down on an axial slice.

  // Set slice axis to point in z (axial slice)
  mri->z_r = 0;
  mri->z_a = 0;
  mri->z_s = 1;

  // Set col axis to point along the vector between v1 and v2
  mri->x_r = +(v2->x - v1->x)/d12;
  mri->x_a = +(v2->y - v1->y)/d12;
  mri->x_s = 0;

  // Set row axis to be orthog to col axis
  mri->y_r = -mri->x_a;
  mri->y_a = mri->x_r;
  mri->y_s = 0;

  double d32r = v3->x - v2->x;
  double d32a = v3->y - v2->y;
  double tmp = mri->y_r * d32r + mri->y_a * d32a;
  if(0 && tmp < 0){
    mri->y_r = -mri->y_a;
    mri->y_a = -mri->y_r;
    printf("Reversing y DC\n");
  }

  MATRIX *vox2ras = mri->get_Vox2RAS(0);

  printf("v1 %d %g %g %g\n",vno1,v1->x,v1->y,v1->z);
  printf("v2 %d %g %g %g\n",vno2,v2->x,v2->y,v2->z);
  printf("x_r %g x_a %g\n",mri->x_r,mri->x_a);
  printf("center %g %g %g\n",mri->c_r,mri->c_a,mri->c_s);
  printf("d12 = %g\n",d12);
  printf("vox2ras = [\n");
  MatrixPrint(stdout,vox2ras);
  printf("];\n");

  //MHT *Hash = MHTcreateVertexTable_Resolution(surf, CURRENT_VERTICES, 16);

  // Go through each point in the mri and find closest vertex in patch
  // Probably a more efficient way to do this with hash, but the data
  // are pretty small. Could probably improve the interp
  MATRIX *vox = MatrixAlloc(4,1,MATRIX_REAL);
  vox->rptr[4][1] = 1;
  MATRIX *ras=NULL;
  int nout = 0;
  for(int c=0; c < mri->width; c++){
    for(int r=0; r < mri->height; r++){
      vox->rptr[1][1] = c;
      vox->rptr[2][1] = r;
      ras = MatrixMultiplyD(vox2ras,vox,ras);
      double dmin = 10e10;
      int nmin = -1;
      for(int n=0; n < vx.size(); n++){
	double dx = vx[n] - ras->rptr[1][1];
	double dy = vy[n] - ras->rptr[2][1];
	double d2 = sqrt(dx*dx + dy*dy);
	if(dmin > d2){
	  dmin = d2;
	  nmin = n;
	}
      }
      if(dmin > 2) nout++;
      if((c==0 && r==0) || (c==mri->width-1 && r==0) || (c==0 && r==mri->height-1) || (c==mri->width-1 && r==mri->height-1)) printf("c=%d r=%d vno=%d\n",c,r,vnolist[nmin]);
      for(int f=0; f < mri->nframes; f++){
	int vno = vnolist[nmin];
	double val = MRIgetVoxVal(ov,vno,0,0,f);
	if(dmin > dthresh) val = 0;
	MRIsetVoxVal(mri,c,r,0,f,val);
      }
    }
  }
  // Could reset the geom so that all subjects have 
  // the same.
  printf("nout %d\n",nout);

  MatrixFree(&vox2ras);
  MatrixFree(&vox);
  MatrixFree(&ras);
  return(mri);
}

class MRISsphere2mri {
public:
  MRIS *sphere=NULL; // just sphere, eg lh.sphere
  MRIS *surfreg=NULL; // sphere.reg of this subject
  LTA *lta = NULL;
  //LABEL *label=NULL;
  double lhfsasph1xyz[3] = {11.40, 15.42, -98.14}; //vno 136722 
  double lhfsasph2xyz[3] = {-4.52, 30.08, -95.26}; //vno 52903
  double lhfsasph3xyz[3] = {2.92, 36.70, -92.98}; //vno 136674
  double rhfsasph1xyz[3] = {-32.1000, 23.9500, -91.6300}; //90206
  double rhfsasph2xyz[3] = {-10.3700, 40.0200, -91.0500}; //89549
  double rhfsasph3xyz[3] = {-17.4600, 42.9400, -88.6100}; //13969
  //90206 89549 13969
  int fsavno1=-1, fsavno2=-1;
  int nvox[2] = {-1,-1};
  double P0[2] = {-180,-90}, P1[2] = {180,90};
  double dthresh = 2;
  MRI *sphere2mri(MRI *ov, int direction=1);
  void xyz2rpt(double x, double y, double z, double *r, double *phi, double *theta){
    *r = sqrt(x*x+y*y+z*z);
    *theta = 90 - acos(z/(*r))*180/M_PI; // -90 to +90
    *phi = atan2(y,x)*180/M_PI;     // -180 to +180
  }
  void rpt2xyz(double r, double phi, double theta, double *x, double *y, double *z){
    double th = (90-theta)*M_PI/180;
    *x = r*sin(th)*cos(phi*M_PI/180);
    *y = r*sin(th)*sin(phi*M_PI/180);
    *z = r*cos(th);
  }
  void rpt2xyz(double r, double phi, double theta, float *x, float *y, float *z){
    double xx,yy,zz;
    this->rpt2xyz(r,phi,theta,&xx,&yy,&zz);
    *x = xx; *y = yy; *z = zz;
  }
};
MRI *MRISsphere2mri::sphere2mri(MRI *ov, int direction)
{
  // direction=1=sphere2mri
  // direction=2=mri2sphere
  if(direction != 1 && direction != 2){
    printf("ERROR: diretion=%d must be 1 or 2\n",direction);
    return(NULL);
  }

  // Three points to define the common coordinate system. Map them to the current surface via sphere.reg=surfreg
  float dmin;
  int vno1, vno2, vno3;
  if(surfreg->hemisphere == 0){
    vno1 = MRISfindClosestVertex(surfreg, lhfsasph1xyz[0], lhfsasph1xyz[1], lhfsasph1xyz[2], &dmin, CURRENT_VERTICES);
    vno2 = MRISfindClosestVertex(surfreg, lhfsasph2xyz[0], lhfsasph2xyz[1], lhfsasph2xyz[2], &dmin, CURRENT_VERTICES);
    vno3 = MRISfindClosestVertex(surfreg, lhfsasph3xyz[0], lhfsasph3xyz[1], lhfsasph3xyz[2], &dmin, CURRENT_VERTICES);
  }
  else{
    vno1 = MRISfindClosestVertex(surfreg, rhfsasph1xyz[0], rhfsasph1xyz[1], rhfsasph1xyz[2], &dmin, CURRENT_VERTICES);
    vno2 = MRISfindClosestVertex(surfreg, rhfsasph2xyz[0], rhfsasph2xyz[1], rhfsasph2xyz[2], &dmin, CURRENT_VERTICES);
    vno3 = MRISfindClosestVertex(surfreg, rhfsasph3xyz[0], rhfsasph3xyz[1], rhfsasph3xyz[2], &dmin, CURRENT_VERTICES);
  }

  // From here on out we use the unregistered ?h.sphere (less distortion)

  // vno1 and vno2 are two vertices that dictate the direction of the columns/width
  VERTEX *v1 = &(sphere->vertices[vno1]);
  VERTEX *v2 = &(sphere->vertices[vno2]);
  VERTEX *v3 = &(sphere->vertices[vno3]);
  //double d12 = sqrt( pow(v1->x-v2->x,2.0) + pow(v1->y-v2->y,2.0) + pow(v1->z-v2->z,2.0));
  //double radius = sqrt( pow(v1->x,2.0) + pow(v1->y,2.0) + pow(v1->z,2.0));
  double radius, phi1, theta1, phi2, theta2, phi3, theta3;
  this->xyz2rpt(v1->x,v1->y,v1->z,&radius,&phi1,&theta1);
  this->xyz2rpt(v2->x,v2->y,v2->z,&radius,&phi2,&theta2);
  this->xyz2rpt(v3->x,v3->y,v3->z,&radius,&phi3,&theta3);
  printf("Prior to any +=360 phi %6.4lf %6.4lf %6.4lf\n",phi1,phi2,phi3);
  // This can happen on rh
  if(phi1 < -90) phi1 += 360;
  if(phi2 < -90) phi2 += 360;
  if(phi3 < -90) phi3 += 360;
  double cphi = (phi1+phi2)/2;
  double ctheta = (theta1+theta2)/2;
  double calpha = atan2((theta2-theta1),(phi2-phi1))*180/M_PI;
  if(surfreg->hemisphere != 0) calpha = 180-calpha;
  double cx,cy,cz;
  this->rpt2xyz(radius,cphi,ctheta,&cx,&cy,&cz);
  int cvno = MRISfindClosestVertex(sphere, cx, cy, cz, &dmin, CURRENT_VERTICES);
  VERTEX *cv = &(sphere->vertices[cvno]);
  printf("vno1=%d (%g,%g,%g) vno2=%d (%g,%g,%g) vno3=%d (%g,%g,%g) \n",vno1,
	 v1->x,v1->y,v1->z, vno2,v2->x,v2->y,v2->z, vno3,v3->x,v3->y,v3->z);
  printf("phi = (%g,%g,%g) theta=(%g,%g,%g)\n",phi1,phi2,phi3,theta1,theta2,theta3);
  printf("center phi=%g theta=%g alpha=%g %4.1f %4.1f %4.1f  cvno=%d  %4.1f %4.1f %4.1f  d=%g\n",
	 cphi,ctheta,calpha,cx,cy,cz,cvno,cv->x,cv->y,cv->z,dmin);
  printf("radius %g\n",radius);
  printf("v1 %d %g %g %g\n",vno1,v1->x,v1->y,v1->z);
  printf("v2 %d %g %g %g\n",vno2,v2->x,v2->y,v2->z);

  double par[12];
  for(int k=0; k<12; k++) par[k]=0;
  par[3] = calpha; // rotation about x
  par[4] = ctheta; // rotation about y
  par[5] = cphi; // rotation about z (applied first)
  par[6] = 1; par[7] = 1; par[8] = 1;
  MATRIX *T = TranformAffineParams2Matrix(par, NULL);

  //VERTEX *v = &(sphere->vertices[cvno]);
  //printf("before %g %g %g\n",v->x,v->y,v->z);
  MRISmatrixMultiply(sphere, T);
  //MRISwrite(sphere,"lh.sphere.T");
  //printf("after  %g %g %g\n",v->x,v->y,v->z);
  MatrixFree(&T);
  MHT *Hash = MHTcreateVertexTable_Resolution(sphere, CURRENT_VERTICES, 16);

  MRI *ovsurf=NULL, *ovslice=NULL;
  if(direction == 1){
    ovsurf = ov;
    ovslice = MRIallocSequence(nvox[0],nvox[1],1,ovsurf->type,ovsurf->nframes);
    MRIcopyPulseParameters(ovsurf,ovslice);
    if(ovsurf->ct) ovslice->ct = CTABdeepCopy(ovsurf->ct);  
    ovslice->valid = 1;
    
    // The coordinates of the MRI structure will be in phi/theta degrees
    ovslice->xsize = fabs((P1[0]-P0[0])/nvox[0]);
    ovslice->ysize = fabs((P1[1]-P0[1])/nvox[1]);
    ovslice->zsize = 1; // flat in this dim

    // The center
    ovslice->c_r  = 0.5;
    ovslice->c_a  = (P1[0]+P0[0])/2;
    ovslice->c_s  = (P1[1]+P0[1])/2;
    // Set col axis to A->P 
    ovslice->x_r = 0;
    ovslice->x_a = 1;
    ovslice->x_s = 0;
    // Set row axis to be orthog to col axis
    ovslice->y_r = 0;
    ovslice->y_a = 0;
    ovslice->y_s = 1;
    // Set slice axis to point in x
    ovslice->z_r = 1;
    ovslice->z_a = 0;
    ovslice->z_s = 0;
  } else {
    ovslice = ov;
    ovsurf = MRIallocSequence(sphere->nvertices,1,1,ovslice->type,ovslice->nframes);
    MRIcopyPulseParameters(ovslice,ovsurf);
    if(ovslice->ct) ovsurf->ct = CTABdeepCopy(ovslice->ct);  
  }
  MATRIX *vox2ras = ovslice->get_Vox2RAS(0);
  printf("x_r %g x_a %g\n",ovslice->x_r,ovslice->x_a);
  printf("center %g %g %g\n",ovslice->c_r,ovslice->c_a,ovslice->c_s);
  //printf("vox2ras = [\n");
  //MatrixPrint(stdout,vox2ras);
  //printf("];\n");

  // Go through each point in the mri and find closest vertex in patch
  // Probably a more efficient way to do this with hash, but the data
  // are pretty small. Could probably improve the interp
  MATRIX *vox = MatrixAlloc(4,1,MATRIX_REAL);
  vox->rptr[4][1] = 1;
  MATRIX *phitheta=NULL;
  for(int c=0; c < ovslice->width; c++){
    for(int r=0; r < ovslice->height; r++){
      if(surfreg->hemisphere == 0) vox->rptr[1][1] = c; //lh 
      if(surfreg->hemisphere == 1) vox->rptr[1][1] = ovslice->width - c - 1; //rh
      vox->rptr[2][1] = r;
      vox->rptr[3][1] = 0;
      phitheta = MatrixMultiplyD(vox2ras,vox,phitheta);
      double p = phitheta->rptr[2][1];
      double t = phitheta->rptr[3][1];
      double x,y,z;
      this->rpt2xyz(radius,p,t,&x,&y,&z);
      int vno = Hash->findClosestVertexNoXYZ(x,y,z, &dmin);
      for(int f=0; f < ovslice->nframes; f++){
	double val;
	if(direction == 1){
	  val = MRIgetVoxVal(ovsurf,vno,0,0,f);
	  MRIsetVoxVal(ovslice,c,r,0,f,val);
	} else {
	  val = MRIgetVoxVal(ovslice,c,r,0,f);
	  // If the slice is 0 here and if the surf already has a
	  // non-zero value here, then don't overwrite it. This is to
	  // try to prevent loss of non-zero voxels at the edge. These
	  // may be lost in the forward (sphere2mri) loop.
	  if(val == 0 && MRIgetVoxVal(ovsurf,vno,0,0,f)!=0) continue;
	  MRIsetVoxVal(ovsurf,vno,0,0,f,val);
	}
      }
    }
  }

  if(direction == 2){ // should only be done with seg
    // Fill holes: Bin. Invert. Cluster. Loop through voxels in clusters 2-N
    MATRIX *ras2vox = MatrixInverse(vox2ras,NULL);
    SURFCLUSTERSUM *scs;
    int NClusters=0;
    MRIScopyMRI(sphere, ovsurf, 0, "val");
    // Set marked by binarizing map as it is
    for(int vno=0; vno < sphere->nvertices; vno++) {
      VERTEX *v = &(sphere->vertices[vno]);
      if(v->val > .001) v->marked = 1;
      else              v->marked = 0;
    }
    // Close using marked by dilating by 2 and then eroding by 1. The
    // goal here is only to make sure that any voids are actual holes
    // and not inlets.
    MRISdilateMarked(sphere,2);
    MRISerodeMarked(sphere,1);
    // Reset val to be the inverse of the closed marked
    for(int vno=0; vno < sphere->nvertices; vno++) {
      VERTEX *v = &(sphere->vertices[vno]);
      v->val = !v->marked;
    }
    // Cluster to find any holes
    scs = sclustMapSurfClusters(sphere,.0001,-1,0,0,&NClusters,NULL,NULL);
    printf("nc = %d\n",NClusters);
    for(int vno=0; vno < sphere->nvertices; vno++){
      VERTEX *v = &(sphere->vertices[vno]);
      int cno = v->undefval;
      // If it is not marked and not in a hole then skip vertex
      if(! v->marked && cno < 2) continue;
      double phi, theta;
      phitheta->rptr[4][1] = 1; // just to make sure
      this->xyz2rpt(v->x,v->y,v->z,&radius,&phi,&theta);
      phitheta->rptr[1][1] = 0;
      phitheta->rptr[2][1] = phi; 
      phitheta->rptr[3][1] = theta;
      vox = MatrixMultiplyD(ras2vox,phitheta,vox);
      int c = round(vox->rptr[1][1]);
      int r = round(vox->rptr[2][1]);
      if(c < 0 || c >= ovslice->width) continue;
      if(r < 0 || r >= ovslice->height) continue;
      if(surfreg->hemisphere == 1) c = ovslice->width - c  - 1; //rh
      double val = MRIgetVoxVal(ovslice,c,r,0,0);
      //printf("cno=%d vno = %d  pt=(%g,%g) cr = (%d %d)  %g\n",cno,vno,phi,theta,c,r,val);
      MRIsetVoxVal(ovsurf,vno,0,0,0,val);
    }
    MatrixFree(&ras2vox);
  }

  MatrixFree(&vox2ras);
  MatrixFree(&vox);
  MatrixFree(&phitheta);
  if(direction == 1) return(ovslice);
  else return(ovsurf);
}

#if 0
class SpatTempCluster {
  // test masking, edge and corner
  // clusters - centroid, pointsets, size
public:
  int topo=1; //1 = volume, 2 = surface
  int nbrtype = 1; //1=face, 2=+edge, 3=+corner
  MRI *binmask=NULL;
  MRI *cnomap =NULL;
  MRIS *surf=NULL;
  std::vector<std::vector<int>> voxlist; //c,r,s,f
  std::vector<std::vector<int>> GetNearestNeighbors(std::vector<int> vox){
    std::vector<std::vector<int>> nbrs;
    if(topo==1){ // volume
      for(int dc=-1; dc <= +1; dc++){
	for(int dr=-1; dr <= +1; dr++){
	  for(int ds=-1; ds <= +1; ds++){
	    if(abs(dc)+abs(dr)+abs(ds) > nbrtype) continue;
	    if(vox[0]+dc < 0 || vox[0]+dc >= binmask->width) continue;
	    if(vox[1]+dr < 0 || vox[1]+dr >= binmask->height) continue;
	    if(vox[2]+ds < 0 || vox[2]+ds >= binmask->depth) continue;
	    for(int df=-1; df <= +1; df++){
	      if(abs(dc)+abs(dr)+abs(ds)+abs(df)==0) continue; // not self
	      // this line forces a face topology across frames
	      if(abs(df) > 0 && (abs(dc)+abs(dr)+abs(ds) > 0) ) continue;
	      if(vox[3]+df < 0 || vox[3]+df >= binmask->nframes) continue;
	      std::vector<int> dcrsf = {vox[0]+dc,vox[1]+dr,vox[2]+ds,vox[3]+df};
	      nbrs.push_back(dcrsf);
	    }
	  }
	}
      }
    } else { // surface 
      int vno = vox[0];
      // Doing frames separately forces face topology across frame
      for(int df=-1; df <= +1; df++){
	if(vox[3]+df < 0 || vox[3]+df >= binmask->nframes) continue;
	std::vector<int> dcrsf = {vox[0],0,0,vox[3]+df};
	nbrs.push_back(dcrsf);
      }
      VERTEX_TOPOLOGY *vtop = &surf->vertices_topology[vno];
      int vnum = vtop->vtotal;
      for(int i = 0; i < vnum; i++) {
	int vnbno = vtop->v[i];
	std::vector<int> dcrsf = {vnbno,0,0,vox[3]};
	nbrs.push_back(dcrsf);
      }
    }
    return(nbrs);
  }
  class Cluster {
  public: 
    int cno = 0;
    int nmembers = 0;
    std::vector<std::vector<int>> crst;
  };
  std::vector<Cluster> ClusterList;
  int GrowOne(std::vector<int> vox, int cno)
  {
    int m;
    m = MRIgetVoxVal(this->binmask,vox[0],vox[1],vox[2],vox[3]);
    if(!m) return(0); // not a set voxel
    m = MRIgetVoxVal(this->cnomap,vox[0],vox[1],vox[2],vox[3]);
    if(m) return(0); // already in a cluster
    // To get here, the vox must be active in the binmask and not in a cluster already
    MRIsetVoxVal(this->cnomap,vox[0],vox[1],vox[2],vox[3], cno+1);
    this->ClusterList[cno].crst.push_back(vox);
    int nhits = 1;
    std::vector<std::vector<int>> nbrs = this->GetNearestNeighbors(vox);
    for(int n=0; n < nbrs.size(); n++){
      nhits += this->GrowOne(nbrs[n],cno);
    }
    this->ClusterList[cno].nmembers += nhits;
    return(nhits);
  }
  int Clusterize(void){
    if(this->cnomap) MRIfree(&this->cnomap);
    this->cnomap = MRIclone(this->binmask,NULL);
    this->ClusterList.clear();
    int pointno = -1;
    int cno = -1;
    while(1){
      pointno ++;
      if(pointno >= this->voxlist.size()) break;
      int done = 0;
      std::vector<int> vox;
      while(1){ // find the next point that has not been marked
	vox = voxlist[pointno];
	int val = MRIgetVoxVal(this->cnomap,vox[0],vox[1],vox[2],vox[3]);
	if(val) { // point already found
	  pointno++;
	  if(pointno >= this->voxlist.size()) {
	    done = 1;
	    break;
	  }
	}
	else break; // found an unmarked point
      }
      if(done) break;
      cno++;
      Cluster cl;
      cl.cno = cno;
      printf("cno=%d pointno=%d  %d %d %d %d =====\n",cno,pointno,vox[0],vox[1],vox[2],vox[3]);
      this->ClusterList.push_back(cl);
      this->GrowOne(vox,cno);
    }
    printf("Found %d clusters\n",(int)this->ClusterList.size());
    this->cnomap->ct = CTABalloc(this->ClusterList.size()+1);
    CTABunique(this->cnomap->ct, 100);
    return(this->ClusterList.size());
  }
  int PrintClusterSum(FILE *fp){
    for(int n=0; n < this->ClusterList.size(); n++){
      Cluster cl = ClusterList[n];
      fprintf(fp,"%2d %4d   %3d %3d %3d  %3d\n",n+1,(int)cl.crst.size(),cl.crst[0][0],cl.crst[0][1],cl.crst[0][2],cl.crst[0][3]);
    }
    fflush(fp);
    return(0);
  }
  int GetBinMask(MRI *ov, double thmin, double thmax, int sign, MRI *mask)
  {
    if(this->binmask) MRIfree(&this->binmask);
    this->binmask = MRIallocSequence(ov->width,ov->height,ov->depth,MRI_INT,ov->nframes);
    MRIcopyHeader(ov,this->binmask);
    this->voxlist.clear();
    //not thread safe
    for(int c=0; c < ov->width; c++){
      for(int r=0; r < ov->height; r++){
	for(int s=0; s < ov->depth; s++){
	  if(mask && MRIgetVoxVal(ov,c,r,s,0)<0.5) continue;
	  for(int f=0; f < ov->nframes; f++){
	    double val = MRIgetVoxVal(ov,c,r,s,f);
	    if(sign ==  0) val = fabs(val);
	    if(sign == -1) val = -val;
	    if(val < thmin || val > thmax) continue;
	    std::vector<int> vox = {c,r,s,f};
	    voxlist.push_back(vox);
	    MRIsetVoxVal(binmask,c,r,s,f,1);
	  }
	}
      }
    }
    printf("thmin=%g thmax=%g sign=%d nhits=%d\n",thmin,thmax,sign,(int)voxlist.size());
    return(voxlist.size());
  }

};
#endif

#if 0
MRI *MRIsegBorder(MRI *seg, int *maxlabel, int topo)
{
  // This is supposed to replicate the border creation as done in
  // MRIvoxelsInLabelWithPartialVolumeEffects(). In that function,
  // the border is only created for a single label, where as here
  // it is created for all labels. That func uses topo=1
  int nnbrs=0;
  if(topo == 1) nnbrs = 6;  // face neighbors
  if(topo == 2) nnbrs = 18; // +edge neighbors
  if(topo == 3) nnbrs = 26; // +corner neighbors
  MRI *border = MRIallocSequence(seg->width,seg->height,seg->depth,MRI_INT,1+nnbrs);
  MRIcopyHeader(seg,border);
  MRIcopyPulseParameters(seg,border);

  *maxlabel = 0;
  int nhits = 0;
  for(int c=0; c < seg->width; c++){
    for(int r=0; r < seg->height; r++){
      for(int s=0; s < seg->depth; s++){
	int voxlabel = MRIgetVoxVal(seg,c,r,s,0);
	if(*maxlabel < voxlabel) *maxlabel = voxlabel;
	int knbr = 0;
	for(int dc = -1; dc < 2; dc++){
	  int cc = c+dc;
	  if(cc < 0 || cc >= seg->width) continue;
	  for(int dr = -1; dr < 2; dr++){
	    int rr = r+dr;
	    if(rr < 0 || rr >= seg->height) continue;
	    for(int ds = -1; ds < 2; ds++){
	      int ss = s+ds;
	      if(ss < 0 || ss >= seg->depth) continue;
	      if(abs(dc)+abs(dr)+abs(ds)>topo) continue;
	      int nbrlabel = MRIgetVoxVal(seg,cc,rr,ss,0);
	      if(voxlabel != nbrlabel) {
		// dont use voxlabel here because label 0 will be 0
		MRIsetVoxVal(border,c,r,s,0,1);
		MRIsetVoxVal(border,c,r,s,1+knbr,nbrlabel);
		knbr++;
		//dc = 2; dr = 2; ds = 2; // break
		nhits++;
	      }
	    }//ds
	  }//dr
	}//dc

      }//s
    }//r
  }//c
  printf("nhits = %d, maxlabel=%d\n",nhits,*maxlabel);
  return(border);
}
#endif

std::vector<MRI*> MRIarraySmooth(std::vector<MRI*>y, std::vector<double> fwhm, MRI *mask, MRIS *surf)
{
  std::vector<MRI*> emptyVec;
  int ny = y.size();
  std::vector<MRI*> yout(ny);
  for(int n=0; n < ny; n++) {
    yout[n] = MRIcopy(y[n], NULL);
    if(!yout[n]) return(emptyVec);
  }
  // Check that fhwm{1,2}=0 with surf
  if(fwhm[0]>0 || fwhm[1]>0 || fwhm[2]>0){
    double gstdc = fwhm[0]/sqrt(log(256.0));
    double gstdr = fwhm[1]/sqrt(log(256.0));
    double gstds = fwhm[2]/sqrt(log(256.0));
    printf("std %g %g %g\n",gstdc,gstdr,gstds);
    for(int n=0; n < ny; n++){
      if(surf==NULL){
	MRI *tmpmri = MRIgaussianSmoothNI(yout[n], gstdc, gstdr, gstds, yout[n]);
	if(!tmpmri) return(emptyVec);
      } else {
	int niters = MRISfwhm2niters(fwhm[0],surf);
	if(n==0) printf("Smoothing input by fwhm=%lf niters=%d \n",fwhm[0],niters);
	MRI *tmpmri = MRISsmoothMRI(surf, yout[n], niters, mask, yout[n]);
	if(!tmpmri) return(emptyVec);
      }
    }
  }
  if(fwhm[3]>0){
    // Create a new volume for each frame putting the array dim as the first dim,
    // then col, row, and slice. Then filter in only the first dim. Then reverse
    // the process. Assumes all inputs have the same number of frames. This applies
    // to both vol and surface because this dimension is outside of both.
    MRI *tmpmri = MRIallocSequence(ny, y[0]->width,y[0]->height, MRI_FLOAT, y[0]->depth);
    tmpmri->xsize = y[0]->tr;
    double gstdc = fwhm[3]/sqrt(log(256.0));
    printf("Array smoothing ny=%d, tr=%g, fwhm=%g, gstd=%g\n",ny,y[0]->tr,fwhm[3],gstdc);
    for(int f=0; f < y[0]->nframes; f++){
      for(int n=0; n < ny; n++) {
	for(int c=0; c < y[0]->width; c++){
	  for(int r=0; r < y[0]->height; r++){
	    for(int s=0; s < y[0]->depth; s++){
	      double val = MRIgetVoxVal(yout[n],c,r,s,f);
	      MRIsetVoxVal(tmpmri,n,c,r,s,val);
	    }
	  }
	}
      }
      MRI *tmpmri2 = MRIgaussianSmoothNI(tmpmri, gstdc, 0, 0, tmpmri);
      if(!tmpmri2) return(emptyVec);
      for(int n=0; n < ny; n++) {
	for(int c=0; c < y[0]->width; c++){
	  for(int r=0; r < y[0]->height; r++){
	    for(int s=0; s < y[0]->depth; s++){
	      double val = MRIgetVoxVal(tmpmri,n,c,r,s);
	      MRIsetVoxVal(yout[n],c,r,s,f,val);
	    }
	  }
	}
      }
    } // frame
  }

  return(yout);
}

#if 0
// Gaussian smoothing in the temporal (frame) dimension. tstd is the
// stddev of the gaussian in units of the src->tr (typically ms). Note
// that the normalization uses normtype=2 (if not specified) to force
// the rows the of the G matrix to sum to 1
MRI *MRItemporalGSmooth(MRI *src, MRI *mask, double tstd, MRI *targ, int normtype)
{
  if(targ == NULL) {
    targ = MRIallocSequence(src->width, src->height, src->depth, MRI_FLOAT, src->nframes);
    if(targ == NULL) {
      printf("ERROR: MRItemporalGSmooth(): could not alloc\n");
      return (NULL);
    }
    MRIcopy(src, targ);
    MRIcopyPulseParameters(src,targ);
  }
  
  printf("tstd = %g, tr=%g, nstd %g, normtype %d\n", tstd,src->tr, tstd/src->tr, normtype);
  MATRIX *G = GaussianMatrix(src->nframes, tstd/src->tr, normtype, NULL);
#ifdef HAVE_OPENMP
#pragma omp parallel for 
#endif
  for(int c = 0; c < src->width; c++){
    for(int r = 0; r < src->height; r++) {
      for(int s = 0; s < src->depth; s++) {
	if(mask && MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;
	MATRIX *v = MatrixAlloc(src->nframes, 1, MATRIX_REAL);
	for(int f=0; f < src->nframes; f++) v->rptr[f+1][1] = MRIgetVoxVal(src, c, r, s, f);
	MATRIX *vg = MatrixMultiply(G, v, NULL);
	for(int f=0; f < src->nframes; f++) MRIsetVoxVal(targ, c, r, s, f, vg->rptr[f+1][1]);
	MatrixFree(&v);
	MatrixFree(&vg);
      } // slice
    } // rows
  }// col
  MatrixFree(&G);
  return(targ);
}
#endif

/*----------------------------------------*/
int main(int argc, char **argv) 
{
  int threads = 20;
#ifdef HAVE_OPENMP
  omp_set_num_threads(threads);
#endif
  //Gdiag_no = 72533;
  MRI *mri;

  // dngtester 1=sphere 2=sphere.reg 3=lta 4=ovsurf 5=ovslice 6=dir <7=ctab>
  int direction = 0;
  sscanf(argv[6],"%d",&direction);
  MRISsphere2mri s2m;
  s2m.sphere = MRISread(argv[1]); // sphere
  s2m.surfreg = MRISread(argv[2]); // sphere.reg
  s2m.lta = LTAread(argv[3]); // lta
  if(s2m.lta->type != REGISTER_DAT){
    printf("Changing LTA to TKR\n");
    LTAchangeType(s2m.lta, REGISTER_DAT);
  }
  LTAmat2RotMat(s2m.lta);
  MRISltaMultiply(s2m.sphere,s2m.lta);
  if(direction == 1){
    s2m.nvox[0] = 400;
    s2m.nvox[1] = 200;
    s2m.P0[0] = -90; s2m.P0[1] = -45;
    s2m.P1[0] = +90; s2m.P1[1] = +45;
    mri = MRIread(argv[4]); //ovsurf is input, ovslice output
  } else {
    mri = MRIread(argv[5]); //ovslice is input, ovsurf output
  }
  if(argc > 7) mri->ct = CTABreadASCII(argv[7]);
  MRI *ss = s2m.sphere2mri(mri,direction);
  if(direction == 1)  MRIwrite(ss,argv[5]);
  else                MRIwrite(ss,argv[4]);

  exit(0);


  mri = MRIread(argv[1]);
  MRI *mrimask = MRIread(argv[2]);
  SpatTempCluster stc;
  stc.nbrtype = 3;
  printf("binmask\n");
  stc.GetBinMask(mri,1.3,1000,0,mrimask);
  printf("cluster\n");
  stc.Clusterize();
  printf("write\n");
  MRIwrite(stc.cnomap,argv[3]);
  FILE *fp = fopen(argv[4],"w");
  stc.PrintClusterSum(fp);
  fclose(fp);

  exit(0);

  mri = MRIread(argv[1]);
  double tstd;
  sscanf(argv[2],"%lf",&tstd);
  MRI *tmri = MRItemporalGSmooth(mri, NULL, tstd, NULL);
  MRIwrite(tmri,argv[3]);
  exit(0);

  std::vector<MRI*> yarray;
  std::vector<double> fwhm(4);
  char *outbase = argv[1];
  for(int n=2; n < 6; n++) {
    sscanf(argv[n],"%lf",&fwhm[n-2]);
    printf("%d %lf\n",n,fwhm[n-2]);
  }
  for(int n=6; n < argc; n++){
    printf("Reading %d %s\n",n,argv[n]);
    MRI *mritmp = MRIread(argv[n]);
    yarray.push_back(mritmp);
  }
  std::vector<MRI*> yout = MRIarraySmooth(yarray, fwhm, NULL, NULL);
  for(int n = 0; n < yout.size(); n++){
    char tmpstr[2000];
    sprintf(tmpstr,"%s.%03d.nii.gz",outbase,n);
    MRIwrite(yout[n],tmpstr);
  }

  exit(0);

  int maxlabel;
  int segid = 17;
  mri = MRIread(argv[1]);
  Timer mytimer;
  MRI *sb = MRIsegBorder(mri,&maxlabel,1);
  int msecBorder = mytimer.milliseconds();
  printf("%d\n",msecBorder);
  MRIwrite(sb,"segborder.mgz");
  MRI *mri_vals = MRIread(argv[2]);
  printf("starting with sb ================\n");
  mytimer.reset();
  MRIvoxelsInLabelWithPartialVolumeEffects(mri,mri_vals,segid,NULL,NULL,maxlabel+1,sb);
  //msecBorder = mytimer.milliseconds();
  //printf("%d\n",msecBorder);
  printf("starting without sb ================\n");
  mytimer.reset();
  MRIvoxelsInLabelWithPartialVolumeEffects(mri,mri_vals,segid,NULL,NULL,maxlabel+1,NULL);
  //msecBorder = mytimer.milliseconds();
  //printf("%d\n",msecBorder);
  exit(0);

  mri = MRIread(argv[1]);
  SpatTempCluster stc2;
  stc.GetBinMask(mri,0.1,std::numeric_limits<double>::infinity(),0,NULL);
  stc.surf = MRISread(argv[2]);
  stc.topo = 2;
  stc.Clusterize();
  MRIwrite(stc.cnomap,argv[3]);
  MRIwrite(stc.binmask,"binmask.mgz");
  stc.PrintClusterSum(stdout);
  exit(0);



  
  MRIS *surf;
  std::vector<int> segids;
  double dist=0;
  LABEL *label2;

  // define an uber label in fsaverage space
  // map label to individual mris_apply_reg
  // convert label to a patch using infalted (label2patch)
  // flatten the patch
  // run this program

  // surf patch label overlay out
  MRISflat2mri f2m;
  f2m.nvox[0] = 100;
  f2m.nvox[1] = 100;
  f2m.voxsize[0] = 0.5;
  f2m.voxsize[1] = 0.5;
  //f2m.vno1 = 95610;
  //f2m.vno2 = 59038;
  f2m.surf = MRISread(argv[1]);
  MRISreadPatch(f2m.surf, argv[2]);
  f2m.surfreg = MRISread(argv[3]);
  //f2m.label = LabelRead(NULL, argv[3]);
  MRI *ov = MRIread(argv[4]);
  MRI *flat = f2m.flat2mri(ov);
  MRIwrite(flat,argv[5]);

  exit(0);
  RapheSeg rs;
  rs.aseg0 = MRIread(argv[1]);
  rs.dasb0 = MRIread(argv[2]);
  rs.dr0 = MRIread(argv[3]);
  rs.mr0 = MRIread(argv[4]);
  MRI *rseg = rs.doit();
  MRIwrite(rseg,argv[5]);
  //MRI *segaq = MRIaseg2aqueduct(mri);
  //MRIsegRaphe(segaq,dasb);

  exit(0);

  surf = MRISread(argv[1]);
  mri = MRIread(argv[2]);
  std::vector<int> segidlist;
  segidlist.push_back(1);
  segidlist.push_back(30);
  LABEL *lab = MRISseg2Label(surf, mri, segidlist);
  LabelWrite(lab, argv[3]);

  exit(0);
  
  // surf label dist outsurf outmask
  surf = MRISread(argv[1]);
  label2 = LabelRead(NULL,argv[2]);
  sscanf(argv[3],"%lf",&dist);
  //mri = SimAtrophy(surf,vno,nhops,navgs);
  mri = SimAtrophyLabel(surf, label2, dist);
  printf("Saving to %s\n",argv[4]);
  MRISwrite(surf,argv[4]);
  MRIwrite(mri,argv[5]);
  exit(0);

  MRISsetNeighborhoodSizeAndDist(surf,2);
  mrisComputeTangentPlanes(surf);

  //vno = 156791; //138268;
  int vno = 1000;
  Curv2(surf, vno);
  //MRIScomputeSecondFundamentalFormAtVertex2(surf, vno);
  exit(0);

  for(vno=0; vno < surf->nvertices; vno++){
    Curv2(surf, vno);
    double k1 = surf->vertices[vno].k1;
    double k2 = surf->vertices[vno].k2;
    MRIScomputeSecondFundamentalFormAtVertex2(surf, vno);
    if(fabs(surf->vertices[vno].k1-k1)> .01 || fabs(surf->vertices[vno].k2-k2)> .001){
      VERTEX_TOPOLOGY *vtop = &(surf->vertices_topology[vno]);
      printf("%d  %6.3f %6.3f   %6.3f %6.3f   %3d\n",vno,surf->vertices[vno].k1,k1,surf->vertices[vno].k2,k2,vtop->vtotal);
    }
  }


  exit(0);

  surf = MRISread(argv[1]);
  int nremove, nstep, nremoved=0,iter=0;
  sscanf(argv[2],"%d",&nremove);
  sscanf(argv[3],"%d",&nstep);
  while(nremoved < nremove){
    iter++;
    int nstep2 = nstep;
    if(nstep > (nremove-nremoved)) nstep2 = nremove-nremoved;
    printf("\niter %3d %4d %4d =========###############==================\n",iter,nremoved,nstep2);
    int nr = remesh(surf, nstep2);
    char tmpstr[2000];
    sprintf(tmpstr,"lh.iter%d",iter);
    MRISwrite(surf,tmpstr);
    nremoved += nr;
  }
  printf("nvertices %d nfaces %d\n",surf->nvertices,surf->nfaces);
  MRISwrite(surf,argv[4]);
  exit(0);

  surf = MRISread(argv[1]);
  int rfaceno[2];
  MRIScollapseEdge(surf, 78029, 78031,rfaceno);
  MRISremoveVertex(surf,150627);
  MRISremoveFace(surf,rfaceno[0]);
  MRISremoveFace(surf,rfaceno[1]);
  MRISwrite(surf,argv[2]);
  exit(0);

  surf = MRISread(argv[1]);
  MRISmarkEdge(surf, NULL, 2, 160, 1);
  mri = MRIcopyMRIS(NULL, surf, 0, "marked");
  MRIwrite(mri,argv[2]);
  MRISremoveHighAngleHinges(surf, 140, 1, 2, 40, NULL);
  MRISwrite(surf,argv[3]);
  mri = MRIcopyMRIS(NULL, surf, 0, "marked");
  MRIwrite(mri,argv[4]);
  exit(0); //------------------------------------

  surf = MRISread(argv[1]);
  mri = MRIread(argv[2]);
  RemoveLeaks(surf, mri, 44, 100);
  MRISwrite(surf,argv[3]);

  int nhits = 0;
  for(int vtxno=0; vtxno < surf->nvertices; vtxno++) {
    VERTEX *v = &surf->vertices[vtxno];
    v->ripflag = !v->marked;
    if(!v->ripflag) nhits++;
  }
  printf("nhits %d\n",nhits);
  INTEGRATION_PARMS parms;
  parms.fill_interior = 0 ;
  parms.projection = NO_PROJECTION ;
  parms.tol = 1e-4 ;
  parms.dt = 0.5f ;
  parms.base_dt = parms.dt ;
  parms.integration_type = INTEGRATE_MOMENTUM ;
  parms.momentum = 0.0 /*0.8*/ ;
  parms.dt_increase = 1.0 /* DT_INCREASE */;
  parms.dt_decrease = 0.50 /* DT_DECREASE*/ ;
  parms.error_ratio = 50.0 /*ERROR_RATIO */;
  parms.niterations = 30;

  //parms.l_location  = 0.500;
  //MRISsaveVertexPositions(surf, TARGET_VERTICES);

  //parms.l_spring_nzr = 10;
  //parms.l_spring_nzr_len = .1;

  parms.l_repulse = 5;


  MRI *involPS = MRIread(argv[4]);
  MRISpositionSurface(surf, involPS, involPS, &parms);
  MRISwrite(surf,argv[5]);


  exit(0);

  surf = MRISread(argv[1]);
  mri = MRIread(argv[2]);
  sscanf(argv[3],"%d",&vno);
  NormalDispersion(surf, mri, vno);
  exit(1);

  mri = MRIread(argv[1]);
  int erode = 1;
  if(erode > 0){
    printf("  Eroding %d\n",erode);
    MRI *mritmp = NULL;
    for(int i=0; i<erode; i++) {
      mritmp = MRIerodeNN(mri,mritmp,NEAREST_NEIGHBOR_FACE); //NEAREST_NEIGHBOR_CORNER
      MRIcopy(mritmp,mri);
    }
    MRIfree(&mritmp);
  }
  MRIwrite(mri,argv[2]);
  exit(0);


  if(1){
  surf = MRISread(argv[1]);
  sscanf(argv[2],"%d",&vno);
  surf->vertices[vno].marked=1;
  //MRISdilateMarked(surf,1);
  MRISnotMarked(surf);
  Gdiag_no = vno;
  MRISsoapBubbleVertexPositions(surf, 100);
  MRISwrite(surf,argv[3]);
  exit(0);
  }


  surf = MRISread(argv[1]);
  MRIScomputeMetricProperties(surf);
  MRISfaceMetric(surf,0);
  mri  = MRIread(argv[2]);
  MarkFaces(surf, mri, 1, -0.9);
  MRI *mriA = MRIcopyMRIS(NULL, surf, 0, "marked");
  MRIwrite(mriA,argv[3]);
  exit(0);




  mri = MRIread(argv[1]);
  surf = MRISread(argv[2]);
  MRI *surfseg = MRIread(argv[3]);
  LTA *reg = LTAread(argv[4]);
  MRI *mriseg = MRISmapAnnot(mri,surf, reg, surfseg, segids, -1, +1, .2, NULL);
  MRIwrite(mriseg,argv[5]);

  exit(0);


  int nvertices, nfaces, nedges,eno;

  surf = MRISread(argv[1]);
  LABEL *label = LabelRead(NULL,argv[2]);
  for(int n=0; n < label->n_points; n++) surf->vertices[label->lv[n].vno].marked=1;
  MRISnotMarked(surf);  // turn off->on and on->off so soap bubble is correct (marked are fixed)
  MRISsoapBubbleVertexPositions(surf, 100);
  MRISwrite(surf,argv[3]);

  exit(0);
  int crsFoV[3];
  MRI *crop;
  double crsCenter[3], rasCenter[3];
  LTA *lta=NULL;
  std::vector<int>iKeep;
  iKeep.push_back(7);

  mri = MRIread(argv[1]);
  sscanf(argv[2],"%lf",&rasCenter[0]);
  sscanf(argv[3],"%lf",&rasCenter[1]);
  sscanf(argv[4],"%lf",&rasCenter[2]);
  sscanf(argv[5],"%d",&crsFoV[0]);
  sscanf(argv[6],"%d",&crsFoV[1]);
  sscanf(argv[7],"%d",&crsFoV[2]);
  lta = LTAread(argv[8]);
  crop = MRIcropAroundRAS(mri,rasCenter,crsFoV,lta,iKeep);
  if(!crop) exit(1);
  MRIwrite(crop,argv[9]);
  exit(0);

  mri = MRIread(argv[1]);
  sscanf(argv[2],"%lf",&crsCenter[0]);
  sscanf(argv[3],"%lf",&crsCenter[1]);
  sscanf(argv[4],"%lf",&crsCenter[2]);
  sscanf(argv[5],"%d",&crsFoV[0]);
  sscanf(argv[6],"%d",&crsFoV[1]);
  sscanf(argv[7],"%d",&crsFoV[2]);
  crop = MRIcropAroundCRS(mri,crsCenter,crsFoV);
  if(!crop) exit(1);
  MRIwrite(crop,argv[8]);
  exit(0);

  printf("argc = %d\n",argc);
  MATRIX *R=NULL;
  MATRIX *scras = NULL;
  MATRIX *tkras = MatrixAlloc(4,1,MATRIX_REAL);
  tkras->rptr[4][1] = 1;
  for(int n=2; n < argc; n++){
    printf("%3d %s\n",n,argv[n]);
    surf = MRISread(argv[n]);
    if(surf == NULL) exit(1);
    if(n==2){
      mri = MRIallocSequence(surf->nvertices,argc-2,1,MRI_FLOAT,3);
      R = surf->vg.get_TkregRAS2RAS();
    }
    for(int vno=0; vno < surf->nvertices; vno++){
      VERTEX *v = &surf->vertices[vno];
      tkras->rptr[1][1] = v->x;
      tkras->rptr[2][1] = v->y;
      tkras->rptr[3][1] = v->z;
      scras = MatrixMultiply(R,tkras,scras);
      MRIsetVoxVal(mri,vno,n-2,0,0,scras->rptr[1][1]);
      MRIsetVoxVal(mri,vno,n-2,0,1,scras->rptr[2][1]);
      MRIsetVoxVal(mri,vno,n-2,0,2,scras->rptr[3][1]);
    }
    MRISfree(&surf);
  }
  MRIwrite(mri,argv[1]);
  exit(0);



  //===============================================================
  mri = MRIread(argv[1]); // 
  MRI *mri2 = MRIread(argv[2]); // 
  MRI *mri3 = MRIread(argv[3]); // 
  double thresh;
  sscanf(argv[4],"%lf",&thresh);
  MRI *emask = MRIgetEditMask(mri, mri2, mri3, thresh, 25);
  MRIwrite(emask,argv[5]);
  exit(0);

  mri = MRIread(argv[1]);
  CountFaces(mri);
  exit(0);

  mri = MRIalloc(3,3,3,MRI_INT);
  MRIsetVoxVal(mri,1,1,1,0,1);
  CountFaces(mri);
  exit(0);

  surf = MRIStessellate(mri,1,0);
  eno = MRIScomputeEulerNumber(surf, &nvertices, &nfaces, &nedges) ;
  printf("v=%d e=%d f=%d   eno=%d\n",nvertices, nedges, nfaces, eno) ;

  MRIsetVoxVal(mri,6,6,5,0,1);
  surf = MRIStessellate(mri,1,0);
  eno = MRIScomputeEulerNumber(surf, &nvertices, &nfaces, &nedges) ;
  printf("v=%d e=%d f=%d   eno=%d\n",nvertices, nedges, nfaces, eno) ;


  exit(0);

  mri = MRIread(argv[1]);
  MRI *out = MRIerodeNN(mri,NULL,NEAREST_NEIGHBOR_FACE,0);
  MRIwrite(out,argv[2]);
  MRI *out2 = MRIerodeNN(mri,NULL,NEAREST_NEIGHBOR_FACE,1);
  MRIwrite(out2,argv[3]);
  exit(0);

  surf = MRISread(argv[1]);
  //SmoothSurf(surf, 1);
  MRISwrite(surf,argv[2]);

  return(0);
}


#if 0

int RapheSeg::AcqueductSeg(void)
{
  if(aseg->type == MRI_UCHAR){
    printf("Changing aseg type from uchar to int\n");
    MRI *tmp = MRISeqchangeType(aseg, MRI_INT, 0.0, 0.999, 1);
    MRIfree(&aseg);
    aseg = tmp;
  }

  if(aseg0) MRIfree(&aseg0); // or keep aseg0 if set?
  aseg0 = MRIcopy(aseg,NULL);

  // Reorient to RAS
  char ostr[5];
  ostr[4] = '\0';
  MRIdircosToOrientationString(aseg,ostr);
  printf("Orientation %s\n",ostr);
  if(strcmp(ostr,"RAS")){
    printf("Reorientating to RAS\n");
    MRI *rastemplate = MRIcopy(aseg,NULL);
    MRIorientationStringToDircos(rastemplate, "RAS");
    MRI *tmp = MRIresample(aseg, rastemplate, SAMPLE_NEAREST);
    MRIfree(&aseg);
    aseg = tmp;
    tmp = MRIresample(t1, rastemplate, SAMPLE_NEAREST);
    MRIfree(&t1);
    t1 = tmp;
    MRIfree(&rastemplate);
  }

  // Update the ctab
  if(!aseg->ct) aseg->ct = CTABreadDefault();
  CTE *cte = aseg->ct->entries[540];
  if(cte == NULL) {
    cte = (CTE *)malloc(sizeof(CTE));
    aseg->ct->entries[540] = cte;
  }
  strcpy(cte->name,"Cerebral-Aqueduct");
  cte->ri=170; cte->gi=85; cte->bi=255; cte->ai=0;
  cte->rf=170/255.0; cte->gf=85/255.0; cte->bf=255/255.0;
  aseg0->ct = CTABdeepCopy(aseg->ct);

  // Get extent limits of 4th vent (and brainstem while there)
  int V4APmin=aseg->height,V4APmax=0,V4SImin=aseg->depth,V4SImax=0; // vox
  int nhits=0, nbs=0, nv4=0;
  double bssum=0, v4sum=0;
  for(int c=0; c < aseg->width; c++){
    for(int r=0; r < aseg->height; r++){
      for(int s=0; s < aseg->depth; s++){
	int segid = MRIgetVoxVal(aseg,c,r,s,0);
	if(segid == 15){
	  nv4++;
	  if(t1) v4sum += MRIgetVoxVal(t1,c,r,s,0);
	}
	if(segid == 16){
	  nbs++;
	  if(t1) bssum += MRIgetVoxVal(t1,c,r,s,0);
	}
	if(segid == 16 || segid == 15){
	  // get vox bounding box for brainstem/4th vent
	  if(bsCorner1[0] > c) bsCorner1[0] = c;
	  if(bsCorner1[1] > r) bsCorner1[1] = r;
	  if(bsCorner1[2] > s) bsCorner1[2] = s;
	  if(bsCorner2[0] < c) bsCorner2[0] = c;
	  if(bsCorner2[1] < r) bsCorner2[1] = r;
	  if(bsCorner2[2] < s) bsCorner2[2] = s;
	}
	if(segid != 15) continue; //15 = 4th vent
	nhits++;
	if(V4APmin > r) V4APmin = r;
	if(V4APmax < r) V4APmax = r;
	if(V4SImin > s) V4SImin = s;
	if(V4SImax < s) V4SImax = s;
      }
    }
  }
  printf("nhits = %d, AP=(%d,%d) SI=(%d,%d)\n",nhits,V4APmin,V4APmax,V4SImin,V4SImax);
  printf("BS (%d,%d,%d) (%d,%d,%d)\n",bsCorner1[0],bsCorner1[1],bsCorner1[2],bsCorner2[0],bsCorner2[1],bsCorner2[2]);
  double bsmean = bssum/nbs;
  double v4mean = v4sum/nv4;
  printf("nbs=%d, bsmean=%g, nv4=%d, v4mean=%g\n",nbs,bsmean,nv4,v4mean);

  // Starting at the most superior slice. Examine each slice to see
  // where the acqueduct is not totally surrounded by brainstem.
  int aqfound=0;
  for(int s=V4SImax; s >= V4SImin; s--){
    double t1mn = 0;
    int nbs = 0;
    for(int c=0; c < aseg->width; c++){
      for(int r=0; r < aseg->height; r++){
	int segid = MRIgetVoxVal(aseg,c,r,s,0);
	if(segid == 16) {
	  if(t1) t1mn += MRIgetVoxVal(t1,c,r,s,0);
	  nbs ++;
	}
      }
    }
    t1mn /= nbs;
    // Get single slice binarizerd to non-brainstem. The 4thvent/aqd
    // will show up as an island
    MRI *slice = MRIalloc(aseg->width,aseg->height,1,MRI_INT);
    int nv4 = 0;
    int V4found=0;
    for(int c=0; c < aseg->width; c++){
      for(int r=0; r < aseg->height; r++){
	int segid = MRIgetVoxVal(aseg,c,r,s,0);
	double t1val = 0;
	if(t1) t1val = MRIgetVoxVal(t1,c,r,s,0);
	int val=1; // not brainstem
	if(segid == 16){ // in brainstem
	  if(0 && t1val < t1mn*t1thresh) val = 1;//dark enough to change label
	  else val=0; // keep as brainstem
	}
	if(val) nv4++;
	MRIsetVoxVal(slice,c,r,0,0,val);
	if(segid == 15) V4found = 1;// should always be true if using V4 slice limits
      }
    }
    if(nv4==0){ // No BS voxels, probably does not happen
      MRIfree(&slice);
      break;
    }
    // Clusterize. 
    int nClusters;
    VOLCLUSTER **ClusterList = clustGetClusters(slice, 0, 0.5, 0, 0, 0, NULL, &nClusters, NULL);
    if(nClusters > 1){
      // Look through all the clusters and relable the 4th vent voxels as acqueduct
      for(int cno = 1; cno < nClusters; cno++){
	VOLCLUSTER *vc = ClusterList[cno];
	for(int n = 0; n < vc->nmembers; n++){
	  int segid = MRIgetVoxVal(aseg,vc->col[n],vc->row[n],s,0);
	  if(segid != 15 && segid != 16) continue;
	  MRIsetVoxVal(aseg, vc->col[n], vc->row[n], s, 0, 540); // acqueduct
	  aqfound = 1;
	}
      }
    }
    if(s==83) {
      //MRIwrite(slice,"slice.mgz");
      //printf(" slice %g\n",MRIgetVoxVal(aseg,131,102,s,0));
    }
    MRIfree(&slice);
    printf("s=%d nv4=%d nc=%d nbs=%d t1mn=%g aqfound=%d V4found=%d\n",s,nv4,nClusters,nbs,t1mn,aqfound,V4found);
    clustFreeClusterList(&ClusterList, nClusters);
    if(aqfound && V4found && nClusters <= 1) break;
  }

  //MRIwrite(aseg,"dng.mgz");

  if(strcmp(ostr,"RAS")){
    printf("Reorientating back to %s\n",ostr);
    MRI *tmp = MRIresample(aseg, aseg0, SAMPLE_NEAREST);
    MRIfree(&aseg);
    aseg = tmp;
    aseg->ct = CTABdeepCopy(aseg0->ct);
  }

  //MRIwrite(aseg,"aq.aseg.mgz");

  return(0);
}
#endif



#if 0   // ok to delete this
  int stindex2index(std::vector<int> stindex){
    int nvox = this->binmask->width * this->binmask->height * this->binmask->depth;
    int index = stindex[0] + stindex[1]*nvox;
    return(index);
  }
  std::vector<int> index2stindex(int index){
    int nvox = this->binmask->width * this->binmask->height * this->binmask->depth;
    int tindex = floor((double)index/nvox);
    int sindex = index - tindex*nvox;
    std::vector<int> stindex = {sindex,tindex};
    return(stindex);
  }
  std::vector<int> sindex2crs(int sindex){
    int w = this->binmask->width;
    int h = this->binmask->height;
    int s = floor((double)sindex/(w*h));
    int r = (double(sindex-s*w*h)/w);
    int c = sindex-r*w-s*w*h;
    std::vector<int> crs = {c,r,s};
    return(crs);
  }
#endif
