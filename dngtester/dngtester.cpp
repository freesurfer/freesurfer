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
#include "surfgrad.h"

#include "romp_support.h"
#undef private

//int MRISfaceNormalFace_AddDeltaVertex = -1;
//long double MRISfaceNormalFace_AddDelta[3]={0,0,0};

double MRISfaceMinEdge(MRIS *surf, int faceno);
DMATRIX *MRISvertex2Matrix(MRIS *surf, int vtxno, DMATRIX *vm);
int MRISmatrix2Vertex(MRIS *surf, int vtxno, DMATRIX *vm);
int MRISscaleVertices(MRIS *surf, double scale);
int OldTestBBRCostFace(BBRFACE *bbrf, BBRPARAMS *bbrpar, int wrtvtxno);

DMATRIX *AllocGradCost(int nvertices);
int FreeGradCost(int nvertices, DMATRIX **pgradCost);
int UpdateVertexPosition(MRIS *surf, double stepsize, DMATRIX *grad);

long double ConjGradBeta(int method, DMATRIX *SteepDir, DMATRIX *PrevSteepDir, DMATRIX *PrevCGDir);
DMATRIX *ConjGradDir(int method, DMATRIX *SteepDir, DMATRIX *PrevSteepDir, DMATRIX *PrevCGDir, DMATRIX *CGDir);
int PlaceSurf(MRI *mri, MRIS *surf);


/*----------------------------------------*/
int main(int argc, char **argv) 
{
  MRIS *surf, *surf2;
  int msec, nvertices; //vtxno=0;
  Timer mytimer;
  MRI *mri, *mriindex, *mri2;
  Geodesics *geod;
  float maxdist;
  double d, dmax;
  int c,r,s;
  LABEL2SURF *l2s;
  LTA *lta;
  DMATRIX *dJ, *gnum;
  VERTEX *v;
  double d1, delta;

#ifdef HAVE_OPENMP
  omp_set_num_threads(10);
#endif

  UnitizeNormalFace = 1;
  //surf = MRISread("/homes/4/greve/subjects/vrfp-mar14-anat/surf/lh.white");
  //surf = MRISread("/homes/4/greve/l/sp1/fsh.github.local/lib/bem/ic1.tri");
  //MRISscaleVertices(surf, 100);
  //MRISwrite(surf,"./lh.tmp0");
  //surf = MRISread(argv[1]);

  //printf("Reading volume\n"); fflush(stdout);
  //mri = MRIread("/homes/4/greve/subjects/vrfp-mar14-anat/mri/orig.mgz");
  //mri = MRIread("/autofs/cluster/fsm/users/greve/subjects/t1fit.fsm010/resample/coreg.0037.T2w_SPC_vNav.conf.mgz");
  mri = MRIread("/autofs/cluster/fsm/users/greve/subjects/fsm010/mri/norm.mgz");

  //printf("Reading surface\n"); fflush(stdout);
  //surf = MRISread("/homes/4/greve/subjects/vrfp-mar14-anat/surf/lh.white");
  surf = MRISread("/autofs/cluster/fsm/users/greve/subjects/fsm010/surf/lh.pial");


  MRISedges(surf);
  MRIScomputeMetricProperties(surf);
  MRISfaceNormalGrad(surf, 0);
  PlaceSurf(mri, surf);
  exit(0);

  BBRFACE *bbrf = BBRFaceAlloc();
  BBRPARAMS *bbrpar = (BBRPARAMS *) calloc(sizeof(BBRPARAMS),1);
  bbrpar->Din  = 0.5;
  bbrpar->Dout = 0.5;
  bbrpar->M = 0.5;
  bbrpar->Q0 = -10.0;
  bbrpar->mri = mri;
  bbrpar->surf = surf;
  BBRPARsras2vox(bbrpar);
  bbrpar->interp = SAMPLE_CUBIC;
  //bbrpar->interp = SAMPLE_TRILINEAR;

  dJ = AllocGradCost(surf->nvertices);
  MRIScomputeMetricProperties(surf);
  MRISfaceNormalGrad(surf, 0);
  d = MRISedgeCost(surf, dJ);
  //d = MRISbbrCost(bbrpar, dJ);
  gnum = DMatrixAlloc(1,3,MATRIX_REAL);
  delta = .0001;
  for(s=0; s < surf->nvertices; s++){
    v = &(surf->vertices[s]);
    for(c=0; c<3; c++){
      double dx = (c==0) ? delta : 0.0;
      double dy = (c==1) ? delta : 0.0;
      double dz = (c==2) ? delta : 0.0;
      MRISsetXYZ(surf,s,
        v->x + dx,
        v->y + dy,
        v->z + dz);
      MRISfaceNormalGrad(surf, 0);
      //d1 = MRISbbrCost(bbrpar, NULL);
      d1 = MRISedgeCost(surf, NULL);
      gnum->rptr[1][c+1] = (d1-d)/delta;
      MRISsetXYZ(surf,s,
        v->x - dx,
        v->y - dy,
        v->z - dz);
    }
    printf("#@# %d ",s);
    for(c=0; c<3; c++) printf("%12.10lf ",100000*(dJ->rptr[s+1][c+1]-gnum->rptr[1][c+1]));
    printf("\n");
    printf("g0   = "); for(c=0; c<3; c++) printf("%12.8lf ",100000*dJ->rptr[s+1][c+1]);  printf("\n");
    printf("gnum = "); for(c=0; c<3; c++) printf("%12.8lf ",100000*gnum->rptr[1][c+1]);  printf("\n");
    fflush(stdout);
  }
  exit(0); //--------------------------------------

  s = 0;
  dmax = 0;
  for(c=0; c < surf->nfaces; c++){
    for(r=0; r < 3; r++){
      d = TestBBRCostFace(bbrpar, c, r, .0001, 0);
      if(dmax < d) dmax = d;
      printf("#@# %6d %d %g %g\n",c,r,d,dmax); fflush(stdout);
      if(d > .05){
	printf("#@# %6d %d %g\n",c,r,d);
	TestBBRCostFace(bbrpar, c, r, .0001, 1);
	fflush(stdout);
	s = s + 1;
      }
    }
  }
  printf("dmax = %g, nhits = %d\n",dmax,s);
  exit(1);

  
  d = MRISbbrCost(bbrpar, dJ);
  for(s=0; s < 10000; s++){
    MRISfaceNormalGrad(surf, 0);
    d = MRISbbrCost(bbrpar, dJ);
    printf("%d %lf ----------\n",s,d); fflush(stdout);
    UpdateVertexPosition(surf, 0.001, dJ);
  }
  MRISwrite(surf,"./lh.tmp");
  exit(0); //--------------------------------------

  for(c=0; c < 3; c++){
    printf("%d ------------------------------\n",c);
    MRIScomputeMetricProperties(bbrpar->surf);
    MRISfaceNormalGrad(bbrpar->surf, 1); // dont compute grad (norm only)
    bbrf = BBRCostFace(82, -1, bbrpar, bbrf); // Compute Cost
    bbrf = BBRCostFace(82, c, bbrpar, bbrf); // Compute grads
    BBRFacePrint(stdout, bbrf, bbrpar);
    TestBBRCostFace(bbrpar, 82, c, 0.001, 0); // wrtvtxno=0 (faceno=1)
  }
  BBRFaceFree(&bbrf);
  DMatrixFree(&bbrpar->sras2vox);
  DMatrixFree(&bbrpar->sras2vox3x3);
  MRIfree(&mri);
  MRISfree(&surf);
  exit(0); //--------------------------------------

  sscanf(argv[1],"%lf",&delta);
  printf("delta = %lf\n",delta);
  surf = MRISread(argv[2]);


  vg_isEqual_Threshold = 10e-4;

  c = L2Stest("fsf03anat");
  exit(c);

  if(0){
    mri  = MRIread(argv[1]);
    surf = MRISread(argv[2]);
    surf2 = MRISread(argv[3]);
    lta = LTAread(argv[4]);
  }else{
    mri  = MRIread("template.nii.gz");
    surf = MRISread("lh.white");
    surf2 = MRISread("rh.white");
    lta = LTAread("register.dof6.lta");
  }
  printf("\n");
  printf("alloc \n");
  l2s = L2Salloc(2, "");
  l2s->mri_template = mri;
  l2s->surfs[0] = surf;
  l2s->surfs[1] = surf2;
  l2s->dmax = 3;
  l2s->hashres = 16;
  l2s->vol2surf = lta;
  l2s->debug = 0;
  //l2s->nhopsmax = 10;
  printf("init \n");
  L2Sinit(l2s);
  printf("nhopsmax %d\n",l2s->nhopsmax);
  printf("loop \n");
  s = 17;
  for(c = 0; c < mri->width; c++){
    //printf("%2d ",c); fflush(stdout);
    for(r = 0; r < mri->height; r++) L2SaddPoint(l2s, c, r, s, 0, 1);
    //printf("\n");
  }
  // remove, erase a few
  //for(c = 10; c < mri->width-10; c++) L2SaddPoint(l2s, c, 15, s, 0, 0);

  // add surface vertices
  for(c = 0; c < 10000; c++) L2SaddPoint(l2s, c, -1, -1, 1, 1);

  // add a surface vertex that is already there
  L2SaddPoint(l2s, 110027, -1, -1, 1, 1);

  // remove some surface vertices
  L2SaddPoint(l2s, 111010, -1, -1, 1, 0);
  L2SaddPoint(l2s,   5000, -1, -1, 1, 0);


  LabelWrite(l2s->labels[0],"./my.label0");
  LabelWrite(l2s->labels[1],"./my.label1");
  MRIwrite(l2s->masks[0],"mask0.mgh");
  MRIwrite(l2s->masks[1],"mask1.mgh");
  L2Sfree(&l2s);

  exit(0);


  // good vertex on the lateral side 108489

  // apply smoothing: 5 args: surf geod input index output
  surf = MRISread(argv[1]);
  printf("reaDing geo\n"); fflush(stdout);
  mytimer.reset() ;
  geod = geodesicsRead(argv[2], &nvertices);
  msec = mytimer.milliseconds() ;
  printf("t = %g min\n",msec/(1000.0*60));
  mri  = MRIread(argv[3]);
  mriindex  = MRIread(argv[4]);
  printf("Smoothing\n");
  mytimer.reset() ;
  mri2 = GeoSmooth(mri, 10, surf, geod, mriindex, NULL);
  msec = mytimer.milliseconds() ;
  printf("t = %g min\n",msec/(1000.0*60));
  fflush(stdout);
  MRIwrite(mri2,argv[5]);
  mri2 = MRIcopyMRIS(NULL,surf,0,"val2bak");
  MRIwrite(mri2,"nnbrs.mgh");

  exit(0);  //-------------------------

  // check distances: surf geod
  surf = MRISread(argv[1]);
  d = MRISsphereDist(surf, &surf->vertices[220], &surf->vertices[3140]);
  geod = geodesicsRead(argv[2], &nvertices);
  geodesicsCheckSphereDist(surf, geod);
  exit(0);  //-------------------------

  // create geod file: 3 args: surf distmax output
  surf = MRISread(argv[1]);
  sscanf(argv[2],"%f",&maxdist);
  mytimer.reset() ;
  geod = computeGeodesics(surf, maxdist);
  msec = mytimer.milliseconds() ;
  printf("done t = %g min\n",msec/(1000.0*60));
  geodesicsWrite(geod, surf->nvertices, argv[3]);
  msec = mytimer.milliseconds() ;
  printf("done write t = %g min\n",msec/(1000.0*60));
  exit(0); //-------------------------


  //GeoDumpVertex("uvtx.108489.dist.dat", geod, 108489);
  mytimer.reset() ;
  geod = geodesicsReadV2(argv[1], &nvertices);
  msec = mytimer.milliseconds() ;
  printf(" read2 t = %g min\n",msec/(1000.0*60));
  exit(0); //-------------------------

  mytimer.reset() ;
  geodesicsWrite(geod, nvertices, "u1.geod");
  msec = mytimer.milliseconds() ;
  printf(" write1 t = %g min\n",msec/(1000.0*60));

  mytimer.reset() ;
  geodesicsWriteV2(geod, nvertices, "u2.geod");
  msec = mytimer.milliseconds() ;
  printf(" write2 t = %g min\n",msec/(1000.0*60));

  exit(0); //-------------------------



}
double MRISfaceMinEdge(MRIS *surf, int faceno)
{
  FACE *f;
  VERTEX *vn, *vm;
  int n, m;
  double dmin,dmax,dv;

  f = &(surf->faces[faceno]);
  dmin = 10e10;
  dmax = 0;
  for(n = 0; n < 3; n++){
    vn = &(surf->vertices[f->v[n]]);
    m = n + 1;
    if(m > 2) m = 0;
    vm = &(surf->vertices[f->v[m]]);
    dv = sqrt((vn->x-vm->x)*(vn->x-vm->x) + (vn->y-vm->y)*(vn->y-vm->y) + (vn->z-vm->z)*(vn->z-vm->z));
    if(dmin > dv) dmin = dv;
    if(dmax < dv) dmax = dv;
  }
  return(dmin/dmax);
}




DMATRIX *MRISvertex2Matrix(MRIS *surf, int vtxno, DMATRIX *vm)
{
  if(vm==NULL) 
    vm = DMatrixAlloc(3,1,MATRIX_REAL);
  if(vm->rows != 3 || vm->cols != 1){
    printf("ERROR: MRISvertex2Matrix(): vm wrong dim %d %d\n",vm->rows,vm->cols);
    return(NULL);
  }
  vm->rptr[1][1] = surf->vertices[vtxno].x;
  vm->rptr[2][1] = surf->vertices[vtxno].y;
  vm->rptr[3][1] = surf->vertices[vtxno].z;
  return(vm);
}

int MRISmatrix2Vertex(MRIS *surf, int vtxno, DMATRIX *vm)
{
  if(vm==NULL){
    printf("ERROR: MRISmatrix2Vertex(): vm is NULL\n");
    return(1);
  }
  if(vm->rows != 3 || vm->cols != 1){
    printf("ERROR: MRISmatrix2Vertex(): vm wrong dim %d %d\n",vm->rows,vm->cols);
    return(1);
  }
  MRISsetXYZ(surf,vtxno,
    vm->rptr[1][1],
    vm->rptr[2][1],
    vm->rptr[3][1]);
  return(0);
}


int MRISscaleVertices(MRIS *surf, double scale)
{
  int vtxno;
  VERTEX *v;
  for(vtxno = 0; vtxno < surf->nvertices; vtxno++){
    v = &(surf->vertices[vtxno]);
    MRISsetXYZ(surf,vtxno,
      v->x * scale,
      v->y * scale,
      v->z * (scale*(1+fabs(v->z))));
  }
  return(0);
}


/*-----------------------------------------------------------------------------*/
DMATRIX *AllocGradCost(int nvertices)
{
  DMATRIX *gradCost;
  gradCost = DMatrixAlloc(nvertices,3,MATRIX_REAL);
  return(gradCost);
}
/*-----------------------------------------------------------------------------*/
int FreeGradCost(int nvertices, DMATRIX **pgradCost)
{
  DMatrixFree(pgradCost);
  return(0);
}

int UpdateVertexPosition(MRIS *surf, double stepsize, DMATRIX *grad)
{
  int vtxno;
  #ifdef HAVE_OPENMP
  #pragma omp parallel for 
  #endif
  for(vtxno = 0; vtxno < surf->nvertices; vtxno++){
    VERTEX *v = &(surf->vertices[vtxno]);
    MRISsetXYZ(surf,vtxno,
      v->x - (stepsize*grad->rptr[vtxno+1][1]),
      v->y - (stepsize*grad->rptr[vtxno+1][2]),
      v->z - (stepsize*grad->rptr[vtxno+1][3]));
  }
  return(0);
}



long double ConjGradBeta(int method, DMATRIX *SteepDir, DMATRIX *PrevSteepDir, DMATRIX *PrevCGDir)
{
  long double beta, ss, ssprev;
  int r,c;

  // Fletcher Reeves
  ss=0;
  for(r=1; r < SteepDir->rows; r++){
    for(c=1; c < SteepDir->cols; c++){
      ss += (SteepDir->rptr[r][c] * SteepDir->rptr[r][c]);
    }
  }
  ssprev=0;
  for(r=1; r < SteepDir->rows; r++){
    for(c=1; c < SteepDir->cols; c++){
      ssprev += (PrevSteepDir->rptr[r][c] * PrevSteepDir->rptr[r][c]);
    }
  }
  beta = ss/(ssprev+10e-10);

  return(beta);
}

DMATRIX *ConjGradDir(int method, DMATRIX *SteepDir, DMATRIX *PrevSteepDir, DMATRIX *PrevCGDir, DMATRIX *CGDir)
{
  long double beta;
  beta = ConjGradBeta(method, SteepDir, PrevSteepDir, PrevCGDir);
  CGDir = DMatrixAddMul(SteepDir, PrevCGDir, 1, beta, CGDir);
  return(CGDir);
}


int PlaceSurf(MRI *mri, MRIS *surf)
{
  DMATRIX *gradBBR, *gradEdge, *SteepDir=NULL;
  double costBBR, costEdge, cost;
  int n;

  BBRPARAMS *bbrpar = (BBRPARAMS *) calloc(sizeof(BBRPARAMS),1);
  bbrpar->Din  = 0.5;
  bbrpar->Dout = 0.5;
  bbrpar->M = 0.5;
  bbrpar->Q0 = -10.0;
  bbrpar->mri = mri;
  bbrpar->surf = surf;
  BBRPARsras2vox(bbrpar);
  //bbrpar->interp = SAMPLE_CUBIC;
  bbrpar->interp = SAMPLE_TRILINEAR;

  MRISfaceNormalGrad(surf, 0);
  gradEdge = AllocGradCost(surf->nvertices);
  gradBBR  = AllocGradCost(surf->nvertices);
  costEdge = MRISedgeCost(surf, gradEdge);
  costBBR = MRISbbrCost(bbrpar, gradBBR);
  cost = costEdge + costBBR;
  for(n=0; n < 1000; n++){
    SteepDir = DMatrixAddMul(gradEdge,gradBBR,1,1,SteepDir);
    UpdateVertexPosition(surf, 1000, SteepDir);
    MRISfaceNormalGrad(surf, 0);
    costEdge = MRISedgeCost(surf, gradEdge);
    costBBR = MRISbbrCost(bbrpar, gradBBR);
    cost = costEdge + costBBR;
    printf("%3d %12.10lf %12.10lf %12.10lf\n",n,costEdge,costBBR,cost);
    fflush(stdout);
    if(n== 20) MRISwrite(surf,"./lh.fsm010.pial.020");
    if(n==100) MRISwrite(surf,"./lh.fsm010.pial.100");
    if(n==200) MRISwrite(surf,"./lh.fsm010.pial.200");
    if(n==400) MRISwrite(surf,"./lh.fsm010.pial.400");
    if(n==700) MRISwrite(surf,"./lh.fsm010.pial.700");
  }
  MRISwrite(surf,"./lh.fsm010.pial");

  return(0);
}

