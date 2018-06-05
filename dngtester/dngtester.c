/**
 * @file  dngtester.c
 * @brief dougs super special test code
 *
 */
/*
 * Original Author: Doug Greve
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2015/04/26 16:15:55 $
 *    $Revision: 1.59 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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

#include "romp_support.h"

//int MRISfaceNormalFace_AddDeltaVertex = -1;
//long double MRISfaceNormalFace_AddDelta[3]={0,0,0};

 
/*!
  \fn int BBRPARsras2vox(BBRPARAMS *bbrpar)
  \brief Creates matrices that map the surface RAS to CRS
*/
int BBRPARsras2vox(BBRPARAMS *bbrpar)
{
  MATRIX *vox2sras, *sras2vox, *sras2vox3x3=NULL;
  vox2sras = MRIxfmCRS2XYZtkreg(bbrpar->mri);
  sras2vox = MatrixInverse(vox2sras,NULL);
  bbrpar->sras2vox = DMatrixCopyFMatrix(sras2vox,bbrpar->sras2vox);
  sras2vox3x3 = MatrixCopyRegion(sras2vox, sras2vox3x3, 1,1,3,3,1,1);
  bbrpar->sras2vox3x3 = DMatrixCopyFMatrix(sras2vox3x3,bbrpar->sras2vox3x3);
  MatrixFree(&vox2sras);
  MatrixFree(&sras2vox);
  MatrixFree(&sras2vox3x3);
  return(0);
}

double MRISfaceMinEdge(MRIS *surf, int faceno);
DMATRIX *MRISvertex2Matrix(MRIS *surf, int vtxno, DMATRIX *vm);
int MRISmatrix2Vertex(MRIS *surf, int vtxno, DMATRIX *vm);
int MRISscaleVertices(MRIS *surf, double scale);
int OldTestBBRCostFace(BBRFACE *bbrf, BBRPARAMS *bbrpar, int wrtvtxno);

double MRISedgeCost(MRIS *surf, DMATRIX **gradCost);
double MRISbbrCost(BBRPARAMS *bbrpar, DMATRIX **gradCost);
DMATRIX **AllocGradCost(int nvertices);
int FreeGradCost(int nvertices, DMATRIX ***pgradCost);
int UpdateVertexPosition(MRIS *surf, double stepsize, DMATRIX **grad);

/*----------------------------------------*/
int main(int argc, char **argv) 
{
  MRIS *surf, *surf2;
  int msec, nvertices; //vtxno=0;
  struct timeb  mytimer;
  MRI *mri, *mriindex, *mri2;
  Geodesics *geod;
  float maxdist;
  double d, dmax;
  int c,r,s;
  LABEL2SURF *l2s;
  LTA *lta;
  DMATRIX **dJ, **dpTprev, *gnum;
  VERTEX *v;
  double d1, delta;

  omp_set_num_threads(10);

  UnitizeNormalFace = 1;
  //surf = MRISread("/homes/4/greve/subjects/vrfp-mar14-anat/surf/lh.white");
  //surf = MRISread("/homes/4/greve/l/sp1/fsh.github.local/lib/bem/ic1.tri");
  //MRISscaleVertices(surf, 100);
  //MRISwrite(surf,"./lh.tmp0");
  //surf = MRISread(argv[1]);

  //printf("Reading volume\n"); fflush(stdout);
  mri = MRIread("/homes/4/greve/subjects/vrfp-mar14-anat/mri/orig.mgz");

  //printf("Reading surface\n"); fflush(stdout);
  surf = MRISread("/homes/4/greve/subjects/vrfp-mar14-anat/surf/lh.orig");
  MRIScomputeMetricProperties(surf);
  MRISfaceNormalGrad(surf, 0);
  MRISedges(surf);
  //MRISedgeGradDot(surf);

  BBRFACE *bbrf = BBRFaceAlloc();
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

  dJ = AllocGradCost(surf->nvertices);
  MRIScomputeMetricProperties(surf);
  MRISfaceNormalGrad(surf, 0);
  d = MRISbbrCost(bbrpar, dJ);
  gnum = DMatrixAlloc(1,3,MATRIX_REAL);
  delta = .0001;
  for(s=0; s < surf->nvertices; s++){
    v = &(surf->vertices[s]);
    for(c=0; c<3; c++){
      if(c==0) v->x += delta;
      if(c==1) v->y += delta;
      if(c==2) v->z += delta;
      //MRIScomputeMetricProperties(surf);
      MRISfaceNormalGrad(surf, 0);
      d1 = MRISbbrCost(bbrpar, NULL);
      gnum->rptr[1][c+1] = (d1-d)/delta;
      if(c==0) v->x -= delta;
      if(c==1) v->y -= delta;
      if(c==2) v->z -= delta;
    }
    printf("#@# %d ",s);
    DMatrixPrintFmt(stdout,"%12.10lf",DMatrixSubtract(dJ[s],gnum,NULL));
    printf("g0   = "); for(c=0; c<3; c++) printf("%12.8lf ",dJ[s]->rptr[1][c+1]);  printf("\n");
    printf("gnum = "); for(c=0; c<3; c++) printf("%12.8lf ",gnum->rptr[1][c+1]);  printf("\n");
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

#if 0
  MRISfaceNormalGrad(surf, 0);
  MRISedgeGradDot(surf);
  dmax = 0;
  s = 0;
  for(c = 0; c < surf->nedges; c++){
    for(r=0; r < 4; r++){
      d = MRISedgeGradDotEdgeVertexTest(surf, c, r, 0.00000001, 0);
      if(d > 0.1){
	printf("#@# %4d %5d %d %12.10lf\n",s,c,r,d);
	MRISedgeGradDotEdgeVertexTest(surf, c, r, 0.00000001, 1);
	s++;
      }
      if(dmax < d) dmax = d;
    }
  }
  printf("%5d %12.10lf\n",s,dmax);
  exit(0);
#endif

  dJ = (DMATRIX **) calloc(surf->nvertices,sizeof(DMATRIX *));
  for(c = 0; c < surf->nvertices; c++){
    dJ[c] = DMatrixAlloc(1,3,MATRIX_REAL);
  }
  dpTprev = (DMATRIX **) calloc(surf->nvertices,sizeof(DMATRIX *));
  for(c = 0; c < surf->nvertices; c++){
    dpTprev[c] = DMatrixAlloc(3,1,MATRIX_REAL);
  }

  long double J;
  DMATRIX *dJev=NULL, *V=NULL, *dp=NULL, *dpT=NULL;
  MRI_EDGE *e;
  int vtxno;
  dJev = DMatrixAlloc(1,3,MATRIX_REAL);
  for(s=0; s < 400000; s++){
    J = 0;
    MRISfaceNormalGrad(surf, 0);
    MRISedgeGradDot(surf);
    for(c = 0; c < surf->nedges; c++){
      for(r=0; r < 4; r++){
	//printf("c = %d  r = %d\n",c,r); fflush(stDout);
	e = &(surf->edges[c]);
	//J += MRISedgeAngleCostEdgeVertex(surf, c, r, &dJev);
	J += pow(1.0-e->dot,2.0);
	//printf("# %d %d %lf\n",c,r,e->gradDot[r]->rptr[1][1]);
	dJev = DMatrixScalarMul(e->gradDot[r],-2*(1.0-e->dot),dJev);
	vtxno = e->vtxno[r];
	dJ[vtxno] = DMatrixAdd(dJ[vtxno],dJev,dJ[vtxno]);
      }
    }
    printf("%3d %Lf\n",s,J);
    fflush(stdout);
    for(c = 0; c < surf->nvertices; c++){
      V = MRISvertex2Matrix(surf, c, V);
      dp = DMatrixScalarMul(dJ[c],-0.001,dp);
      dpT = DMatrixTranspose(dp,dpT);
      if(c>0) dpT = DMatrixAddMul(dpT,dpTprev[c],1,-0.9,dpT);
      V = DMatrixAdd(V,dpT,V);
      MRISmatrix2Vertex(surf, c, V);
      dpTprev[c] = DMatrixCopy(dpT,dpTprev[c]);
    }
  }
  MRISwrite(surf,"./lh.tmp");

  MRISfree(&surf);
  exit(0);

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
  TimerStart(&mytimer) ;
  geod = geodesicsRead(argv[2], &nvertices);
  msec = TimerStop(&mytimer) ;
  printf("t = %g min\n",msec/(1000.0*60));
  mri  = MRIread(argv[3]);
  mriindex  = MRIread(argv[4]);
  printf("Smoothing\n");
  TimerStart(&mytimer) ;
  mri2 = GeoSmooth(mri, 10, surf, geod, mriindex, NULL);
  msec = TimerStop(&mytimer) ;
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
  TimerStart(&mytimer) ;
  geod = computeGeodesics(surf, maxdist);
  msec = TimerStop(&mytimer) ;
  printf("done t = %g min\n",msec/(1000.0*60));
  geodesicsWrite(geod, surf->nvertices, argv[3]);
  msec = TimerStop(&mytimer) ;
  printf("done write t = %g min\n",msec/(1000.0*60));
  exit(0); //-------------------------


  //GeoDumpVertex("uvtx.108489.dist.dat", geod, 108489);
  TimerStart(&mytimer) ;
  geod = geodesicsReadV2(argv[1], &nvertices);
  msec = TimerStop(&mytimer) ;
  printf(" read2 t = %g min\n",msec/(1000.0*60));
  exit(0); //-------------------------

  TimerStart(&mytimer) ;
  geodesicsWrite(geod, nvertices, "u1.geod");
  msec = TimerStop(&mytimer) ;
  printf(" write1 t = %g min\n",msec/(1000.0*60));

  TimerStart(&mytimer) ;
  geodesicsWriteV2(geod, nvertices, "u2.geod");
  msec = TimerStop(&mytimer) ;
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
  surf->vertices[vtxno].x = vm->rptr[1][1];
  surf->vertices[vtxno].y = vm->rptr[2][1];
  surf->vertices[vtxno].z = vm->rptr[3][1];
  return(0);
}


int MRISscaleVertices(MRIS *surf, double scale)
{
  int vtxno;
  VERTEX *v;
  for(vtxno = 0; vtxno < surf->nvertices; vtxno++){
    v = &(surf->vertices[vtxno]);
    v->x *= scale;
    v->y *= scale;
    v->z *= (scale*(1+fabs(v->z)));
  }
  return(0);
}


double MRISbbrCost(BBRPARAMS *bbrpar, DMATRIX **gradCost)
{
  double cost = 0;
  int faceno,nthreads,threadno;
  BBRFACE **bbrfth;
  DMATRIX ***gradCostth=NULL;

  // This must already have been run
  //MRISfaceNormalGrad(surf, 0);
  // Not sure if the gradient part can be easily parallelized

  // Get number of threads
  nthreads = 1;
  #ifdef HAVE_OPENMP
  nthreads = omp_get_max_threads();  // using max should be ok
  #endif
  bbrfth = (BBRFACE **) calloc(sizeof(BBRFACE*),nthreads);
  for(threadno = 0; threadno < nthreads; threadno++){
    bbrfth[threadno] = BBRFaceAlloc();
  }
  if(gradCostth==NULL){
    gradCostth = (DMATRIX ***) calloc(sizeof(DMATRIX **),nthreads);
    for(threadno = 0; threadno < nthreads; threadno++){
      if(gradCost) gradCostth[threadno] = AllocGradCost(bbrpar->surf->nvertices);
    }
    // if static, then these have to be zeroed
  }

  cost = 0;
  #ifdef HAVE_OPENMP
  #pragma omp parallel for reduction(+ : cost)
  #endif
  for(faceno = 0; faceno < bbrpar->surf->nfaces; faceno++){
    BBRFACE *bbrf;

    int thno=0;
    #ifdef HAVE_OPENMP
    thno = omp_get_thread_num();
    #endif

    bbrf = bbrfth[thno];
    bbrf = BBRCostFace(faceno, -1, bbrpar, bbrf); // Compute Cost
    cost += bbrf->cost;
    if(gradCost){
      DMATRIX **thisGradCost;
      FACE *f;
      int  wrtvtxno, svtxno;
      thisGradCost = gradCostth[thno];
      f = &(bbrpar->surf->faces[faceno]);
      for(wrtvtxno=0; wrtvtxno < 3; wrtvtxno++){
	bbrf = BBRCostFace(faceno, wrtvtxno, bbrpar, bbrf); // Compute grads
	svtxno = f->v[wrtvtxno];
	thisGradCost[svtxno] = DMatrixAdd(thisGradCost[svtxno],bbrf->gradCost[wrtvtxno],thisGradCost[svtxno]);
      }
    }
  }

  // Merge threads
  if(gradCost){
    int  svtxno;
    //#ifdef HAVE_OPENMP
    //#pragma omp parallel for 
    //#endif
    for(svtxno=0; svtxno < bbrpar->surf->nvertices; svtxno++){
      DMATRIX *g;
      int threadno;
      g = gradCost[svtxno];
      DMatrixConstVal(0.0, 1, 3, g);
      // threadno here refers to threadno above
      for(threadno = 0; threadno < nthreads; threadno++){
	g = DMatrixAdd(g,gradCostth[threadno][svtxno],g);
      }
      // don't divde by number of threads
    }
  }

  for(threadno = 0; threadno < nthreads; threadno++){
    BBRFaceFree(&bbrfth[threadno]);
    if(gradCost) FreeGradCost(bbrpar->surf->nvertices, &gradCostth[threadno]);
  }

  free(bbrfth);
  free(gradCostth);

  //return(cost/bbrpar->surf->nvertices);
  return(cost);
}

/*-----------------------------------------------------------------------------*/
DMATRIX **AllocGradCost(int nvertices)
{
  DMATRIX **gradCost;
  int vtxno;

  gradCost = (DMATRIX **) calloc(nvertices,sizeof(DMATRIX *));
  for(vtxno = 0; vtxno < nvertices; vtxno++){
    gradCost[vtxno] = DMatrixAlloc(1,3,MATRIX_REAL);
  }
  return(gradCost);
}
/*-----------------------------------------------------------------------------*/
int FreeGradCost(int nvertices, DMATRIX ***pgradCost)
{
  int vtxno;
  DMATRIX **gradCost = *pgradCost;
  for(vtxno = 0; vtxno < nvertices; vtxno++){
    if(gradCost[vtxno]) DMatrixFree(&gradCost[vtxno]);
  }
  free(gradCost);
  return(0);
}

int UpdateVertexPosition(MRIS *surf, double stepsize, DMATRIX **grad)
{
  int vtxno;

  #ifdef HAVE_OPENMP
  #pragma omp parallel for 
  #endif
  for(vtxno = 0; vtxno < surf->nvertices; vtxno++){
    DMATRIX *V=NULL,*dp=NULL,*dpT=NULL;
    dp = DMatrixScalarMul(grad[vtxno],-stepsize,dp);
    dpT = DMatrixTranspose(dp,dpT);
    V = MRISvertex2Matrix(surf, vtxno, V);
    V = DMatrixAdd(V,dpT,V);
    MRISmatrix2Vertex(surf, vtxno, V);
    DMatrixFree(&V);
    DMatrixFree(&dp);
    DMatrixFree(&dpT);
  }
  return(0);
}



double MRISedgeCost(MRIS *surf, DMATRIX **gradCost)
{
  double cost = 0;
  int edgeno, wrtvtxno, svtxno;
  MRI_EDGE *e;
  DMATRIX *gradCostEV;

  // These must have been already run
  //MRISfaceNormalGrad(surf, 0);
  //MRISedgeGradDot(surf);

  for(edgeno = 0; edgeno < surf->nedges; edgeno++){
    e = &(surf->edges[edgeno]);
    // easy enough to compute actual cost here
    cost += pow(1.0-e->dot, 2.0);
    if(gradCost){
      for(wrtvtxno=0; wrtvtxno < 4; wrtvtxno++){
	gradCostEV = DMatrixScalarMul(e->gradDot[wrtvtxno],-2*(1.0-e->dot),gradCostEV);
	svtxno = e->vtxno[wrtvtxno];
	gradCost[svtxno] = DMatrixAdd(gradCost[svtxno],gradCostEV,gradCost[svtxno]);
      }
    }
  }
  return(cost);

}
