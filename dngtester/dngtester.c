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

#include "romp_support.h"

//int MRISfaceNormalFace_AddDeltaVertex = -1;
//long double MRISfaceNormalFace_AddDelta[3]={0,0,0};

 
typedef struct {
  double Din, Dout; // abs distance (mm) to project in and out
  double M, Q0; // BBR slope and offset
  MRI *mri;
  MRIS *surf;
  DMATRIX *sras2vox, *sras2vox3x3; // maps surface RAS to mri CRS 
} 
BBRPARAMS;

typedef struct {
  int faceno;
  DMATRIX *pc, *pin, *pout; // xyz of center, inward point and outward point
  DMATRIX *vin, *vout; // col,row,slice of inward point and outward point
  double Iin,Iout; // intensity of inward point and outward point
  double Q,cost; // Q=relative contrast, cost=bbr cost
  DMATRIX *gradPin[3], *gradVin[3], *gradIin[3], *gradPout[3], *gradVout[3], *gradIout[3];
  DMATRIX *gradInterpIn[3], *gradInterpOut[3];
  DMATRIX *gradQ[3], *gradCost[3];
} 
BBRFACE;


double MRISfaceMinEdge(MRIS *surf, int faceno);
DMATRIX *MRISvertex2Matrix(MRIS *surf, int vtxno, DMATRIX *vm);
int MRISmatrix2Vertex(MRIS *surf, int vtxno, DMATRIX *vm);
int MRISscaleVertices(MRIS *surf, double scale);
BBRFACE *BBRFaceAlloc(void);
int BBRFaceFree(BBRFACE **pbbrf);
BBRFACE *BBRCostFace(MRIS *surf, int faceno, int wrtvtxno, BBRPARAMS *bbrp, BBRFACE *bbrf);
int BBRFacePrint(FILE *fp, BBRFACE *bbrf, BBRPARAMS *bbrpar);
int TestBBRCostFace(BBRFACE *bbrf, BBRPARAMS *bbrpar, int wrtvtxno);

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
  long double delta;
  DMATRIX **dJ, **dpTprev;
  MATRIX *vox2sras=NULL, *sras2vox=NULL, *sras2vox3x3=NULL;

  UnitizeNormalFace = 1;
  //surf = MRISread("/homes/4/greve/subjects/vrfp-mar14-anat/surf/lh.white");
  //surf = MRISread("/homes/4/greve/l/sp1/fsh.github.local/lib/bem/ic1.tri");
  //MRISscaleVertices(surf, 100);
  //MRISwrite(surf,"./lh.tmp0");
  //surf = MRISread(argv[1]);

  printf("Reading surface\n"); fflush(stdout);
  surf = MRISread("/homes/4/greve/subjects/vrfp-mar14-anat/surf/lh.sphere");
  MRISfaceNormalGrad(surf, 0);
  MRISedgeGradDot(surf);
  MRISedges(surf);

  printf("Reading volume\n"); fflush(stdout);
  mri = MRIread("/homes/4/greve/subjects/vrfp-mar14-anat/mri/orig.mgz");

  BBRFACE *bbrf = BBRFaceAlloc();
  BBRPARAMS *bbrpar = (BBRPARAMS *) calloc(sizeof(BBRPARAMS),1);
  bbrpar->Din  = 0.5;
  bbrpar->Dout = 0.5;
  bbrpar->M = 0.01;
  bbrpar->Q0 = 0.0;
  bbrpar->mri = mri;
  bbrpar->surf = surf;

  vox2sras = MRIxfmCRS2XYZtkreg(mri);
  sras2vox = MatrixInverse(vox2sras,NULL);
  //sras2vox->rptr[1][4] += 0.0;
  bbrpar->sras2vox = DMatrixCopyFMatrix(sras2vox,bbrpar->sras2vox);
  sras2vox3x3 = MatrixCopyRegion(sras2vox, sras2vox3x3, 1,1,3,3,1,1);
  bbrpar->sras2vox3x3 = DMatrixCopyFMatrix(sras2vox3x3,bbrpar->sras2vox3x3);
  MatrixFree(&vox2sras);
  MatrixFree(&sras2vox);

  for(c=0; c < 3; c++){
    printf("------------------------------\n");
    bbrf = BBRCostFace(surf, 10, -1, bbrpar, bbrf); // Compute Cost
    bbrf = BBRCostFace(surf, 10, c, bbrpar, bbrf); // Compute grads
    BBRFacePrint(stdout, bbrf, bbrpar);
    TestBBRCostFace(bbrf, bbrpar, c); // wrtvtxno=0 (faceno=1)
  }
  BBRFaceFree(&bbrf);
  DMatrixFree(&bbrpar->sras2vox);
  DMatrixFree(&bbrpar->sras2vox3x3);
  MRIfree(&mri);
  MRISfree(&surf);
  exit(0);

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

  sscanf(argv[1],"%Lf",&delta);
  printf("delta = %Lf\n",delta);
  surf = MRISread(argv[2]);

  s = 0;
  for(c=1; c < surf->nfaces; c++){
    for(r=0; r < 3; r++){
      d = MRISfaceNormalGradFaceTest(surf, c, r, delta,0);
      if(dmax < d) dmax = d;
      if(d > .01){
	printf("%6d %d %g\n",c,r,d);
	MRISfaceNormalGradFaceTest(surf, c, r, delta,1);
	s = s + 1;
      }
    }
  }
  printf("dmax = %g, nhits = %d\n",dmax,s);
  exit(0);

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

BBRFACE *BBRFaceAlloc(void)
{
  BBRFACE *bbrf;
  int n;
  bbrf = (BBRFACE *) calloc(sizeof(BBRFACE),1);
  bbrf->pc = DMatrixAlloc(3,1,MATRIX_REAL);
  bbrf->pin = DMatrixAlloc(4,1,MATRIX_REAL);
  bbrf->vin = DMatrixAlloc(4,1,MATRIX_REAL);
  bbrf->pout = DMatrixAlloc(4,1,MATRIX_REAL);
  bbrf->vout = DMatrixAlloc(4,1,MATRIX_REAL);
  bbrf->pin->rptr[4][1] = 1.0;
  bbrf->vin->rptr[4][1] = 1.0;
  bbrf->pout->rptr[4][1] = 1.0;
  bbrf->vout->rptr[4][1] = 1.0;
  for(n=0; n < 3; n++){
    bbrf->gradPin[n] = DMatrixAlloc(3,3,MATRIX_REAL);
    bbrf->gradVin[n] = DMatrixAlloc(3,3,MATRIX_REAL);
    bbrf->gradInterpIn[n] = DMatrixAlloc(1,3,MATRIX_REAL);
    bbrf->gradIin[n] = DMatrixAlloc(1,3,MATRIX_REAL);
    bbrf->gradPout[n] = DMatrixAlloc(3,3,MATRIX_REAL);
    bbrf->gradVout[n] = DMatrixAlloc(3,3,MATRIX_REAL);
    bbrf->gradIout[n] = DMatrixAlloc(1,3,MATRIX_REAL);
    bbrf->gradInterpOut[n] = DMatrixAlloc(1,3,MATRIX_REAL);
    bbrf->gradQ[n] = DMatrixAlloc(1,3,MATRIX_REAL);
    bbrf->gradCost[n] = DMatrixAlloc(1,3,MATRIX_REAL);
  }
  return(bbrf);
}
int BBRFaceFree(BBRFACE **pbbrf)
{
  BBRFACE *bbrf = *pbbrf;
  int n;
  if(bbrf->pc) DMatrixFree(&bbrf->pc);
  if(bbrf->pin) DMatrixFree(&bbrf->pin);
  if(bbrf->vin) DMatrixFree(&bbrf->vin);
  if(bbrf->pout) DMatrixFree(&bbrf->pout);
  if(bbrf->vout) DMatrixFree(&bbrf->pout);
  for(n=0; n < 3; n++){
    DMatrixFree(&bbrf->gradPin[n]);
    DMatrixFree(&bbrf->gradVin[n]);
    DMatrixFree(&bbrf->gradInterpIn[n]);
    DMatrixFree(&bbrf->gradIin[n]);
    DMatrixFree(&bbrf->gradPout[n]);
    DMatrixFree(&bbrf->gradVout[n]);
    DMatrixFree(&bbrf->gradInterpOut[n]);
    DMatrixFree(&bbrf->gradIout[n]);
    DMatrixFree(&bbrf->gradQ[n]);
    DMatrixFree(&bbrf->gradCost[n]);
  }
  // don't free ras2vox or MRI
  free(bbrf);
  *pbbrf = NULL;
  return(0);
}

int BBRFacePrint(FILE *fp, BBRFACE *bbrf, BBRPARAMS *bbrpar)
{
  FACE *f;
  int n;
  f = &(bbrpar->surf->faces[bbrf->faceno]);
  fprintf(fp,"faceno %d\n",bbrf->faceno);
  fprintf(fp,"area %12.10f\n",f->area);
  fprintf(fp,"pin ");   for(n=0; n<3; n++)   fprintf(fp,"%6.4lf ",bbrf->pin->rptr[n+1][1]);
  fprintf(fp,"\n");
  fprintf(fp,"vin ");   for(n=0; n<3; n++)   fprintf(fp,"%6.4lf ",bbrf->vin->rptr[n+1][1]);
  fprintf(fp,"\n");
  fprintf(fp,"pout ");   for(n=0; n<3; n++)   fprintf(fp,"%6.4lf ",bbrf->pout->rptr[n+1][1]);
  fprintf(fp,"\n");
  fprintf(fp,"vout ");   for(n=0; n<3; n++)   fprintf(fp,"%6.4lf ",bbrf->vout->rptr[n+1][1]);
  fprintf(fp,"\n");
  fprintf(fp,"Iin  %lf\n",bbrf->Iin);
  fprintf(fp,"Iout %lf\n",bbrf->Iout);
  fprintf(fp,"Q %lf\n",bbrf->Q);
  fprintf(fp,"cost %12.10lf\n",bbrf->cost);
  fflush(fp);
  return(0);
}


BBRFACE *BBRCostFace(MRIS *surf, int faceno, int wrtvtxno, BBRPARAMS *bbrpar, BBRFACE *bbrf)
{
  FACE *f;
  VERTEX *v;
  int n;
  double c,r,s;
  double OutInSum, OutInDiff, vtanh;
  DMATRIX *tmp1, *tmp2;

  if(bbrf==NULL) {
    printf("ERROR: BBRCostFace(): bbrf is NULL\n");
    return(NULL);
  }

  bbrf->faceno = faceno;
  f = &(surf->faces[faceno]);

  if(wrtvtxno < 0){
    // Compute the center of the triangle
    bbrf->pc = DMatrixConstVal(0, 4, 1, bbrf->pc);
    for(n=0; n<3; n++){
      v = &(surf->vertices[f->v[n]]);
      bbrf->pc->rptr[1][1] += (v->x);
      bbrf->pc->rptr[2][1] += (v->y);
      bbrf->pc->rptr[3][1] += (v->z);
      // Add to vertex here if testing
    }
    DMatrixScalarMul(bbrf->pc,1/3.0,bbrf->pc);
    
    // Compute the XYZ at the interior point = pc - norm*Din
    // Note: pc and norm are 3x1 whereas pin is 4x1
    for(n=0; n<3; n++)
      bbrf->pin->rptr[n+1][1] = bbrf->pc->rptr[n+1][1] - (bbrpar->Din*f->norm->rptr[n+1][1]);
    
    // Compute the CRS at the interior point = ras2vox*pin
    bbrf->vin = DMatrixMultiply(bbrpar->sras2vox,bbrf->pin,bbrf->vin);
    
    // Compute the intensity at the interior point
    c=bbrf->vin->rptr[1][1]; r=bbrf->vin->rptr[1][2]; s=bbrf->vin->rptr[1][3];
    MRIsampleVolume(bbrpar->mri, c,r,s, &bbrf->Iin);
    
    // Compute the XYZ at the exterior point = pc + norm*Dout
    // Note: pc and norm are 3x1 whereas pin is 4x1
    for(n=0; n<3; n++)
      bbrf->pout->rptr[n+1][1] = bbrf->pc->rptr[n+1][1] + (bbrpar->Dout*f->norm->rptr[n+1][1]);
    
    // Compute the CRS at the exterior point = ras2vox*pout
    bbrf->vout = DMatrixMultiply(bbrpar->sras2vox,bbrf->pout,bbrf->vout);
    
    // Compute the intensity at the exterior point
    c=bbrf->vout->rptr[1][1]; r=bbrf->vout->rptr[1][2]; s=bbrf->vout->rptr[1][3];
    MRIsampleVolume(bbrpar->mri, c,r,s, &bbrf->Iout);

    // Finally, compute the cost
    OutInSum  = bbrf->Iout + bbrf->Iin;
    OutInDiff = bbrf->Iout - bbrf->Iin;
    bbrf->Q = 100*(OutInDiff)/(0.5*(OutInSum));
    vtanh = tanh(bbrpar->M*(bbrf->Q - bbrpar->Q0));
    bbrf->cost = 1 + vtanh;

    return(bbrf);
  } // if(wrtvtxno < 0)

  // Now compute the gradients

  // Recompute these
  OutInSum  = bbrf->Iout + bbrf->Iin;
  OutInDiff = bbrf->Iout - bbrf->Iin;
  vtanh = tanh(bbrpar->M*(bbrf->Q - bbrpar->Q0));

  // gradPin = gradPc - Din*gradNorm, gradPc = eye(3)/3
  bbrf->gradPin[wrtvtxno] = DMatrixScalarMul(f->gradNorm[wrtvtxno],-bbrpar->Din,bbrf->gradPin[wrtvtxno]);
  for(n=0; n<3; n++) bbrf->gradPin[wrtvtxno]->rptr[n+1][n+1] += 1/3.0;

  // gradVin = ras2vox*gradPin
  bbrf->gradVin[wrtvtxno] = DMatrixMultiply(bbrpar->sras2vox3x3,bbrf->gradPin[wrtvtxno],bbrf->gradVin[wrtvtxno]);

  // Compute gradient of the sampled intensity WRT a change in col, row, or slice
  c=bbrf->vin->rptr[1][1]; r=bbrf->vin->rptr[1][2]; s=bbrf->vin->rptr[1][3];
  bbrf->gradInterpIn[wrtvtxno] = MRIgradTrilinInterp(bbrpar->mri, c,r,s, bbrf->gradInterpIn[wrtvtxno]);

  // Compute gradient of the sampled intensity WRT a change vertex position
  bbrf->gradIin[wrtvtxno] = 
    DMatrixMultiply(bbrf->gradInterpIn[wrtvtxno],bbrf->gradVin[wrtvtxno],bbrf->gradIin[wrtvtxno]);


  // gradPout = gradPc + Dout*gradNorm, gradPc = eye(3)/3
  bbrf->gradPout[wrtvtxno] = DMatrixScalarMul(f->gradNorm[wrtvtxno],+bbrpar->Dout,bbrf->gradPout[wrtvtxno]);
  for(n=0; n<3; n++) bbrf->gradPout[wrtvtxno]->rptr[n+1][n+1] += 1/3.0;

  // gradVout = ras2vox*gradPout
  bbrf->gradVout[wrtvtxno] = DMatrixMultiply(bbrpar->sras2vox3x3,bbrf->gradPout[wrtvtxno],bbrf->gradVout[wrtvtxno]);

  // Compute gradient of the sampled intensity WRT a change in col, row, or slice
  c=bbrf->vout->rptr[1][1]; r=bbrf->vout->rptr[1][2]; s=bbrf->vout->rptr[1][3];
  bbrf->gradInterpOut[wrtvtxno] = MRIgradTrilinInterp(bbrpar->mri, c,r,s, bbrf->gradInterpOut[wrtvtxno]);

  // Compute gradient of the sampled intensity WRT a change vertex position
  bbrf->gradIout[wrtvtxno] = 
    DMatrixMultiply(bbrf->gradInterpOut[wrtvtxno],bbrf->gradVout[wrtvtxno],bbrf->gradIout[wrtvtxno]);

  // Compute the gradient of Q = 200*((gradOut-gradIn)/(Iout+Iin) - (gradOut+gradIn)/(Iout-Iin)/((Iout-Iin)^2));
  tmp1 = DMatrixSubtract(bbrf->gradIout[wrtvtxno],bbrf->gradIin[wrtvtxno],NULL);
  tmp1 = DMatrixScalarMul(tmp1,1.0/OutInSum,tmp1);
  tmp2 = DMatrixAdd(bbrf->gradIout[wrtvtxno],bbrf->gradIin[wrtvtxno],NULL);
  tmp2 = DMatrixScalarMul(tmp2,-OutInDiff/(OutInSum*OutInSum),tmp2);
  bbrf->gradQ[wrtvtxno] = DMatrixAdd(tmp1,tmp2, bbrf->gradQ[wrtvtxno]);
  bbrf->gradQ[wrtvtxno] = DMatrixScalarMul(bbrf->gradQ[wrtvtxno],200,bbrf->gradQ[wrtvtxno]);
  DMatrixFree(&tmp1);
  DMatrixFree(&tmp2);

  // Finally, compute the gadient of the cost = M*(1-tanh^2)*gradQ
  bbrf->gradCost[wrtvtxno] = DMatrixScalarMul(bbrf->gradQ[wrtvtxno],bbrpar->M*(1.0 - vtanh*vtanh),bbrf->gradCost[wrtvtxno]);

  return(bbrf);
}

int TestBBRCostFace(BBRFACE *bbrf,   BBRPARAMS *bbrpar, int wrtvtxno)
{
  int wrtdimno,svtxno;
  FACE *f;
  VERTEX *v;
  DMATRIX *Pin0=NULL, *Pin1=NULL, *g0, *gnum, *d=NULL;
  double delta = 0.001, Iin0, Iin1;

  Pin0 = DMatrixCopy(bbrf->vin,NULL);
  //printf("P0 ---------------\n");
  //DMatrixPrintFmt(stdout,"%12.10lf",Pin0);
  Iin0 = bbrf->cost;
  printf("Iin0 = %g\n",Iin0);

  g0 = bbrf->gradCost[wrtvtxno];

  f = &(bbrpar->surf->faces[bbrf->faceno]);
  printf("fnorm0 ---------------\n");
  DMatrixPrintFmt(stdout,"%12.10lf",f->norm);
  svtxno = f->v[wrtvtxno];
  v = &(bbrpar->surf->vertices[svtxno]);
  gnum = DMatrixAlloc(1,3,MATRIX_REAL);
  for(wrtdimno=0; wrtdimno<3; wrtdimno++){
    printf("wrtdim %d ============================\n",wrtdimno);
    if(wrtdimno==0) v->x += delta;
    if(wrtdimno==1) v->y += delta;
    if(wrtdimno==2) v->z += delta;
    MRISfaceNormalGrad(bbrpar->surf, 1); // dont compute grad (norm only)
    //printf("fnorm%d ---------------\n",wrtdimno);
    //DMatrixPrintFmt(stdout,"%12.10lf",f->norm);
    bbrf = BBRCostFace(bbrpar->surf, bbrf->faceno, -1, bbrpar, bbrf); // dont compute grad
    Pin1 = DMatrixCopy(bbrf->vin,NULL);
    Iin1 = bbrf->cost;
    printf("Iin1 = %g\n",Iin0);
    //printf("P1 ---------------\n");
    //DMatrixPrintFmt(stdout,"%12.10lf",Pin1);
    d = DMatrixSubtract(Pin1,Pin0,d);
    d = DMatrixScalarMul(d,1.0/delta,d);
    gnum->rptr[1][wrtdimno+1] = (Iin1-Iin0)/delta;
    if(wrtdimno==0) v->x -= delta;
    if(wrtdimno==1) v->y -= delta;
    if(wrtdimno==2) v->z -= delta;
  }
  printf("============================\n");
  printf("g0 ---------------\n");
  DMatrixPrintFmt(stdout,"%12.10lf",g0);
  printf("gnum ---------------\n");
  DMatrixPrintFmt(stdout,"%12.10lf",gnum);

  DMatrixFree(&Pin0);
  DMatrixFree(&Pin1);
  DMatrixFree(&gnum);
  DMatrixFree(&d);

  return(0);
}


