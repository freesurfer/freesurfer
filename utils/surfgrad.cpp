/**
 * @brief Utilities to compute gradients on the surface
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

#include <ctype.h>
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "diag.h"
#include "error.h"
#include "matrix.h"
#include "dmatrix.h"
#include "mrisurf.h"
#include "romp_support.h"
#undef private

#define _SURFGRAD_SRC
#include "surfgrad.h"


/*!
  \fn int MRISfaceMetric(MRIS *surf, int DoGrad)
  \brief Computes face metrics for all faces. The only
  metric it computes are the face normals. 
 */
int MRISfaceMetric(MRIS *surf, int DoGrad)
{
  int faceno;

  ROMP_PF_begin
  #ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(assume_reproducible) 
  #endif
  for(faceno=0; faceno < surf->nfaces; faceno++){
    ROMP_PFLB_begin
    if(DoGrad)
      MRISfaceNormalGradFace(surf, faceno);
    else
      MRISfaceNormalFace(surf, faceno, NULL, NULL);
    ROMP_PFLB_end
  }
  ROMP_PF_end
  return(0);
}

/*!
  \fn int MRISedgeMetric(MRIS *surf)
  \brief Computes the edge metric for all edges.
 */
int MRISedgeMetric(MRIS *surf, int DoGrad)
{
  int edgeno;
  ROMP_PF_begin
  #ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(assume_reproducible) 
  #endif
  for(edgeno = 0; edgeno < surf->nedges; edgeno++){
    ROMP_PFLB_begin
    MRISedgeMetricEdge(surf, edgeno, DoGrad);
    ROMP_PFLB_end
  }
  ROMP_PF_end
  return(0);
}

/*!
  \fn double MRIScornerMetric(MRIS *surf, const int DoGrad)
  \brief Computes the metrics (dot, angle) of all corners. If
  requested, the gradient of the dot is also computed WRT each vertex
  in all corners.
*/
int MRIScornerMetric(MRIS *surf, const int DoGrad)
{
  int cornerno;
  if(surf->corners == NULL) MRIScorners(surf);

  ROMP_PF_begin
  #ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(assume_reproducible) 
  #endif
  for(cornerno = 0; cornerno < surf->ncorners; cornerno++){
    ROMP_PFLB_begin
    MRIScornerMetricCorner(surf, cornerno, DoGrad);
    ROMP_PFLB_end
  }
  ROMP_PF_end

  return(0);
}



/*!
  \fn int MRISfaceNormalFace(MRIS *surf, int faceno, DMATRIX **pc, double *pcL)
  \brief Compute the normal to the given face. If pc and pcL are not NULL,
  then the un-normalized cross product (c) and lenght of c (cL) are returned.
  There are two global variables: MRISfaceNormalFace_AddDeltaVertex and
  MRISfaceNormalFace_AddDelta[3] which are used for testing. If AddDeltaVertex
  is less than 0, then these have no effect.
 */
int MRISfaceNormalFace(MRIS *surf, int faceno, DMATRIX **pc, double *pcL)
{
  extern int MRISfaceNormalFace_AddDeltaVertex;
  extern long double MRISfaceNormalFace_AddDelta[3];
  VERTEX *v0, *v1, *v2;
  FACE *f;
  long double ax, ay, az, bx, by, bz, cx, cy, cz, cL;

  if(faceno < 0 || faceno >= surf->nfaces){
    printf("MRISfaceNormalFace(): faceno = %d (0,%d)\n",faceno,surf->nfaces);
    return(1);
  }

  f = &(surf->faces[faceno]);

  v0 = &(surf->vertices[f->v[0]]);
  v1 = &(surf->vertices[f->v[1]]);
  v2 = &(surf->vertices[f->v[2]]);

  // Vector a points from v0 to v1
  ax = (long double)v1->x - v0->x;
  ay = (long double)v1->y - v0->y;
  az = (long double)v1->z - v0->z;

  // Vector b points from v1 to v2
  bx = (long double)v2->x - v1->x;
  by = (long double)v2->y - v1->y;
  bz = (long double)v2->z - v1->z;

  if(MRISfaceNormalFace_AddDeltaVertex > -1){
    /* Add the given amount to the given vertex. Good for testing
       theoretical derivative computation.  Set AddDeltaVertex to the
       vertex with respect to the numertical derivative will be
       computed (0, 1, 2); also set AddDelta[3] to the appropriate values.
       It is done in this way to preserve the precision. If a value is just
       added to the given xyz in the vertex structure, very small values
       needed to compute accurate derivatives are lost. */
    //printf("Adding to vertex %d %Lf %Lf %Lf\n",MRISfaceNormalFace_AddDeltaVertex,
    //	   MRISfaceNormalFace_AddDelta[0],MRISfaceNormalFace_AddDelta[1],MRISfaceNormalFace_AddDelta[2]);
    switch(MRISfaceNormalFace_AddDeltaVertex){
    case 0:
      // Add delta onto vertex 0 (use neg because a = v1-v0)
      ax -= MRISfaceNormalFace_AddDelta[0];
      ay -= MRISfaceNormalFace_AddDelta[1];
      az -= MRISfaceNormalFace_AddDelta[2];
      break;
    case 1:
      // Add delta onto vertex 1 (add to a, subtract from b)
      ax += MRISfaceNormalFace_AddDelta[0];
      ay += MRISfaceNormalFace_AddDelta[1];
      az += MRISfaceNormalFace_AddDelta[2];
      bx -= MRISfaceNormalFace_AddDelta[0];
      by -= MRISfaceNormalFace_AddDelta[1];
      bz -= MRISfaceNormalFace_AddDelta[2];
      break;
    case 2:
      // Add delta onto vertex 2 (add to b)
      bx += MRISfaceNormalFace_AddDelta[0];
      by += MRISfaceNormalFace_AddDelta[1];
      bz += MRISfaceNormalFace_AddDelta[2];
      break;
    }
  }

  // Compute cross product
  cx = ay*bz - az*by;
  cy = az*bx - ax*bz;
  cz = ax*by - ay*bx;
  cL = sqrt(cx*cx + cy*cy + cz*cz);
  if(pcL != NULL) *pcL = cL;

  // Store c if needed
  if(pc != NULL){
    if(*pc==NULL) {
      *pc = DMatrixAlloc(3,1,MATRIX_REAL);
    }
    (*pc)->rptr[1][1] = cx;
    (*pc)->rptr[2][1] = cy;
    (*pc)->rptr[3][1] = cz;
  }

  // Store in the face norm matrix
  if(f->norm == NULL)
    f->norm = DMatrixAlloc(3,1,MATRIX_REAL);
  f->norm->rptr[1][1] = cx/cL;
  f->norm->rptr[2][1] = cy/cL;
  f->norm->rptr[3][1] = cz/cL;

  //printf("a %Lf %Lf %Lf\n",ax,ay,az);
  //printf("b %Lf %Lf %Lf\n",bx,by,bz);
  //printf("c %Lf %Lf %Lf   %Lf\n",cx,cy,cz,cL);
  //DMatrixPrint(stdout,f->norm);

  return(0);
}

/*!
  \fn int MRISfaceNormalGradFace(MRIS *surf, int faceno)

  \brief Computes the gradient of the face normal with respect to changes in
the the position of each vertex.  The matrices are stored in
face->gradNorm[refvtxno] (also computes the norm itself).
The gradient is a 3x3 matrix of the form
   dnx/dvx dnx/dvy dnx/dvz  
   dny/dvx dny/dvy dny/dvz  
   dnz/dvx dnz/dvy dnz/dvz  
Where nx is the x component of the normal and vx is the x component of the vertex.
Since there are three vertices, there are three such matrices.
The norm is stored in face->norm
These equations themselves have been verified emperically, though there could be
some issues with accuracy in certain cases because the data structures use float.
See MRISfaceNormalGradTest().
*/
int MRISfaceNormalGradFace(MRIS *surf, int faceno)
{
  VERTEX *vSrc, *vWRT, *vC;
  FACE *f;
  DMATRIX *c=NULL, *cT=NULL, *gradc=NULL, *gradcL=NULL, *tmp1=NULL, *tmp2=NULL;
  long double ax, ay, az, bx, by, bz;
  double cL;
  int wrtvtxno,k,m;

  if(faceno < 0 || faceno >= surf->nfaces){
    printf("MRISfaceNormalGrad(): faceno = %d (0,%d)\n",faceno,surf->nfaces);
    return(1);
  }

  MRISfaceNormalFace(surf, faceno, &c, &cL);
  cT = DMatrixTranspose(c,cT);

  f = &(surf->faces[faceno]);

  // Compute gradient with respect to (wrt) each vertex in the triangle
  for(wrtvtxno = 0; wrtvtxno < 3; wrtvtxno++){

    // This is the vertex we are computing the gradient for
    vWRT = &(surf->vertices[f->v[wrtvtxno]]);

    // This is the "source" vertex. Vectors used in the cross product will 
    // point away from this vertex
    k = wrtvtxno-1;
    if(k<0) k = 2;
    vSrc = &(surf->vertices[f->v[k]]);

    // This is just the 3rd vertex
    m = wrtvtxno+1;
    if(m>2) m = 0;
    vC = &(surf->vertices[f->v[m]]);

    // Vector a must go from source vertex to WRT vertex
    ax = (long double)vWRT->x - vSrc->x;
    ay = (long double)vWRT->y - vSrc->y;
    az = (long double)vWRT->z - vSrc->z;

    // Vector b must go from source to the other vertex (and
    // specifically not include the WRT vertex)
    bx = (long double)vC->x - vSrc->x;
    by = (long double)vC->y - vSrc->y;
    bz = (long double)vC->z - vSrc->z;

    // Gradient of the unnormalized cross 
    if(gradc == NULL) gradc = DMatrixAlloc(3,3,MATRIX_REAL);
    gradc->rptr[1][2] =  bz;
    gradc->rptr[1][3] = -by;
    gradc->rptr[2][1] = -bz;
    gradc->rptr[2][3] =  bx;
    gradc->rptr[3][1] =  by;
    gradc->rptr[3][2] = -bx;

    // Gradient of the length of c = c'*gradc/cL
    gradcL = DMatrixMultiply(cT,gradc,gradcL);
    gradcL = DMatrixScalarMul(gradcL, 1.0/cL, gradcL); // scale
    
    // Gradient of normalized cross
    // gradNorm = gradc/cL - [c*gradcL]/(cL^2) = tmp1-tmp2
    tmp1 = DMatrixScalarMul(gradc, 1.0/cL, tmp1); // scale
    tmp2 = DMatrixMultiply(c,gradcL,tmp2);
    tmp2 = DMatrixScalarMul(tmp2, 1.0/(cL*cL), tmp2); // scale
    f->gradNorm[wrtvtxno] = DMatrixSubtract(tmp1,tmp2,f->gradNorm[wrtvtxno]); // subtract
    
    if(0){
      printf("a = [%Lf %Lf %Lf]';",ax,ay,az);
      printf("b = [%Lf %Lf %Lf]';",bx,by,bz);
      printf("c = [");
      DMatrixPrint(stdout,c);
      printf("];\n");
      printf("gradc = [");
      DMatrixPrint(stdout,gradc);
      printf("];\n");
      printf("gradcL = [");
      DMatrixPrint(stdout,gradcL);
      printf("];\n");
      printf("grad = [");
      DMatrixPrint(stdout,f->gradNorm[wrtvtxno]);
      printf("];\n");
      printf("[g0 gradc0 gradcL0]= gradnormcross(a,b);\n");
      printf("mar(g0-grad)\n");
    }
  }

  DMatrixFree(&c); 
  DMatrixFree(&gradc); 
  DMatrixFree(&cT); 
  DMatrixFree(&gradcL); 
  DMatrixFree(&tmp1); 
  DMatrixFree(&tmp2);

  return(0);
}


double MRISfaceNormalGradFaceTest(MRIS *surf, int faceno, int wrtvtxno, long double delta, int verbose)
{
  FACE *f;
  VERTEX *v;
  DMATRIX *n0=NULL, *g0, *n1=NULL, *g1=NULL, *gnum;
  int wrtdimno, c;
  double diff, maxdiff, maxrdiff, gmax;
  extern int MRISfaceNormalFace_AddDeltaVertex;
  extern long double MRISfaceNormalFace_AddDelta[3];

  f = &(surf->faces[faceno]);

  // Compute the normal at baseline position
  MRISfaceNormalFace_AddDeltaVertex = -1;
  MRISfaceNormalFace(surf, faceno, NULL, NULL);
  n0 = DMatrixCopy(f->norm,n0);
  if(verbose){
    printf("n0 = [");
    DMatrixPrint(stdout,n1);
    printf("];\n");
  }

  // Compute the theoretical gradient, does all 3 wrt vertices
  MRISfaceNormalGradFace(surf, faceno);
  g0 = f->gradNorm[wrtvtxno];
  gmax = DMatrixMaxAbs(g0);

  // Compute the numerical gradient
  MRISfaceNormalFace_AddDeltaVertex = wrtvtxno;
  gnum = DMatrixAlloc(3,3,MATRIX_REAL);
  maxdiff = 0;
  // Go through each dim in the WRT vertex position
  for(wrtdimno=0; wrtdimno < 3; wrtdimno++){
    // Set the delta for this dim
    for(c=0; c < 3; c++) MRISfaceNormalFace_AddDelta[c] = 0;
    MRISfaceNormalFace_AddDelta[wrtdimno] = delta;
    // Recompute the normal
    MRISfaceNormalFace(surf, faceno, NULL, NULL);
    n1 = DMatrixCopy(surf->faces[faceno].norm,n1);
    if(verbose){
      printf("n1c%d = [",wrtdimno);
      DMatrixPrint(stdout,n1);
      printf("];\n");
    }
    // numerical grad = diff in the normal div by delta
    g1 = DMatrixSubtract(n1,n0,g1);
    g1 = DMatrixScalarMul(g1,1.0/delta,g1);
    // Fill in the matrix and compute difference with theoretical
    for(c=0; c < 3; c++){
      gnum->rptr[c+1][wrtdimno+1] = g1->rptr[c+1][1];
      diff = fabs(gnum->rptr[c+1][wrtdimno+1] - g0->rptr[c+1][wrtdimno+1]);
      if(maxdiff < diff) maxdiff = diff;
    }
  }

  // The final metric is the maxabsdiff divided by the max component in the 
  // theoretical to make sure scaling is correct
  maxrdiff = maxdiff/gmax;

  // Undo everything
  MRISfaceNormalFace_AddDeltaVertex = -1;
  for(c=0; c < 3; c++) MRISfaceNormalFace_AddDelta[c] = 0;
  MRISfaceNormalFace(surf, wrtvtxno, NULL, NULL);

  if(verbose){
    printf("area = %g\n",f->area);
    for(c=0; c < 3; c++){
      v = &(surf->vertices[f->v[0]]);
      printf("vtx%d = [%g %g %g]';\n",c,v->x,v->y,v->z);
    }
    printf("gmax = %lf\n",gmax);
    printf("g0 = [");
    DMatrixPrintFmt(stdout,"%12.10f",g0);
    printf("];\n");
    printf("gnum = [");
    DMatrixPrintFmt(stdout,"%12.10f",gnum);
    printf("];\n");
  }

  DMatrixFree(&n0);
  DMatrixFree(&n1);
  DMatrixFree(&g1);
  DMatrixFree(&gnum);

  return(maxrdiff);
}


/*!
  \fn int MRISedgeMetricEdge(MRIS *surf, int edgeno, int DoGrad)
  \brief Computes edge-related metrics including the length, the dot
  product of the angle of the faces across the edge, and the angle
  (def).
 */
int MRISedgeMetricEdge(MRIS *surf, int edgeno, int DoGrad)
{
  MRI_EDGE *e;
  FACE *f0, *f1;
  VERTEX *v0, *v1;
  double v0xyz[3],v1xyz[3];
  DMATRIX *gradU=NULL;
  e = &(surf->edges[edgeno]);
  v0 = &(surf->vertices[e->vtxno[0]]);
  v1 = &(surf->vertices[e->vtxno[1]]);
  StuffVertexCoords(surf, e->vtxno[0], &v0xyz[0]);
  StuffVertexCoords(surf, e->vtxno[1], &v1xyz[0]);
  if(DoGrad == 1 || DoGrad == 3){
    if(e->gradU == NULL) e->gradU = DMatrixAlloc(3,3,MATRIX_REAL);
    gradU = e->gradU;
  }
  e->len = SurfGradUnitVector(v0xyz, v1xyz, e->u, gradU);
  if(DoGrad == 2 || DoGrad == 3){
    int wrtvtx; // compute d->gradDot[4]
    for(wrtvtx = 0; wrtvtx < 4; wrtvtx ++){
      MRISedgeGradDotEdgeVertex(surf, edgeno, wrtvtx);
    }
  }
  // Note f->norm must have been computed with MRISfaceNormalGrad(), grad not needed
  f0 = &(surf->faces[e->faceno[0]]);
  f1 = &(surf->faces[e->faceno[1]]);
  e->area[0] = f0->area;
  e->area[1] = f1->area;
  if(e->area[0] > e->area[1]) e->maxarea = e->area[0];
  else                        e->maxarea = e->area[1];
  e->dot = DVectorDot(f0->norm,f1->norm);
  if(e->dot > 1.0)  e->dot = +1.0; // might be a slight overflow
  if(e->dot < -1.0) e->dot = -1.0; // might be a slight overflow
  e->angle = acos(e->dot)*180/M_PI;
  return(0);
}
/*!
\fn int MRISedgeCompare(const void *a, const void *b, void *psorttype)
\brief Compare two edges in a way suitable for qsort_r(). psorttype is an int with:
  case 1: // Length then Max Area, Descending (high to low)
  case 2: // Length then Max Area, Ascending  (low to high)
  case 3: // Max Area then Length, Descending (high to low)
  case 4: // Max Area then Length, Ascending  (low to high)
 */
#ifdef QSORT_R_THUNK_FIRST
int MRISedgeCompare(void *psorttype, const void *a, const void *b)
#else
int MRISedgeCompare(const void *a, const void *b, void *psorttype)
#endif
{
  MRI_EDGE *e1 = (MRI_EDGE*)a;
  MRI_EDGE *e2 = (MRI_EDGE*)b;
  int *SortType = (int*)psorttype;
  switch(*SortType){
  case 1: // Length then Max Area, Descending
    if(e1->len > e2->len) return(-1);
    if(e1->len < e2->len) return(+1);
    if(e1->maxarea > e2->maxarea) return(-1);
    if(e1->maxarea < e2->maxarea) return(+1);
    break;
  case 2: // Length then Max Area, Ascending
    if(e1->len < e2->len) return(-1);
    if(e1->len > e2->len) return(+1);
    if(e1->maxarea < e2->maxarea) return(-1);
    if(e1->maxarea > e2->maxarea) return(+1);
    break;
  case 3: // Max Area then Length, Descending
    if(e1->maxarea > e2->maxarea) return(-1);
    if(e1->maxarea < e2->maxarea) return(+1);
    if(e1->len > e2->len) return(-1);
    if(e1->len < e2->len) return(+1);
    break;
  case 4: // Max Area then Length, Ascending
    if(e1->maxarea < e2->maxarea) return(-1);
    if(e1->maxarea > e2->maxarea) return(+1);
    if(e1->len < e2->len) return(-1);
    if(e1->len > e2->len) return(+1);
    break;
  default:
    printf("ERROR: MRISEdgeCompare() SorType = %d\n",*SortType);
    fflush(stdout);
    return(0);
    break;
  }
  return(0); // should never get here
}

/*!
\fn MRI_EDGE *MRISedgeSort(MRIS *surf, const int SortType, MRI_EDGE *edges)
\brief Sort the edges in surf->edges (returns a new array). Could be done
in place, but probably not a good idea.
SortType is
  case 1: // Length then Max Area, Descending (high to low)
  case 2: // Length then Max Area, Ascending  (low to high)
  case 3: // Max Area then Length, Descending (high to low)
  case 4: // Max Area then Length, Ascending  (low to high)
 */
MRI_EDGE *MRISedgeSort(MRIS *surf, const int SortType, MRI_EDGE *edges)
{
  // Make sure the edges are built and metrics computed
  MRISedges(surf);
  MRISfaceMetric(surf,0);
  MRISedgeMetric(surf,0);
  if(edges == NULL) edges = (MRI_EDGE*)calloc(surf->nedges,sizeof(MRI_EDGE));
  memcpy(edges,surf->edges,surf->nedges*sizeof(MRI_EDGE));
  portable_qsort_r(edges, surf->nedges, sizeof(MRI_EDGE), MRISedgeCompare, (void*)&SortType);
  return(edges);
}

/*!
  \fn double MRISedgeAngleCostEdgeVertex(MRIS *surf, int edgeno, int wrtvtxno, DMATRIX **pgrad)
  \brief Computes the cost of an edge as well as its gradient with
  respect to the given edge vertex (0-3). cost = (1-dot).^2 all
  faces. If pgrad=NULL, then the gradient is not computed. If *pgrad
  is NULL, then the 1x3 grad is allocated.
 */
double MRISedgeAngleCostEdgeVertex(MRIS *surf, int edgeno, int wrtvtxno, DMATRIX **pgrad)
{
  MRI_EDGE *e;
  double J;
  DMATRIX *grad=NULL;

  if(pgrad) grad = *pgrad;

  if(wrtvtxno < 0 || wrtvtxno > 3){
    printf("ERROR: MRISedgeAngleCostEdgeVertex(): wrtvtxno=%d, must be 0-3\n",wrtvtxno);
    if(grad) grad = NULL;
    return(-1);
  }

  e = &(surf->edges[edgeno]);
  J = pow(1.0-e->dot,2.0);
  if(grad != NULL){
    grad = DMatrixScalarMul(e->gradDot[wrtvtxno], -2*(1.0-e->dot) , grad);
  }
  return(J);
}

int MRISedgeAngleCostEdgeVertexTest(MRIS *surf, int edgeno, int wrtvtxno, long double delta)
{
  MRI_EDGE *e;
  DMATRIX *g0, *gnum;
  double J0,J1;
  int wrtdimno, c;

  MRISfaceNormalGrad(surf, 0);
  //MRISedgeGradDot(surf);
  
  g0 = DMatrixAlloc(1,3,MATRIX_REAL);
  J0 = MRISedgeAngleCostEdgeVertex(surf, edgeno, wrtvtxno, &g0);
  e = &(surf->edges[edgeno]);
  
  int      const surfvtxno = e->vtxno[wrtvtxno];
  VERTEX * const v         = &surf->vertices[surfvtxno];
  
  //MRISedgePrint(surf, edgeno, stdout);
  //printf("%3d %g %g %g  %g\n",surfvtxno,v->x,v->y,v->z,J0);

  MRISfreeDistsButNotOrig(surf);
    // MRISsetXYZ will invalidate all of these,
    // so make sure they are recomputed before being used again!

  gnum = DMatrixAlloc(1,3,MATRIX_REAL);  
  for(wrtdimno=0; wrtdimno < 3; wrtdimno++){

    float const old_x = v->x, old_y = v->y, old_z = v->z;
    float x = old_x, y = old_y, z = old_z;
      
    if(wrtdimno == 0) x += delta;
    if(wrtdimno == 1) y += delta;
    if(wrtdimno == 2) z += delta;
    MRISsetXYZ(surf, surfvtxno, x,y,z);
    
    MRISfaceNormalGrad(surf, 0);
    //MRISedgeGradDot(surf);
    J1 = MRISedgeAngleCostEdgeVertex(surf, edgeno, wrtvtxno, NULL);
    //printf("  %g %g %g    %g\n",v->x,v->y,v->z,J1);
    gnum->rptr[1][wrtdimno+1] = (J1-J0)/delta;

    MRISsetXYZ(surf, surfvtxno, old_x,old_y,old_z);
  }
  if(0){
  printf("g0    = [");
  for(c=0; c<3; c++) printf("%12.8f ",g0->rptr[1][c+1]);
  printf("];\n");
  printf("gnum  = [");
  for(c=0; c<3; c++) printf("%12.8f ",gnum->rptr[1][c+1]);
  printf("];\n");
  printf("gdiff = [");
  for(c=0; c<3; c++) printf("%12.8f ",gnum->rptr[1][c+1]-g0->rptr[1][c+1]);
  printf("];\n");
  }

  DMatrixFree(&g0);
  DMatrixFree(&gnum);

  // BUG forgets to recompute the metric properties
  
  return(0);
}

/*!
  \fn int MRISfaceNormalGrad(MRIS *surf, int NormOnly)
  \brief Computes the face normal and gradient of the face normal for
  all faces. If NormOnly=1, then the gradient is not computed.
 */
int MRISfaceNormalGrad(MRIS *surf, int NormOnly)
{
  int r = MRISfaceMetric(surf, !NormOnly);
  return(r);
}

/*!
  \fn double MRISedgeGradDotEdgeVertexTest(MRIS *surf, int edgeno, int wrtvtxno, long double delta, int verbose)
  \brief Test of the gradient of the edge dot product for the given
  edge and vertex. Tests by computing a numerical gradient by changing
  a vertex position by delta (mm). Returns the max relative difference
  between the theoretical and numerical.
 */
double MRISedgeGradDotEdgeVertexTest(MRIS *surf, int edgeno, int wrtvtxno, long double delta, int verbose)
{
  int c, wrtdimno,vnosurf,nthface;
  MRI_EDGE *e;
  VERTEX *v;
  double dot0,dot1;
  DMATRIX *g0, *gnum=NULL;
  double diff, maxdiff, maxrdiff, gmax;
  extern int MRISfaceNormalFace_AddDeltaVertex; // cant run test in parallel
  extern long double MRISfaceNormalFace_AddDelta[3];

  e = &(surf->edges[edgeno]);
  vnosurf = e->vtxno[wrtvtxno]; // vertex no wrt the surface
  v = &(surf->vertices[vnosurf]);

  // Compute the dot at baseline
  MRISfaceNormalFace_AddDeltaVertex = -1;
  MRISfaceNormalGradFace(surf, e->faceno[0]);
  MRISfaceNormalGradFace(surf, e->faceno[1]);
  MRISedgeMetricEdge(surf, edgeno, 3); // does all wrtvtxno
  dot0 = e->dot;
  //MRISedgeGradDotEdgeVertex(surf, edgeno, wrtvtxno);
  g0 = e->gradDot[wrtvtxno];
  gmax = DMatrixMaxAbs(g0);

  // Compute the numerical gradient
  gnum = DMatrixAlloc(1,3,MATRIX_REAL);
  maxdiff = 0;
  // Go through each dim in the WRT vertex position
  for(wrtdimno=0; wrtdimno < 3; wrtdimno++){
    // Set the delta for this dim
    if(wrtvtxno==0 || wrtvtxno==1){
      for(nthface = 0; nthface < 2; nthface++){
	MRISfaceNormalFace_AddDeltaVertex = e->corner[wrtvtxno][nthface];
	for(c=0; c < 3; c++) MRISfaceNormalFace_AddDelta[c] = 0;
	MRISfaceNormalFace_AddDelta[wrtdimno] = delta;
	MRISfaceNormalGradFace(surf, e->faceno[nthface]);
	for(c=0; c < 3; c++) MRISfaceNormalFace_AddDelta[c] = 0;
      }
    }
    if(wrtvtxno==2){
      nthface = 0;
      MRISfaceNormalFace_AddDeltaVertex = e->corner[wrtvtxno][nthface];
      for(c=0; c < 3; c++) MRISfaceNormalFace_AddDelta[c] = 0;
      MRISfaceNormalFace_AddDelta[wrtdimno] = delta;
      MRISfaceNormalGradFace(surf, e->faceno[nthface]);
      for(c=0; c < 3; c++) MRISfaceNormalFace_AddDelta[c] = 0;
      MRISfaceNormalFace_AddDeltaVertex = -1;
      MRISfaceNormalGradFace(surf, e->faceno[1]);
    }
    if(wrtvtxno==3){
      nthface = 1;
      MRISfaceNormalFace_AddDeltaVertex = e->corner[wrtvtxno][nthface];
      for(c=0; c < 3; c++) MRISfaceNormalFace_AddDelta[c] = 0;
      MRISfaceNormalFace_AddDelta[wrtdimno] = delta;
      MRISfaceNormalGradFace(surf, e->faceno[nthface]);
      for(c=0; c < 3; c++) MRISfaceNormalFace_AddDelta[c] = 0;
      MRISfaceNormalFace_AddDeltaVertex = -1;
      MRISfaceNormalGradFace(surf, e->faceno[0]);
    }

    MRISedgeMetricEdge(surf, edgeno, 0);
    dot1 = e->dot;
    gnum->rptr[1][wrtdimno+1] = (dot1-dot0)/delta;
    diff = fabs(gnum->rptr[1][wrtdimno+1] - g0->rptr[1][wrtdimno+1]);
    if(maxdiff < diff) maxdiff = diff;
  }

  // The final metric is the maxabsdiff divided by the max component in the 
  // theoretical to make sure scaling is correct
  maxrdiff = maxdiff/gmax;

  // Undo everything
  MRISfaceNormalFace_AddDeltaVertex = -1;
  for(c=0; c < 3; c++) MRISfaceNormalFace_AddDelta[c] = 0;
  for(nthface = 0; nthface < 2; nthface++)
    MRISfaceNormalGradFace(surf, e->faceno[nthface]);

  if(verbose){
    printf("gmax = %lf\n",gmax);
    printf("g0 = [");
    DMatrixPrintFmt(stdout,"%12.10f",g0);
    printf("];\n");
    printf("gnum = [");
    DMatrixPrintFmt(stdout,"%12.10f",gnum);
    printf("];\n");
  }

  DMatrixFree(&gnum);

  return(maxrdiff);
}

#if 0
/*!
  \fn int MRISedgeGradDot(MRIS *surf)
  \brief Computes gradient of the dot product of the normals of the
  two faces adjacent to the all edges with respect to (wrt) all
  four adjacent vertices.
 */
int MRISedgeGradDot(MRIS *surf)
{
  int edgeno;
  #ifdef HAVE_OPENMP
  #pragma omp parallel for 
  #endif
  for(edgeno=0; edgeno < surf->nedges; edgeno++){
    MRI_EDGE *e;
    int vtxno;
    e = &(surf->edges[edgeno]);
    for(vtxno=0; vtxno < 4; vtxno++){
      MRISedgeGradDotEdgeVertex(surf, edgeno, vtxno);
    }
  }
  return(0);
}
#endif

/*!
  \fn int MRISedgeGradDotEdgeVertex(MRIS *surf, int edgeno, int wrtvtxno)
  \brief Computes gradient of the dot product of the normals of the
  two faces adjacent to the given edge with respect to (wrt) the given
  vertex; wrtvtxno = 0, 1, 2, 3.
 */
int MRISedgeGradDotEdgeVertex(MRIS *surf, int edgeno, int wrtvtxno)
{
  MRI_EDGE *e;
  FACE *f0, *f1;
  int face0wrtvtxno,face1wrtvtxno;
  DMATRIX *n0=NULL, *n1=NULL, *gradN0=NULL, *gradN1=NULL;
  DMATRIX *n0T=NULL, *gradN0T=NULL, *tmp1=NULL, *tmp1T=NULL, *tmp2=NULL;

  e = &(surf->edges[edgeno]);

  f0 = &(surf->faces[e->faceno[0]]);
  n0 = f0->norm;
  n0T = DMatrixTranspose(n0,NULL);

  f1 = &(surf->faces[e->faceno[1]]);
  n1 = f1->norm;

  e->dot = DVectorDot(f0->norm,f1->norm);

  face0wrtvtxno = e->corner[wrtvtxno][0];
  if(face0wrtvtxno < 3){
    gradN0 = f0->gradNorm[face0wrtvtxno];
    gradN0T = DMatrixTranspose(gradN0,NULL);
    tmp1 = DMatrixMultiply(gradN0T,n1,NULL);
    tmp1T = DMatrixTranspose(tmp1,NULL);
  }
  face1wrtvtxno = e->corner[wrtvtxno][1];
  if(face1wrtvtxno < 3){
    gradN1 = f1->gradNorm[face1wrtvtxno];
    tmp2 = DMatrixMultiply(n0T,gradN1,NULL);
  }

  // gradDot = (gradN0T)*n1 + n0T*(gradN1) = 1x3
  if(tmp1T && tmp2)
    e->gradDot[wrtvtxno] = DMatrixAdd(tmp1T,tmp2,e->gradDot[wrtvtxno]);
  else if(tmp1T) 
    e->gradDot[wrtvtxno] = DMatrixCopy(tmp1T,e->gradDot[wrtvtxno]);
  else
    e->gradDot[wrtvtxno] = DMatrixCopy(tmp2,e->gradDot[wrtvtxno]);

  // Don't free n0, n1, gradN0, gradN1
  DMatrixFree(&n0T);
  if(gradN0T) DMatrixFree(&gradN0T);
  DMatrixFree(&tmp1);
  DMatrixFree(&tmp1T);
  DMatrixFree(&tmp2);

  return(0);
}

/*!
  \fn int MRISedgePrint(MRIS *surf, int edgeno, FILE *fp)
  \brief Prints some info about an edge
 */
int MRISedgePrint(MRIS *surf, int edgeno, FILE *fp)
{
  MRI_EDGE *e;
  int c;
  VERTEX *v;
  e = &(surf->edges[edgeno]);
  fprintf(fp,"edgeno %d\n",edgeno);
  fprintf(fp,"len %lf\n",e->len);
  fprintf(fp,"dot %lf\n",e->dot);
  fprintf(fp,"angle %lf\n",e->angle);
  fprintf(fp,"faces %d %d\n",e->faceno[0],e->faceno[1]);
  fprintf(fp,"vertices %d %d %d %d\n",e->vtxno[0],e->vtxno[1],e->vtxno[2],e->vtxno[3]);
  for(c=0; c<4; c++){
    v = &(surf->vertices[e->vtxno[c]]);
    printf("v%d = [%f %f %f];\n",c,v->x,v->y,v->z);
  }
  for(c=0; c<4; c++)
    fprintf(fp,"corner %d %d %d\n",c,e->corner[c][0],e->corner[c][1]);
  MRISfacePrint(surf, e->faceno[0], fp);
  MRISfacePrint(surf, e->faceno[1], fp);
  return(0);
}

/*!
  \fn int MRISfacePrint(MRIS *surf, int faceno, FILE *fp)
  \brief Prints some info about a face.
 */
int MRISfacePrint(MRIS *surf, int faceno, FILE *fp)
{
  FACE *f;
  VERTEX *v;
  int c;
  f = &(surf->faces[faceno]);
  fprintf(fp,"faceno %d\n",faceno);
  fprintf(fp,"area %f\n",f->area);
  if(f->norm)
    fprintf(fp,"norm %f %f %f\n",f->norm->rptr[1][1],f->norm->rptr[2][1],f->norm->rptr[3][1]);
  fprintf(fp,"vertices %d %d %d\n",f->v[0],f->v[1],f->v[2]);
  for(c=0; c<3; c++){
    v = &(surf->vertices[f->v[c]]);
    printf("v%d = [%f %f %f];\n",c,v->x,v->y,v->z);
  }
  return(0);
}


/*!
  \fn BBRFACE *BBRFaceAlloc(void)
  \brief Allocates the BBRFACE structure
 */
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
/*!
  \fn int BBRFaceFree(BBRFACE **pbbrf)
  \brief Frees the BBRFACE structure
 */
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
  free(bbrf);
  *pbbrf = NULL;
  return(0);
}
/*!
  \fn BBRFACE *BBRFaceDiff(const BBRFACE *bbrf1, const BBRFACE *bbrf2, const double delta, BBRFACE *bbrfout)
  \brief Compute the difference between two BBRFACE structures. This can be useful when computing numerical
  gradients by setting delta to the amount a vertex position was changed between computing bbrf1 and bbrf2.
 */
BBRFACE *BBRFaceDiff(const BBRFACE *bbrf1, const BBRFACE *bbrf2, const double delta, BBRFACE *bbrfout)
{
  int n;
  double f = 1.0/delta;
  if(bbrfout == NULL) bbrfout = BBRFaceAlloc();
  bbrfout->faceno = bbrf1->faceno;

  DMatrixAddMul(bbrf1->pc,  bbrf2->pc,f,-f,bbrfout->pc);
  DMatrixAddMul(bbrf1->pin, bbrf2->pin,f,-f,bbrfout->pin);
  DMatrixAddMul(bbrf1->pout,bbrf2->pout,f,-f,bbrfout->pout);
  DMatrixAddMul(bbrf1->vin, bbrf2->vin,f,-f,bbrfout->vin);
  DMatrixAddMul(bbrf1->vout,bbrf2->vout,f,-f,bbrfout->vout);
  for(n=0; n < 3; n++){
    DMatrixAddMul(bbrf1->gradPin[n],      bbrf2->gradPin[n],f,-f,bbrfout->gradPin[n]);
    DMatrixAddMul(bbrf1->gradPout[n],     bbrf2->gradPout[n],f,-f,bbrfout->gradPout[n]);
    DMatrixAddMul(bbrf1->gradVin[n],      bbrf2->gradVin[n],f,-f,bbrfout->gradVin[n]);
    DMatrixAddMul(bbrf1->gradVout[n],     bbrf2->gradVout[n],f,-f,bbrfout->gradVout[n]);
    DMatrixAddMul(bbrf1->gradInterpIn[n], bbrf2->gradInterpIn[n],f,-f,bbrfout->gradInterpIn[n]);
    DMatrixAddMul(bbrf1->gradInterpOut[n],bbrf2->gradInterpOut[n],f,-f,bbrfout->gradInterpOut[n]);
    DMatrixAddMul(bbrf1->gradIin[n],      bbrf2->gradIin[n],f,-f,bbrfout->gradIin[n]);
    DMatrixAddMul(bbrf1->gradIout[n],     bbrf2->gradIout[n],f,-f,bbrfout->gradIout[n]);
    DMatrixAddMul(bbrf1->gradQ[n],        bbrf2->gradQ[n],f,-f,bbrfout->gradQ[n]);
    DMatrixAddMul(bbrf1->gradCost[n],     bbrf2->gradCost[n],f,-f,bbrfout->gradCost[n]);
  }
  bbrfout->Iin  = (bbrf1->Iin  - bbrf2->Iin)/delta;
  bbrfout->Iout = (bbrf1->Iout - bbrf2->Iout)/delta;
  bbrfout->Q    = (bbrf1->Q    - bbrf2->Q)/delta;
  bbrfout->cost = (bbrf1->cost - bbrf2->cost)/delta;
  return(bbrfout);
}

/*!
  \fn int BBRFacePrint(FILE *fp, BBRFACE *bbrf, BBRPARAMS *bbrpar)
  \brief Prints info about the BBRFACE structure used to manage the BBR
  costs and gradients.
 */
int BBRFacePrint(FILE *fp, BBRFACE *bbrf, BBRPARAMS *bbrpar)
{
  FACE *f;
  VERTEX *v;
  int n;
  f = &(bbrpar->surf->faces[bbrf->faceno]);
  fprintf(fp,"faceno %d\n",bbrf->faceno);
  fprintf(fp,"area %12.10f\n",f->area);
  fprintf(fp,"interp %d\n",bbrpar->interp);
  fprintf(fp,"vtxno ");   for(n=0; n<3; n++)   fprintf(fp,"%d ",f->v[n]);  fprintf(fp,"\n");
  for(n=0; n<3; n++){
    v = &(bbrpar->surf->vertices[f->v[n]]);
    fprintf(fp,"vtx%d %g %g %g\n",n,v->x,v->y,v->z);
  }
  fprintf(fp,"fnorm ");for(n=0; n<3; n++)  fprintf(fp,"%12.10lf ",f->norm->rptr[n+1][1]); fprintf(fp,"\n");
  fprintf(fp,"pc ");   for(n=0; n<3; n++)  fprintf(fp,"%6.4lf ",bbrf->pc->rptr[n+1][1]); fprintf(fp,"\n");
  fprintf(fp,"pin ");   for(n=0; n<3; n++) fprintf(fp,"%6.4lf ",bbrf->pin->rptr[n+1][1]); fprintf(fp,"\n");
  fprintf(fp,"pout ");   for(n=0; n<3; n++)fprintf(fp,"%6.4lf ",bbrf->pout->rptr[n+1][1]); fprintf(fp,"\n");
  fprintf(fp,"vin ");   for(n=0; n<3; n++) fprintf(fp,"%6.4lf ",bbrf->vin->rptr[n+1][1]); fprintf(fp,"\n");
  fprintf(fp,"vout ");   for(n=0; n<3; n++)fprintf(fp,"%6.4lf ",bbrf->vout->rptr[n+1][1]); fprintf(fp,"\n");
  fprintf(fp,"Iin  %lf\n",bbrf->Iin);
  fprintf(fp,"Iout %lf\n",bbrf->Iout);
  fprintf(fp,"Q %lf\n",bbrf->Q);
  fprintf(fp,"cost %12.10lf\n",bbrf->cost);
  fflush(fp);
  return(0);
}

/*!
  \fn BBRFACE *BBRCostFace(int faceno, int wrtvtxno, BBRPARAMS *bbrpar, BBRFACE *bbrf)
  \brief Computes the BBR cost at the given face, and, optionally, its
  gradient with respect to (wrt) the given vertex. If wrtvtxno<0, then
  only the cost is computed.  If wrtvtxno>0, then only the gradient of
  cost is computed. wrtvtxno should be 0, 1, or 2, corresponding to
  the three vertices of the triangle. Intermediate values and
  gradients are computed. All values are stored in the BBRFACE
  structure.
 */
BBRFACE *BBRCostFace(int faceno, int wrtvtxno, BBRPARAMS *bbrpar, BBRFACE *bbrf)
{
  FACE *f;
  VERTEX *v;
  int n,oob=0;
  double c,r,s;
  double OutInSum, OutInDiff, vtanh;
  DMATRIX *tmp1, *tmp2;
  MRIS *surf = bbrpar->surf;

  if(bbrf==NULL)
    bbrf = BBRFaceAlloc();

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
    if(bbrpar->interp == SAMPLE_TRILINEAR)
      oob = MRIsampleVolume(bbrpar->mri, c,r,s, &bbrf->Iin);
    if(bbrpar->interp == SAMPLE_CUBIC)
      oob = MRIcubicSampleVolume(bbrpar->mri, c,r,s, &bbrf->Iin);
    if(oob){
      printf("ERROR: BBRCostFace(): in-point out of volume\n");fflush(stdout);
      BBRFacePrint(stdout, bbrf, bbrpar);
      return(NULL);
    }
    
    // Compute the XYZ at the exterior point = pc + norm*Dout
    // Note: pc and norm are 3x1 whereas pin is 4x1
    for(n=0; n<3; n++)
      bbrf->pout->rptr[n+1][1] = bbrf->pc->rptr[n+1][1] + (bbrpar->Dout*f->norm->rptr[n+1][1]);
    
    // Compute the CRS at the exterior point = ras2vox*pout
    bbrf->vout = DMatrixMultiply(bbrpar->sras2vox,bbrf->pout,bbrf->vout);
    
    // Compute the intensity at the exterior point
    c=bbrf->vout->rptr[1][1]; r=bbrf->vout->rptr[1][2]; s=bbrf->vout->rptr[1][3];
    if(bbrpar->interp == SAMPLE_TRILINEAR)
      oob = MRIsampleVolume(bbrpar->mri, c,r,s, &bbrf->Iout);
    if(bbrpar->interp == SAMPLE_CUBIC)
      oob = MRIcubicSampleVolume(bbrpar->mri, c,r,s, &bbrf->Iout);
    if(oob){
      printf("ERROR: BBRCostFace(): out-point out of volume\n");fflush(stdout);
      BBRFacePrint(stdout, bbrf, bbrpar);
      return(NULL);
    }

    // Finally, compute the cost
    OutInSum  = bbrf->Iout + bbrf->Iin + 10e-10;
    OutInDiff = bbrf->Iout - bbrf->Iin;
    bbrf->Q = 100*(OutInDiff)/(0.5*(OutInSum));
    vtanh = tanh(bbrpar->M*(bbrf->Q - bbrpar->Q0));
    bbrf->cost = 1 + vtanh;

    return(bbrf);
  } // if(wrtvtxno < 0)

  // Now compute the gradients

  // Recompute these
  OutInSum  = bbrf->Iout + bbrf->Iin + 10e-10;
  OutInDiff = bbrf->Iout - bbrf->Iin;
  vtanh = tanh(bbrpar->M*(bbrf->Q - bbrpar->Q0));

  // gradPin = gradPc - Din*gradNorm, gradPc = eye(3)/3
  bbrf->gradPin[wrtvtxno] = DMatrixScalarMul(f->gradNorm[wrtvtxno],-bbrpar->Din,bbrf->gradPin[wrtvtxno]);
  for(n=0; n<3; n++) bbrf->gradPin[wrtvtxno]->rptr[n+1][n+1] += 1/3.0;

  // gradVin = ras2vox*gradPin
  bbrf->gradVin[wrtvtxno] = DMatrixMultiply(bbrpar->sras2vox3x3,bbrf->gradPin[wrtvtxno],bbrf->gradVin[wrtvtxno]);

  // Compute gradient of the sampled intensity WRT a change in col, row, or slice
  c=bbrf->vin->rptr[1][1]; r=bbrf->vin->rptr[1][2]; s=bbrf->vin->rptr[1][3];
  if(bbrpar->interp == SAMPLE_TRILINEAR)
    bbrf->gradInterpIn[wrtvtxno] = MRIgradTrilinInterp(bbrpar->mri, c,r,s, bbrf->gradInterpIn[wrtvtxno]);
  if(bbrpar->interp == SAMPLE_CUBIC)
    bbrf->gradInterpIn[wrtvtxno] = MRIgradCubicInterp(bbrpar->mri, c,r,s, bbrf->gradInterpIn[wrtvtxno]);
  if(bbrf->gradInterpIn[wrtvtxno] == NULL){
    printf("ERROR: BBRCostFace(): grad %d in-point out of volume\n",wrtvtxno);fflush(stdout);
    BBRFacePrint(stdout, bbrf, bbrpar);
    return(NULL);
  }

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
  if(bbrpar->interp == SAMPLE_TRILINEAR)
    bbrf->gradInterpOut[wrtvtxno] = MRIgradTrilinInterp(bbrpar->mri, c,r,s, bbrf->gradInterpOut[wrtvtxno]);
  if(bbrpar->interp == SAMPLE_CUBIC)
    bbrf->gradInterpOut[wrtvtxno] = MRIgradCubicInterp(bbrpar->mri, c,r,s, bbrf->gradInterpOut[wrtvtxno]);
  if(bbrf->gradInterpOut[wrtvtxno] == NULL){
    printf("ERROR: BBRCostFace(): grad %d out-point out of volume\n",wrtvtxno);fflush(stdout);
    BBRFacePrint(stdout, bbrf, bbrpar);
    return(NULL);
  }

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

/*!
  \fn double TestBBRCostFace(BBRPARAMS *bbrpar, int faceno, int wrtvtxno, double delta, int verbose)
  \brief Runs a test to numerically measure the gradient of the BBR
  cost and compares it against the theortical for the given face and
  vertex. Returns the maximum relative difference.
*/
double TestBBRCostFace(BBRPARAMS * const bbrpar, int const faceno, int const wrtvtxno, double const delta, int const verbose)
{
  // Compute baseline cost and gradient
  MRISfaceNormalGrad(bbrpar->surf, 0);

  BBRFACE *bbrf1;
  bbrf1 = BBRCostFace(faceno, -1,       bbrpar, NULL);
  bbrf1 = BBRCostFace(faceno, wrtvtxno, bbrpar, bbrf1);

  if (verbose) {
    printf("wrtvtxno = %d, delta = %lf \n",wrtvtxno,delta);
    printf("BBRF1 ----------\n");
    BBRFacePrint(stdout,bbrf1,bbrpar);
  }

  FACE*   const f      = &bbrpar->surf->faces[faceno];
  int     const svtxno = f->v[wrtvtxno];
  MRIS*   const mris   = bbrpar->surf;
  VERTEX* const v      = &mris->vertices[svtxno];
  
  double maxgrad = 0, maxdiff = 0;
  
  MRISfreeDistsButNotOrig(mris);
    // MRISsetXYZ will invalidate all of these,
    // so make sure they are recomputed before being used again!

  // Go through each dimension in the WRT vertex
  //
  BBRFACE* bbrfgrad = NULL;     // these are reused around the loop
  BBRFACE* bbrf2    = NULL; 

  int wrtdimno;
  for (wrtdimno = 0; wrtdimno < 3; wrtdimno++) {

    // Doing the delta here is convenient but sometimes inaccurate
    // because xyz are floats so limited on how small delta can be.
    //
    float const old_x = v->x, old_y = v->y, old_z = v->z;
    
    float x = old_x, y = old_y, z = old_z;
    if (wrtdimno==0) x += delta;
    if (wrtdimno==1) y += delta;
    if (wrtdimno==2) z += delta;
    MRISsetXYZ(mris,svtxno,x,y,z);
    
    // Compute face normal (but not grad of norm). 
    // Runs for all faces, inefficient, but only a test
    //
    MRISfaceNormalGrad(bbrpar->surf, 1); 

    // Compute BBRFace cost, etc, not grads
    bbrf2 = BBRCostFace(faceno, -1, bbrpar, bbrf2); 

    // Numerically compute the grad by taking the diff and dividing by delta.  
    bbrfgrad = BBRFaceDiff(bbrf2, bbrf1, delta, bbrfgrad);

    // This will work for all the parameters, but just tesing cost
    // gradient -- if that is correct, then everything upstream must
    // be correct. But upstream grads could be computed in the same way.
    double const diff = bbrf1->gradCost[wrtvtxno]->rptr[1][wrtdimno+1]-bbrfgrad->cost;

    // Keep track of the max gradient and the max diff bet theoretical and numerical
    if (maxgrad < fabs(bbrf1->gradCost[wrtvtxno]->rptr[1][wrtdimno+1]))
        maxgrad = fabs(bbrf1->gradCost[wrtvtxno]->rptr[1][wrtdimno+1]);
        
    if (maxdiff < fabs(diff)) 
        maxdiff = fabs(diff) ;

    if (verbose) {
      int n;
    
      printf("wrtdim %d =================================================\n",wrtdimno);
      printf("%5d %5d %d %d %12.10lf %12.10lf %12.10lf\n",faceno,svtxno,wrtvtxno,wrtdimno,
	     bbrf1->gradCost[wrtvtxno]->rptr[1][wrtdimno+1],bbrfgrad->cost,diff);
      printf("BBRF2 ----------\n");
      BBRFacePrint(stdout,bbrf2,bbrpar);
      printf("BBRF Theoretical Grad ----------\n");
      printf("gradPin  = "); for(n=0;n<3;n++) printf("%lf ",bbrf1->gradPin[wrtvtxno]->rptr[n+1][wrtdimno+1]); printf("\n");
      printf("gradPout = "); for(n=0;n<3;n++) printf("%lf ",bbrf1->gradPout[wrtvtxno]->rptr[n+1][wrtdimno+1]); printf("\n");
      printf("gradVin  = "); for(n=0;n<3;n++) printf("%lf ",bbrf1->gradVin[wrtvtxno]->rptr[n+1][wrtdimno+1]); printf("\n");
      printf("gradVout = "); for(n=0;n<3;n++) printf("%lf ",bbrf1->gradVout[wrtvtxno]->rptr[n+1][wrtdimno+1]); printf("\n");
      printf("gradIin  = %lf\n",bbrf1->gradIin[wrtvtxno]->rptr[1][wrtdimno+1]);
      printf("gradIout = %lf\n",bbrf1->gradIout[wrtvtxno]->rptr[1][wrtdimno+1]);
      printf("gradQ    = %lf\n",bbrf1->gradQ[wrtvtxno]->rptr[1][wrtdimno+1]);
      printf("gradCost = %lf\n",bbrf1->gradCost[wrtvtxno]->rptr[1][wrtdimno+1]);
      printf("BBRF Numerical Grad ----------\n");
      BBRFacePrint(stdout,bbrfgrad,bbrpar);
    }

    // restore
    // the old code subtracted delta, but floating point arithmetic does not require that (f+d)-f == f
    //
    MRISsetXYZ(mris,svtxno,old_x,old_y,old_z);
  }
    
  // Compute the maximum diff relative to the maximum gradient
  double const maxrdiff = maxdiff/(maxgrad+10e-10);

  if(verbose){
    printf("maxgrad = %lf maxdiff = %lf  maxrdiff = %lf\n",maxgrad,maxdiff,maxrdiff);
    fflush(stdout);
  }

  // Clean up
  BBRFaceFree(&bbrf1);
  BBRFaceFree(&bbrf2);
  BBRFaceFree(&bbrfgrad);

  // Recompute baseline cost and gradient
  MRIScomputeMetricProperties(bbrpar->surf);
  MRISfaceNormalGrad(bbrpar->surf, 0);

  return(maxrdiff);
}

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

long double MRISbbrCost(BBRPARAMS *bbrpar, DMATRIX *gradCost)
{
  long double cost = 0;
  int faceno,nthreads,threadno;
  BBRFACE **bbrfth;
  static DMATRIX **gradCostTh=NULL;
  static int nthreadsalloc=0;
  extern int MRISbbrCostFree;

  // This must already have been run
  //MRISfaceNormalGrad(surf, 0);

  // Get number of threads
  nthreads = 1;
  #ifdef HAVE_OPENMP
  nthreads = omp_get_max_threads();  // using max should be ok
  #endif

  if(nthreadsalloc > 0 && (nthreadsalloc < nthreads || MRISbbrCostFree) ){
    // Free the static gradCostTh if
    // threads have been alloced AND either:
    //   The number of threads that have been alloced is less than the number of threads
    //   The caller wants to free the data
    for(threadno = 0; threadno < nthreads; threadno++)
      DMatrixFree(&gradCostTh[threadno]);
    free(gradCostTh);
    gradCostTh = NULL;
    nthreadsalloc = 0;
    if(MRISbbrCostFree) return(0);
  }

  if(gradCost){
    if(gradCostTh==NULL){
      nthreadsalloc = nthreads;
      gradCostTh = (DMATRIX **) calloc(sizeof(DMATRIX *),nthreads);
      for(threadno = 0; threadno < nthreads; threadno++){
	gradCostTh[threadno] = DMatrixAlloc(bbrpar->surf->nvertices,3,MATRIX_REAL);
      }
    }
    else{
      // The static gradCostTh have already been alloced, so zero all the threads
      for(threadno = 0; threadno < nthreads; threadno++)
	DMatrixZero(0,0,gradCostTh[threadno]);
    }
  }

  bbrfth = (BBRFACE **) calloc(sizeof(BBRFACE*),nthreads);
  for(threadno = 0; threadno < nthreads; threadno++){
    bbrfth[threadno] = BBRFaceAlloc();
  }

  cost = 0;
  #ifdef HAVE_OPENMP
  #pragma omp parallel for private(threadno) reduction(+ : cost)
  #endif
  for(faceno = 0; faceno < bbrpar->surf->nfaces; faceno++){
    BBRFACE *bbrf;

    threadno=0;
    #ifdef HAVE_OPENMP
    threadno = omp_get_thread_num();
    #endif

    bbrf = bbrfth[threadno];
    bbrf = BBRCostFace(faceno, -1, bbrpar, bbrf); // Compute cost only, no grad
    cost += bbrf->cost;
    if(gradCost){
      FACE *f = &(bbrpar->surf->faces[faceno]);
      int  wrtvtxno, svtxno, n;
      for(wrtvtxno=0; wrtvtxno < 3; wrtvtxno++){
	bbrf = BBRCostFace(faceno, wrtvtxno, bbrpar, bbrf); // Compute grads
	svtxno = f->v[wrtvtxno];
	for(n=0; n<3; n++)
	  gradCostTh[threadno]->rptr[svtxno+1][n+1] += bbrf->gradCost[wrtvtxno]->rptr[1][n+1];
      }
    }
  }
  // Change cost to average cost over number of faces
  cost = cost/bbrpar->surf->nfaces;

  // Merge gradients from different threads
  if(gradCost){
    DMatrixZero(0,0,gradCost);
    for(threadno = 0; threadno < nthreads; threadno++)
      gradCost = DMatrixAdd(gradCost,gradCostTh[threadno],gradCost);
    DMatrixScalarMul(gradCost,1.0/bbrpar->surf->nfaces,gradCost);
  }

  // Free the BBRFACE structs
  for(threadno = 0; threadno < nthreads; threadno++)
    BBRFaceFree(&bbrfth[threadno]);
  free(bbrfth);

  return(cost);
}

double MRISedgeCost(MRIS *surf, DMATRIX *gradCost)
{
  long double cost;
  int edgeno, wrtvtxno, svtxno,n;
  MRI_EDGE *e;
  DMATRIX *gradCostEV=NULL;

  cost = 0;
  for(edgeno = 0; edgeno < surf->nedges; edgeno++){
    e = &(surf->edges[edgeno]);
    // easy enough to compute actual cost here
    cost += pow(1.0-e->dot, 2.0);
    if(gradCost){
      for(wrtvtxno=0; wrtvtxno < 4; wrtvtxno++){
	gradCostEV = DMatrixScalarMul(e->gradDot[wrtvtxno],-2*(1.0-e->dot),gradCostEV);
	svtxno = e->vtxno[wrtvtxno];
	for(n=0; n<3; n++)
	  gradCost->rptr[svtxno+1][n+1] += gradCostEV->rptr[1][n+1];
      }
    }
  }
  cost /= surf->nedges;

  if(gradCost)
    DMatrixScalarMul(gradCost,1.0/surf->nedges,gradCost);

  return(cost);
}

/*!
  \fn long double MRISedgeLengthCostEdge(MRIS *surf, int edgeno, double L0, DMATRIX **pgradv0, DMATRIX **pgradv1)
  \brief Compute the cost and cost gradient of the given edge assuming
  its desired length is L0.  This is a spring cost with a
  (potentially) non-zero resting length of L0. If *pgradv{0,1} is
  non-NULL, then the gradient of the cost with respect to a change in
  vertex{0,1} position is computed. *pgradvX should be a 3x1 DMATRIX.
 */
long double MRISedgeLengthCostEdge(MRIS *surf, int edgeno, double L0, DMATRIX **pgradv0, DMATRIX **pgradv1)
{
  long double cost, dL, a;
  MRI_EDGE *e;
  VERTEX *v0, *v1;

  e = &(surf->edges[edgeno]);

  // diff between the actual and desired edge length
  dL = e->len - L0; 
  cost = dL*dL;
  if(*pgradv0==NULL && *pgradv1==NULL) return(cost); // done

  // Compute the gradient of the cost wrt the vertex position
  a = 2.0*dL/e->len; // compute constant for speed

  v0 = &(surf->vertices[e->vtxno[0]]);
  v1 = &(surf->vertices[e->vtxno[1]]);

  if(*pgradv0){
    // change in the cost wrt a change in vertex0 position
    (*pgradv0)->rptr[1][1] = a*(v0->x - v1->x);
    (*pgradv0)->rptr[2][1] = a*(v0->y - v1->y);
    (*pgradv0)->rptr[3][1] = a*(v0->z - v1->z);
  }

  if(*pgradv1){
    // change in the cost wrt a change in vertex1 position
    (*pgradv1)->rptr[1][1] = a*(v1->x - v0->x);
    (*pgradv1)->rptr[2][1] = a*(v1->y - v0->y);
    (*pgradv1)->rptr[3][1] = a*(v1->z - v0->z);
  }

  return(cost);
}

/*!
  \fn double MRISedgeLengthCost(MRIS *surf, double L0, int DoGrad)
  \brief Compute the total edge length cost over all edges assuming
  the desired length of each edge is L0. This is a spring cost with a
  (potentially) non-zero resting length of L0. If DoGrad==1, then the
  gradient of the total cost wrt each vertex position is computed and
  stored in v->t{xyz}. If weight>=0, then the total cost and gradients
  are multiplied by weight and the gradients are stored in
  v->d{xyz}. The total cost and gradients are normalized by the number
  of edges. Skips edges where both vertices are ripped.
*/
double MRISedgeLengthCost(MRIS *surf, double L0, double weight, int DoGrad)
{
  int eno, vno, nhits;
  long double totcost,ecost;
  DMATRIX *gradv0=NULL, *gradv1=NULL;
  MRI_EDGE *e;
  VERTEX *v;
  double *tx,*ty,*tz;

  if(weight == 0.0) return(0);

  if(DoGrad){
    gradv0 = DMatrixAlloc(3,1,MATRIX_REAL);
    gradv1 = DMatrixAlloc(3,1,MATRIX_REAL);
    tx = (double*)calloc(sizeof(double),surf->nvertices);
    ty = (double*)calloc(sizeof(double),surf->nvertices);
    tz = (double*)calloc(sizeof(double),surf->nvertices);
  }

  // To parallelized this, need txyz for each thread, then combine them 
  totcost = 0;
  nhits = 0;
  for(eno = 0; eno < surf->nedges; eno++){
    e = &(surf->edges[eno]);
    // Skip this edge if both vertices are ripped
    // Allowing one to be ripped makes a "soap bubble" possible
    VERTEX *v0 = &(surf->vertices[e->vtxno[0]]);
    VERTEX *v1 = &(surf->vertices[e->vtxno[1]]);
    if(v0->ripflag && v1->ripflag)  continue;

    nhits++;
    ecost = MRISedgeLengthCostEdge(surf, eno, L0, &gradv0, &gradv1);
    totcost += ecost;
    if(!DoGrad) continue;
    tx[e->vtxno[0]] += gradv0->rptr[1][1];
    ty[e->vtxno[0]] += gradv0->rptr[2][1];
    tz[e->vtxno[0]] += gradv0->rptr[3][1];
    tx[e->vtxno[1]] += gradv1->rptr[1][1];
    ty[e->vtxno[1]] += gradv1->rptr[2][1];
    tz[e->vtxno[1]] += gradv1->rptr[3][1];
  }
  // normalize to the number of edges
  totcost /= nhits;
  if(weight > 0.0) totcost *= weight;

  // This can be easily parallelized
  if(DoGrad){
    double a = weight/nhits;
    for(vno = 0; vno < surf->nvertices; vno++){
      v = &(surf->vertices[vno]);
      if (v->ripflag)  continue;
      if(weight > 0.0){
	// This is for when called from the surface placement optimizer
	v->dx += -a*tx[vno];
	v->dy += -a*ty[vno];
	v->dz += -a*tz[vno];
      }
      else {
	v->tx = tx[vno]/nhits;
	v->ty = ty[vno]/nhits;
	v->tz = tz[vno]/nhits;
      }
    }
    DMatrixFree(&gradv0);
    DMatrixFree(&gradv1);
    free(tx);
    free(ty);
    free(tz);
  }

  return(totcost);
}

/*!
  \fn int MRISedgeLengthCostTest(MRIS *surf, double delta, double thresh, double L0)
  \brief Tests the theoretical gradient calculation by comparing it to
  an emperical value esimated by changing a vertex position by delta,
  recomputing the cost, and dividing the cost difference by
  delta. This is done for x, y, and z. The angle between the
  theoretical and emperical gradients are then computed.  An "error"
  is generated when this angle is greather than delta. The number of
  errors is returned.  At this point, I expect that differences
  between the theoretical and emperical are due to errors in the
  emperical.  In principle, one should be able to reduce these error
  to arbitrary amounts by making delta very small, but there are
  limits due to the floating point precision used in the surface
  structure. In theory, this would only need to be rerun if the
  underlying code is changed (ie, it should not matter what the
  surface is). Modifies vertex t{xyz} (theoretical gradients), t2{xyz}
  (emperical gradients), and val2bak (angle).
 */
int MRISedgeLengthCostTest(MRIS *surf, double delta, double thresh, double L0)
{
  int vno, k, nerrs, vmpeak0, vmpeak1;
  long double cost0, cost, dcost[3], expdcost[3], d2sum, expd2sum;
  double dot, angle, angleavg, anglemax, expdcostnorm[3], dcostnorm[3];
  double x,y,z, meanlen;
  VERTEX *v0;
  int edgeno;

  // Set up the surface
  MRISedges(surf);
  MRIScomputeMetricProperties(surf);
  MRISedgeMetric(surf,1); //1=DoGrad= only do gradU needed for this cost

  meanlen = 0;
  for(edgeno = 0; edgeno < surf->nedges; edgeno++){
    meanlen += surf->edges[edgeno].len;
  }
  meanlen /= surf->nedges;

  printf("Starting MRISedgeLengthCostTest()\n");
  printf("delta = %g, thresh = %g, L0=%g\n",delta,thresh,L0);
  printf("nvertices %d, nedges = %d, meanlen = %g\n",surf->nvertices,surf->nedges,meanlen);
  if(L0 < 0){
    printf("Setting L0 to mean length\n");
    L0 = meanlen;
  }

  // Compute the cost with this configuration of vertices. Also
  // compute the gradient at each vertex.
  cost0 = MRISedgeLengthCost(surf, L0, -1.0, 1);
  printf("cost0 = %Lf\n",cost0);

  vmpeak0 = GetVmPeak();
  printf("#VMPC# VmPeak  %d\n",vmpeak0);

  // Go thru each vertex
  nerrs = 0;
  anglemax = 0;
  angleavg = 0;
  for(vno = 0; vno < surf->nvertices; vno++){
    v0 = &(surf->vertices[vno]);
    // Expected (theoretical) gradient at this vertx
    expdcost[0] = v0->tx;
    expdcost[1] = v0->ty;
    expdcost[2] = v0->tz;
    // Get the original position so can restore
    x = v0->x;
    y = v0->y;
    z = v0->z;
    // perturb x, y, and z separately
    expd2sum = 0;
    d2sum = 0;
    for(k=0; k < 3; k++){
      if(k==0) v0->x += delta;
      if(k==1) v0->y += delta;
      if(k==2) v0->z += delta;
      // This step is slow. Could just do the edges connected to this vertex, but
      // probably should use the full function
      MRISedgeMetric(surf,0);
      // In theory, this could just be done just for the connected edges, but, again, 
      // probably should use the full function
      cost = MRISedgeLengthCost(surf, L0, -1.0, 0);
      // Compute the emperical gradient
      dcost[k] = (cost-cost0)/delta;
      // Keep track of the sum^2 for dot product calc below
      d2sum    += (dcost[k]*dcost[k]);
      expd2sum += (expdcost[k]*expdcost[k]);
      // restore for next cycle
      v0->x = x;
      v0->y = y;
      v0->z = z;
    } //k

    // Store emperical gradient
    v0->t2x = dcost[0];
    v0->t2y = dcost[1];
    v0->t2z = dcost[2];

    // Compute the dot product between the theoretical and emperical gradients
    dot = 0;
    for(k=0; k < 3; k++){
      expdcostnorm[k] = expdcost[k]/sqrt(expd2sum);
      dcostnorm[k]    = dcost[k]/sqrt(d2sum);
      dot += (expdcostnorm[k]*dcostnorm[k]);
    }
    if(dot > +1.0) dot = +1.0; // might be a slight overflow
    if(dot < -1.0) dot = -1.0; // might be a slight overflow
    // Now compute the angle (ideally, this would be 0)
    angle = acos(dot)*180/M_PI;
    v0->val2bak = angle;
    angleavg += angle;
    if(anglemax < angle) anglemax = angle;
    if(angle>thresh){
      nerrs++;
      printf("%5d %4d  %6.4f %6.4f    ",nerrs,vno,angle,anglemax);
      for(k=0; k < 3; k++) printf("%14.9Lf  %14.9Lf   ",expdcost[k],dcost[k]);
      printf("\n");
      fflush(stdout);
    }
  }
  angleavg /= surf->nvertices;

  vmpeak1 = GetVmPeak();
  printf("#VMPC# VmPeak  %d\n",vmpeak1);
  if(vmpeak1 > vmpeak0)  printf("WARNING: there might be a memory leak\n");
  else                   printf("No memory leak detected\n");

  printf("error rate %g, anglemax = %g, angleavg = %g\n",
	 (double)nerrs/(3.0*surf->nvertices),anglemax,angleavg);
  fflush(stdout);
  return(nerrs);
}

/*!
  \fn long double MRISedgeAngleCost(MRIS *surf, double weight, int DoGrad)
  \brief Computes the total angle cost of all non-ripped vertices. If
  DoGrad=1, then computes the gradient of the cost with respect to
  each vertex and saves in v->t{xyz}. If weight>0, then the cost and
  gradients are multipled by weight, and the negative of the weighted
  gradient is added to v->d{xyz} to make it compatible with the
  surface placement optimization code. edgeAngle=hinge
 */
long double MRISedgeAngleCost(MRIS *surf, double weight, int DoGrad)
{
  int eno, wrtvtxno, vno, skip, nhits;
  long double totcost,ecost;
  DMATRIX *gradv=NULL;
  MRI_EDGE *e;
  VERTEX *v;
  double *tx,*ty,*tz;

  if(weight == 0.0) return(0);

  if(DoGrad){
    gradv = DMatrixAlloc(1,3,MATRIX_REAL); // row vector
    tx = (double*)calloc(sizeof(double),surf->nvertices);
    ty = (double*)calloc(sizeof(double),surf->nvertices);
    tz = (double*)calloc(sizeof(double),surf->nvertices);
  }

  nhits = 0;
  totcost = 0;
  for(eno = 0; eno < surf->nedges; eno++){
    e = &(surf->edges[eno]);

    skip = 0;
    for(wrtvtxno = 0; wrtvtxno < 4; wrtvtxno++){
      v = &(surf->vertices[e->vtxno[wrtvtxno]]);
      if (v->ripflag) {
	skip = 1;
	break;
      }
    }
    if(skip) continue;
    nhits ++;

    for(wrtvtxno = 0; wrtvtxno < 4; wrtvtxno++){
      ecost = MRISedgeAngleCostEdgeVertex(surf, eno, wrtvtxno, &gradv);
      // same cost each time, only count cost once. 
      if(wrtvtxno == 0) totcost += ecost; 
      if(!DoGrad) continue;
      // gradv is a row vector
      tx[e->vtxno[wrtvtxno]] += gradv->rptr[1][1];
      ty[e->vtxno[wrtvtxno]] += gradv->rptr[1][2];
      tz[e->vtxno[wrtvtxno]] += gradv->rptr[1][3];
      if(e->vtxno[wrtvtxno] == 195){
	//printf("#@& eno= %4d v=239 wrt=%d  %8.7f %8.7f %8.7f\n",eno,wrtvtxno,gradv->rptr[1][1],gradv->rptr[1][2],gradv->rptr[1][3]);
      }
    }
  }

  // normalize to the number of edges
  totcost /= nhits;
  if(weight > 0.0) totcost *= weight;

  if(DoGrad){
    // Stuff gradients a
    double a = 1.0/nhits;
    if(weight > 0.0) a = weight/nhits;
    for(vno = 0; vno < surf->nvertices; vno++){
      v = &(surf->vertices[vno]);
      if (v->ripflag)  continue;
      if(weight > 0.0){
	// This is for when called from the surface placement optimizer
	v->dx += -a*tx[vno];
	v->dy += -a*ty[vno];
	v->dz += -a*tz[vno];
      }
      else {
	v->tx = a*tx[vno];
	v->ty = a*ty[vno];
	v->tz = a*tz[vno];
      }
    }
    DMatrixFree(&gradv);
    free(tx);
    free(ty);
    free(tz);
  }

  return(totcost);
}

int MRISedgeAngleCostTest(MRIS *surf, double delta, double anglethresh, double magthresh)
{
  int vno, k, nerrs, vmpeak0, vmpeak1;
  long double cost0, cost, dcost[3], expdcost[3], d2sum, expd2sum;
  double dot, angle, angleavg, anglemax, expdcostnorm[3], dcostnorm[3];
  double x,y,z, meanlen;
  VERTEX *v0;
  int edgeno;
  double expmag, empmag, magerr, magerrmax;

  printf("Starting MRISedgeAngleCostTest()\n");

  // Set up the surface
  MRISedges(surf);
  MRIScomputeMetricProperties(surf);
  MRISedgeMetric(surf,1);
  MRISfaceNormalGrad(surf, 0);
  //MRISedgeGradDot(surf);

  meanlen = 0;
  for(edgeno = 0; edgeno < surf->nedges; edgeno++){
    meanlen += surf->edges[edgeno].len;
  }
  meanlen /= surf->nedges;

  printf("delta = %g, anglethresh = %g, magthresh = %g\n",delta,anglethresh,magthresh);
  printf("nvertices %d, nedges = %d, meanlen = %g\n",surf->nvertices,surf->nedges,meanlen);

  // Compute the cost with this configuration of vertices. Also
  // compute the gradient at each vertex.
  cost0 =  MRISedgeAngleCost(surf, -1.0, 1);
  printf("cost0 = %Lf\n",cost0);

  vmpeak0 = GetVmPeak();
  printf("#VMPC# VmPeak  %d\n",vmpeak0);

  // Go thru each vertex
  nerrs = 0;
  magerrmax = 0;
  anglemax = 0;
  angleavg = 0;
  for(vno = 0; vno < surf->nvertices; vno++){
    v0 = &(surf->vertices[vno]);
    // Expected (theoretical) gradient at this vertx
    expdcost[0] = v0->tx;
    expdcost[1] = v0->ty;
    expdcost[2] = v0->tz;
    // Get the original position so can restore
    x = v0->x;
    y = v0->y;
    z = v0->z;
    // perturb x, y, and z separately
    expd2sum = 0;
    d2sum = 0;
    for(k=0; k < 3; k++){
      if(k==0) v0->x += delta;
      if(k==1) v0->y += delta;
      if(k==2) v0->z += delta;
      // This step is slow. Could just do the edges connected to this vertex, but
      // probably should use the full function
      MRIScomputeMetricProperties(surf);
      MRISedgeMetric(surf,0);
      MRISfaceNormalGrad(surf, 0);
      //MRISedgeGradDot(surf);
      // In theory, this could just be done just for the connected edges, but, again, 
      // probably should use the full function
      cost = MRISedgeAngleCost(surf, -1.0, 0);
      // Compute the emperical gradient
      dcost[k] = (cost-cost0)/delta;
      // Keep track of the sum^2 for dot product calc below
      d2sum    += (dcost[k]*dcost[k]);
      expd2sum += (expdcost[k]*expdcost[k]);
      // restore for next cycle
      v0->x = x;
      v0->y = y;
      v0->z = z;
    } //k

    // Store emperical gradient in vertex t2{xyz}
    v0->t2x = dcost[0];
    v0->t2y = dcost[1];
    v0->t2z = dcost[2];

    // Compute the gradient magnitude
    expmag = sqrt(expd2sum); // expected  magnitude
    empmag = sqrt(d2sum);    // emperical magnitude
    magerr = fabs(expmag-empmag)/(expmag+FLT_MIN);
    if(magerrmax < magerr) magerrmax = magerr;

    // Compute the dot product between the theoretical and emperical gradients
    dot = 0;
    for(k=0; k < 3; k++){
      expdcostnorm[k] = expdcost[k]/expmag;
      dcostnorm[k]    = dcost[k]/empmag;
      dot += (expdcostnorm[k]*dcostnorm[k]);
    }
    if(dot > +1.0) dot = +1.0; // might be a slight overflow
    if(dot < -1.0) dot = -1.0; // might be a slight overflow
    // Now compute the angle (ideally, this would be 0)
    angle = 180*acos(dot)/M_PI;
    v0->val2bak = angle;
    angleavg += angle;
    if(anglemax < angle) anglemax = angle;
    if(angle>anglethresh || magerr > magthresh){
      nerrs++;
      printf("%5d %4d  %6.4f %6.4f   %8.8f %8.8f %8.8f ",nerrs,vno,angle,anglemax,expmag,empmag,magerr);
      for(k=0; k < 3; k++) printf("%14.9Lf ",expdcost[k]);
      printf("  ");
      for(k=0; k < 3; k++) printf("%14.9Lf ",dcost[k]);
      printf("\n");
      fflush(stdout);
    }
    //printf("d2sum = %20.19Lf, expd2sum = %20.19Lf\n",d2sum,expd2sum);

  } // loop over vertices
  angleavg /= surf->nvertices;

  vmpeak1 = GetVmPeak();
  printf("#VMPC# VmPeak  %d\n",vmpeak1);
  if(vmpeak1 > vmpeak0)  printf("WARNING: there might be a memory leak\n");
  else                   printf("No memory leak detected\n");

  printf("error rate %g, anglemax = %g, angleavg = %g, magerrmax = %g\n",
	 (double)nerrs/(3.0*surf->nvertices),anglemax,angleavg,magerrmax);
  fflush(stdout);
  return(nerrs);
}


/*!
  \fn long double MRISedgeLengthCostEdge(MRIS *surf, int edgeno, double L0, DMATRIX **pgradv0, DMATRIX **pgradv1)
  \brief Compute the cost and cost gradient of the given edge assuming
  its desired length is L0.  This is a spring cost with a
  (potentially) non-zero resting length of L0. If *pgradv{0,1} is
  non-NULL, then the gradient of the cost with respect to a change in
  vertex{0,1} position is computed. *pgradvX should be a 3x1 DMATRIX.
 */
long double MRIStargetCostVertex(const MRIS *surf, const int vno, long double *dc)
{
  long double cost, dx, dy, dz;
  VERTEX *v;

  v = &(surf->vertices[vno]);
  dx = v->targx - v->x;
  dy = v->targy - v->y;
  dz = v->targz - v->z;

  cost = dx*dx + dy*dy + dz*dz;
  if(dc==NULL) return(cost); // done

  // change in the cost wrt a change in vertex position
  dc[0] = -2*dx;
  dc[1] = -2*dy;
  dc[2] = -2*dz;

  return(cost);
}

/*!
  \fn long double MRIStargetCost(MRIS *surf, const double weight, int DoGrad)
  \brief Compute the total target vertex cost over all vertices. The
  cost for a single vertex is the the distance from xyz to targ{xyz}.
  If DoGrad==1, then the gradient of the total cost wrt each vertex
  position is computed and stored in v->t{xyz}. If weight>0, then the
  total cost and gradients are multiplied by weight and the gradients
  are stored in v->d{xyz}. The total cost and gradients are normalized
  by the number of edges. Uses OpenMP.
*/
#if GCC_VERSION > 80000
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
long double MRIStargetCost(MRIS *surf, const double weight, int DoGrad)
{
  int vno, nhits;
  long double totcost,vcost;
  double *tx,*ty,*tz;

  if(weight == 0.0) return(0);

  if(DoGrad){
    tx = (double*)calloc(sizeof(double),surf->nvertices);
    ty = (double*)calloc(sizeof(double),surf->nvertices);
    tz = (double*)calloc(sizeof(double),surf->nvertices);
  }

  // This particular cost function is easy to parallelize
  totcost = 0;
  nhits = 0;
  ROMP_PF_begin
  #ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(assume_reproducible) reduction(+ : nhits, totcost)
  #endif
  for(vno = 0; vno < surf->nvertices; vno++){
    ROMP_PFLB_begin
    VERTEX  *v;
    v = &(surf->vertices[vno]);
    if (v->ripflag)  continue;
    nhits++;
    long double dc[3];
    vcost = MRIStargetCostVertex(surf, vno, &dc[0]);
    totcost += vcost;
    if(!DoGrad) continue;
    // This is thread safe
    tx[vno] = dc[0];
    ty[vno] = dc[1];
    tz[vno] = dc[2];
    ROMP_PFLB_end
  }
  ROMP_PF_end
  // normalize to the number of vertices
  totcost /= nhits;
  if(weight > 0.0) totcost *= weight;

  if(! DoGrad) return(totcost);

  double a;
  if(weight > 0.0) a = weight/nhits;
  else           a = 1.0/nhits;

  ROMP_PF_begin
  #ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(assume_reproducible) 
  #endif
  for(vno = 0; vno < surf->nvertices; vno++){
    ROMP_PFLB_begin
    VERTEX *v;
    v = &(surf->vertices[vno]);
    if (v->ripflag)  continue;
    if(weight > 0.0){
      // This is for when called from the surface placement optimizer
      v->dx += -a*tx[vno];
      v->dy += -a*ty[vno];
      v->dz += -a*tz[vno];
    }
    else {
      v->tx = a*tx[vno];
      v->ty = a*ty[vno];
      v->tz = a*tz[vno];
    }
    ROMP_PFLB_end
  }
  ROMP_PF_end

  free(tx);
  free(ty);
  free(tz);

  return(totcost);
}
#if GCC_VERSION > 80000
#pragma GCC diagnostic pop
#endif


int MRIStargetCostTest(MRIS *surf, const double delta, const double anglethresh, const double magthresh)
{
  int vno, k, nerrs, vmpeak0, vmpeak1;
  long double cost0, cost, dcost[3], expdcost[3], d2sum, expd2sum;
  double dot, angle, angleavg, anglemax, expdcostnorm[3], dcostnorm[3];
  double x,y,z;
  VERTEX *v0;
  double expmag, empmag, magerr, magerrmax;

  printf("Starting MRIStargetCostTest()\n");

  // Set up the surface
  MRISedges(surf);
  MRIScomputeMetricProperties(surf);
  printf("delta = %g, anglethresh = %g, magthresh = %g\n",delta,anglethresh,magthresh);
  printf("nvertices %d, \n",surf->nvertices);

  // Compute the cost with this configuration of vertices. Also
  // compute the gradient at each vertex.
  Timer timer;
  cost0 =  MRIStargetCost(surf, -1.0, 1);
  printf("cost0 = %Lf\n",cost0);
  printf("cost computation time %9.6f sec\n",timer.seconds());

  vmpeak0 = GetVmPeak();
  printf("#VMPC# VmPeak  %d\n",vmpeak0);

  // Go thru each vertex
  nerrs = 0;
  magerrmax = 0;
  anglemax = 0;
  angleavg = 0;
  for(vno = 0; vno < surf->nvertices; vno++){
    v0 = &(surf->vertices[vno]);
    // Expected (theoretical) gradient at this vertx
    expdcost[0] = v0->tx;
    expdcost[1] = v0->ty;
    expdcost[2] = v0->tz;
    // Get the original position so can restore
    x = v0->x;
    y = v0->y;
    z = v0->z;
    // perturb x, y, and z separately
    expd2sum = 0;
    d2sum = 0;
    for(k=0; k < 3; k++){
      if(k==0) v0->x += delta;
      if(k==1) v0->y += delta;
      if(k==2) v0->z += delta;
      // This step is slow. Could just do the edges connected to this vertex, but
      // probably should use the full function
      MRIScomputeMetricProperties(surf);
      // In theory, this could just be done just for the connected edges, but, again, 
      // probably should use the full function
      cost = MRIStargetCost(surf, -1.0, 0);
      // Compute the emperical gradient
      dcost[k] = (cost-cost0)/delta;
      // Keep track of the sum^2 for dot product calc below
      d2sum    += (dcost[k]*dcost[k]);
      expd2sum += (expdcost[k]*expdcost[k]);
      // restore for next cycle
      v0->x = x;
      v0->y = y;
      v0->z = z;
    } //k

    // Store emperical gradient in vertex t2{xyz}
    v0->t2x = dcost[0];
    v0->t2y = dcost[1];
    v0->t2z = dcost[2];

    // Compute the gradient magnitude
    expmag = sqrt(expd2sum); // expected  magnitude
    empmag = sqrt(d2sum);    // emperical magnitude
    magerr = fabs(expmag-empmag)/(expmag+FLT_MIN);
    if(magerrmax < magerr) magerrmax = magerr;

    // Compute the dot product between the theoretical and emperical gradients
    dot = 0;
    for(k=0; k < 3; k++){
      expdcostnorm[k] = expdcost[k]/expmag;
      dcostnorm[k]    = dcost[k]/empmag;
      dot += (expdcostnorm[k]*dcostnorm[k]);
    }
    if(dot > 1.0)  dot = +1.0; // might be a slight overflow
    if(dot < -1.0) dot = -1.0; // might be a slight overflow
    // Now compute the angle (ideally, this would be 0)
    angle = 180*acos(dot)/M_PI;
    v0->val2bak = angle;
    angleavg += angle;
    if(anglemax < angle) anglemax = angle;
    if(angle>anglethresh || magerr > magthresh){
      nerrs++;
      printf("%5d %4d  %6.4f %6.4f   %8.8f %8.8f %8.8f ",nerrs,vno,angle,anglemax,expmag,empmag,magerr);
      for(k=0; k < 3; k++) printf("%14.9Lf ",expdcost[k]);
      printf("  ");
      for(k=0; k < 3; k++) printf("%14.9Lf ",dcost[k]);
      printf("\n");
      fflush(stdout);
    }
    //printf("d2sum = %20.19Lf, expd2sum = %20.19Lf\n",d2sum,expd2sum);

  } // loop over vertices
  angleavg /= surf->nvertices;

  vmpeak1 = GetVmPeak();
  printf("#VMPC# VmPeak  %d\n",vmpeak1);
  if(vmpeak1 > vmpeak0)  printf("WARNING: there might be a memory leak\n");
  else                   printf("No memory leak detected\n");

  printf("error rate %g, anglemax = %g, angleavg = %g, magerrmax = %g\n",
	 (double)nerrs/(3.0*surf->nvertices),anglemax,angleavg,magerrmax);
  printf("Test time %4.1f sec\n",timer.seconds());
  fflush(stdout);
  return(nerrs);
}

/*
  \fn double SurfGradUnitVector(const double v0[3], const double v1[3], double u[3], DMATRIX *grad)
  \brief Computes the unit vector from v0 to v1. Returns the distance
  between v0 and v1. If grad is non-NULL, then it computes the
  gradient of u with respect to v0 which equals that with respect to
  v1 (ie, you do not need to change the sign of grad).  grad should be
  a 3x3 matrix.  grad[m][k] = du[m]/dv0[k]. See also SurfGradUnitVectorTest().
 */
long double SurfGradUnitVector(const double v0[3], const double v1[3], double u[3], DMATRIX *grad)
{
  int k,m;
  long double w[3],L=0;

  // It is inefficient to compute length here. It should
  // be done as part of the edge metric computation
  for(k=0;k<3;k++) {
    w[k] = (long double)v1[k]-v0[k];
    L += (w[k]*w[k]);
  }

  L = sqrt(L);
  if(L==0){
    for(k=0;k<3;k++) u[k] = 0;
    if(grad){
      for(k=0;k<3;k++) {
        for(m=0;m<3;m++) {
          grad->rptr[k+1][m+1] = 0.0;
        }
      }
    }
    return(L);
  }
  for(k=0;k<3;k++) u[k] = w[k]/L;
  if(grad==NULL) return(L);

  // grad = (w*w')/(L^3) - I/L;
  // Note that w*w' is invariant to the sign of w, ie, if v0 and v1
  // are swapped w*w' does not change. This makes the gradient the
  // same if v0 and v1 are swapped. Therefore, do not apply a sign
  // flip to the gradient.
  long double iL3 = 1.0/(L*L*L);
  long double iL = 1.0/L;
  long double A;
  for(k=0;k<3;k++) {
    for(m=k;m<3;m++) {
      A = (w[k]*w[m])*iL3;
      if(k==m) A -= iL;
      // note: grad is double, not long double
      grad->rptr[k+1][m+1] = A;
      // Matrix is symmetric
      if(k!=m) grad->rptr[m+1][k+1] = A;
    }
  }
  return(L);
}

/*
  \fn int SurfGradUnitVectorTest(long double delta, double thresh, int ntrials)
  \brief Test for the gradient part of SurfGradUnitVector. The test is
  done by emperically computing the gradient by moving each dimension
  of v0 by delta then computing the change in the unit vector. This is
  repeated ntrials. If an error measure is greater than thresh, then
  an error is logged. Returns the number of errors. There may be some
  errors due to round-off. The error measure should be interpreted as
  a number between 0 and 1. Eg, this invocation currently has no errors:
  SurfGradUnitVectorTest(.0001, .01, 1000);
 */
int SurfGradUnitVectorTest(long double delta, double thresh, int ntrials)
{
  int k,m,n,nerrs;
  double v0[3], v1[3],  u[3], L, v0b[3], ub[3];
  DMATRIX *grad, *gradnum;
  double e,emax,g,gmax, err, errmax;

  grad = DMatrixAlloc(3,3,MATRIX_REAL);
  gradnum = DMatrixAlloc(3,3,MATRIX_REAL);

  nerrs = 0;
  errmax = 0;
  for(n=0; n<ntrials; n++){
    for(k=0;k<3;k++) {
      v0[k] = drand48();
      v1[k] = drand48();
    }
    
    L = SurfGradUnitVector(v0,v1, u, grad);
    //printf("L=%g\n",L);
 
    for(k=0;k<3;k++) {
      for(m=0;m<3;m++) v0b[m] = v0[m];
      v0b[k] = v0[k] + delta;
      SurfGradUnitVector(v0b,v1, ub, NULL);
      for(m=0; m<3; m++) gradnum->rptr[m+1][k+1] = (ub[m]-u[m])/delta;
    }

    emax = 0;
    gmax = 0;
    for(k=0;k<3;k++) {
      for(m=0; m<3; m++) {
	e = fabs(grad->rptr[m+1][k+1]-gradnum->rptr[m+1][k+1]);
	if(emax < e) emax = e;
	g = fabs(grad->rptr[m+1][k+1]);
	if(gmax < g) gmax = g;
      }
    }
    err = emax/gmax;
 
    if(err > thresh){
      printf("n = %d\n",n);
      printf("err = %g, thresh=%g, emax=%g, gmax=%g\n",err,thresh,emax,gmax);
      printf("v0 = ["); for(k=0;k<3;k++) printf(" %g ",v0[k]);printf("]';\n");
      printf("v1 = ["); for(k=0;k<3;k++) printf(" %g ",v1[k]);printf("]';\n");
      printf("grad = [\n");
      DMatrixPrint(stdout,grad);
      printf("];\n");
      printf("gradnum = [\n");
      DMatrixPrint(stdout,gradnum);
      printf("];\n");
      nerrs++;
    }
    if(errmax < err) errmax = err;
  }

  printf("nerrs = %d/%d, thresh=%g, delta=%Lf, errmax = %g\n",nerrs,ntrials,thresh,delta,errmax);

  return(nerrs);
}

/*!
  \fn double MRIScornerMetricCorner(MRIS *surf, const int cornerno, const int DoGrad)
  \brief Computes the metrics (dot, angle) of the corner. If requested, the gradient
  of the dot is also computed WRT each vertex in the corner.
*/
double MRIScornerMetricCorner(MRIS *surf, const int cornerno, const int DoGrad)
{
  int k;
  MRI_CORNER *c;
  MRI_EDGE *e1, *e2;

  c  = &(surf->corners[cornerno]);
  e1 = &(surf->edges[c->edgeno[0]]);
  e2 = &(surf->edges[c->edgeno[1]]);

  // Compute the dot product, apply appropriate direction sign
  c->dot = 0;
  for(k=0; k<3; k++) c->dot += (c->edgedir[0]*e1->u[k])*(c->edgedir[1]*e2->u[k]);
  // Compute the angle of the corner, might not be needed
  c->angle = 180*acos(c->dot)/M_PI; 

  // If not computing the gradient, then return now
  if(!DoGrad) return(c->dot);

  // Now compute the gradient of the dot WRT each vertex position

  // Put u1 and u2 into matrices, transpose along the way and apply the sign
  // Note that the edge direction sign does not need to be applyed to gradU
  DMATRIX *u1mat, *u2mat;
  u1mat = DMatrixAlloc(1,3,MATRIX_REAL);
  for(k=0; k<3; k++) u1mat->rptr[1][k+1] = (c->edgedir[0]*e1->u[k]);
  u2mat = DMatrixAlloc(1,3,MATRIX_REAL);
  for(k=0; k<3; k++) u2mat->rptr[1][k+1] = (c->edgedir[1]*e2->u[k]);

  //gradv1Dot = -u2'*gradv0u1; (do negation below)
  c->gradDot[1] = DMatrixMultiply(u2mat,e1->gradU,c->gradDot[1]);

  //gradv2Dot = -u1'*gradv0u2; (do negation below)
  c->gradDot[2] = DMatrixMultiply(u1mat,e2->gradU,c->gradDot[2]);

  //gradv0Dot = (u2'*gradv0u1 + u1'*gradv0u2);
  c->gradDot[0] = DMatrixAdd(c->gradDot[1],c->gradDot[2],c->gradDot[0]);

  // Now negate 1 and 2
  DMatrixScalarMul(c->gradDot[1], -1.0, c->gradDot[1]);
  DMatrixScalarMul(c->gradDot[2], -1.0, c->gradDot[2]);

  DMatrixFree(&u1mat);
  DMatrixFree(&u2mat);

  return(c->dot);
}


/*!
  \fn double *MRIScornerStats(MRIS *surf, int metricid, MRI *mask, double *stats)
  \brief Computes stats (ncorners, mean, stddev, min, max) over
  all corners that do not have a ripped vertex. If mask is non-null,
  then the mask value of a vertex must be greater than 0.5 to
  be included in the list. metricid: 0=dot, 1=angle.
  Will create the corner structure if not already there. Runs
  MRIScomputeMetricProperties() and MRIScornerMetric(surf).
 */
double *MRIScornerStats(MRIS *surf, int metricid, MRI *mask, double *stats)
{
  int cornerno, ncorners, nthv;
  MRI_CORNER *c;
  double *metric;

  if(surf->corners == NULL)  MRIScorners(surf);
  MRIScomputeMetricProperties(surf);
  MRIScornerMetric(surf,0);

  metric = (double*)calloc(sizeof(double),surf->ncorners);
  ncorners = 0;
  for(cornerno = 0; cornerno < surf->ncorners; cornerno++){
    c = &(surf->corners[cornerno]);
    int skip = 0;
    for(nthv=0; nthv < 3; nthv++){
      int vno = c->vtxno[nthv];
      VERTEX  * const v = &(surf->vertices[vno]);
      if(v->ripflag) skip = 1;
      if(mask && MRIgetVoxVal(mask,vno,0,0,0) < 0.5) skip = 1;
    }
    if(skip) continue;
    switch(metricid){
    case 0: metric[ncorners] = c->dot; break;
    case 1: metric[ncorners] = c->angle; break;
    default:
      printf("ERROR: MRIScornerStats() metricid %d unrecognized\n",metricid);
      return(NULL);
    }
    //printf("%lf\n",metric[ncorners]);
    ncorners++;
  }

  stats = DListStats(metric, ncorners, stats);
  free(metric);
  return(stats);
}
/*!
  \fn int MRIScornersPrint(FILE *fp, MRIS *surf)
  \brief Print corner topology and metrics to a stream. Runs
  MRIScomputeMetricProperties() and MRIScornerMetric()
 */
int MRIScornersPrint(FILE *fp, MRIS *surf)
{
  int cornerno;
  MRI_CORNER *c;
  double cost;
  MRIScomputeMetricProperties(surf);
  MRIScornerMetric(surf,0);
  for(cornerno = 0; cornerno < surf->ncorners; cornerno++){
    c = &(surf->corners[cornerno]);
    cost = (c->dot - 0.5)*(c->dot - 0.5); // cost for an equalateral triangle
    fprintf(fp,"%6d %6d %6d %6d %5d %10.8f %8.4f %11.9f\n",cornerno,
	    c->vtxno[0],c->vtxno[1],c->vtxno[2],c->faceno,c->dot,c->angle,cost);
  }
  fflush(fp);
  return(0);
}

/*!
  \fn int MRIScornersWrite(char *filename, MRIS *surf)
  \brief Print corner topology and metrics to a file. Runs
  MRIScomputeMetricProperties() and MRIScornersMetric()
  via MRIScornersPrint(). 
 */
int MRIScornersWrite(char *filename, MRIS *surf)
{
  FILE *fp;
  if(surf->corners == NULL) MRIScorners(surf);
  fp = fopen(filename,"w");
  if(fp == NULL) return(1);
  MRIScornersPrint(fp,surf);
  return(0);
}
int MRIScornerPrint(FILE *fp, const MRIS *surf, const int cornerno)
{
  double cost;
  int k;
  VERTEX *v;
  MRI_CORNER *c;
  MRI_EDGE *e;
  VERTEX_TOPOLOGY *vt;

  c = &(surf->corners[cornerno]);
  vt = &(surf->vertices_topology[c->vtxno[0]]);
  v = &(surf->vertices[c->vtxno[0]]);
  cost = (c->dot - 0.5)*(c->dot - 0.5); // cost for an equalateral triangle
  fprintf(fp,"%6d %6d %6d %6d %5d %5d %5d %10.8f %8.4f %11.9f\n",cornerno,
	  c->vtxno[0],c->vtxno[1],c->vtxno[2],c->faceno,c->edgeno[0],c->edgeno[1],c->dot,c->angle,cost);
  fprintf(fp,"Edge dir: %d %d\n",c->edgedir[0],c->edgedir[1]);
  for(k=0; k<2; k++){
    e = &(surf->edges[c->edgeno[k]]);
    fprintf(fp,"Edge %d  %6d %6d\n",k+1,e->vtxno[0],e->vtxno[1]);
  }
  for(k=0; k<vt->vtotal; k++){
    fprintf(fp,"Vnbr %d  %d  %d\n",k,vt->v[k],vt->e[k]);
  }
  for(k=0; k<3; k++){
    v = &(surf->vertices[c->vtxno[k]]);
    fprintf(fp,"v%d = [%g %g %g]';\n",k,v->x,v->y,v->z);
  }
  for(k=0; k<3; k++){
    fprintf(fp,"gradDot%d = [\n",k);
    DMatrixPrint(fp,c->gradDot[k]);
    fprintf(fp,"];\n");
  }
  fflush(fp);

  return(0);
}

/*!
  \fn int MRISedgeGradU(MRIS *surf)
  \brief Computes gradient of the unit vector pointing
  from v0 to v1 WRT v0.
 */
int MRISedgeGradU(MRIS *surf)
{
  int edgeno;

  for(edgeno=0; edgeno < surf->nedges; edgeno++){
    MRI_EDGE *e;
    double v0[3], v1[3];
    e = &(surf->edges[edgeno]);
    StuffVertexCoords(surf, e->vtxno[0], &v0[0]);
    StuffVertexCoords(surf, e->vtxno[1], &v1[0]);
    if(e->gradU == NULL) e->gradU = DMatrixAlloc(3,3,MATRIX_REAL);
    e->len = SurfGradUnitVector(v0, v1, e->u, e->gradU);
  }
  return(0);
}

/*!
  \fn long double MRIScornerDotCost(MRIS *surf, double weight, int DoGrad)
  \brief Computes the total corner dot cost of all non-ripped
  vertices. If DoGrad=1, then computes the gradient of the cost with
  respect to each vertex and saves in v->t{xyz}. If weight>0, then the
  cost and gradients are multipled by weight, and the negative of the
  weighted gradient is added to v->d{xyz} to make it compatible with
  the surface placement optimization code, and v->t{xyz} is not set.
  See MRIScornerDotCostTest() for how to test.
 */
long double MRIScornerDotCost(MRIS *surf, double weight, int DoGrad)
{
  int cno, wrtvtxno, vno, skip, nhits;
  long double totcost,ccost,A;
  MRI_CORNER *c;
  VERTEX *v;
  double *tx,*ty,*tz;

  if(weight == 0.0) return(0);

  if(DoGrad){
    tx = (double*)calloc(sizeof(double),surf->nvertices);
    ty = (double*)calloc(sizeof(double),surf->nvertices);
    tz = (double*)calloc(sizeof(double),surf->nvertices);
  }

  nhits = 0;
  totcost = 0;
  for(cno = 0; cno < surf->ncorners; cno++){
    c = &(surf->corners[cno]);

    // Determine whether to skip this corner
    skip = 0;
    for(wrtvtxno = 0; wrtvtxno < 3; wrtvtxno++){
      v = &(surf->vertices[c->vtxno[wrtvtxno]]);
      if (v->ripflag) {
	skip = 1;
	break;
      }
    }
    if(skip) continue;
    nhits ++;

    MRIScornerMetricCorner(surf, cno, DoGrad);
    A = (c->dot-0.5);
    ccost = A*A;
    totcost += ccost;

    if(DoGrad){
      // gradCcost = 2*(dot-dot0)*gradDot;
      int vtxno;
      for(wrtvtxno = 0; wrtvtxno < 3; wrtvtxno++){
	// gradDot is a row vector
	vtxno = c->vtxno[wrtvtxno];
	tx[vtxno] += (2*A*c->gradDot[wrtvtxno]->rptr[1][1]);
	ty[vtxno] += (2*A*c->gradDot[wrtvtxno]->rptr[1][2]);
	tz[vtxno] += (2*A*c->gradDot[wrtvtxno]->rptr[1][3]);
      }
    }
  }

  // normalize to the number of corners
  totcost /= nhits;
  if(weight > 0.0) totcost *= weight;

  if(DoGrad){
    // Stuff gradients intp proper surface field
    double a = 1.0/nhits;
    if(weight > 0.0) a = weight/nhits;
    for(vno = 0; vno < surf->nvertices; vno++){
      v = &(surf->vertices[vno]);
      if (v->ripflag)  continue;
      if(weight > 0.0){
	// This is for when called from the surface placement optimizer
	v->dx += -a*tx[vno];
	v->dy += -a*ty[vno];
	v->dz += -a*tz[vno];
      }
      else {
	v->tx = a*tx[vno];
	v->ty = a*ty[vno];
	v->tz = a*tz[vno];
      }
    }
    free(tx);
    free(ty);
    free(tz);
  }

  return(totcost);
}

/*!
  \fn int MRIScornerDotCostTest(MRIS *surf, const double delta, const double anglethresh, const double magthresh)
  \brief Test of gradients computed by MRIScornerDot(). Works by
  changing the x,y,z of a each vertex by delta and then measuring the
  change in the cost (ie, emperically computing gradient). There are
  two performance metrics: (1) the angle between the theoretical and
  empirical gradient directions, and (2) a comparison of the
  length/magnitude of the theoretical and empirical gradient vectors
  computed as a fractional change from the theoretical. A "failure"
  occurs if the angle>anglethresh or the relative magnitude difference
  > magthresh.  Eg, MRIScornerDotCostTest(surf, .001, 1, .05);
  Failures may occur due to round-off errors in the emperical
  calcluation.  In theory, this can be controlled by making delta
  arbitrarily small, but that can also create rouding errors. This is
  just testing the theoretical computation, so oen would expect it to
  fail at virtually every vertex if the theoretical calculation is
  off.  In other words, a few failures due to round-off error are
  nothing to worry about.
 */
int MRIScornerDotCostTest(MRIS *surf, const double delta, const double anglethresh, const double magthresh)
{
  int vno, k, nerrs, vmpeak0, vmpeak1;
  long double cost0, cost, dcost[3], expdcost[3], d2sum, expd2sum;
  double dot, angle, angleavg, anglemax, expdcostnorm[3], dcostnorm[3];
  double x,y,z;
  VERTEX *v0;
  double expmag, empmag, magerr, magerrmax;

  printf("Starting MRIScornerDotCostTest()\n");

  // Set up the surface
  MRISedges(surf);
  MRIScomputeMetricProperties(surf);
  MRISedgeMetric(surf,1); // gradU needed for corner
  //MRISedgeGradU(surf);
  MRIScornerMetric(surf, 1);
  printf("delta = %g, anglethresh = %g, magthresh = %g\n",delta,anglethresh,magthresh);
  printf("nvertices %d, nfaces %d, ncorners %d\n",surf->nvertices,surf->nfaces,surf->ncorners);

  // Time without gradient computation
  Timer timer0;
  cost0 =  MRIScornerDotCost(surf, -1.0, 0);
  printf("cost computation time without gradient %9.6f sec\n",timer0.seconds());
  printf("cost0 = %Lf\n",cost0);

  // Compute the cost with this configuration of vertices. Also
  // computes the gradient at each vertex.
  Timer timer;
  cost0 =  MRIScornerDotCost(surf, -1.0, 1);
  printf("cost computation time with gradient %9.6f sec\n",timer.seconds());
  printf("cost0 = %Lf\n",cost0);

  // print out the first corner info
  MRIScornerPrint(stdout, surf, 0);

  vmpeak0 = GetVmPeak(); // For checking for memory leaks

  // Go thru each vertex
  nerrs = 0;
  magerrmax = 0;
  anglemax = 0;
  angleavg = 0;
  for(vno = 0; vno < surf->nvertices; vno++){
    v0 = &(surf->vertices[vno]);
    // Expected (theoretical) gradient at this vertx computed above
    expdcost[0] = v0->tx;
    expdcost[1] = v0->ty;
    expdcost[2] = v0->tz;
    // Get the original position so can restore
    x = v0->x;
    y = v0->y;
    z = v0->z;
    // perturb x, y, and z separately
    expd2sum = 0;
    d2sum = 0;
    for(k=0; k < 3; k++){
      if(k==0) v0->x += delta;
      if(k==1) v0->y += delta;
      if(k==2) v0->z += delta;
      // This step is slow. Could just do the edges connected to this vertex, but
      // probably should use the full function to be sure
      MRIScomputeMetricProperties(surf);
      MRISedgeMetric(surf,0);
      MRIScornerMetric(surf, 0);
      cost = MRIScornerDotCost(surf, -1.0, 0);
      // Compute the empirical gradient
      dcost[k] = (cost-cost0)/delta;
      // Keep track of the sum^2 for dot product calc below
      d2sum    += (dcost[k]*dcost[k]);
      expd2sum += (expdcost[k]*expdcost[k]);
      // restore for next cycle
      v0->x = x;
      v0->y = y;
      v0->z = z;
    } //k

    // Store emperical gradient in vertex t2{xyz}
    v0->t2x = dcost[0];
    v0->t2y = dcost[1];
    v0->t2z = dcost[2];

    // Compute the gradient magnitude
    expmag = sqrt(expd2sum); // expected  magnitude
    empmag = sqrt(d2sum);    // emperical magnitude
    magerr = fabs(expmag-empmag)/(expmag+FLT_MIN);
    if(magerrmax < magerr) magerrmax = magerr;

    // Compute the dot product between the theoretical and emperical gradients
    dot = 0;
    for(k=0; k < 3; k++){
      expdcostnorm[k] = expdcost[k]/expmag;
      dcostnorm[k]    = dcost[k]/empmag;
      dot += (expdcostnorm[k]*dcostnorm[k]);
    }
    if(dot > +1.0) dot = +1.0; // might be a slight overflow
    if(dot < -1.0) dot = -1.0; // might be a slight overflow
    // Now compute the angle (ideally, this would be 0)
    angle = 180*acos(dot)/M_PI;
    v0->val2bak = angle;
    angleavg += angle;
    if(anglemax < angle) anglemax = angle;
    if(angle>anglethresh || magerr > magthresh){
      nerrs++;
      printf("%5d   %6d  %6.4f %6.4f   %8.8f %8.8f %8.8f ",nerrs,vno,angle,anglemax,expmag,empmag,magerr);
      for(k=0; k < 3; k++) printf("%14.9Lf ",expdcost[k]);
      printf("  ");
      for(k=0; k < 3; k++) printf("%14.9Lf ",dcost[k]);
      printf("\n");
      fflush(stdout);
    }
  } // loop over vertices
  angleavg /= surf->nvertices;

  vmpeak1 = GetVmPeak();
  printf("#VMPC# VmPeak  %d %d\n",vmpeak1,vmpeak0);
  if(vmpeak1 > vmpeak0)  printf("WARNING: there might be a memory leak (%d,%d)\n",vmpeak1,vmpeak0);
  else                   printf("No memory leak detected\n");

  printf("error rate %g, anglemax = %g, angleavg = %g, magerrmax = %g\n",
	 (double)nerrs/(3.0*surf->nvertices),anglemax,angleavg,magerrmax);
  printf("Test time %4.1f sec\n",timer.seconds());
  fflush(stdout);
  return(nerrs);
}

