/**
 * @file  surfgrad.c
 * @brief Utilities to compute gradients on the surface
 *
 */
/*
 * Original Author: Doug Greve
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2016/07/08 21:14:24 $
 *    $Revision: 1.153 $
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
Computes the gradient of the face normal with respect to changes 
in the the position of each vertex (also computes the norm itself).
The gradient is a 3x3 matrix of the form
   dnx/dvx dnx/dvy dnx/dvz  
   dny/dvx dny/dvy dny/dvz  
   dnz/dvx dnz/dvy dnz/dvz  
Where nx is the x component of the normal and vx is the x component of the vertex.
Since there are three vertices, there are three such matrices.
The matrices are stored in face->gradNorm[refvtxno]
The norm is stored in face->norm
If NormOnly != 0, then only the norm is computed
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
  MRISedgeGradDot(surf);
  
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
    MRISedgeGradDot(surf);
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
  int faceno;

  #ifdef HAVE_OPENMP
  #pragma omp parallel for 
  #endif
  for(faceno=0; faceno < surf->nfaces; faceno++){
    if(NormOnly)
      MRISfaceNormalFace(surf, faceno, NULL, NULL);
    else
      MRISfaceNormalGradFace(surf, faceno);
  }
  return(0);
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
  MRISedgeMetricEdge(surf, edgeno);
  dot0 = e->dot;
  MRISedgeGradDotEdgeVertex(surf, edgeno, wrtvtxno);
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

    MRISedgeMetricEdge(surf, edgeno);
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
  \fn int MRISedgeMetricEdge(MRIS *surf, int edgeno)
  \brief Computes edge-related metrics including the length, the dot
  product of the angle of the faces across the edge, and the angle
  (def).
 */
int MRISedgeMetricEdge(MRIS *surf, int edgeno)
{
  MRI_EDGE *e;
  FACE *f0, *f1;
  VERTEX *v0, *v1;
  int nthface,faceno;
  e = &(surf->edges[edgeno]);
  v0 = &(surf->vertices[e->vtxno[0]]);
  v1 = &(surf->vertices[e->vtxno[1]]);
  e->len = sqrt((v0->x-v1->x)*(v0->x-v1->x) + (v0->y-v1->y)*(v0->y-v1->y) + (v0->z-v1->z)*(v0->z-v1->z));
  for(nthface = 0; nthface < 2; nthface++){
    faceno = e->faceno[nthface];
    f0 = &(surf->faces[faceno]);
    if(!f0->norm) MRISfaceNormalFace(surf, faceno, NULL, NULL);
  }
  f0 = &(surf->faces[e->faceno[0]]);
  f1 = &(surf->faces[e->faceno[1]]);
  e->dot = DVectorDot(f0->norm,f1->norm);
  e->angle = acos(e->dot)*180/M_PI;
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

  // These must have been already run
  //MRISfaceNormalGrad(surf, 0);
  //MRISedgeGradDot(surf);

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
