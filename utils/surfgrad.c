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

#include "surfgrad.h"

int MRISfaceNormalFace_AddDeltaVertex = -1;
long double MRISfaceNormalFace_AddDelta[3]={0,0,0};

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
  int surfvtxno, wrtdimno, c;
  VERTEX *v;

  MRISfaceNormalGrad(surf, 0);
  MRISedgeGradDot(surf);
  g0 = DMatrixAlloc(1,3,MATRIX_REAL);
  J0 = MRISedgeAngleCostEdgeVertex(surf, edgeno, wrtvtxno, &g0);
  e = &(surf->edges[edgeno]);
  surfvtxno = e->vtxno[wrtvtxno];
  v = &(surf->vertices[surfvtxno]);
  //MRISedgePrint(surf, edgeno, stdout);
  //printf("%3d %g %g %g  %g\n",surfvtxno,v->x,v->y,v->z,J0);

  gnum = DMatrixAlloc(1,3,MATRIX_REAL);  
  for(wrtdimno=0; wrtdimno < 3; wrtdimno++){
    if(wrtdimno == 0) v->x += delta;
    if(wrtdimno == 1) v->y += delta;
    if(wrtdimno == 2) v->z += delta;
    MRISfaceNormalGrad(surf, 0);
    MRISedgeGradDot(surf);
    J1 = MRISedgeAngleCostEdgeVertex(surf, edgeno, wrtvtxno, NULL);
    //printf("  %g %g %g    %g\n",v->x,v->y,v->z,J1);
    gnum->rptr[1][wrtdimno+1] = (J1-J0)/delta;
    if(wrtdimno == 0) v->x -= delta;
    if(wrtdimno == 1) v->y -= delta;
    if(wrtdimno == 2) v->z -= delta;
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

  return(0);
}

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
