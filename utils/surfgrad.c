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

#if 0
/*!
  \fn int MRISfaceNormalGradTest(MRIS *surf, char *surfpath)

  This function tests MRISfaceNormalGrad() by numerically computing
  the gradient and comparing it to the theoretical. The idea is
  simple, but the application is more complicated for two reasons. (1)
  If a triangle is very small or has a very short edge, then the
  perturbation needs to be very small. (2) Very small perturbations
  can get overwhelmed by precision noise because the underlying FS and
  matrix code uses float. As a result, you can sometimes get fairly
  large differences between the theoretical and numerical. There are
  not too many mismatches, but there are enough and they are big
  enough that it makes it hard to set a good threshold.

  Returns the number of errors as indicated by an arbitrary threshold
  of triangles with areas above another threshold.

  If surf is passed, then it uses that surface as input. If surfpath
  is passed, then it loads that surface (and deallocs at the end).
  If both are NULL, then it uses ic1.tri.

 */
int MRISfaceNormalGradTest(MRIS *surf, char *surfpath, double delta)
{
  int faceno, wrtvtx, c,k, Dealloc, nhits;
  VERTEX *v;
  FACE *f;
  DMATRIX *norm1=NULL, *norm2[3]={NULL,NULL,NULL}, *dnorm=NULL, *numgrad=NULL;
  double maxabsnumgrad, absdiff, maxabsdiff, rdiff, maxrdiff;
  float vx0, vy0, vz0;
  char tmpstr[2000];

  Dealloc=1;
  if(surf) Dealloc=0;
  else {
    if(surfpath == NULL) {
      sprintf(tmpstr,"%s/lib/bem/ic1.tri",getenv("FREESURFER_HOME"));
      surfpath = tmpstr;
    }
    fprintf(stderr,"Reading in surf %s\n",surfpath);
    surf = MRISread(surfpath);
    if(surf==NULL) exit(1);
  }

  maxrdiff = -1;
  numgrad = DMatrixAlloc(3,3,MATRIX_REAL);
  nhits = 0;
  // go through each face
  for(faceno=0; faceno < surf->nfaces; faceno++){
    f = &(surf->faces[faceno]);
    // Compute the normal and the gradients for this face
    MRISfaceNormalGradFace(surf, faceno);
    norm1 = DMatrixCopy(f->norm,norm1); // make a copy of the norm

    // Empirically compute gradient WRT each vertex in the triangle
    DMatrixConstVal(3, 3, 0, numgrad);
    for(wrtvtx = 0; wrtvtx < 3; wrtvtx++){
      v = &(surf->vertices[f->v[wrtvtx]]);
      vx0 = v->x; vy0 = v->y;  vz0 = v->z; // store orig pos

      // Go thru each WRT vertex coordinate
      for(c=0; c<3; c++){
	// Change coord by a small amount
	if(c==0) v->x += delta;
	if(c==1) v->y += delta;
	if(c==2) v->z += delta;
	// Compute only the normal
	MRISfaceNormalGrad(surf, faceno, 1);
	norm2[c] = DMatrixCopy(f->norm,norm2[c]); // make a copy
	dnorm = DMatrixSubtract(norm2[c],norm1,dnorm); 
	// Numerically compute the gradient
	for(k=0; k<3; k++) 
	  numgrad->rptr[k+1][c+1] = (dnorm->rptr[k+1][1]/delta);
	v->x = vx0; v->y = vy0;	v->z = vz0; //restore
      } // coordinate
      // Compute the maximum abs gradient and maximum abs difference
      maxabsnumgrad = -1;
      maxabsdiff = -1;
      for(k=0; k<3; k++) {
	for(c=0; c<3; c++){
	  if(maxabsnumgrad < fabs(numgrad->rptr[k+1][c+1])) maxabsnumgrad = fabs(numgrad->rptr[k+1][c+1]);
	  absdiff = fabs(numgrad->rptr[k+1][c+1] - f->gradNorm[wrtvtx]->rptr[k+1][c+1]);
	  if(maxabsdiff < absdiff) maxabsdiff = absdiff;
	}
      }

      // Compute the difference relative to the maximum
      rdiff = maxabsdiff/maxabsnumgrad;
      if(maxrdiff < rdiff && f->area > 0.1) maxrdiff = rdiff;

      // Log an error if rdiff is over threshold and the face area is big enough
      // This is not perfect as the area may be ok but one edge may be really short
      if(rdiff > 0.02 && f->area > 0.1){
	printf("#@# %4d %4d %d %6.4lf  %6.4lf  %6.4lf %10.8f\n",nhits,faceno,wrtvtx,maxabsdiff,maxabsnumgrad,rdiff,f->area);
	fflush(stdout);
	nhits ++;
	if(0){
	  printf("delta = %lf\n",delta);
	  for(k = 0; k < 3; k++){
	    v = &(surf->vertices[f->v[k]]);
	    printf("vtx%d = [%g,%g,%g]\n",k,v->x,v->y,v->z);
	  }
	  printf("norm1 = [\n");
	  DMatrixPrintFmt(stdout,"8.5lf",norm1);
	  printf("]\n");
	  for(c=0; c<3; c++){
	    printf("norm2 = [\n");
	    DMatrixPrintFmt(stdout,"8.5lf",norm2[c]);
	    printf("]\n");
	  }
	  printf("dnorm = [\n");
	  DMatrixPrintFmt(stdout,"8.5lf",dnorm);
	  printf("]\n");
	  printf("gradNorm = [\n");
	  DMatrixPrintFmt(stdout,"8.5lf",f->gradNorm[wrtvtx]);
	  printf("]\n");
	  printf("numgrad = [\n");
	  DMatrixPrintFmt(stdout,"8.5lf",numgrad);
	  printf("]\n");
	  printf("mar(gradNorm-numgrad)\n");
	  //exit(1);
	}
      }
    } // corner vertex
  } // face
  printf("%4d %4d %lf\n",surf->nfaces,nhits,maxrdiff);

  if(Dealloc) MRISfree(&surf);

  DMatrixFree(&norm1);
  for(k = 0; k < 3; k++) DMatrixFree(&norm2[k]);
  DMatrixFree(&dnorm);
  DMatrixFree(&numgrad);

  return(nhits);
}
#endif

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

