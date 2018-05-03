/**
 * @file  surfgrad.c
 * @brief Utilities to compute gradients on the surface
 *
 */
/*
 * Original Author: Doug Greve
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2015/03/31 20:27:34 $
 *    $Revision: 1.90 $
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


#ifndef SURFGRAD_H
#define SURFGRAD_H

#if defined(__cplusplus)
extern "C" {
#endif

#ifdef X
#undef X
#endif

#include "mrisurf.h"

#ifdef _SURFGRAD_SRC
int MRISfaceNormalFace_AddDeltaVertex = -1;
long double MRISfaceNormalFace_AddDelta[3]={0,0,0};
#else
extern int MRISfaceNormalFace_AddDeltaVertex;
extern long double MRISfaceNormalFace_AddDelta[3];
#endif

int MRISfaceNormalGradFace(MRIS *surf, int faceno);
int MRISfaceNormalGradTest(MRIS *surf, char *surfpath, double delta);
int MRISfaceNormalFace(MRIS *surf, int faceno, DMATRIX **pc, double *pcL);
double MRISfaceNormalGradFaceTest(MRIS *surf, int faceno, int wrtvtxno, long double delta, int verbose);

  double MRISedgeAngleCostEdgeVertex(MRIS *surf, int edgeno, int wrtvtxno, DMATRIX **pgrad);
  int MRISedgeGradDotEdgeVertex(MRIS *surf, int edgeno, int wrtvtxno);
  int MRISedgeGradDot(MRIS *surf);
  double MRISedgeGradDotEdgeVertexTest(MRIS *surf, int edgeno, int wrtvtxno, long double delta, int verbose);
  int MRISfaceNormalGrad(MRIS *surf, int NormOnly);
  int MRISedgeAngleCostEdgeVertexTest(MRIS *surf, int edgeno, int wrtvtxno, long double delta);
  int MRISedgePrint(MRIS *surf, int edgeno, FILE *fp);
  int MRISedgeMetricEdge(MRIS *surf, int edgeno);
  int MRISfacePrint(MRIS *surf, int faceno, FILE *fp);

#if defined(__cplusplus)
};
#endif


#endif
