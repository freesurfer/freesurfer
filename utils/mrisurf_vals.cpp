#define COMPILING_MRISURF_TOPOLOGY_FRIEND_CHECKED
/*
 *
 */
/*
 * surfaces Author: Bruce Fischl, extracted from mrisurf.c by Bevin Brett
 *
 * $ Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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
#include "mrisurf_vals.h"

#include "mrisurf_base.h"


// Vals are scalar properties of vertexs or faces
// that are independent of 
//
// vertices
//              marked  marked2     annotations     flags   marked3     val     mean    cropped
//  count       x
//  set all     x
//  set         x                                   x
//  clear       x       x           x               x
//  invert      x
//  negate                                                              x
//  dilate      x
//  erode       x
//      two-operand
//          copy    
//              marked2=marked      marked3=marked
//              marked=marked2      marked=marked3
//              anotation=val
//              val=mean            val=mean_imag       val=std_error
//
// faces
//              marks
//  clear       x
//

//=============================================================================
// vertex
//
/*!
  \fn int MRIScountAllMarked(MRIS *mris)
  \brief Returns the total number of vertices have v->marked > 0
 */
int MRIScountAllMarked(MRIS *mris)
{
  int nmarked=0, vno;
  for (vno = 0 ; vno < mris->nvertices ; vno++){
    if(mris->vertices[vno].marked>0) nmarked++;
  }
  return(nmarked);
}
/*!
  \fn int MRIScountMarked(MRIS *mris, int mark_threshold)
  \brief Returns the total number of non-ripped vertices 
  that have v->marked >= threshold
 */
int MRIScountMarked(MRIS *mris, int mark_threshold)
{
  int vno, total_marked;
  VERTEX *v;

  for (total_marked = vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (v->marked >= mark_threshold) {
      total_marked++;
    }
  }
  return (total_marked);
}


int MRISmarkedVertices(MRIS *mris)
{
  int vno, nvertices, nmarked;

  nvertices = mris->nvertices;
  for (vno = nmarked = 0; vno < nvertices; vno++)
    if (!mris->vertices[vno].ripflag && mris->vertices[vno].marked > 0) {
      nmarked++;
    }

  return (nmarked);
}


// set all the marks to a user specified value INCLUDING RIPPED VERTICES!
int MRISsetAllMarks(MRIS *mris, int mark)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    v->marked = mark;
  }
  return (NO_ERROR);
}


int MRISsetMarks(MRIS *mris, int mark)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->marked = mark;
  }
  return (NO_ERROR);
}

int MRISclearMarks(MRIS *mris)
{
  return MRISsetMarks(mris, 0);
}

int MRISinvertMarks(MRIS *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->marked = !v->marked;
  }
  return (NO_ERROR);
}


int MRISnotMarked(MRIS *mris)
{
  return MRISinvertMarks(mris);
}


int MRISthresholdValIntoMarked(MRI_SURFACE *mris, float thresh)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->marked = v->val >= thresh;
  }

  return (NO_ERROR);
}


/* assume that the mark is 1 */
int MRISexpandMarked(MRIS *mris)
{
  int vno, n;

  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->marked == 1) {
      for (n = 0; n < vt->vnum; n++) {
        VERTEX * const vn = &mris->vertices[vt->v[n]];
        if (vn->marked == 0) {
          vn->marked = 2;
        }
      }
    }
  }
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX * const v = &mris->vertices[vno];
    if (v->marked == 2) {
      v->marked = 1;
    }
  }
  return (NO_ERROR);
}

/*! -----------------------------------------------------
  \fn int MRISdilateMarked(MRIS *mris, int ndil)
  \brief Dilates the marked vertices by marking a vertex
  if any of its non-ripped neighbors is ripped.
  ------------------------------------------------------*/
int MRISdilateMarked(MRIS *mris, int ndil)
{
  int vno, i, n, mx;

  // Loop through each dilation
  for (i = 0; i < ndil; i++) {

    // Set v->tx to 0 for unripped vertices
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      if(v->ripflag) continue;
      v->tx = 0;
    }

    // Loop through vertices (skip ripped)
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if(v->ripflag) continue;
      // set v->tx=1 if this vertex or any of its neightbors is marked
      mx = v->marked;
      for (n = 0; n < vt->vnum; n++) {
        VERTEX const * const vn = &mris->vertices[vt->v[n]];
        mx = MAX(vn->marked, mx);
      }
      v->tx = mx;
    }

    // Now copy tx into marked
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      if (v->ripflag) continue;
      v->marked = (int)v->tx;
    }

  }// end loop over dilations
  return (NO_ERROR);
}

int MRISerodeMarked(MRIS *mris, int num)
{
  int vno, i, n, mn;

  for (i = 0; i < num; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }
      v->tx = 0;
    }

    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) {
        continue;
      }
      mn = v->marked;
      for (n = 0; n < vt->vnum; n++) {
        VERTEX * const vn = &mris->vertices[vt->v[n]];
        mn = MIN(vn->marked, mn);
      }
      v->tx = mn;
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }
      v->marked = (int)v->tx;
    }
  }
  return (NO_ERROR);
}


int MRIScloseMarked(MRIS *mris, int order)
{
  MRISdilateMarked(mris, order);
  MRISerodeMarked(mris, order);
  return (NO_ERROR);
}


int MRISopenMarked(MRIS *mris, int order)
{
  MRISerodeMarked(mris, order);
  MRISdilateMarked(mris, order);
  return (NO_ERROR);
}



// marks2
//
int MRISclearMark2s(MRIS *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->marked2 = 0;
  }
  return (NO_ERROR);
}


// annotations
//
int MRISclearAnnotations(MRIS *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->annotation = 0;
  }
  return (NO_ERROR);
}


int MRISsetCroppedToZero(MRI_SURFACE *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) continue;
    v->cropped = 0;
  }
  return (NO_ERROR);
}


// flags
//
int MRISsetFlags(MRIS *mris, int flags)
{
  int vno;

  for (vno = 0; vno < mris->nvertices; vno++) {
    mris->vertices[vno].flags |= flags;
  }
  return (NO_ERROR);
}


int MRISclearFlags(MRIS *mris, int flags)
{
  int vno;

  for (vno = 0; vno < mris->nvertices; vno++) {
    mris->vertices[vno].flags &= (~flags);
  }
  return (NO_ERROR);
}


// val
//

void mrisSetVal(MRIS *mris, float val)
{
  int n;
  for (n = 0; n < mris->nvertices; n++) 
    mris->vertices[n].val = val;
}

void mrisSetValAndClearVal2(MRIS *mris, float val) {
  int n;
  for (n = 0; n < mris->nvertices; n++) { 
    mris->vertices[n].val  = val;
    mris->vertices[n].val2 = 0.0f;
  }
}

int mrisClearGradient(MRI_SURFACE *mris)
{
  int vno, nvertices;
  VERTEX *v;

  nvertices = mris->nvertices;
  for (vno = 0; vno < nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) continue;
    v->dx = 0;
    v->dy = 0;
    v->dz = 0;
  }
  return (NO_ERROR);
}


#if 0
int mrisClearExtraGradient(MRI_SURFACE *mris)
{
  int     vno, nvertices ;
  VERTEX  *v ;

  nvertices = mris->nvertices ;
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (mris->dx2)
      mris->dx2[vno] = mris->dy2[vno] = mris->dz2[vno] = 0 ;
  }
  return(NO_ERROR) ;
}
#endif



int MRISaddToValues(MRI_SURFACE *mris, float val)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->val += val;
  }
  return (NO_ERROR);
}


int MRISnegateValues(MRI_SURFACE *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->val *= -1.0;
  }
  return (NO_ERROR);
}


int MRISsetVals(MRI_SURFACE *mris, float val)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->val = val;
  }
  return (NO_ERROR);
}


int MRISscaleVals(MRI_SURFACE *mris, float scale)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->val *= scale;
  }
  return (NO_ERROR);
}


int MRISmulVal(MRI_SURFACE *mris, float mul)
{
  int vno, nvertices;
  VERTEX *v;

  nvertices = mris->nvertices;
  for (vno = 0; vno < nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->val *= mul;
  }
  return (NO_ERROR);
}


int MRISsqrtVal(MRI_SURFACE *mris)
{
  int vno, nvertices;
  VERTEX *v;

  nvertices = mris->nvertices;
  for (vno = 0; vno < nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (v->val > 0) {
      v->val = sqrt(v->val);
    }
  }
  return (NO_ERROR);
}


int MRISlogOdds(MRI_SURFACE *mris, LABEL *area, double slope)
{
  int vno;
  VERTEX *v;
  double p;

  MRISdistanceTransform(mris, area, DTRANS_MODE_SIGNED);
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    p = v->val;
    if (p < 0) {
      p = 0;
    }
    p = exp(-p * slope);
    v->val = p;
  }
  return (NO_ERROR);
}



int MRISsetVal2(MRI_SURFACE *mris, float val)
{
  VERTEX *v;
  int vno;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->val2 = val;
  }
  return (NO_ERROR);
}


int MRISsetValues(MRI_SURFACE *mris, float val)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->val = val;
  }
  return (NO_ERROR);
}


int MRISclearFixedValFlags(MRI_SURFACE *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->fixedval = FALSE;
  }
  return (NO_ERROR);
}


float* MRISexportCurv(MRIS* mris) {
  float* p = (float*)malloc(mris->nvertices * sizeof(float));
  MRISextractCurvatureVector(mris, p);
  return p;
}

int MRISextractCurvatureVector(MRIS *mris, float *curvs)
{
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    curvs[vno] = mris->vertices[vno].curv;
  }

  return (NO_ERROR);
}


int MRISextractCurvatureVectorDouble(MRI_SURFACE *mris, double *curvs, int offset)
{
  int vno;

  for (vno = 0; vno < mris->nvertices; vno++) {
    curvs[vno + offset] = (double)mris->vertices[vno].curv;
  }

  return (NO_ERROR);
}


int MRISextractCurvatureDoubleVector(MRI_SURFACE *mris, double *curvs)
{
  int vno;

  for (vno = 0; vno < mris->nvertices; vno++) {
    curvs[vno] = (double)mris->vertices[vno].curv;
  }

  return (NO_ERROR);
}


int MRISimportCurvatureVector(MRI_SURFACE *mris, float *curvs)
{
  int vno;

  for (vno = 0; vno < mris->nvertices; vno++) {
    mris->vertices[vno].curv = curvs[vno];
  }

  return (NO_ERROR);
}


int MRISimportValVector(MRI_SURFACE *mris, float *vals)
{
  int vno;

  for (vno = 0; vno < mris->nvertices; vno++) {
    mris->vertices[vno].val = vals[vno];
  }

  return (NO_ERROR);
}


int MRISimportValFromMatrixColumn(MRI_SURFACE *mris, MATRIX *m, int col)
{
  int vno;

  for (vno = 0; vno < mris->nvertices; vno++) {
    mris->vertices[vno].val = *MATRIX_RELT(m, vno + 1, col);
  }

  return (NO_ERROR);
}


int MRISexportValVector(MRI_SURFACE *mris, float *vals)
{
  int vno;

  for (vno = 0; vno < mris->nvertices; vno++) {
    vals[vno] = mris->vertices[vno].val;
  }

  return (NO_ERROR);
}


int MRISexportValVectorDouble(MRI_SURFACE *mris, double *vals, int offset)
{
  int vno;

  for (vno = 0; vno < mris->nvertices; vno++) {
    vals[vno + offset] = (double)mris->vertices[vno].val;
  }

  return (NO_ERROR);
}


// grad
//
int MRISclearGradient(MRI_SURFACE *mris)
{
  int vno, nvertices;
  VERTEX *v;

  nvertices = mris->nvertices;
  for (vno = 0; vno < nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->dx = v->dy = v->dz = 0;
  }
  return (NO_ERROR);
}


int mrisClearExtraGradient(MRI_SURFACE *mris)
{
  int vno, nvertices;
  VERTEX *v;

  nvertices = mris->nvertices;
  for (vno = 0; vno < nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (mris->dx2) {
      mris->dx2[vno] = mris->dy2[vno] = mris->dz2[vno] = 0;
    }
  }
  return (NO_ERROR);
}


int mrisClearMomentum(MRI_SURFACE *mris)
{
  int vno, nvertices;
  VERTEX *v;

  nvertices = mris->nvertices;
  for (vno = 0; vno < nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->odx = v->ody = v->odz = 0;
  }
  return (NO_ERROR);
}

int MRISabsVals(MRI_SURFACE *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->val = fabs(v->val);
  }
  return (NO_ERROR);
}

int MRISabsCurvature(MRI_SURFACE *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->curv = fabs(v->curv);
  }
  return (NO_ERROR);
}


// Two-operand
//
int MRIScopyMarkedToMarked2(MRIS *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    v->marked2 = v->marked;
  }
  return (NO_ERROR);
}


int MRIScopyMarked2ToMarked(MRIS *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    v->marked = v->marked2;
  }
  return (NO_ERROR);
}


int MRIScopyMarkedToMarked3(MRIS *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    v->marked3 = v->marked;
  }
  return (NO_ERROR);
}


int MRIScopyMarked3ToMarked(MRIS *mris)
{
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX * const v = &mris->vertices[vno];
    v->marked = v->marked3;
  }
  return (NO_ERROR);
}


int MRIScopyValsToAnnotations(MRIS *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    v->annotation = v->val;
  }
  return (NO_ERROR);
}


int MRIScopyMeansToValues(MRI_SURFACE *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->val = v->mean;
  }
  return (NO_ERROR);
}


int MRIScopyImaginaryMeansToValues(MRI_SURFACE *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->val = v->mean_imag;
  }
  return (NO_ERROR);
}


int MRIScopyStandardErrorsToValues(MRI_SURFACE *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->val = v->std_error;
  }
  return (NO_ERROR);
}

int MRIScopyStatsToValues(MRI_SURFACE *mris)
{
  int vno, nvertices;
  VERTEX *v;

  nvertices = mris->nvertices;
  for (vno = 0; vno < nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->val = v->stat;
  }
  return (NO_ERROR);
}


int MRIScopyStatsFromValues(MRI_SURFACE *mris)
{
  int vno, nvertices;
  VERTEX *v;

  nvertices = mris->nvertices;
  for (vno = 0; vno < nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->stat = v->val ;
  }
  return (NO_ERROR);
}

int MRIScopyValuesToImagValues(MRI_SURFACE *mris)
{
  int vno, nvertices;
  VERTEX *v;

  nvertices = mris->nvertices;
  for (vno = 0; vno < nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->imag_val = v->val;
  }
  return (NO_ERROR);
}


int MRIScopyImagValuesToValues(MRI_SURFACE *mris)
{
  int vno, nvertices;
  VERTEX *v;

  nvertices = mris->nvertices;
  for (vno = 0; vno < nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->val = v->imag_val;
  }
  return (NO_ERROR);
}

int MRIScopyVal2ToVal(MRI_SURFACE *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->val = v->val2;
  }
  return (NO_ERROR);
}


int MRIScopyCurvatureToValues(MRI_SURFACE *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->val = v->curv;
  }
  return (NO_ERROR);
}


int MRIScopyVal2BakToVal(MRI_SURFACE *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->val = v->val2bak;
  }
  return (NO_ERROR);
}


int MRIScopyCurvatureToImagValues(MRI_SURFACE *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->imag_val = v->curv;
  }
  return (NO_ERROR);
}


int MRIScopyValuesToCurvature(MRI_SURFACE *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->curv = v->val;
  }
  return (NO_ERROR);
}


int MRIScopyImagValuesToCurvature(MRI_SURFACE *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->curv = v->imag_val;
  }
  return (NO_ERROR);
}


int MRIScopyCurvatureFromValues(MRI_SURFACE *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->curv = v->val;
  }
  return (NO_ERROR);
}


int MRIScopyCurvatureFromImagValues(MRI_SURFACE *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->curv = v->imag_val;
  }
  return (NO_ERROR);
}


int MRIScopyValToVal2(MRI_SURFACE *mris)
{
  int vno, nvertices;
  VERTEX *v;

  nvertices = mris->nvertices;
  for (vno = 0; vno < nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->val2 = v->val;
  }
  return (NO_ERROR);
}


int MRIScopyValuesToVal2Bak(MRI_SURFACE *mris)
{
  int vno, nvertices;
  VERTEX *v;

  nvertices = mris->nvertices;
  for (vno = 0; vno < nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->val2bak = v->val;
  }
  return (NO_ERROR);
}


int MRIScopyValToValBak(MRI_SURFACE *mris)
{
  int vno, nvertices;
  VERTEX *v;

  nvertices = mris->nvertices;
  for (vno = 0; vno < nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->valbak = v->val;
  }
  return (NO_ERROR);
}


int MRIScopyValToVal2Bak(MRI_SURFACE *mris)
{
  int vno, nvertices;
  VERTEX *v;

  nvertices = mris->nvertices;
  for (vno = 0; vno < nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->val2bak = v->val;
  }
  return (NO_ERROR);
}


int MRIScopyFromCropped(MRI_SURFACE *mris, int which)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) continue;

    switch (which) {
      case VERTEX_VALS:
        v->val = v->cropped;
        break;
      case VERTEX_MARKS:
        v->marked = v->cropped;
        break;
      default:
        ErrorExit(ERROR_UNSUPPORTED, "MRIScopyFromCropped: %d unsupported target", which);
    }
  }
  return (NO_ERROR);
}


int MRIScopyToCropped(MRI_SURFACE *mris, int which)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) continue;

    switch (which) {
      case VERTEX_VALS:
        v->cropped = v->val;
        break;
      case VERTEX_MARKS:
        v->cropped = v->marked;
        break;
      default:
        ErrorExit(ERROR_UNSUPPORTED, "MRIScopyFromCropped: %d unsupported target", which);
    }
  }
  return (NO_ERROR);
}


int MRISorigAreaToCurv(MRI_SURFACE *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->curv = v->origarea;
  }
  return (NO_ERROR);
}

int MRISareaToCurv(MRI_SURFACE *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->curv = v->area;
  }
  return (NO_ERROR);
}

int MRISmarkedToCurv(MRI_SURFACE *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->curv = (float)v->marked;
  }
  return (NO_ERROR);
}

int MRIScurvToMarked(MRI_SURFACE *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->marked = (int)v->curv;
  }
  return (NO_ERROR);
}

int MRISdToCurv(MRI_SURFACE *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->curv = (float)v->d;
  }
  return (NO_ERROR);
}

int MRIScurvToD(MRI_SURFACE *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->d = v->curv;
  }
  return (NO_ERROR);
}

int MRIScopyFixedValFlagsToMarks(MRI_SURFACE *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->marked = v->fixedval;
  }
  return (NO_ERROR);
}


//======
// faces
//
int MRISclearFaceMarks(MRIS *mris)
{
  int fno;
  FACE *f;

  for (fno = 0; fno < mris->nfaces; fno++) {
    f = &mris->faces[fno];
    f->marked = 0;
  }
  return (NO_ERROR);
}


//======
// ripped
//
int MRISripLabel(MRIS *mris, LABEL *area)
{
  int i;
  VERTEX *v;

  for (i = 0; i < area->n_points; i++) {
    v = &mris->vertices[area->lv[i].vno];
    v->ripflag = 1;
  }
  return (NO_ERROR);
}


int MRISripNotLabel(MRIS *mris, LABEL *area)
{
  int i, vno;
  VERTEX *v;

  for (i = 0; i < area->n_points; i++) {
    v = &mris->vertices[area->lv[i].vno];
    v->marked = 1;
  }
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->marked) {
      continue;
    }
    v->ripflag = 1;
  }
  for (i = 0; i < area->n_points; i++) {
    v = &mris->vertices[area->lv[i].vno];
    v->marked = 0;
  }
  return (NO_ERROR);
}

/*!
  \fn int MRISripMarked(MRIS *mris)
  \brief Sets v->ripflag=1 if v->marked==1.
  Note: does not unrip any vertices.
*/
int MRISripMarked(MRIS *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->marked) {
      v->ripflag = 1;
    }
  }
  return (NO_ERROR);
}


int MRISripUnmarked(MRIS *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->marked == 0) {
      v->ripflag = 1;
    }
  }
  return (NO_ERROR);
}


/*-------------------------------------------------------*/
/*!
  \fn int MRISripZeros(MRIS *surf, MRI *mri)
  \brief Sets ripflag=1 for vertices where the mri value is 0
    (actually less than 1e-5). If mri is null, then uses the
    val field. No change to a vertex if ripflag already = 1.
*/
int MRISripZeros(MRIS *surf, MRI *mri)
{
  int k;
  double v;

  if (mri) {
    if (mri->width != surf->nvertices) {
      printf("ERROR: MRISripZeros(): dimension mismatch\n");
      return (1);
    }
  }

  for (k = 0; k < surf->nvertices; k++) {
    if (surf->vertices[k].ripflag) continue;
    if (mri)
      v = MRIgetVoxVal(mri, k, 0, 0, 0);
    else
      v = surf->vertices[k].val;
    if (fabs(v) < 1e-5) surf->vertices[k].ripflag = 1;
  }
  return (0);
}


/*!
  \fn int MRISripUnknown(MRIS *surf)
  \brief Sets the ripflag = 1 in places where the annotation is unknown
 */
int MRISripUnknown(MRIS *surf)
{
  int nripped = 0, vtxno, annot, annotid;

  for (vtxno = 0; vtxno < surf->nvertices; vtxno++) {
    annot = surf->vertices[vtxno].annotation;
    CTABfindAnnotation(surf->ct, annot, &annotid);
    if (annotid == 0 || annotid == -1) {
      surf->vertices[vtxno].ripflag = 1;
      nripped++;
    }
  }
  return (0);
}


int MRISripZeroThicknessRegions(MRI_SURFACE *mris)
{
  VERTEX *v;
  int vno;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (FZERO(v->whitex - v->pialx) && FZERO(v->whitey - v->pialy) && FZERO(v->whitez - v->pialz)) {
      v->ripflag = 1;
    }
  }
  return (NO_ERROR);
}


int MRISdilateRipped(MRIS *mris, int ndil)
{
  int vno, i, n;

  for (i = 0; i < ndil; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      v->tx = v->ripflag;
    }

    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      
      if (v->ripflag == 0) {
        continue;
      }

      // turn on ripflag of all neighbors of this (ripped) vertex
      for (n = 0; n < vt->vnum; n++) {
        VERTEX * const vn = &mris->vertices[vt->v[n]];
        if (vn->ripflag == 0) {
          vn->tx = 1;
        }
      }
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      v->ripflag = (int)v->tx;
    }
  }
  MRISsetRipInFacesWithRippedVertices(mris);
  return (NO_ERROR);
}

int MRISerodeRipped(MRIS *mris, int ndil)
{
  int vno, i, n, mn;

  for (i = 0; i < ndil; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      v->tx = 0;
    }

    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      mn = v->ripflag;
      for (n = 0; n < vt->vnum; n++) {
        VERTEX * const vn = &mris->vertices[vt->v[n]];
        mn = MIN(vn->ripflag, mn);
      }
      v->tx = mn;
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      v->ripflag = (int)v->tx;
    }
  }
  MRISsetRipInFacesWithRippedVertices(mris);
  return (NO_ERROR);
}

//=========================================================NOT YET REFURBISHED==================================================


/*---------------------------------------------------------------------
  MRI *MRISdilateConfined() - dilates surface mask niters iterations.
  If annotidmask >= 0, then dilation is confined to the annotidmask
  annotation. If newid >= 0, then the surface annot field of vertices
  in the dilated mask is set to the annot corresponding to newid.
  -------------------------------------------------------------------*/
MRI *MRISdilateConfined(MRIS *surf, MRI *mask, int annotidmask, int niters, int newid)
{
  int vtxno, annot, annotid, nnbrs, nbrvtxno, nthnbr, nthiter, new_annot;
  MRI *mri1, *mri2;
  float val;

  mri1 = MRIcopy(mask, NULL);
  mri2 = MRIcopy(mask, NULL);

  for (nthiter = 0; nthiter < niters; nthiter++) {
    // printf("iter %d\n",nthiter);

    for (vtxno = 0; vtxno < surf->nvertices; vtxno++) {

      /* Set to 0 if not in annotidmask (ie, dont dilate outside
      of annotidmask (if it is set) */
      if (annotidmask > -1) {
        annot = surf->vertices[vtxno].annotation;
        CTABfindAnnotation(surf->ct, annot, &annotid);
        if (annotid != annotidmask) {
          MRIsetVoxVal(mri2, vtxno, 0, 0, 0, 0);
          continue;
        }
      }

      // Check whether this vertex has been set
      val = MRIgetVoxVal(mri1, vtxno, 0, 0, 0);
      if (val) {
        MRIsetVoxVal(mri2, vtxno, 0, 0, 0, 1);
        continue;
      }

      // If it gets here, the vtx is in the annot and has not been set
      nnbrs = surf->vertices_topology[vtxno].vnum;
      for (nthnbr = 0; nthnbr < nnbrs; nthnbr++) {
        nbrvtxno = surf->vertices_topology[vtxno].v[nthnbr];
        if (surf->vertices[nbrvtxno].ripflag) {
          continue;  // skip ripped vtxs
        }
        val = MRIgetVoxVal(mri1, nbrvtxno, 0, 0, 0);
        if (val) {
          MRIsetVoxVal(mri2, vtxno, 0, 0, 0, 1);
          continue;
        }
      }
    }
    MRIcopy(mri2, mri1);
  }

  MRIfree(&mri2);
  // MRIwrite(mri1,"mymask2.mgh");

  if (newid > 0) {
    // Set annots in this mask to the given annot
    CTABannotationAtIndex(surf->ct, newid, &new_annot);
    for (vtxno = 0; vtxno < surf->nvertices; vtxno++) {
      if (MRIgetVoxVal(mri1, vtxno, 0, 0, 0)) {
        surf->vertices[vtxno].annotation = new_annot;
      }
    }
  }

  return (mri1);
}


int mrisMarkBadEdgeVertices(MRIS *mris, int mark)
{
  int vno, n, nfaces, m, vno2, nmarked;

  for (nmarked = vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    for (n = 0; n < vt->vnum; n++) {
      vno2 = vt->v[n];
      if (vno2 < vno) {
        continue;
      }
      for (nfaces = m = 0; m < vt->vnum; m++) {
        if (vt->v[m] == vno2) {
          continue;
        }
        if (vertexNeighbor(mris, vno2, vt->v[m]) && isFace(mris, vno, vno2, vt->v[m])) {
          nfaces++;
        }
      }
      if (nfaces != 2) {
        v->marked = mark;
        mris->vertices[vno2].marked = mark;
        nmarked += 2;
        break;
      }
    }
  }
  return (nmarked);
}


int MRISmarkRandomVertices(MRIS *mris, float prob_marked)
{
  int vno;
  VERTEX *v;
  float r;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    r = randomNumber(0.0, 1.0);
    if (r < prob_marked) {
      v->marked = 1;
    }
  }
  return (NO_ERROR);
}

bool mrisAnyVertexOfFaceMarked(MRIS *mris, int fno)
{
  FACE const * const f = &mris->faces[fno];

  int marked = 0;
  int n;
  for (n = 0; n < VERTICES_PER_FACE; n++) {     // check all, to avoid branch mispredicts
    marked |= mris->vertices[f->v[n]].marked;   // 
  }

  return !!marked;
}


bool triangleMarked(MRI_SURFACE *mris, int fno)
{
  return mrisAnyVertexOfFaceMarked(mris, fno);
}



int findNonMarkedFace(MRIS *mris, int vno, int vn1)
{
  cheapAssertValidVno(mris,vno);
  cheapAssert        (vno != vn1);
  VERTEX_TOPOLOGY const * const v = &mris->vertices_topology[vno];
  
  int nf;
  for (nf = 0; nf < v->num; nf++) {
    int const          fn = v->f[nf];
    FACE const * const f  = &mris->faces[fn];
    if (f->marked) continue;

    int i;
    for (i = 0; i < VERTICES_PER_FACE; i++) {
      if (f->v[i] == vn1) {
        return fn;
      }
    }
  }

  return -1;
}


void MRISsetCurvaturesToValues(MRIS *mris, int fno)
{
  int n;
  VERTEX *v;
  for (n = 0; n < mris->nvertices; n++) {
    v = &mris->vertices[n];
    if (v->ripflag) {
      continue;
    }
    ((VALS_VP *)v->vp)->vals[fno] = v->curv;
  }
}

void MRISsetCurvaturesToOrigValues(MRIS *mris, int fno)
{
  int n;
  VERTEX *v;
  for (n = 0; n < mris->nvertices; n++) {
    v = &mris->vertices[n];
    if (v->ripflag) {
      continue;
    }
    ((VALS_VP *)v->vp)->orig_vals[fno] = v->curv;
  }
}

void MRISsetOrigValuesToCurvatures(MRIS *mris, int fno)
{
  int n;
  VERTEX *v;
  for (n = 0; n < mris->nvertices; n++) {
    v = &mris->vertices[n];
    if (v->ripflag) {
      continue;
    }
    v->curv = ((VALS_VP *)v->vp)->orig_vals[fno];
  }
}

void MRISsetOrigValuesToValues(MRIS *mris, int fno)
{
  int n;
  VERTEX *v;
  for (n = 0; n < mris->nvertices; n++) {
    v = &mris->vertices[n];
    if (v->ripflag) {
      continue;
    }
    ((VALS_VP *)v->vp)->vals[fno] = ((VALS_VP *)v->vp)->orig_vals[fno];
  }
}


void MRISsetValuesToCurvatures(MRIS *mris, int fno)
{
  int n;
  VERTEX *v;
  for (n = 0; n < mris->nvertices; n++) {
    v = &mris->vertices[n];
    if (v->ripflag) {
      continue;
    }
    v->curv = ((VALS_VP *)v->vp)->vals[fno];
  }
}


int MRISmarkVerticesWithValOverThresh(MRI_SURFACE *mris, float thresh)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) 
      continue;
    v->marked = (int)(v->val > thresh) ;
  }
  return (NO_ERROR);
}


int MRISminFilterCurvatures(MRI_SURFACE *mris, int niter)
{
  int i, vno, vnb, vnum;
  float curv;

  for (i = 0; i < niter; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) {
        continue;
      }
      curv = v->curv;
      int const * pnb  = vt->v;
      vnum = vt->vnum;
      for (vnb = 0; vnb < vnum; vnb++) {
        VERTEX const * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
        if (vn->ripflag) {
          continue;
        }
        if (vn->curv < curv) {
          curv = vn->curv;
        }
      }
      v->tdx = curv;
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }
      v->curv = v->tdx;
    }
  }
  mrisComputeCurvatureMinMax(mris);
  return (NO_ERROR);
}


int MRISmaxFilterCurvatures(MRI_SURFACE *mris, int niter)
{
  int i, vno, vnb, vnum;
  float curv;

  for (i = 0; i < niter; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) {
        continue;
      }
      curv = v->curv;
      int const * pnb  = vt->v;
      vnum = vt->vnum;
      for (vnb = 0; vnb < vnum; vnb++) {
        VERTEX const * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
        if (vn->ripflag) {
          continue;
        }
        if (vn->curv > curv) {
          curv = vn->curv;
        }
      }
      v->tdx = curv;
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }
      v->curv = v->tdx;
    }
  }
  mrisComputeCurvatureMinMax(mris);
  return (NO_ERROR);
}

/*-----------------------------------------------------
int MRISaverageCurvatures(MRI_SURFACE *mris, int navgs)
Performs navgs steps of iterative spatial smoothing on curv.
------------------------------------------------------*/
int MRISaverageCurvatures(MRI_SURFACE *mris, int navgs)
{
  int i, vno, vnb, vnum;
  float curv, num;

  for (i = 0; i < navgs; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) {
        continue;
      }
      curv = v->curv;
      int const * pnb  = vt->v;
      vnum = vt->vnum;
      for (num = 0.0f, vnb = 0; vnb < vnum; vnb++) {
        VERTEX const * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
        if (vn->ripflag) {
          continue;
        }
        num++;
        curv += vn->curv;
      }
      num++; /*  account for central vertex */
      v->tdx = curv / num;
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }
      v->curv = v->tdx;
    }
  }
  mrisComputeCurvatureMinMax(mris);
  return (NO_ERROR);
}

int MRISaverageMarkedCurvatures(MRI_SURFACE *mris, int navgs)
{
  int i, vno, vnb, vnum;
  float curv, num;

  for (i = 0; i < navgs; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag || !v->marked) {
        continue;
      }
      curv = v->curv;
      int const * pnb  = vt->v;
      vnum = vt->vnum;
      for (num = 0.0f, vnb = 0; vnb < vnum; vnb++) {
        VERTEX const * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
        if (vn->ripflag) {
          continue;
        }
        num++;
        curv += vn->curv;
      }
      num++; /*  account for central vertex */
      v->tdx = curv / num;
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      if (v->ripflag || v->marked == 0) {
        continue;
      }
      v->curv = v->tdx;
    }
  }
  mrisComputeCurvatureMinMax(mris);
  return (NO_ERROR);
}

int MRISaverageMarkedVals(MRI_SURFACE *mris, int navgs)
{
  int i, vno, vnb, vnum;
  float val, num;

  for (i = 0; i < navgs; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag || v->marked == 0) {
        continue;
      }
      val = v->val;
      int const * pnb  = vt->v;
      vnum = vt->vnum;
      for (num = 0.0f, vnb = 0; vnb < vnum; vnb++) {
        VERTEX const * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
        if (vn->ripflag || vn->marked == 0) {
          continue;
        }
        num++;
        val += vn->val;
      }
      num++; /*  account for central vertex */
      v->tdx = val / num;
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      if (v->ripflag || v->marked == 0) {
        continue;
      }
      v->val = v->tdx;
    }
  }
  return (NO_ERROR);
}

int MRISaverageVals(MRI_SURFACE *mris, int navgs)
{
  int i, vno;

  for (i = 0; i < navgs; i++) {
    ROMP_PF_begin
#ifdef HAVE_OPENMP
    #pragma omp parallel for if_ROMP(experimental) shared(mris, i) schedule(static, 1)
#endif
    for (vno = 0; vno < mris->nvertices; vno++) {
      ROMP_PFLB_begin
      
      int vnb, vnum;
      float val, num;

      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) ROMP_PF_continue;

      val = v->val;
      int const * pnb  = vt->v;
      vnum = vt->vnum;
      for (num = 0.0f, vnb = 0; vnb < vnum; vnb++) {
        VERTEX const * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
        if (vn->ripflag) continue;

        num++;
        val += vn->val;
      }
      num++; /*  account for central vertex */
      v->tdx = val / num;
      ROMP_PFLB_end
    }
    ROMP_PF_end

    ROMP_PF_begin
#ifdef HAVE_OPENMP
    #pragma omp parallel for if_ROMP(experimental) shared(mris, i) schedule(static, 1)
#endif
    for (vno = 0; vno < mris->nvertices; vno++) {
      ROMP_PFLB_begin
      VERTEX * const v = &mris->vertices[vno];
      if (v->ripflag) ROMP_PF_continue;

      v->val = v->tdx;
      ROMP_PFLB_end
    }
    ROMP_PF_end
  }

  return (NO_ERROR);
}

// spatially average the v->d field
int MRISaverageD(MRI_SURFACE *mris, int navgs)
{
  int i, vno, vnb, vnum;
  float val, num;

  for (i = 0; i < navgs; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) {
        continue;
      }
      val = v->d;
      int const * pnb  = vt->v;
      vnum = vt->vnum;
      for (num = 0.0f, vnb = 0; vnb < vnum; vnb++) {
        VERTEX const * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
        if (vn->ripflag) {
          continue;
        }
        num++;
        val += vn->d;
      }
      num++; /*  account for central vertex */
      v->tdx = val / num;
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }
      v->d = v->tdx;
    }
  }
  return (NO_ERROR);
}

int MRISmedianFilterVals(MRI_SURFACE *mris, int nmedians)
{
  int i, vno, vnb, vnum, num;
  float val_list[MAX_NEIGHBORS];

  for (i = 0; i < nmedians; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) {
        continue;
      }
      int const * pnb = vt->v;
      //      vnum = vt->vtotal ;   BRF - used to be vnum
      vnum = vt->vnum;
      val_list[0] = v->val;
      for (num = 1, vnb = 0; vnb < vnum; vnb++) {
        VERTEX const * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
        if (vn->ripflag) {
          continue;
        }

        val_list[num++] = vn->val;
      }
      qsort(val_list, num, sizeof(val_list[0]), mris_sort_compare_float);
      if (ISODD(num)) {
        v->tdx = val_list[(num - 1) / 2];
      }
      else {
        v->tdx = (val_list[num / 2] + val_list[num / 2 - 1]) / 2;
      }
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }
      v->val = v->tdx;
    }
  }
  return (NO_ERROR);
}


#if 0
int MRISmedianFilterVerexPositions(MRI_SURFACE *mris, int nmedians)
{
  int i, vno, vnb, vnum, num;
  float val_list[MAX_NEIGHBORS];

  for (i = 0; i < nmedians; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) {
        continue;
      }
      int const * pnb  = vt->v;
      vnum = vt->vnum;
      val_list[0] = v->x;
      for (num = 1, vnb = 0; vnb < vnum; vnb++) {
        VERTEX const * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
        if (vn->ripflag) {
          continue;
        }

        val_list[num++] = vn->x;
      }
      qsort(val_list, num, sizeof(val_list[0]), mris_sort_compare_float);
      if (ISODD(num)) {
        v->tdx = val_list[(num - 1) / 2];
      }
      else {
        v->tdx = (val_list[num / 2] + val_list[num / 2 - 1]) / 2;
      }

      pnb  = vt->v;
      vnum = vt->vnum;
      val_list[0] = v->y;
      for (num = 1, vnb = 0; vnb < vnum; vnb++) {
        VERTEX const * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
        if (vn->ripflag) {
          continue;
        }

        val_list[num++] = vn->y;
      }
      qsort(val_list, num, sizeof(val_list[0]), mris_sort_compare_float);
      if (ISODD(num)) {
        v->tdy = val_list[(num - 1) / 2];
      }
      else {
        v->tdy = (val_list[num / 2] + val_list[num / 2 - 1]) / 2;
      }

      pnb  = vt->v;
      vnum = vt->vnum;
      val_list[0] = v->z;
      for (num = 1, vnb = 0; vnb < vnum; vnb++) {
        VERTEX const * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
        if (vn->ripflag) {
          continue;
        }

        val_list[num++] = vn->z;
      }
      qsort(val_list, num, sizeof(val_list[0]), mris_sort_compare_float);
      if (ISODD(num)) {
        v->tdz = val_list[(num - 1) / 2];
      }
      else {
        v->tdz = (val_list[num / 2] + val_list[num / 2 - 1]) / 2;
      }
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }
      v->x = v->tdx;
      v->y = v->tdy;
      v->z = v->tdz;
    }
  }
  return (NO_ERROR);
}
#endif


int MRISgaussianFilterD(MRI_SURFACE *mris, double wt)
{
  int vno, vnb, *pnb, vnum;
  double nbr_wt, val, norm;

  nbr_wt = 1 - wt;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    pnb  = vt->v;
    vnum = vt->vnum;

    val = wt * v->d;

    for (norm = wt, vnb = 0; vnb < vnum; vnb++) {
      VERTEX const * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
      if (vn->ripflag) {
        continue;
      }

      val += nbr_wt * vn->d;
      norm += nbr_wt;
    }
    v->tdx = val / norm;
  }
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX * const v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->d = v->tdx;
  }
  return (NO_ERROR);
}

int MRISmedianFilterD(MRI_SURFACE *mris, int nmedians, int vtotal)
{
  int i, vno, vnb, vnum, num;
  float val_list[MAX_NEIGHBORS];

  for (i = 0; i < nmedians; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) {
        continue;
      }
      int const * pnb = vt->v;
      if (vtotal) {
        vnum = vt->vtotal;
      }
      else {
        vnum = vt->vnum;
      }
      val_list[0] = v->d;
      for (num = 1, vnb = 0; vnb < vnum; vnb++) {
        VERTEX const * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
        if (vn->ripflag) {
          continue;
        }

        val_list[num++] = vn->d;
      }
      qsort(val_list, num, sizeof(val_list[0]), mris_sort_compare_float);
      if (ISODD(num)) {
        v->tdx = val_list[(num - 1) / 2];
      }
      else {
        v->tdx = (val_list[num / 2] + val_list[num / 2 - 1]) / 2;
      }
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }
      v->d = v->tdx;
    }
  }
  return (NO_ERROR);
}


int MRISmedianFilterCurvature(MRI_SURFACE *mris, int nmedians)
{
  int i, vno, vnb, vnum, num;
  float val_list[MAX_NEIGHBORS];

  for (i = 0; i < nmedians; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) {
        continue;
      }
      int const * pnb  = vt->v;
      vnum = vt->vnum;
      val_list[0] = v->curv;
      for (num = 1, vnb = 0; vnb < vnum; vnb++) {
        VERTEX const * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
        if (vn->ripflag) {
          continue;
        }

        val_list[num++] = vn->curv;
      }
      qsort(val_list, num, sizeof(val_list[0]), mris_sort_compare_float);
      if (ISODD(num)) {
        v->tdx = val_list[(num - 1) / 2];
      }
      else {
        v->tdx = (val_list[num / 2] + val_list[num / 2 - 1]) / 2;
      }
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }
      v->curv = v->tdx;
    }
  }
  return (NO_ERROR);
}


int MRISmedianFilterVal2s(MRI_SURFACE *mris, int nmedians)
{
  int i, vno, vnb, vnum, num;
  float val_list[MAX_NEIGHBORS];

  for (i = 0; i < nmedians; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) {
        continue;
      }
      int const * pnb  = vt->v;
      vnum = vt->vnum;
      val_list[0] = v->val2;
      for (num = 1, vnb = 0; vnb < vnum; vnb++) {
        VERTEX const * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
        if (vn->ripflag) {
          continue;
        }

        val_list[num++] = vn->val2;
      }
      qsort(val_list, num, sizeof(val_list[0]), mris_sort_compare_float);
      if (ISODD(num)) {
        v->tdx = val_list[(num - 1) / 2];
      }
      else {
        v->tdx = (val_list[num / 2] + val_list[num / 2 - 1]) / 2;
      }
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }
      v->val2 = v->tdx;
    }
  }
  return (NO_ERROR);
}


int MRISmedianFilterVal2baks(MRI_SURFACE *mris, int nmedians)
{
  int i, vno, vnb, vnum, num;
  float val_list[MAX_NEIGHBORS];

  for (i = 0; i < nmedians; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) {
        continue;
      }
      int const * pnb  = vt->v;
      vnum = vt->vnum;
      val_list[0] = v->val2bak;
      for (num = 1, vnb = 0; vnb < vnum; vnb++) {
        VERTEX const * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
        if (vn->ripflag) {
          continue;
        }

        val_list[num++] = vn->val2bak;
      }
      qsort(val_list, num, sizeof(val_list[0]), mris_sort_compare_float);
      if (ISODD(num)) {
        v->tdx = val_list[(num - 1) / 2];
      }
      else {
        v->tdx = (val_list[num / 2] + val_list[num / 2 - 1]) / 2;
      }
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }
      v->val2bak = v->tdx;
    }
  }
  return (NO_ERROR);
}


int MRISaverageVal2s(MRI_SURFACE *mris, int navgs)
{
  int i, vno, vnb, vnum;
  float val, num;

  for (i = 0; i < navgs; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) {
        continue;
      }
      val  = v->val2;
      int const * pnb  = vt->v;
      vnum = vt->vnum;
      for (num = 0.0f, vnb = 0; vnb < vnum; vnb++) {
        VERTEX const * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
        if (vn->ripflag) {
          continue;
        }
        num++;
        val += vn->val2;
      }
      num++; /*  account for central vertex */
      v->tdx = val / num;
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }
      v->val2 = v->tdx;
    }
  }
  return (NO_ERROR);
}


int MRISaccumulateMeansOnSurface(MRI_SURFACE *mris, int total_dof, int new_dof)
{
  int vno, ndof;
  VERTEX *v;

  ndof = total_dof + new_dof;
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }

    v->mean = (v->mean * total_dof + v->val * new_dof) / ndof;
  }

  return (NO_ERROR);
}


int MRISaccumulateImaginaryMeansOnSurface(MRI_SURFACE *mris, int total_dof, int new_dof)
{
  int vno, ndof;
  VERTEX *v;

  ndof = total_dof + new_dof;
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }

    v->mean_imag = (v->mean_imag * total_dof + v->val * new_dof) / ndof;
  }

  return (NO_ERROR);
}


int MRISaccumulateStandardErrorsOnSurface(MRI_SURFACE *mris, int total_dof, int new_dof)
{
  int vno, ndof;
  VERTEX *v;
  double var;

  ndof = total_dof + new_dof;
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }

    var = (v->std_error * SQR(total_dof) + v->val * SQR(new_dof)) / SQR(ndof);
    v->std_error = var;
  }

  return (NO_ERROR);
}


int MRIScomputeAverageCircularPhaseGradient(MRI_SURFACE *mris, LABEL *area, float *pdx, float *pdy, float *pdz)
{
  int N, vno, n, i;
  VECTOR *vdf, *vfz;
  MATRIX *mz, *mzt, *mztz, *mztz_inv, *mztz_inv_zt;
  double x0, y0, z0, f0, x1, y1, z1, f1, dx, dy, dz, df;

  dx = dy = dz = 0.0;
  for (i = 0; i < area->n_points; i++) {
    vno = area->lv[i].vno;
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    x0 = v->x;
    y0 = v->y;
    z0 = v->z;
    f0 = atan2(v->imag_val, v->val);

    /* first count # of valid neighbors */
    for (N = n = 0; n < vt->vnum; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      if (vn->ripflag) {
        continue;
      }
      N++;
    }
    /* now allocate vectors and matrix */
    vdf = VectorAlloc(N, MATRIX_REAL); /* function deltas */
    if (mris->patch) {
      mz = MatrixAlloc(N, 2, MATRIX_REAL); /* vertex spacing deltas */
    }
    else {
      mz = MatrixAlloc(N, 3, MATRIX_REAL); /* vertex spacing deltas */
    }

    /* now fill in matrix and vector entries */
    for (n = 0; n < vt->vnum; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      if (vn->ripflag) {
        continue;
      }
      x1 = vn->x;
      y1 = vn->y;
      z1 = vn->z;
      f1 = atan2(vn->imag_val, vn->val);
      df = f1 - f0;
      while (df > M_PI) {
        df -= 2 * M_PI;
      }
      while (df < -M_PI) {
        df += 2 * M_PI;
      }
      VECTOR_ELT(vdf, n + 1) = df;
      *MATRIX_RELT(mz, n + 1, 1) = x1 - x0;
      *MATRIX_RELT(mz, n + 1, 2) = y1 - y0;
      if (!mris->patch) {
        *MATRIX_RELT(mz, n + 1, 3) = z1 - z0;
      }
    }
    mzt = MatrixTranspose(mz, NULL);
    mztz = MatrixMultiply(mzt, mz, NULL);
    mztz_inv = MatrixSVDInverse(mztz, NULL);
    if (mztz_inv) {
      mztz_inv_zt = MatrixMultiply(mztz_inv, mzt, NULL);
      vfz = MatrixMultiply(mztz_inv_zt, vdf, NULL);
      v->dx = VECTOR_ELT(vfz, 1);
      v->dy = VECTOR_ELT(vfz, 2);
      if (!mris->patch) {
        v->dz = VECTOR_ELT(vfz, 3);
      }
      else {
        v->dz = 0.0f;
      }
      dx += v->dx;
      dy += v->dy;
      dz += v->dz;
    }

    VectorFree(&vdf);
    MatrixFree(&mz);
    MatrixFree(&mzt);
    MatrixFree(&mztz);
    if (mztz_inv) {
      VectorFree(&vfz);
      MatrixFree(&mztz_inv);
      MatrixFree(&mztz_inv_zt);
    }
  }

  *pdx = dx /= (float)area->n_points;
  *pdy = dy /= (float)area->n_points;
  *pdz = dz /= (float)area->n_points;
  return (NO_ERROR);
}


int MRISaverageVal2baks(MRI_SURFACE *mris, int navgs)
{
  int i, vno, vnb, vnum;
  float val, num;

  for (i = 0; i < navgs; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) {
        continue;
      }
      val  = v->val2bak;
      int const * pnb  = vt->v;
      vnum = vt->vnum;
      for (num = 0.0f, vnb = 0; vnb < vnum; vnb++) {
        VERTEX const * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
        if (vn->ripflag) {
          continue;
        }
        num++;
        val += vn->val2bak;
      }
      num++; /*  account for central vertex */
      v->tdx = val / num;
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }
      v->val2bak = v->tdx;
    }
  }
  return (NO_ERROR);
}


int MRISmodeFilterVals(MRI_SURFACE *mris, int niter)
{
  int *histo;
  int i, n, vno, ino, index, max_histo, max_index, nchanged, nzero, max_val;

  for (max_val = vno = 0; vno < mris->nvertices; vno++) {
    VERTEX const * const v = &mris->vertices[vno];
    if (v->val > max_val) max_val = v->val;
  }
  histo = (int *)calloc(max_val + 1, sizeof(int));
  if (histo == NULL)
    ErrorExit(ERROR_NOMEMORY, "MRISmodeFilterVals: could not allocate histo array of %d ints", max_val + 1);

  MRISclearMarks(mris); /* v->marked = 0 means it hasn't converged yet */
  for (ino = 0; ino < niter; ino++) {
    nzero = nchanged = 0;
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (vno == Gdiag_no) DiagBreak();

      if (v->ripflag || v->marked) continue;

      if (vno == Gdiag_no) DiagBreak();

      if (nint(v->val) == 0) nzero++;

      // initialize
      memset(histo, 0, (max_val + 1) * sizeof(*histo));
      // create histogram
      for (n = 0; n < vt->vnum; n++) {
        VERTEX const * const vn = &mris->vertices[vt->v[n]];
        index = (int)nint(vn->val);
        if (index < 0) continue;

        histo[index]++;
      }
      max_histo = histo[0];
      max_index = 0;
      for (i = 1; i <= max_val; i++) {
        if (histo[i] > max_histo) {
          max_histo = histo[i];
          max_index = i;
        }
      }
      v->valbak = max_index;
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      if (v->ripflag) continue;

      if (vno == Gdiag_no) DiagBreak();

      if (v->val != v->valbak) {
        v->marked = 0; /* process it again */
        nchanged++;
      }
      else
        v->marked = 1; /* didn't change */

      v->val = v->valbak;
    }

    /* unmark all nbrs of unmarked vertices */
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX          const * const v  = &mris->vertices         [vno];
      if (v->ripflag || v->marked == 1) continue;

      if (vno == Gdiag_no) DiagBreak();
      for (n = 0; n < vt->vnum; n++) {
        VERTEX * const vn = &mris->vertices[vt->v[n]];
        vn->marked = 0; /* process it again */
      }
    }
    printf("iter %d: %d changed, %d zero\n", ino, nchanged, nzero);
    if (!nchanged) break;
  }
  free(histo);
  MRISclearMarks(mris);
  return (NO_ERROR);
}

int MRISmodeFilterZeroVals(MRI_SURFACE *mris)
{
  int *histo, i, n, vno, ino, max_val, index, max_histo, max_index, nchanged, nzero;

  for (max_val = vno = 0; vno < mris->nvertices; vno++) {
    VERTEX const * const v = &mris->vertices[vno];
    if (v->val > max_val) max_val = v->val;
  }
  histo = (int *)calloc(max_val + 1, sizeof(int));
  if (histo == NULL)
    ErrorExit(ERROR_NOMEMORY, "MRISmodeFilterVals: could not allocate histo array of %d ints", max_val + 1);

  MRISclearMarks(mris); /* v->marked = 0 means it hasn't converged yet */
  ino = 0;
  do {
    nzero = nchanged = 0;
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (vno == Gdiag_no) DiagBreak();

      if (v->ripflag || v->marked) continue;

      if (nint(v->val) == 0)
        nzero++;
      else {
        v->valbak = v->val;
        continue;  // only process vertices that have v->val == 0
      }

      memset(histo, 0, (max_val + 1) * sizeof(*histo));
      for (n = 0; n < vt->vnum; n++) {
        VERTEX const * const vn = &mris->vertices[vt->v[n]];
        index = (int)nint(vn->val);
        if (index < 0) continue;

        histo[index]++;
      }
      max_histo = 0;
      max_index = 0;
      for (i = 1; i <= max_val; i++) {
        if (histo[i] > max_histo) {
          max_histo = histo[i];
          max_index = i;
        }
      }
      v->valbak = max_index;
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      if (v->ripflag) continue;

      if (vno == Gdiag_no) DiagBreak();

      if (v->val != v->valbak) {
        v->marked = 0; /* process it again */
        nchanged++;
      }
      else
        v->marked = 1; /* didn't change */

      v->val = v->valbak;
    }

    /* unmark all nbrs of unmarked vertices */
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX          const * const v  = &mris->vertices         [vno];
      if (v->ripflag || v->marked == 1) continue;

      if (vno == Gdiag_no) DiagBreak();

      for (n = 0; n < vt->vnum; n++) {
        VERTEX * const vn = &mris->vertices[vt->v[n]];
        vn->marked = 0; /* process it again */
      }
    }
    printf("iter %d: %d changed, %d zero\n", ino++, nchanged, nzero);
    if (!nchanged) {
      break;
    }
  } while (nchanged > 0 && nzero > 0);
  MRISclearMarks(mris);
  free(histo);
  return (NO_ERROR);
}

int MRISmodeFilterAnnotations(MRI_SURFACE *mris, int niter)
{
  int *histo, i, n, vno, ino, index, max_histo, max_index, max_annotation, *annotations, nchanged = 0;

  for (max_index = vno = 0; vno < mris->nvertices; vno++) {
    VERTEX const * const v = &mris->vertices[vno];
    index = annotation_to_index(v->annotation);
    if (index > max_index) max_index = index;
  }
  histo = (int *)calloc(max_index + 1, sizeof(int));
  if (histo == NULL)
    ErrorExit(ERROR_NOMEMORY, "MRISmodeFilterVals: could not allocate histo array of %d ints", max_index + 1);
  annotations = (int *)calloc(max_index + 1, sizeof(int));
  if (annotations == NULL)
    ErrorExit(ERROR_NOMEMORY, "MRISmodeFilterVals: could not allocate annotation array of %d ints", max_index + 1);

  // reset the annotation table using the surface's own
  // colortable when it's available
  if (mris->ct != NULL) set_atable_from_ctable(mris->ct);

  for (ino = 0; ino < niter; ino++) {
    for (nchanged = vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) continue;

      if (vno == Gdiag_no) DiagBreak();

      memset(histo, 0, (max_index + 1) * sizeof(*histo));
      memset(annotations, 0, (max_index + 1) * sizeof(*annotations));
      for (n = 0; n < vt->vtotal; n++) {
        VERTEX const * const vn = &mris->vertices[vt->v[n]];
        index = annotation_to_index(vn->annotation);
        if (index < 0) continue;

        histo[index]++;
        annotations[index] = vn->annotation;
      }
      index = annotation_to_index(v->annotation);
      if (index >= 0) {
        annotations[index] = v->annotation;
        histo[index]++;
        max_histo = histo[index];
        max_annotation = v->annotation;
      }
      else
        max_histo = max_annotation = 0;

      for (i = 1; i <= max_index; i++) {
        if (histo[i] > max_histo) {
          max_histo = histo[i];
          max_annotation = annotations[i];
        }
      }
      v->undefval = max_annotation;
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      if (v->ripflag) continue;

      if (vno == Gdiag_no) DiagBreak();

      if (v->annotation != v->undefval) nchanged++;

      v->annotation = v->undefval;
    }
    if (nchanged == 0) break;
  }
  printf("%d filter iterations complete (%d requested, %d changed)\n", ino, niter, nchanged);

  free(histo);
  free(annotations);
  return (NO_ERROR);
}


int MRISsoapBubbleVals(MRI_SURFACE *mris, int navgs)
{
  int vno, n, i, cmpt, nmarked;
  double mean;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "performing soap bubble smoothing of vals for %d iterations.\n", navgs);
  for (i = 0; i < navgs; i++) {
    for (nmarked = vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag || (v->marked == 1)) {
        continue;
      }

      /* compute average of self and neighbors */
      mean = 0.0;
      cmpt = 0;
      for (n = 0; n < vt->vnum; n++) {
        VERTEX const * const vn = &mris->vertices[vt->v[n]];
        if (vn->marked) {
          mean += vn->val;
          cmpt++;
        }
      }
      if (cmpt > 0) /* had some neighbors with real values */
      {
        v->val = mean / (double)(cmpt);
        if (!v->marked) /* has never been computed before */
        {
          nmarked++;
        }
        v->marked = 2; /* has a real value, but is not fixed */
      }
    }
#if 0
    if (!nmarked)
    {
      break ;
    }
#endif
  }

  /*  fprintf(stdout, "\n") ;*/
  return (NO_ERROR);
}


/*-----------------------------------------------------
    apply a soap bubble to the vertex->d field
     (used for target distances frequently)
  ------------------------------------------------------*/
int MRISsoapBubbleD(MRI_SURFACE *mris, int navgs)
{
  int vno, n, i, cmpt, nmarked;

  double mean;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "performing soap bubble smoothing of D vals for %d iterations.\n", navgs);
  for (i = 0; i < navgs; i++) {
    for (nmarked = vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag || (v->marked == 1)) {
        continue;
      }

      /* compute average of self and neighbors */
      mean = 0.0;
      cmpt = 0;
      for (n = 0; n < vt->vnum; n++) {
        VERTEX const * const vn = &mris->vertices[vt->v[n]];
        if (vn->marked) {
          mean += vn->d;
          cmpt++;
        }
      }
      if (cmpt > 0) /* had some neighbors with real values */
      {
        v->d = mean / (double)(cmpt);
        if (!v->marked) /* has never been computed before */
        {
          nmarked++;
        }
        v->marked = 2; /* has a real value, but is not fixed */
      }
    }
#if 0
    if (!nmarked)
    {
      break ;
    }
#endif
  }

  /*  fprintf(stdout, "\n") ;*/
  return (NO_ERROR);
}

int detectContrast(MRIS *mris)
{
  int n, nv;
  double wm, gm;
  VERTEX *v;

  fprintf(WHICH_OUTPUT, "Detecting contrast direction:");
  wm = gm = 0.0;
  nv = 0;
  for (n = 0; n < mris->nvertices; n++) {
    v = &mris->vertices[n];
    if (v->ripflag) {
      continue;
    }
    nv++;
    wm += v->val2;
    gm += v->val2bak;
  }
  wm /= nv;
  gm /= nv;

  if (wm > gm) {
    fprintf(WHICH_OUTPUT, " right direction [ GM (%3.2f) < WM (%3.2f) ]\n", gm, wm);
    return 1;
  }
  else {
    fprintf(WHICH_OUTPUT, " inverted [ WM (%3.2f) < GM (%3.2f) ]\n", wm, gm);
    return -1;
  }
}


int SetHop(int CenterVtx, MRI_SURFACE *Surf, int HopNo, int MaxHopNo)
{
  int nbr, nbr_vtx, nbr_hopno, nbr_has_been_center;

  if (HopNo >= MaxHopNo) return (0);

  Surf->vertices[CenterVtx].valbak = 1;  // has been a center

  for (nbr = 0; nbr < Surf->vertices_topology[CenterVtx].vnum; nbr++) {
    nbr_vtx = Surf->vertices_topology[CenterVtx].v[nbr];
    nbr_hopno = Surf->vertices[nbr_vtx].undefval;
    if (nbr_hopno != 0 && nbr_hopno < HopNo + 1) continue;
    Surf->vertices[nbr_vtx].undefval = HopNo + 1;
  }

  for (nbr = 0; nbr < Surf->vertices_topology[CenterVtx].vnum; nbr++) {
    nbr_vtx = Surf->vertices_topology[CenterVtx].v[nbr];
    nbr_has_been_center = Surf->vertices[nbr_vtx].valbak;
    if (nbr_has_been_center) continue;
    SetHop(nbr_vtx, Surf, HopNo + 1, MaxHopNo);
  }

  return (0);
}

SURFHOPLIST *SetSurfHopListAlloc(MRI_SURFACE *Surf, int nHops)
{
  SURFHOPLIST *shl;
  int nthhop;
  shl = (SURFHOPLIST *)calloc(sizeof(SURFHOPLIST), 1);
  shl->nhops = nHops;
  shl->hit = (char *)calloc(sizeof(char), Surf->nvertices);
  shl->nperhop = (int *)calloc(sizeof(int), nHops);
  shl->nperhop_alloced = (int *)calloc(sizeof(int), nHops);
  shl->vtxlist = (int **)calloc(sizeof(int *), nHops);
  for (nthhop = 0; nthhop < nHops; nthhop++) {
    shl->nperhop_alloced[nthhop] = 100;
    shl->vtxlist[nthhop] = (int *)calloc(sizeof(int), 100);
  }
  return (shl);
}

int SurfHopListFree(SURFHOPLIST **shl0)
{
  SURFHOPLIST *shl = *shl0;
  int nthhop;
  for (nthhop = 0; nthhop < shl->nhops; nthhop++) free(shl->vtxlist[nthhop]);
  free(shl->vtxlist);
  free(shl->nperhop_alloced);
  free(shl->nperhop);
  free(shl->hit);
  free(shl);
  *shl0 = NULL;
  return (0);
}

/*!
  \fn SURFHOPLIST *SetSurfHopList(int CenterVtx, MRI_SURFACE *Surf, int nHops)
  \brief Fills in a SURFHOPLIST structure. This is a structure that indicates
  which vertices are a given number of links (hops) away. This can be used to
  set neighborhoods or compute the spatial autocorrelation function.
*/

SURFHOPLIST *SetSurfHopList(int CenterVtx, MRI_SURFACE *Surf, int nHops)
{
  SURFHOPLIST *shl;
  int nthhop, nhits, nthvtx, vtxno, nthnbr, nbr_vtxno;

  shl = SetSurfHopListAlloc(Surf, nHops);
  shl->cvtx = CenterVtx;

  // 0 hop is the center vertex
  shl->nperhop[0] = 1;
  shl->vtxlist[0][0] = CenterVtx;
  shl->hit[CenterVtx] = 1;

  // go through all hops
  for (nthhop = 1; nthhop < nHops; nthhop++) {
    nhits = 0;
    // go through each vertex of the previous hop
    for (nthvtx = 0; nthvtx < shl->nperhop[nthhop - 1]; nthvtx++) {
      vtxno = shl->vtxlist[nthhop - 1][nthvtx];
      // go through the neighbors of this vertex
      for (nthnbr = 0; nthnbr < Surf->vertices_topology[vtxno].vnum; nthnbr++) {
        nbr_vtxno = Surf->vertices_topology[vtxno].v[nthnbr];
        if (shl->hit[nbr_vtxno]) continue;  // ignore if it has been hit already
        shl->hit[nbr_vtxno] = 1;
        if (nhits >= shl->nperhop_alloced[nthhop]) {
          // realloc if need to
          shl->nperhop_alloced[nthhop] += 100;
          shl->vtxlist[nthhop] = (int *)realloc(shl->vtxlist[nthhop], sizeof(int) * shl->nperhop_alloced[nthhop]);
        }
        // assign the vertex
        shl->vtxlist[nthhop][nhits] = nbr_vtxno;
        nhits++;
      }
    }
    shl->nperhop[nthhop] = nhits;
  }

  // This assigns the hop number to the undefval, good for debugging
  for (nthhop = 0; nthhop < nHops; nthhop++) {
    // printf("nper hop %d %6d\n",nthhop,shl->nperhop[nthhop]);
    for (nthvtx = 0; nthvtx < shl->nperhop[nthhop]; nthvtx++) {
      vtxno = shl->vtxlist[nthhop][nthvtx];
      Surf->vertices[vtxno].undefval = nthhop;
    }
  }

  return (shl);
}


// average the stats on the surface, and propagate marks
// outwards to new non-zero locations
int MRISaverageMarkedStats(MRI_SURFACE *mris, int navgs)
{
  int i, vno, vnb, *pnb, vnum;
  float val, num;

  for (i = 0; i < navgs; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (vno == Gdiag_no) {
        DiagBreak();
      }
      if (v->ripflag) {
        continue;
      }
      val = v->stat;
      pnb = vt->v;
      vnum = vt->vnum;
      for (num = 0.0f, vnb = 0; vnb < vnum; vnb++) {
        VERTEX const * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
        if (vn->ripflag) {
          continue;
        }
        num++;
        val += vn->stat;
      }
      num++; /*  account for central vertex */
      v->tdx = val / num;
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v  = &mris->vertices         [vno];
      if (vno == Gdiag_no) {
        DiagBreak();
      }
      if (v->ripflag) {
        continue;
      }
      v->stat = v->tdx;
      if (v->marked == 0 && !FZERO(v->stat)) {
        v->marked = 2;
      }
      if (FZERO(v->stat)) {
        DiagBreak();
      }
    }
  }
  return (NO_ERROR);
}


LABEL *MRISannotation_to_label(MRI_SURFACE *mris, int annot_index)
{
  int vno, nv, annot;
  VERTEX *v;
  LABEL *area;

  CTABannotationAtIndex(mris->ct, annot_index, &annot);
  for (nv = vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (v->annotation == annot) {
      nv++;
    }
  }
  if (nv == 0) {
    return (NULL);
  }
  area = LabelAlloc(nv, NULL, mris->ct->entries[annot_index]->name);
  for (nv = vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (v->annotation == annot) {
      area->lv[nv].vno = vno;
      area->n_points++;
      nv++;
    }
  }
  return (area);
}


#include "histo.h"
HISTOGRAM *MRISgetHistogram(MRI_SURFACE *mris, int nbins, int field)
{
  int vno, bin;
  VERTEX *v;
  double fmin, fmax, bin_size, val;
  HISTOGRAM *h;

  fmin = 1e10;
  fmax = -fmin;
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
#if GCC_VERSION > 80000
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#endif
    switch (field) {
      default:
        ErrorExit(ERROR_BADPARM, "MRISgetHistogram: unknown field %d", field);
      case SCLV_VAL:
        val = (v)->val;
        break;
      case SCLV_VAL2:
        val = v->val2;
        break;
      case SCLV_VALBAK:
        val = v->valbak;
        break;
      case SCLV_VAL2BAK:
        val = v->val2bak;
        break;
      case SCLV_VALSTAT:
        val = v->stat;
        break;
      case SCLV_IMAG_VAL:
        val = v->imag_val;
        break;
      case SCLV_MEAN:
        val = v->mean;
        break;
      case SCLV_MEAN_IMAG:
        val = v->mean_imag;
        break;
      case SCLV_STD_ERROR:
        val = v->std_error;
        break;
    }
#if GCC_VERSION > 80000
#pragma GCC diagnostic pop
#endif
    if (val < fmin) {
      fmin = val;
    }
    if (val > fmax) {
      fmax = val;
    }
  }
  if (nbins < 0) {
    nbins = 1000;
  }
  h = HISTOalloc(nbins);
  if (fmax == fmin) {
    bin_size = 1;
  }
  else {
    bin_size = (fmax - fmin) / ((float)h->nbins - 1);
  }
  h->bin_size = bin_size;
  h->min = fmin;
  h->max = fmax;
  for (bin = 0; bin < h->nbins; bin++) {
    h->bins[bin] = bin * h->bin_size + h->min;
  }
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
#if GCC_VERSION > 80000
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#endif
    switch (field) {
      default:
        ErrorExit(ERROR_BADPARM, "MRISgetHistogram: unknown field %d", field);
      case SCLV_VAL:
        val = (v)->val;
        break;
      case SCLV_VAL2:
        val = v->val2;
        break;
      case SCLV_VALBAK:
        val = v->valbak;
        break;
      case SCLV_VAL2BAK:
        val = v->val2bak;
        break;
      case SCLV_VALSTAT:
        val = v->stat;
        break;
      case SCLV_IMAG_VAL:
        val = v->imag_val;
        break;
      case SCLV_MEAN:
        val = v->mean;
        break;
      case SCLV_MEAN_IMAG:
        val = v->mean_imag;
        break;
      case SCLV_STD_ERROR:
        val = v->std_error;
        break;
    }
#if GCC_VERSION > 80000
#pragma GCC diagnostic pop
#endif
    bin = HISTOvalToBinDirect(h, val);
    h->counts[bin]++;
  }
  return (h);
}


void MRISclearD(MRIS *mris)
{
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    mris->vertices[vno].d = 0.0f;
  }
}


int MRISsampleAtEachDistance(MRI_SURFACE *mris, int nbhd_size, int nbrs_per_distance)
{
  int n, nbrs_array[MAX_NBHD_SIZE];

  if (!nbhd_size) {
    return (NO_ERROR);
  }

  if (Gdiag & (DIAG_HEARTBEAT | DIAG_SHOW)) {
    fprintf(stdout, "sampling long-range distances");
  }
  for (n = 0; n <= nbhd_size; n++) {
    nbrs_array[n] = nbrs_per_distance;
  }
  return MRISsampleDistances(mris, nbrs_array, nbhd_size);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Expand the list of neighbors of each vertex, 
  reallocating the v->v array to hold the expanded list.
  
  It can expand the list way beyond the normal nsizeMax of 3
  but the ones beyond may be a subset rather than all neighbors at that distance.
  
  This also leaves v->vtotal beyond v->nsize's v->v#num
  It also computes the approx avg_nbrs
    
  ------------------------------------------------------*/
#define MAX_V 5000                                         /* max for any one node, actually way too big */
#define TRIANGLE_DISTANCE_CORRECTION 1.09f                 /*1.1f*/
/*1.066f*/ /*1.12578*/ /* 1.13105f*/                       /*1.1501f  (1.1364f)*/
#define QUADRANGLE_DISTANCE_CORRECTION ((1 + sqrt(2)) / 2) /* 1.2071  */

static int MRISsampleDistances_new(MRI_SURFACE *mris, int *nbrs, int max_nbhd, FILE* trace);

int MRISsampleDistances(MRI_SURFACE *mris, int *nbrs, int max_nbhd) {

  static int  traceCount;
  
  const char* traceFnm = 
    (++traceCount > 0) 
    ? NULL
    : "./MRISsampleDistances_new.txt";
            
  FILE* trace = traceFnm ? fopen(traceFnm, "w") : NULL;
  if (trace) fprintf(trace, "MRISsampleDistances %d max_nbhd:%d\n", traceCount, max_nbhd);
  
  int result = MRISsampleDistances_new(mris, nbrs, max_nbhd, trace);
  
  if (trace) {
    mrisDumpShape(trace, mris);
    fclose(trace);
    fprintf(stdout, "%s:%d IMMEDIATE EXIT\n",__FILE__,__LINE__);
    exit(1);
  }
  
  return result;
}

static int MRISsampleDistances_new(MRI_SURFACE *mris, int *nbrs, int max_nbhd, FILE* trace)
{
  if (trace) fprintf(trace, "MRISsampleDistances initial mris.nsize:%d\n", mris->nsize);

  mrisCheckVertexFaceTopology(mris);
  
  int diag_vno1 = -1, diag_vno2 = -1;
  {
    const char* cp;
    if ((cp = getenv("VDIAG1")) != NULL) diag_vno1 = atoi(cp);
    if ((cp = getenv("VDIAG2")) != NULL) diag_vno2 = atoi(cp);
    if (diag_vno1 >= 0) {
      printf("\nlooking for vertex pair %d, %d\n", diag_vno1, diag_vno2);
    }
  }
  
  /* adjustment for Manhattan distance */
  float const dist_scale =
    IS_QUADRANGULAR(mris) ? (1.0 + sqrt(2.0)) / 2.0f : TRIANGLE_DISTANCE_CORRECTION;

  // TBD
  //
  float c[100];
  int  nc[100];
  if (max_nbhd >= 100) {
    ErrorExit(ERROR_OUT_OF_BOUNDS,
      "MRISsampleDistances: temporaries are too small");
  }
  memset(c,  0, 100 * sizeof(float));
  memset(nc, 0, 100 * sizeof(int));
  
  // Temps needed below
  //
  VECTOR* v1 = VectorAlloc(3, MATRIX_REAL);
  VECTOR* v2 = VectorAlloc(3, MATRIX_REAL);

  // Holds the extended set that will be sampled from
  //
  int * const vnbrs = (int *)calloc(MAX_NBHD_VERTICES, sizeof(int));
  int * const vall  = (int *)calloc(MAX_NBHD_VERTICES, sizeof(int));
  int * const vnb   = (int *)calloc(MAX_NBHD_VERTICES, sizeof(int));
  if (!vnbrs || !vall || !vnb)
    ErrorExit(ERROR_NOMEMORY, "MRISsampleDistances: could not allocated v buffer");
  
  // Estimate the needs and show some stats
  //
  int max_possible = 0;     // The requested number upto and beyond them
  {
    int vtotal     = 0;     // The requested number beyond the already known neighbors
    { int n;
      for (n = 1; n <= max_nbhd; n++) {
        max_possible += nbrs[n];
        if (n > mris->nsize) {
          vtotal += nbrs[n];
        }
      }
    }
  
    if (Gdiag & DIAG_HEARTBEAT)
      fprintf(stdout,
            "\nsampling %d dists/vertex (%2.1f at each dist) = %2.1fMB\n",
            vtotal,
            (float)vtotal / ((float)max_nbhd - (float)mris->nsize),
            (float)vtotal * MRISvalidVertices(mris) * sizeof(float) * 3.0f / (1024.0f * 1024.0f));
  }
  
  //int i, n, vno, vnum, old_vnum, total_nbrs, max_possible, max_v, vtotal;
  //int found, n2, vnbrs_num, vall_num, nbhd_size, done, checks = 0;
  //float xd, yd, zd, min_dist, dist;
  //float min_angle, angle;
  //char *cp;

  bool adjusted_mris_nsize = false;
  
  int total_nbrs = 0;
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
  
    if ((Gdiag & DIAG_HEARTBEAT) && (!(vno % (mris->nvertices / 10))))
      fprintf(stdout, "%%%1.0f done\n", 100.0f * (float)vno / (float)mris->nvertices);

    VERTEX_TOPOLOGY * const vt = &mris->vertices_topology[vno];    
    VERTEX          * const v  = &mris->vertices         [vno];
    
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    if (v->ripflag) {
      continue;
    }

    /* small neighborhood is always fixed, don't overwrite them */
    
    // Expand each vertex to its nsizeMax.
    // This is assumed later.
    //
    cheapAssert(vt->nsizeMax > 0);
    MRIS_setNsizeCur(mris, vno, vt->nsizeMax);

    if (!adjusted_mris_nsize) {
      adjusted_mris_nsize = true;
      mris->nsize = vt->nsizeCur;
    } else {
      if (mris->nsize != vt->nsizeCur) {
        ErrorExit(ERROR_OUT_OF_BOUNDS,
          "MRISsampleDistances: different vertexs had different nsizeMax");
      }
    }
   
    // Expand v, dist, and dist_orig to accomodate the predicted extras
    // Zero the expansion to prevent weird variations
    //
    int const vCapacity = vt->vtotal + max_possible;
    if (trace && vno == 0) fprintf(trace, "vCapacity:%d\n", vCapacity);
    
    {
      vt->v = (int *)realloc(vt->v, vCapacity*sizeof(int));
      MRISgrowDist    (mris,vno,vCapacity);
      MRISgrowDistOrig(mris,vno,vCapacity);
      int i;
      for (i = vt->vtotal; i < vCapacity; i++) vt->v[i] = v->dist[i] = v->dist_orig[i] = 0;
    }

    /*
     find all the neighbors at each extent (i.e. 1-neighbors, then
     2-neighbors, etc..., marking their corrected edge-length distances
     as you go.
    */
    
    // Start with the vertex itself
    //
    vall[0] = vno;
    int old_vnum = 0;   // the start of the known set
    int vall_num = 1;   // the end of the known set
    v->marked = 1;      // mark it as in the set
    
    // Add each next ring of neighbors until max_nbhd is reached
    //
    int nbhd_size;
    for (nbhd_size = 1; nbhd_size <= max_nbhd; nbhd_size++) {

      /* expand neighborhood outward by a ring of vertices */
      int vnbrs_num = 0;        /* will count neighbors in this ring */
      int vnum      = vall_num; // where we are adding onto the set
      int found     = 0;        // how many got added
      
      int n;
      for (n = old_vnum; n < vall_num; n++) {
        VERTEX_TOPOLOGY const * const vnt = &mris->vertices_topology[vall[n]];
        VERTEX          const * const vn  = &mris->vertices         [vall[n]];
        
        if (vn->ripflag) {
          continue;
        }

        /* search through vn's neighbors to find an unmarked vertex */
        int n2;
        for (n2 = 0; n2 < vnt->vnum; n2++) {
          VERTEX * const vn2 = &mris->vertices[vnt->v[n2]];
          if (vn2->ripflag || vn2->marked) {
            continue;
          }

          /* found one, mark it and put it in the vall list */
          
          if (vnum >= MAX_NBHD_VERTICES) 
            ErrorExit(ERROR_OUT_OF_BOUNDS,
              "MRISsampleDistances: too many neighbors");
          
          found++;
          vn2->marked = nbhd_size;
          vnb [vnum] = vall[n];
          vall[vnum] = vnt->v[n2];
          vnum++;
          
          if (nbrs[nbhd_size] > 0) /* want to store this distance */
          {
            cheapAssert(vnbrs_num < vnum);
            vnbrs[vnbrs_num++] = vnt->v[n2];
          }
        }
      } /* done with all neighbors at previous distance */

      /* found all neighbors at this extent - calculate distances */
      
      // describe the added ring
      old_vnum = vall_num; /* old_vnum is index of 1st nbr at this distance*/
      vall_num += found;   /* vall_num is total # of nbrs */
      cheapAssert(vall_num == vnum);
      
      // process the added ring
      for (n = old_vnum; n < vall_num; n++) {
        VERTEX_TOPOLOGY const * const vnt = &mris->vertices_topology[vall[n]];
        VERTEX                * const vn  = &mris->vertices         [vall[n]];
        if (vn->ripflag) {
          continue;
        }
        
        float min_dist = UNFOUND_DIST;
        int n2;
        for (n2 = 0; n2 < vnt->vnum; n2++) {
          VERTEX const * const vn2 = &mris->vertices[vnt->v[n2]];
          if (vn2->ripflag) {
            continue;
          }
          if (!vn2->marked || vn2->marked == nbhd_size) {
            continue;
          }
          
          // BUG - THESE SHOULD BE double!  NOT CHANGED TO AVOID RESULTS CHANGING.
          //
          float xd = vn2->x - vn->x;
          float yd = vn2->y - vn->y;
          float zd = vn2->z - vn->z;
          float dist = sqrt(xd * xd + yd * yd + zd * zd);
#if !MULTI_DIST_SCALING
          if (nbhd_size > 1) {
            dist /= dist_scale;
          }
          dist += vn2->d;
#else
          dist = (dist + vn2->d * distance_scale[vn2->marked]) / distance_scale[vn2->marked + 1];
#endif
          if (min_dist > dist) {
              min_dist = dist;
          }
        }
        
        vn->d = min_dist;
        if (nbhd_size <= 2) {
          // BUG - THESE SHOULD BE double!  NOT CHANGED TO AVOID RESULTS CHANGING.
          //
          float xd = vn->x - v->x;
          float yd = vn->y - v->y;
          float zd = vn->z - v->z;
          float dist = sqrt(xd * xd + yd * yd + zd * zd);
          vn->d = dist;
        }
        
        if (vn->d >= UNFOUND_DIST / 2) {
          printf(
              "***** WARNING - surface distance not found at "
              "vno %d, vall[%d] = %d (vnb[%d] = %d ******",
              vno,
              n,
              vall[n],
              n,
              vnb[n]);
          DiagBreak();
          exit(1);
        }
        
        if ((vall[n] == diag_vno1 && vno == diag_vno2) || (vall[n] == diag_vno2 && vno == diag_vno1))
          printf("vn %d is %2.3f mm from v %d at ring %d\n", vall[n], vn->d, vno, nbhd_size);
      }

      /*
       now check each to see if a neighbor at the same 'distance'
       is actually closer than neighbors which are 'nearer' (i.e. maybe
       the path through a 3-neighbor is shorter than that through any
       of the 2-neighbors.
      */
      for (n = old_vnum; n < vall_num; n++) {
        VERTEX_TOPOLOGY const * const vnt = &mris->vertices_topology[vall[n]];
        VERTEX                * const vn  = &mris->vertices         [vall[n]];
        if (vn->ripflag) {
          continue;
        }

        float min_dist = vn->d;
        
        int n2;
        for (n2 = 0; n2 < vnt->vnum; n2++) {
          VERTEX const * const vn2 = &mris->vertices[vnt->v[n2]];
          if (!vn2->marked || vn2->marked != nbhd_size || vn2->ripflag) {
            continue;
          }

          // BUG - THESE SHOULD BE double!  NOT CHANGED TO AVOID RESULTS CHANGING.
          //
          float xd = vn2->x - vn->x;
          float yd = vn2->y - vn->y;
          float zd = vn2->z - vn->z;
          float dist = sqrt(xd * xd + yd * yd + zd * zd);
#if !MULTI_DIST_SCALING
          if (nbhd_size > 1) {
            dist /= dist_scale;
          }
          if (vn2->d + dist < min_dist) {
            min_dist = vn2->d + dist;
          }
#else
          dist = (dist + vn2->d * distance_scale[vn2->marked]) / distance_scale[vn2->marked + 1];
          if (dist < min_dist) {
            min_dist = dist;
          }
#endif
        }
        vn->d = min_dist;
        if ((min_dist > 5 * nbhd_size) || min_dist > 60) {
          DiagBreak();
        }
        {
          // BUG - THESE SHOULD BE double!  NOT CHANGED TO AVOID RESULTS CHANGING.
          //
          float xd = vn->x - v->x;
          float yd = vn->y - v->y;
          float zd = vn->z - v->z;
          float dist = sqrt(xd * xd + yd * yd + zd * zd);
          c [nbhd_size] += vn->d / dist;
          nc[nbhd_size]++;
        }
      }

      /* if this set of neighbors are to be stored, sample from them */
      if (nbrs[nbhd_size] <= 0) {
        continue;
      }

      /* make sure the points are not too close together */
      float min_angle = 0.9 * 2.0 * M_PI / (float)nbrs[nbhd_size];

      /*
       at this point, the vall list contains all the
       neighbors currently found
       at ALL distances, while the vnbrs list contains ONLY the
       nbhd_size-neighbors.
      */
      if (trace && vno == 0) fprintf(trace, "nbhd_size:%d found:%d\n", nbhd_size, found);

      if (found <= nbrs[nbhd_size]) /* just copy them all in */
      {
        // THIS IS WEIRD BECAUSE IT IS PUTTING REPEATS OF THE 1..v->nsize rings
        // ONTO THE END OF THE LIST
        //
        for (n = 0; n < found; n++, vt->vtotal++) {
          cheapAssert(vt->vtotal < vCapacity);
          vt->v       [vt->vtotal] =                vnbrs[n];
          v->dist_orig[vt->vtotal] = mris->vertices[vnbrs[n]].d;
          if (v->dist_orig[vt->vtotal] > 60) {
            DiagBreak();
          }
          if (trace && vno == 0) fprintf(trace, "copy dist_orig[%d]:%f\n", vt->vtotal, v->dist_orig[vt->vtotal]);
        }
      }
      else /* randomly sample from them */
      {
        int vstart = vt->vtotal;
        for (n = 0; n < nbrs[nbhd_size]; n++, vt->vtotal++) {
          int niter = 0;
          int done;
          int i;
          do {
            do {
              i = nint(randomNumber(0.0, (double)found - 1));
            } while (vnbrs[i] < 0);
            /*
             now check to make sure that the angle between this
             point and the others already selected is not too
             small to make sure the points are not bunched.
            */
            VERTEX const * const vn = &mris->vertices[vnbrs[i]];
            VECTOR_LOAD(v1, vn->x - v->x, vn->y - v->y, vn->z - v->z);
            done = 1;
            int j;
            for (j = vstart; done && j < vt->vtotal; j++) {
              VERTEX const * const vn2 = &mris->vertices[vt->v[j]];
              VECTOR_LOAD(v2, vn2->x - v->x, vn2->y - v->y, vn2->z - v->z);
              float angle = Vector3Angle(v1, v2);
              if (angle < min_angle) {
                done = 0;
              }
            }
            if (++niter > found) /* couldn't find enough at this difference */
            {
              min_angle *= 0.75f; /* be more liberal */
              niter = 0;
            }
          } while (!done && !FZERO(min_angle));
          
          VERTEX const * const vn = &mris->vertices[vnbrs[i]];
          cheapAssert(vt->vtotal < vCapacity);
          vt->v       [vt->vtotal] = vnbrs[i];
          v->dist_orig[vt->vtotal] = vn->d;
          if (v->dist_orig[vt->vtotal] > 60) {
            DiagBreak();
          }
          if (FZERO(vn->d)) {
            DiagBreak();
          }
          vnbrs[i] = -1;
          if (trace && vno == 0) fprintf(trace, "rand dist_orig[%d]:%f\n", vt->vtotal, v->dist_orig[vt->vtotal]);
        }
      }
    }

    if ((Gdiag_no == vno) && DIAG_VERBOSE_ON) {
      FILE *fp;
      char fname[STRLEN];
      int n;
      
      sprintf(fname, "v%d", vno);
      fp = fopen(fname, "w");
      fprintf(fp, "%d\n", vall_num);
      for (n = 0; n < vall_num; n++) {
        fprintf(fp, "%d\n", vall[n]);
      }
      fclose(fp);

      sprintf(fname, "vn%d", vno);
      fp = fopen(fname, "w");
      fprintf(fp, "%d\n", vt->vtotal);
      for (n = 0; n < vt->vtotal; n++) {
        fprintf(fp, "%d\n", vt->v[n]);
      }
      fclose(fp);
      
      for (n = 0; n < mris->nvertices; n++) {
        VERTEX * const vn = &mris->vertices[n];
#if 0
        if (vn->ripflag)
        {
          continue ;
        }
#endif
        vn->curv = vn->d;
      }
      
      sprintf(fname, "%s.dist", mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh");
      MRISwriteCurvature(mris, fname);
    }

    /*
     done building arrays - allocate distance vectors and
     sample from the found neighbors list.
    */
    /* now unmark them all */
    int n;
    for (n = 0; n < vall_num; n++) {
      mris->vertices[vall[n]].marked = 0;
      mris->vertices[vall[n]].d = 0.0;
    }

    total_nbrs += vt->vtotal;
  }

  /* now fill in immediate neighborhood(Euclidean) distances */
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];    
    VERTEX          const * const v  = &mris->vertices         [vno];

    if (vno == Gdiag_no) {
      DiagBreak();
    }
    if (v->ripflag) {
      continue;
    }

    // The above code made nsizeCur == nsizeMax, but vtotal may be beyond these
    // so must determine the nsizeMax v#num
    //
    int immediateCount = 0;
    if (vt->nsizeMax == 3) {
      immediateCount = vt->v3num;
    }
    else if (vt->nsizeMax == 2) {
      immediateCount = vt->v2num;
    }
    else {
      immediateCount = vt->vnum;
    }
    
    int n;
    for (n = 0; n < immediateCount; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      if (vn->ripflag) {
        continue;
      }
      
      // BUG - THESE SHOULD BE double!  NOT CHANGED TO AVOID RESULTS CHANGING.
      //
      float xd = v->x - vn->x;
      float yd = v->y - vn->y;
      float zd = v->z - vn->z;

      // WEIRD - WHY IS DIST NOT CHANGED?
      //
      v->dist_orig[n] = sqrt(xd * xd + yd * yd + zd * zd);
      if (trace && vno == 0) fprintf(trace, "eucliod dist_orig[%d]:%f\n", n, v->dist_orig[n]);
    }
  }

  mris->vtotalsMightBeTooBig = 1;   // TODO find the right place to turn this off again

  mrisCheckVertexFaceTopology(mris);
  
  // make sure distances are symmetric
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];    
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    int n;
    for (n = 0; n < vt->vtotal; n++) {
      VERTEX_TOPOLOGY const * const vnt = &mris->vertices_topology[vt->v[n]];
      VERTEX                * const vn  = &mris->vertices         [vt->v[n]];
      if (vn->ripflag) {
        continue;
      }
      int i;
      for (i = 0; i < vnt->vtotal; i++) {
        if (vnt->v[i] == vno)  // distance in both lists - make it the average
        {
          double dist = (vn->dist_orig[i] + v->dist_orig[n]) / 2;
          vn->dist_orig[i] = dist;
          v ->dist_orig[n] = dist;
          break;
        }
      }
    }
  }

  /* check reasonableness of distances */
  size_t zeroesCount    = 0;
  size_t nonZeroesCount = 0;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];    
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    int n;
    for (n = 0; n < vt->vtotal; n++) {
      if (!DZERO(v->dist_orig[n])) {
        if (!nonZeroesCount && vno != 0) fprintf(stderr, "first zero distance at v %d, n %d (vn = %d)\n", vno, n, vt->v[n]);
        nonZeroesCount++;
      } else {
        zeroesCount++;
        if      (zeroesCount <  50) fprintf(stderr, "zero distance at v %d, n %d (vn = %d)\n", vno, n, vt->v[n]);
        else if (zeroesCount == 50) fprintf(stderr, "further zeroes elided\n");
      }
    }
  }
  if (zeroesCount > 0)
    fprintf(stderr, "Of %ld neighbours, %ld are at distance zero\n", zeroesCount + nonZeroesCount, zeroesCount);

  if (Gdiag_no >= 0) {
    char fname[STRLEN];

    sprintf(fname, "v%d.log", Gdiag_no);
    FILE *fp = fopen(fname, "w");

    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[Gdiag_no];    
    VERTEX          const * const v  = &mris->vertices         [Gdiag_no];

    int i;
    for (i = 0; i < vt->vtotal; i++) {
      fprintf(fp, "%d: %d, %f\n", i, vt->v[i], v->dist_orig[i]);
    }
    fclose(fp);
  }

  mris->avg_nbrs = (float)total_nbrs / (float)MRISvalidVertices(mris);

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "avg_nbrs = %2.1f\n", mris->avg_nbrs);
  }

#if MULTI_DIST_SCALING
  if (Gdiag & DIAG_SHOW) {
    int n;
    for (n = 0; n <= max_nbhd; n++) {
      if (nc[n]) {
        c[n] /= (float)nc[n];
      }
      fprintf(stdout, "c[%d] = %2.5f (%d samples)\n", n, c[n], nc[n]);
    }
    fprintf(stdout, "c[] = { ");
    for (n = 0; n <= max_nbhd; n++) {
      fprintf(stdout, "%2.5f", c[n]);
      if (n < max_nbhd) {
        fprintf(stdout, ", ");
      }
    }
  }
#endif
  free(vnbrs);
  free(vall);
  free(vnb);
  VectorFree(&v1);
  VectorFree(&v2);
  if (Gdiag & DIAG_HEARTBEAT) {
    fprintf(stdout, " done.\n");
  }

  return (NO_ERROR);
}


// Some general purpose operations
//
void MRISclearWhichAndVal2(MRIS *mris, int which)
{
  switch (which) {
    case VERTEX_AREA:
      MRISclearOrigAreaAndVal2(mris);
      return;
    case VERTEX_CURV:
      MRISclearCurvAndVal2(mris);
      return;
    case VERTEX_LOGODDS:
    case VERTEX_VAL:
      mrisSetValAndClearVal2(mris, 0.0f);
      return;
  }

  int vno;

  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX *v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    switch (which) {
      default:
        cheapAssert(false);
    }
    v->val2 = 0;            // surprise!
  }
}
