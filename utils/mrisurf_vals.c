/*
 * @file utilities operating on Original
 *
 */
/*
 * surfaces Author: Bruce Fischl, extracted from mrisurf.c by Bevin Brett
 *
 * $ © copyright-2014,2018 The General Hospital Corporation (Boston, MA) "MGH"
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

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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
    v->dx = 0;
    v->dy = 0;
    v->dz = 0;
  }
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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

/*! ----------------------------------------------------
  \fn int mrisClearMomentum(MRI_SURFACE *mris)
  \brief sets v->od{xyz}=0 for unripped vertices.
  ------------------------------------------------------*/
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
    v->odx = 0;
    v->ody = 0;
    v->odz = 0;
  }
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/

int
MRISmarkVerticesWithValOverThresh(MRI_SURFACE *mris, float thresh)
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


#if 0
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static int
mrisCheck(MRI_SURFACE *mris)
{
  int       vno ;
  VERTEX    *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    if (fabs(v->curv) > 10000.0)
    {
      DiagBreak() ;
      return(ERROR_BADPARM) ;
    }
  }
  return(NO_ERROR) ;
}
#endif
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISminFilterCurvatures(MRI_SURFACE *mris, int niter)
{
  int i, vno, vnb, *pnb, vnum;
  float curv;

  for (i = 0; i < niter; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) {
        continue;
      }
      curv = v->curv;
      pnb  = vt->v;
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISmaxFilterCurvatures(MRI_SURFACE *mris, int niter)
{
  int i, vno, vnb, *pnb, vnum;
  float curv;

  for (i = 0; i < niter; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) {
        continue;
      }
      curv = v->curv;
      pnb  = vt->v;
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
  int i, vno, vnb, *pnb, vnum;
  float curv, num;

  for (i = 0; i < navgs; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) {
        continue;
      }
      curv = v->curv;
      pnb  = vt->v;
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISaverageMarkedCurvatures(MRI_SURFACE *mris, int navgs)
{
  int i, vno, vnb, *pnb, vnum;
  float curv, num;

  for (i = 0; i < navgs; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag || !v->marked) {
        continue;
      }
      curv = v->curv;
      pnb  = vt->v;
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISaverageMarkedVals(MRI_SURFACE *mris, int navgs)
{
  int i, vno, vnb, *pnb, vnum;
  float val, num;

  for (i = 0; i < navgs; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag || v->marked == 0) {
        continue;
      }
      val = v->val;
      pnb  = vt->v;
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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
      
      int vnb, *pnb, vnum;
      float val, num;

      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) ROMP_PF_continue;

      val = v->val;
      pnb  = vt->v;
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  spatially average the v->d field
  ------------------------------------------------------*/
int MRISaverageD(MRI_SURFACE *mris, int navgs)
{
  int i, vno, vnb, *pnb, vnum;
  float val, num;

  for (i = 0; i < navgs; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) {
        continue;
      }
      val = v->d;
      pnb  = vt->v;
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

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISmedianFilterVals(MRI_SURFACE *mris, int nmedians)
{
  int i, vno, vnb, *pnb, vnum, num;
  float val_list[MAX_NEIGHBORS];

  for (i = 0; i < nmedians; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) {
        continue;
      }
      pnb = vt->v;
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
      qsort(val_list, num, sizeof(val_list[0]), compare_sort_vals);
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISmedianFilterVerexPositions(MRI_SURFACE *mris, int nmedians)
{
  int i, vno, vnb, *pnb, vnum, num;
  float val_list[MAX_NEIGHBORS];

  for (i = 0; i < nmedians; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) {
        continue;
      }
      pnb  = vt->v;
      vnum = vt->vnum;
      val_list[0] = v->x;
      for (num = 1, vnb = 0; vnb < vnum; vnb++) {
        VERTEX const * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
        if (vn->ripflag) {
          continue;
        }

        val_list[num++] = vn->x;
      }
      qsort(val_list, num, sizeof(val_list[0]), compare_sort_vals);
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
      qsort(val_list, num, sizeof(val_list[0]), compare_sort_vals);
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
      qsort(val_list, num, sizeof(val_list[0]), compare_sort_vals);
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISmedianFilterD(MRI_SURFACE *mris, int nmedians, int vtotal)
{
  int i, vno, vnb, *pnb, vnum, num;
  float val_list[MAX_NEIGHBORS];

  for (i = 0; i < nmedians; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) {
        continue;
      }
      pnb = vt->v;
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
      qsort(val_list, num, sizeof(val_list[0]), compare_sort_vals);
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISmedianFilterCurvature(MRI_SURFACE *mris, int nmedians)
{
  int i, vno, vnb, *pnb, vnum, num;
  float val_list[MAX_NEIGHBORS];

  for (i = 0; i < nmedians; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) {
        continue;
      }
      pnb  = vt->v;
      vnum = vt->vnum;
      val_list[0] = v->curv;
      for (num = 1, vnb = 0; vnb < vnum; vnb++) {
        VERTEX const * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
        if (vn->ripflag) {
          continue;
        }

        val_list[num++] = vn->curv;
      }
      qsort(val_list, num, sizeof(val_list[0]), compare_sort_vals);
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISmedianFilterVal2s(MRI_SURFACE *mris, int nmedians)
{
  int i, vno, vnb, *pnb, vnum, num;
  float val_list[MAX_NEIGHBORS];

  for (i = 0; i < nmedians; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) {
        continue;
      }
      pnb  = vt->v;
      vnum = vt->vnum;
      val_list[0] = v->val2;
      for (num = 1, vnb = 0; vnb < vnum; vnb++) {
        VERTEX const * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
        if (vn->ripflag) {
          continue;
        }

        val_list[num++] = vn->val2;
      }
      qsort(val_list, num, sizeof(val_list[0]), compare_sort_vals);
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISmedianFilterVal2baks(MRI_SURFACE *mris, int nmedians)
{
  int i, vno, vnb, *pnb, vnum, num;
  float val_list[MAX_NEIGHBORS];

  for (i = 0; i < nmedians; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) {
        continue;
      }
      pnb  = vt->v;
      vnum = vt->vnum;
      val_list[0] = v->val2bak;
      for (num = 1, vnb = 0; vnb < vnum; vnb++) {
        VERTEX const * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
        if (vn->ripflag) {
          continue;
        }

        val_list[num++] = vn->val2bak;
      }
      qsort(val_list, num, sizeof(val_list[0]), compare_sort_vals);
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISaverageVal2s(MRI_SURFACE *mris, int navgs)
{
  int i, vno, vnb, *pnb, vnum;
  float val, num;

  for (i = 0; i < navgs; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) {
        continue;
      }
      val  = v->val2;
      pnb  = vt->v;
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  actually these are squared standard errors
  ------------------------------------------------------*/
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


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  actually these are squared standard errors
  ------------------------------------------------------*/
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

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISaverageVal2baks(MRI_SURFACE *mris, int navgs)
{
  int i, vno, vnb, *pnb, vnum;
  float val, num;

  for (i = 0; i < navgs; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) {
        continue;
      }
      val  = v->val2bak;
      pnb  = vt->v;
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


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
   copy the v->val2bak field into the v->val for every vertex
  ------------------------------------------------------*/
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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
  Parameters:

  Returns value:

  Description
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

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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



/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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
/*-----------------------------------------------------

Parameters:

Returns value:

Description
------------------------------------------------------*/
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
/*-----------------------------------------------------

Parameters:

Returns value:

Description
------------------------------------------------------*/
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

int MRISextractCurvatureVector(MRI_SURFACE *mris, float *curvs)
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
    bin = HISTOvalToBinDirect(h, val);
    h->counts[bin]++;
  }
  return (h);
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

/*-----------------------------------------------------
  Parameters:

    Returns value:

    Description
    ------------------------------------------------------*/
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRIScanonicalToWorld(MRI_SURFACE *mris, double phi, double theta, double *pxw, double *pyw, double *pzw)
{
  double x, y, z, radius;

  radius = mris->radius;
  *pxw = x = radius * sin(phi) * cos(theta);
  *pyw = y = radius * sin(phi) * sin(theta);
  *pzw = z = radius * cos(phi);
  return (NO_ERROR);
}




/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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
  Expand the list of neighbors of each vertex, reallocating
  the v->v array to hold the expanded list.
  ------------------------------------------------------*/
#define MAX_V 5000                                         /* max for any one node, actually way too big */
#define TRIANGLE_DISTANCE_CORRECTION 1.09f                 /*1.1f*/
/*1.066f*/ /*1.12578*/ /* 1.13105f*/                       /*1.1501f  (1.1364f)*/
#define QUADRANGLE_DISTANCE_CORRECTION ((1 + sqrt(2)) / 2) /* 1.2071  */

static int MRISsampleDistances_old(MRI_SURFACE *mris, int *nbrs, int max_nbhd);
static int MRISsampleDistances_new(MRI_SURFACE *mris, int *nbrs, int max_nbhd);

int MRISsampleDistances(MRI_SURFACE *mris, int *nbrs, int max_nbhd) {
    static bool later_time, use_old;
    if (!later_time) {
      later_time = false;
      use_old = getenv("USE_OLD_MRISsampleDistances");
    }
   if (use_old) return MRISsampleDistances_old(mris, nbrs, max_nbhd);
   else         return MRISsampleDistances_new(mris, nbrs, max_nbhd);
}

static int MRISsampleDistances_old(MRI_SURFACE *mris, int *nbrs, int max_nbhd)
{
  int i, n, vno, vnum, old_vnum, total_nbrs, max_possible, max_v, vtotal;
  int *vnbrs, *vall, *vnb, found, n2, vnbrs_num, vall_num, nbhd_size, done, checks = 0;
  float xd, yd, zd, min_dist, dist, dist_scale, *old_dist;
  float min_angle, angle;
  VECTOR *v1, *v2;
  float c[100], *orig_old_dist;
  int nc[100], *old_v;
  int diag_vno1, diag_vno2;
  char *cp;

  orig_old_dist = old_dist = (float *)calloc(MAX_V, sizeof(float));
  old_v = (int *)calloc(MAX_V, sizeof(int));
  if (old_dist == NULL || old_v == NULL)
    ErrorExit(ERROR_NOMEMORY, "MRISsampleDistances: could not allocated old_dist buffer");

  if ((cp = getenv("VDIAG1")) != NULL) {
    diag_vno1 = atoi(cp);
  }
  else {
    diag_vno1 = -1;
  }
  if ((cp = getenv("VDIAG2")) != NULL) {
    diag_vno2 = atoi(cp);
  }
  else {
    diag_vno2 = -1;
  }
  if (diag_vno1 >= 0) {
    printf("\nlooking for vertex pair %d, %d\n", diag_vno1, diag_vno2);
  }

  memset(c, 0, 100 * sizeof(float));
  memset(nc, 0, 100 * sizeof(int));
  v1 = VectorAlloc(3, MATRIX_REAL);
  v2 = VectorAlloc(3, MATRIX_REAL);

  /* adjust for Manhattan distance */
  if (IS_QUADRANGULAR(mris)) {
    dist_scale = (1.0 + sqrt(2.0)) / 2.0f;
  }
  else {
    dist_scale = TRIANGLE_DISTANCE_CORRECTION;
  }

  vnbrs = (int *)calloc(MAX_NBHD_VERTICES, sizeof(int));
  vall = (int *)calloc(MAX_NBHD_VERTICES, sizeof(int));
  vnb = (int *)calloc(MAX_NBHD_VERTICES, sizeof(int));
  vtotal = total_nbrs = 0;
  for (vtotal = max_possible = 0, n = 1; n <= max_nbhd; n++) {
    max_possible += nbrs[n];
    if (n > mris->nsize) {
      vtotal += nbrs[n];
    }
  }
  if (Gdiag & DIAG_HEARTBEAT)
    fprintf(stdout,
            "\nsampling %d dists/vertex (%2.1f at each dist) = %2.1fMB\n",
            vtotal,
            (float)vtotal / ((float)max_nbhd - (float)mris->nsize),
            (float)vtotal * MRISvalidVertices(mris) * sizeof(float) * 3.0f / (1024.0f * 1024.0f));

  for (vno = 0; vno < mris->nvertices; vno++) {
  
    if ((Gdiag & DIAG_HEARTBEAT) && (!(vno % (mris->nvertices / 10))))
      fprintf(stdout, "%%%1.0f done\n", 100.0f * (float)vno / (float)mris->nvertices);
    if ((vno > 139000 || (!(vno % 100))) && 0) {
      if (checks++ == 0) {
        printf("checking surface at vno %d\n", vno);
      }
    }
    VERTEX_TOPOLOGY * const vt = &mris->vertices_topology[vno];    
    VERTEX          * const v  = &mris->vertices         [vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    if (v->ripflag) {
      continue;
    }

    /* small neighborhood is always fixed, don't overwrite them */
    vtotal = vt->vtotal;
    if (vt->nsizeMax == 3) {
      vt->vtotal = vt->v3num;
    }
    else if (vt->nsizeMax == 2) {
      vt->vtotal = vt->v2num;
    }
    else {
      vt->vtotal = vt->vnum;
    }

    max_v = vt->vtotal + max_possible;
    if (vtotal < max_v) /* won't fit in current allocation,
                         reallocate stuff */
    {
      /* save and restore neighbor list */
      memmove(old_v, vt->v, vt->vtotal * sizeof(vt->v[0]));
      free(vt->v);
      vt->v = (int *)calloc(max_v, sizeof(int));
      if (!vt->v)
        ErrorExit(ERROR_NO_MEMORY,
                  "MRISsampleDistances: could not allocate list of %d "
                  "nbrs at v=%d",
                  max_v,
                  vno);
      memmove(vt->v, old_v, vt->vtotal * sizeof(vt->v[0]));

      /* save and restore distance vector */
      if (vt->vtotal >= MAX_V || old_dist != orig_old_dist)
        printf("!!!!!!!!!!!!! v %d has too many (%d) vertices to save distances for !!!!!!!!!!!!!!!!!\n", vno, vtotal);
      memmove(old_dist, v->dist, vt->vtotal * sizeof(v->dist[0]));
      free(v->dist);
      v->dist = (float *)calloc(max_v, sizeof(float));
      if (!v->dist)
        ErrorExit(ERROR_NOMEMORY,
                  "MRISsampleDistances: could not allocate list of %d "
                  "dists at v=%d",
                  max_v,
                  vno);
      memmove(v->dist, old_dist, vt->vtotal * sizeof(v->dist[0]));

      /* save and restore original distance vector */
      {
        int k;
        for (k = 0; k < vt->vtotal; k++) {
          old_dist[k] = 0;
          v->dist_orig[k] *= 1;
          old_dist[k] = v->dist_orig[k];
        }
        if (vno == Gdiag_no) DiagBreak();
      }
      memmove(old_dist, v->dist_orig, vt->vtotal * sizeof(v->dist_orig[0]));
      free(v->dist_orig);
      v->dist_orig = (float *)calloc(max_v, sizeof(float));
      if (!v->dist_orig)
        ErrorExit(ERROR_NOMEMORY,
                  "MRISsampleDistances: could not allocate list of %d "
                  "dists at v=%d",
                  max_v,
                  vno);
      memmove(v->dist_orig, old_dist, vt->vtotal * sizeof(v->dist_orig[0]));
    }

    if ((vno > 139000 || !(vno % 100)) && 0) {
      if (checks++ == 0) {
        printf("checking surface at vno %d\n", vno);
      }
    }
    /*
     find all the neighbors at each extent (i.e. 1-neighbors, then
     2-neighbors, etc..., marking their corrected edge-length distances
     as you go.
    */
    vall[0] = vno;
    vall_num = 1;
    old_vnum = 0;
    v->marked = 1; /* a hack - it is a zero neighbor */
    for (nbhd_size = 1; vall_num < MAX_NBHD_VERTICES && nbhd_size <= max_nbhd; nbhd_size++) {

      /* expand neighborhood outward by a ring of vertices */
      vnbrs_num = 0; /* will count neighbors in this ring */
      vnum = vall_num;
      for (found = 0, n = old_vnum; vall_num < MAX_NBHD_VERTICES && n < vall_num; n++) {
        VERTEX_TOPOLOGY const * const vnt = &mris->vertices_topology[vall[n]];
        VERTEX          const * const vn  = &mris->vertices         [vall[n]];
        if (vn->ripflag) {
          continue;
        }

        /* search through vn's neighbors to find an unmarked vertex */
        for (n2 = 0; n2 < vnt->vnum; n2++) {
          VERTEX * const vn2 = &mris->vertices[vnt->v[n2]];
          if (vn2->ripflag || vn2->marked) {
            continue;
          }

          /* found one, mark it and put it in the vall list */
          found++;
          vn2->marked = nbhd_size;
          vnb[vnum] = vall[n];
          vall[vnum++] = vnt->v[n2];
          if (nbrs[nbhd_size] > 0) /* want to store this distance */
          {
            vnbrs[vnbrs_num++] = vnt->v[n2];
          }
        }
      } /* done with all neighbors at previous distance */

      /* found all neighbors at this extent - calculate distances */
      old_vnum = vall_num; /* old_vnum is index of 1st
                          nbr at this distance*/
      vall_num += found;   /* vall_num is total # of nbrs */
      for (n = old_vnum; n < vall_num; n++) {
        VERTEX_TOPOLOGY const * const vnt = &mris->vertices_topology[vall[n]];
        VERTEX                * const vn  = &mris->vertices         [vall[n]];
        if (vn->ripflag) {
          continue;
        }
        for (min_dist = UNFOUND_DIST, n2 = 0; n2 < vnt->vnum; n2++) {
          VERTEX const * const vn2 = &mris->vertices[vnt->v[n2]];
          if (vn2->ripflag) {
            continue;
          }
          if (!vn2->marked || vn2->marked == nbhd_size) {
            continue;
          }
          xd = vn2->x - vn->x;
          yd = vn2->y - vn->y;
          zd = vn2->z - vn->z;
          dist = sqrt(xd * xd + yd * yd + zd * zd);
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
        if (nbhd_size <= 2) {
          xd = vn->x - v->x;
          yd = vn->y - v->y;
          zd = vn->z - v->z;
          dist = sqrt(xd * xd + yd * yd + zd * zd);
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
        min_dist = vn->d;
        for (n2 = 0; n2 < vnt->vnum; n2++) {
          VERTEX const * const vn2 = &mris->vertices[vnt->v[n2]];
          if (!vn2->marked || vn2->marked != nbhd_size || vn2->ripflag) {
            continue;
          }
          xd = vn2->x - vn->x;
          yd = vn2->y - vn->y;
          zd = vn2->z - vn->z;
          dist = sqrt(xd * xd + yd * yd + zd * zd);
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
          xd = vn->x - v->x;
          yd = vn->y - v->y;
          zd = vn->z - v->z;
          dist = sqrt(xd * xd + yd * yd + zd * zd);
          c[nbhd_size] += vn->d / dist;
          nc[nbhd_size]++;
        }
      }

      /* if this set of neighbors are to be stored, sample from them */
      if (nbrs[nbhd_size] <= 0) {
        continue;
      }

      /* make sure the points are not too close together */
      min_angle = 0.9 * 2.0 * M_PI / (float)nbrs[nbhd_size];

      /*
       at this point, the vall list contains all the
       neighbors currently found
       at ALL distances, while the vnbrs list contains ONLY the
       nbhd_size-neighbors.
      */
      if (found <= nbrs[nbhd_size]) /* just copy them all in */
      {
        for (n = 0; n < found; n++, vt->vtotal++) {
          vt->v[vt->vtotal] = vnbrs[n];
          v->dist_orig[vt->vtotal] = mris->vertices[vnbrs[n]].d;
          if (v->dist_orig[vt->vtotal] > 60) {
            DiagBreak();
          }
        }
      }
      else /* randomly sample from them */
      {
        int vstart = vt->vtotal;
        for (n = 0; n < nbrs[nbhd_size]; n++, vt->vtotal++) {
          int j, niter = 0;
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
            for (j = vstart; done && j < vt->vtotal; j++) {
              VERTEX const * const vn2 = &mris->vertices[vt->v[j]];
              VECTOR_LOAD(v2, vn2->x - v->x, vn2->y - v->y, vn2->z - v->z);
              angle = Vector3Angle(v1, v2);
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
          vt->v[vt->vtotal] = vnbrs[i];
          v->dist_orig[vt->vtotal] = vn->d;
          if (v->dist_orig[vt->vtotal] > 60) {
            DiagBreak();
          }
          if (FZERO(vn->d)) {
            DiagBreak();
          }
          vnbrs[i] = -1;
        }
      }
    }

    if ((vno > 9.0 * mris->nvertices / 10.0) && 0) {
      if (checks++ == 0) {
        printf("checking surface at vno %d\n", vno);
      }
    }

    if ((Gdiag_no == vno) && DIAG_VERBOSE_ON) {
      FILE *fp;
      char fname[STRLEN];

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
    if (vt->nsizeMax == 3) {
      vtotal = vt->v3num;
    }
    else if (vt->nsizeMax == 2) {
      vtotal = vt->v2num;
    }
    else {
      vtotal = vt->vnum;
    }
    for (n = 0; n < vtotal; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      if (vn->ripflag) {
        continue;
      }
      xd = v->x - vn->x;
      yd = v->y - vn->y;
      zd = v->z - vn->z;
      v->dist_orig[n] = sqrt(xd * xd + yd * yd + zd * zd);
    }
  }

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
    for (n = 0; n < vt->vtotal; n++) {
      VERTEX_TOPOLOGY const * const vnt = &mris->vertices_topology[vt->v[n]];
      VERTEX                * const vn  = &mris->vertices         [vt->v[n]];
      if (vn->ripflag) {
        continue;
      }
      for (i = 0; i < vnt->vtotal; i++) {
        if (vnt->v[i] == vno)  // distance in both lists - make it the average
        {
          double dist;
          dist = (vn->dist_orig[i] + v->dist_orig[n]) / 2;
          vn->dist_orig[i] = dist;
          v->dist_orig[n] = dist;
          break;
        }
      }
    }
  }

  /* check reasonableness of distances */
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];    
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    for (n = 0; n < vt->vtotal; n++) {
      if (DZERO(v->dist_orig[n])) fprintf(stderr, "zero distance at v %d, n %d (vn = %d)\n", vno, n, vt->v[n]);
    }
  }

  if (Gdiag_no >= 0) {
    FILE *fp;
    char fname[STRLEN];
    int i;

    sprintf(fname, "v%d.log", Gdiag_no);
    fp = fopen(fname, "w");

    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[Gdiag_no];    
    VERTEX          const * const v  = &mris->vertices         [Gdiag_no];
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
  free(old_dist);
  free(old_v);
  VectorFree(&v1);
  VectorFree(&v2);
  if (Gdiag & DIAG_HEARTBEAT) {
    fprintf(stdout, " done.\n");
  }

  return (NO_ERROR);
}

static int MRISsampleDistances_new(MRI_SURFACE *mris, int *nbrs, int max_nbhd)
{
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
    if (vt->nsizeMax == 3) {
      vt->vtotal = vt->v3num;
    }
    else if (vt->nsizeMax == 2) {
      vt->vtotal = vt->v2num;
    }
    else { cheapAssert(vt->nsizeMax == 1);
      vt->vtotal = vt->vnum;
    }
    vt->nsizeCur = vt->nsizeMax;

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
    {
      vt->v        = (int  *)realloc(vt->v,        vCapacity*sizeof(int  ));
      v->dist      = (float*)realloc(v->dist,      vCapacity*sizeof(float));
      v->dist_orig = (float*)realloc(v->dist_orig, vCapacity*sizeof(float));
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

    // The above code made nsizeCur == nsizeMax
    //
    // cheapAssert(vt->nsizeCur == vt->nsizeMax);
    int n;
    for (n = 0; n < vt->vtotal; n++) {
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
    }
  }

  // mris->vtotalsMightBeTooBig = 1;   // TODO find the right place to turn this off again

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
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];    
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    int n;
    for (n = 0; n < vt->vtotal; n++) {
      if (DZERO(v->dist_orig[n])) fprintf(stderr, "zero distance at v %d, n %d (vn = %d)\n", vno, n, vt->v[n]);
    }
  }

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
