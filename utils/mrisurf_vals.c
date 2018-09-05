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

