#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "volume_io.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "mrisurf.h"
#include "label.h"
#include "mri.h"

static Transform *labelLoadTransform(char *subject_name, char *sdir,
                                     General_transform *transform) ;

#define MAX_VERTICES 500000
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
LABEL *
LabelRead(char *subject_name, char *label_name)
{
  LABEL  *area ;
  char   fname[STRLEN], *cp, line[STRLEN], subjects_dir[STRLEN], lname[STRLEN];
  FILE   *fp ;
  int    vno, nlines ;
  float  x, y, z, stat ;

  area = (LABEL *)calloc(1, sizeof(LABEL)) ;
  if (!area)
    ErrorExit(ERROR_NOMEMORY,"%s: could not allocate LABEL struct.",Progname);

  if (subject_name && !strchr(label_name, '/'))
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, 
                "%s: no subject's directory specified in environment "
                "(SUBJECTS_DIR)", Progname) ;
    strcpy(subjects_dir, cp) ;
    strcpy(lname, label_name) ;
    cp = strstr(lname, ".label") ;
    if (cp)
      *cp = 0 ;
    cp = strrchr(lname, '/') ;
    if (cp)
      label_name = cp+1 ;
    else
      label_name = lname ;
    sprintf(fname, "%s/%s/label/%s.label", subjects_dir,subject_name,
            label_name);
    strcpy(area->subject_name, subject_name) ;
  }
  else
  {
    strcpy(fname, label_name) ;
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      strcpy(subjects_dir, ".") ;
    else
      strcpy(subjects_dir, cp) ;
    strcpy(lname, label_name) ;
    cp = strstr(lname, ".label") ;
    if (cp == NULL)
      sprintf(fname, "%s.label", lname);
    else
      strcpy(fname, label_name) ;
  }

  strcpy(area->name, label_name) ;

  /* read in the file */
  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorReturn(NULL, (ERROR_NOFILE, "%s: could not open label file %s",
                       Progname, fname)) ;

  cp = fgetl(line, 199, fp) ;
  if (!cp)
    ErrorReturn(NULL,
                (ERROR_BADFILE, "%s: empty label file %s", Progname, fname)) ;
  if (!sscanf(cp, "%d", &area->n_points))
    ErrorReturn(NULL,
                (ERROR_BADFILE, "%s: could not scan # of lines from %s",
                 Progname, fname)) ;
  area->max_points = area->n_points ;
  area->lv = (LABEL_VERTEX *)calloc(area->n_points, sizeof(LABEL_VERTEX)) ;
  if (!area->lv)
    ErrorExit(ERROR_NOMEMORY, 
              "%s: LabelRead(%s) could not allocate %d-sized vector",
              Progname, label_name, sizeof(LV)*area->n_points) ;
  nlines = 0 ;
  while ((cp = fgetl(line, 199, fp)) != NULL)
  {
    if (sscanf(cp, "%d %f %f %f %f", &vno, &x, &y, &z, &stat) != 5)
      ErrorReturn(NULL, (ERROR_BADFILE, "%s: could not parse %dth line in %s",
                Progname, area->n_points, fname)) ;
    area->lv[nlines].x = x ;
    area->lv[nlines].y = y ;
    area->lv[nlines].z = z ;
    area->lv[nlines].stat = stat ;
    area->lv[nlines].vno = vno ;
    nlines++ ;
  }

  fclose(fp) ;
  if (!nlines)
    ErrorReturn(NULL,
                (ERROR_BADFILE, 
                 "%s: no data in label file %s", Progname, fname));
  if (subject_name)
  {
    area->linear_transform = 
      labelLoadTransform(subject_name, subjects_dir, &area->transform) ;
    area->inverse_linear_transform = 
      get_inverse_linear_transform_ptr(&area->transform) ;
  }
  return(area) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
LabelDump(FILE *fp, LABEL *area)
{
  int  n ;

  fprintf(fp, "label %s, from subject %s\n", area->name, area->subject_name) ;
  for (n = 0 ; n < area->n_points ; n++)
    fprintf(fp, "%d  %2.3f  %2.3f  %2.3f\n", area->lv[n].vno, area->lv[n].x, 
            area->lv[n].y, area->lv[n].z) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
LabelFree(LABEL **parea)
{
  LABEL *area ;

  area = *parea ;
  *parea = NULL ;
  free(area->lv) ; 
  free(area) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
LabelToCanonical(LABEL *area, MRI_SURFACE *mris)
{
  int     n, vno ;
  VERTEX  *v ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    vno = area->lv[n].vno ;
    v = &mris->vertices[vno] ;
    area->lv[n].vno = -1 ;      /* not associated with a vertex anymore */
    area->lv[n].x = v->cx ;
    area->lv[n].y = v->cy ;
    area->lv[n].z = v->cz ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
LabelFromCanonical(LABEL *area, MRI_SURFACE *mris)
{
  int     n, vno, ui, vi, vlist[500000], nvertices ;
  VERTEX  *v ;
  MRI_SP  *mrisp ;
  LV     *lv ;

  mrisp = MRISPalloc(0.0f, 1) ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    MRISPcoordinate(mrisp, area->lv[n].x,area->lv[n].y,area->lv[n].z, &ui,&vi);
    MRISPvox(mrisp, ui, vi) = 1 ;
  }

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    ImageWrite(mrisp->Ip, "label.tif") ;

  /* now go through the surface and build a list of the vertices and
     keep track of those which fall into volume regions marked as within
     the area.
     */
  for (nvertices = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    MRISPcoordinate(mrisp, v->cx, v->cy, v->cz, &ui, &vi) ;
    if (nvertices+1 >= MAX_VERTICES)
      break ;
    if (MRISPvox(mrisp, ui, vi) > 0.0f)
      vlist[nvertices++] = vno ;
  }

  if (area->max_points < nvertices)
  {
    free(area->lv) ;
    area->max_points = nvertices ;
    area->lv = (LABEL_VERTEX *)calloc(nvertices, sizeof(LABEL_VERTEX)) ;
    if (!area->lv)
      ErrorExit(ERROR_NOMEMORY, 
                "%s: LabelFromCanonical could not allocate %d-sized vector",
                Progname, sizeof(LV)*nvertices) ;
  }
  area->n_points = nvertices ;
  for (n = 0 ; n < nvertices ; n++)
  {
    lv = &area->lv[n] ;
    lv->vno = vlist[n] ;
    v = &mris->vertices[vlist[n]] ;
    lv->x = v->cx ; lv->y = v->cy ; lv->z = v->cz ;
  }
  MRISPfree(&mrisp) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#if 0
#define VRES   2.0f
#define VSIZE  256
#define VDIM   (VSIZE/VRES)
int
LabelFromTalairach(LABEL *area, MRI_SURFACE *mris)
{
  int     n, vno, vlist[MAX_VERTICES], nvertices, xv, yv, zv ;
  VERTEX  *v ;
  MRI     *mri ;
  Real    xw, yw, zw ;
  LV     *lv ;

  /* first write it into a volume */
  mri = MRIalloc(VDIM,VDIM,VDIM, MRI_UCHAR) ;
  mri->inverse_linear_transform = mris->inverse_linear_transform ;
  mri->linear_transform = mris->linear_transform ;
  MRIsetResolution(mri, VRES, VRES, VRES) ;
  for (n = 0 ; n < area->n_points ; n++)
  {
    MRItalairachToVoxel(mri, area->lv[n].x,area->lv[n].y,area->lv[n].z,
                        &xw, &yw, &zw) ;
    xv = nint(xw) ; yv = nint(yw) ; zv = nint(zw) ;
    if (xv < 0) xv = 0 ; if (yv < 0) yv = 0 ; if (zv < 0) zv = 0 ;
    if (xv >= mri->width)  xv = mri->width-1 ;
    if (yv >= mri->height) yv = mri->height-1 ;
    if (zv >= mri->depth)  zv = mri->depth-1 ;
    MRIvox(mri, xv, yv, zv) = 1 ;
  }

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri, "label.mnc") ;
  /* now go through the surface and build a list of the vertices and
     keep track of those which fall into volume regions marked as within
     the area.
     */
  for (nvertices = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    MRISorigVertexToVoxel(v, mri, &xw, &yw, &zw) ;
    xv = nint(xw) ; yv = nint(yw) ; zv = nint(zw) ;
    if (xv < 0) xv = 0 ; if (yv < 0) yv = 0 ; if (zv < 0) zv = 0 ;
    if (xv >= mri->width)  xv = mri->width-1 ;
    if (yv >= mri->height) yv = mri->height-1 ;
    if (zv >= mri->depth)  zv = mri->depth-1 ;
    if (nvertices+1 >= MAX_VERTICES)
      break ;
    if (MRIvox(mri, xv, yv, zv) > 0)
      vlist[nvertices++] = vno ;
  }

  if (area->max_points < nvertices)
  {
    free(area->lv) ;
    area->max_points = nvertices ;
    area->lv = (LABEL_VERTEX *)calloc(nvertices, sizeof(LABEL_VERTEX)) ;
    if (!area->lv)
      ErrorExit(ERROR_NOMEMORY, 
                "%s: LabelFromTalairach could not allocate %d-sized vector",
                Progname, sizeof(LV)*nvertices) ;
  }
  area->n_points = nvertices ;
  for (n = 0 ; n < nvertices ; n++)
  {
    lv = &area->lv[n] ;
    lv->vno = vlist[n] ;
    v = &mris->vertices[vlist[n]] ;
    lv->x = v->origx ; lv->y = v->origy ; lv->z = v->origz ;
  }
  mri->inverse_linear_transform = NULL ;
  mri->linear_transform = NULL ;
  MRIfree(&mri) ;
  return(NO_ERROR) ;
}
#else
int
LabelFromTalairach(LABEL *area, MRI_SURFACE *mris)
{
  int     n, vno ;
  VERTEX  *v ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    DiagHeartbeat((float)(n+1) / (float)area->n_points) ;
    vno = MRIStalairachToVertex(mris,area->lv[n].x,area->lv[n].y,area->lv[n].z);
    v = &mris->vertices[vno] ;
    area->lv[n].vno = vno ;
    area->lv[n].x = v->origx ;
    area->lv[n].y = v->origy ;
    area->lv[n].z = v->origz ;
  }
  return(NO_ERROR) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
LabelToFlat(LABEL *area, MRI_SURFACE *mris)
{
  int     n, vno ;
  VERTEX  *v ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    vno = area->lv[n].vno ;
    if (vno >= 0)   /* already have associated vertex */
    {
      v = &mris->vertices[vno] ;
      area->lv[n].x = v->x ;
      area->lv[n].y = v->y ;
      area->lv[n].z = v->z ;
    }
    else    /* in canonical coordinate system - find closest vertex */
    {
      vno = 
        MRISfindClosestVertex(mris, area->lv[n].x,area->lv[n].y,area->lv[n].z);
      v = &mris->vertices[vno] ;
      area->lv[n].vno = vno ;
      area->lv[n].x = v->x ;
      area->lv[n].y = v->y ;
      area->lv[n].z = v->z ;
    }
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
LabelWrite(LABEL *area, char *label_name)
{
  FILE   *fp ;
  int    n, num ;
  char   fname[STRLEN], *cp, subjects_dir[STRLEN], lname[STRLEN] ;

  strcpy(lname, label_name) ;
  cp = strstr(lname, ".label") ;
  if (cp)
    *cp = 0 ;
#if 0
  cp = strrchr(lname, '/') ;
  if (cp)
    label_name = cp+1 ;
  else
    label_name = lname ;
#else
  label_name = lname ;
#endif
  cp = strrchr(lname, '/') ;
  if ((cp == NULL) && strlen(area->subject_name) > 0)
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, 
                "%s: no subject's directory specified in environment "
                "(SUBJECTS_DIR)", Progname) ;
    strcpy(subjects_dir, cp) ;
    sprintf(fname, "%s/%s/label/%s.label", subjects_dir,area->subject_name,
            label_name);
  }
  else
  {
    cp = strstr(lname, ".label") ;
    if (cp)
      strcpy(fname, label_name) ;
    else
      sprintf(fname, "%s.label", label_name) ;
  }

  for (num = n = 0 ; n < area->n_points ; n++)
    if (!area->lv[n].deleted)
      num++ ;

  fp = fopen(fname, "w") ;
  if (!fp)
    ErrorReturn(ERROR_NOFILE, (ERROR_NO_FILE,
                               "%s: could not open label file %s",
                               Progname, fname)) ;

#if 1
  fprintf(fp, "#!ascii label %s, from subject %s\n", 
          area->name, area->subject_name);
#endif
  fprintf(fp, "%d\n", num) ;
  for (n = 0 ; n < area->n_points ; n++)
    if (!area->lv[n].deleted)
      fprintf(fp, "%d  %2.3f  %2.3f  %2.3f %f\n", 
              area->lv[n].vno, area->lv[n].x, 
              area->lv[n].y, area->lv[n].z, area->lv[n].stat) ;
  fclose(fp) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
LabelRipRestOfSurface(LABEL *area, MRI_SURFACE *mris)
{
  int    vno, n ;
  VERTEX *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    v->ripflag = 1 ;
  }

  for (n = 0 ; n < area->n_points ; n++)
  {
    vno = area->lv[n].vno ;
    v = &mris->vertices[vno] ;
    v->ripflag = 0 ;
  }
  MRISripFaces(mris) ;
  MRISremoveRipped(mris) ;
  return(NO_ERROR) ;
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
LabelRemoveOverlap(LABEL *area1, LABEL *area2)
{
  int   n1, n2, vno ;

  for (n1 = 0 ; n1 < area1->n_points ; n1++)
  {
    vno = area1->lv[n1].vno ;
    for (n2 = 0 ; n2 < area2->n_points ; n2++)
    {
      if (vno == area2->lv[n2].vno)
      {
        area1->lv[n1].deleted = 1 ;
        break ;
      }
    }
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
LABEL *
LabelAlloc(int max_points, char *subject_name, char *label_name)
{
  LABEL  *area ;
  char   *cp, subjects_dir[100] ;

  area = (LABEL *)calloc(1, sizeof(LABEL)) ;
  if (!area)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate LABEL struct.",Progname);
  if (subject_name)
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, 
                "%s: no subject's directory specified in environment "
                "(SUBJECTS_DIR)", Progname) ;
    strcpy(subjects_dir, cp) ;
    strcpy(area->subject_name, subject_name) ;
  }
  else

  if (label_name)
    strcpy(area->name, label_name) ;

  area->n_points = 0 ;
  area->max_points = max_points ;
  area->lv = (LABEL_VERTEX *)calloc(area->max_points, sizeof(LABEL_VERTEX)) ;
  if (!area->lv)
    ErrorExit(ERROR_NOMEMORY, 
              "%s: LabelAlloc(%s) could not allocate %d-sized vector",
              Progname, label_name ? label_name : "", 
              sizeof(LV)*area->n_points) ;
  return(area) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
LabelCurvFill(LABEL *area, int *vertex_list, int nvertices, 
              int max_vertices, MRI_SURFACE *mris)
{
  int    n, nfilled, nv ;
  float  min_curv, max_curv, curv_thresh ;
  VERTEX *v, *vn ;
  LV     *lv, *lvn ;

  if (!max_vertices)
    max_vertices = area->max_points ;

  MRISclearMarks(mris) ;
  max_curv = min_curv = 0 ;
  for (n = 0 ; n < nvertices ; n++)
  {
    v = &mris->vertices[vertex_list[n]] ;
    if (v->ripflag)
      continue ;
    v->marked = 1 ;
    lv = &area->lv[n] ;
    lv->vno = vertex_list[n] ;
    lv->x = v->x ; lv->y = v->y ; lv->z = v->z ;
    area->n_points++ ;
    if (v->curv > max_curv)
      max_curv = v->curv ;
    if (v->curv < min_curv)
      min_curv = v->curv ;
  }
  if (-min_curv > max_curv)
    curv_thresh = min_curv ;
  else
    curv_thresh = max_curv ;
  curv_thresh *= 0.1 ;   /* 10% of max */
  fprintf(stderr, "starting fill with curvature threshold = %2.3f\n", 
          curv_thresh) ;
  do
  {
    nfilled = 0 ;

    for (n = 0; n < area->n_points && 
           area->n_points+nfilled < area->max_points ; n++)
    {
      lv = &area->lv[n] ;
      v = &mris->vertices[lv->vno] ;
      if (v->ripflag)
        continue ;
      for (nv = 0 ; nv < v->vnum ; nv++)   /* go through neighbors */
      {
        vn = &mris->vertices[v->v[nv]] ;
        if (vn->ripflag || vn->marked)
          continue ;
        if (((curv_thresh > 0) && (vn->curv > curv_thresh)) ||
            ((curv_thresh < 0) && (vn->curv < curv_thresh)))
        {
          vn->marked = 1 ;
          lvn = &area->lv[area->n_points+nfilled] ;
          nfilled++ ;
          lvn->vno = v->v[nv] ;
          lvn->x = vn->x ; lvn->y = vn->y ; lvn->z = vn->z ;
        }
        if (area->n_points+nfilled >= area->max_points)
          break ;
      }
      if (area->n_points+nfilled >= area->max_points)
        break ;
    }
    fprintf(stderr, "%d vertices added.\n", nfilled) ;
    area->n_points += nfilled ;
  } while ((nfilled > 0) && (area->n_points < max_vertices)) ;
  MRISclearMarks(mris) ;
  fprintf(stderr, "%d vertices in label %s\n",area->n_points, area->name);
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
LabelFillMarked(LABEL *area, MRI_SURFACE *mris)
{
  int    n, nfilled, nv ;
  VERTEX *v, *vn ;
  LV     *lv, *lvn ;


  do
  {
    nfilled = 0 ;

    for (n = 0; n < area->n_points && 
           area->n_points+nfilled < area->max_points ; n++)
    {
      lv = &area->lv[n] ;
      v = &mris->vertices[lv->vno] ; v->marked = 2 ;
      if (v->ripflag)
        continue ;
      for (nv = 0 ; nv < v->vnum ; nv++)   /* go through neighbors */
      {
        vn = &mris->vertices[v->v[nv]] ;
        if (vn->ripflag)
          continue ;
        if (vn->marked == 1)  /* add it to the label */
        {
          vn->marked = 2 ;
          lvn = &area->lv[area->n_points+nfilled] ;
          nfilled++ ;
          lvn->vno = v->v[nv] ;
          lvn->x = vn->x ; lvn->y = vn->y ; lvn->z = vn->z ;
        }
        if (area->n_points+nfilled >= area->max_points)
          break ;
      }
      if (area->n_points+nfilled >= area->max_points)
        break ;
    }
    /*    fprintf(stderr, "%d vertices added.\n", nfilled) ;*/
    area->n_points += nfilled ;
  } while ((nfilled > 0) && (area->n_points < area->max_points)) ;
#if 0
  fprintf(stderr, "%d vertices in label %s\n",area->n_points, area->name);
#endif
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
LabelFillAnnotated(LABEL *area, MRI_SURFACE *mris)
{
  int    n, nfilled, nv, annotation ;
  VERTEX *v, *vn ;
  LV     *lv, *lvn ;


  annotation = mris->vertices[area->lv[0].vno].annotation ;
  do
  {
    nfilled = 0 ;

    for (n = 0; n < area->n_points && 
           area->n_points+nfilled < area->max_points ; n++)
    {
      lv = &area->lv[n] ;
      v = &mris->vertices[lv->vno] ; v->marked = 1 ;
      if (v->ripflag)
        continue ;
      for (nv = 0 ; nv < v->vnum ; nv++)   /* go through neighbors */
      {
        vn = &mris->vertices[v->v[nv]] ;
        if (vn->ripflag || vn->marked)
          continue ;
        if (vn->annotation == annotation)  /* add it to the label */
        {
          vn->marked = 1 ;
          lvn = &area->lv[area->n_points+nfilled] ;
          nfilled++ ;
          lvn->vno = v->v[nv] ;
          lvn->x = vn->x ; lvn->y = vn->y ; lvn->z = vn->z ;
        }
        if (area->n_points+nfilled >= area->max_points)
          break ;
      }
      if (area->n_points+nfilled >= area->max_points)
        break ;
    }
    /*    fprintf(stderr, "%d vertices added.\n", nfilled) ;*/
    area->n_points += nfilled ;
  } while ((nfilled > 0) && (area->n_points < area->max_points)) ;
#if 0
  fprintf(stderr, "%d vertices in label %s\n",area->n_points, area->name);
#endif
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
LabelFillAll(LABEL *area, int *vertex_list, int nvertices, 
              int max_vertices, MRI_SURFACE *mris)
{
  int    n, nfilled, nv ;
  VERTEX *v, *vn ;
  LV     *lv, *lvn ;

  if (!max_vertices)
    max_vertices = area->max_points ;

  for (n = 0 ; n < nvertices ; n++)
  {
    v = &mris->vertices[vertex_list[n]] ;
    if (v->ripflag)
      continue ;
    v->marked = 1 ;
    lv = &area->lv[n] ;
    lv->vno = vertex_list[n] ;
    lv->x = v->x ; lv->y = v->y ; lv->z = v->z ;
    area->n_points++ ;
  }
  MRISclearMarks(mris) ;
  do
  {
    nfilled = 0 ;

    for (n = 0; n < area->n_points && 
           area->n_points+nfilled < area->max_points ; n++)
    {
      lv = &area->lv[n] ;
      v = &mris->vertices[lv->vno] ;
      if (v->ripflag)
        continue ;
      for (nv = 0 ; nv < v->vnum ; nv++)   /* go through neighbors */
      {
        vn = &mris->vertices[v->v[nv]] ;
        if (vn->ripflag || vn->marked)
          continue ;
        vn->marked = 1 ;
        lvn = &area->lv[area->n_points+nfilled] ;
        nfilled++ ;
        lvn->vno = v->v[nv] ;
        lvn->x = vn->x ; lvn->y = vn->y ; lvn->z = vn->z ;
        if (area->n_points+nfilled >= area->max_points)
          break ;
      }
      if (area->n_points+nfilled >= area->max_points)
        break ;
    }
    fprintf(stderr, "%d vertices added.\n", nfilled) ;
    area->n_points += nfilled ;
  } while ((nfilled > 0) && (area->n_points < max_vertices)) ;
  MRISclearMarks(mris) ;
  fprintf(stderr, "%d vertices in label %s\n",area->n_points, area->name);
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static Transform *
labelLoadTransform(char *subject_name, char *sdir,General_transform *transform)
{
  char xform_fname[STRLEN] ;

  sprintf(xform_fname, "%s/%s/mri/transforms/talairach.xfm",
          sdir, subject_name) ;
  if (input_transform_file(xform_fname, transform) != OK)
    ErrorReturn(NULL,
                (ERROR_NOFILE, "%s: could not load transform file '%s'", 
                 Progname, xform_fname)) ;

  return(get_linear_transform_ptr(transform)) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
LabelTalairachTransform(LABEL *area, MRI_SURFACE *mris)
{
  int    n ;
  LV     *lv ;
  Real   x, y, z, xt, yt, zt ;
  VERTEX *v ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    lv = &area->lv[n] ;
    v = &mris->vertices[lv->vno] ;
    x = v->origx ; y = v->origy ; z = v->origz ;
    transform_point(area->linear_transform, x, y, z, &xt, &yt, &zt) ;
    lv->x = xt ; lv->y = yt ; lv->z = zt ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
LabelSphericalTransform(LABEL *area, MRI_SURFACE *mris)
{
  int     n ;
  LV      *lv ;
  VERTEX  *v ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    lv = &area->lv[n] ;
    v = &mris->vertices[lv->vno] ;
    lv->x = v->cx ; lv->y = v->cy ; lv->z = v->cz ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MATRIX *
LabelCovarianceMatrix(LABEL *area, MATRIX *mCov)
{
  MATRIX   *mInputs ;
  int      i ;
  LV       *lv ;

  mInputs = MatrixAlloc(area->n_points, 3, MATRIX_REAL) ;
  for (i = 0 ; i < area->n_points ; i++)
  {
    lv = &area->lv[i] ;
    *MATRIX_RELT(mInputs,i+1, 1) = lv->x ;
    *MATRIX_RELT(mInputs,i+1, 2) = lv->y ;
    *MATRIX_RELT(mInputs,i+1, 3) = lv->z ;
  }

  mCov = MatrixCovariance(mInputs, mCov, NULL) ;

  return(mCov) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
LABEL *
LabelCombine(LABEL *asrc, LABEL *adst)
{
  LABEL  *atmp ;
  int    n ;
  LV     *vsrc, *vdst ;

  if (!adst)
    adst = LabelClone(asrc) ;
  if (adst->max_points < asrc->n_points+adst->n_points)/* won't fit - expand */
  {
    atmp = LabelAlloc(2*(asrc->n_points+adst->n_points),
                      asrc->subject_name, asrc->name) ;
    LabelCopy(adst, atmp) ;
    LabelFree(&adst) ;
    adst = atmp ;
  }
  
  for (n = 0 ; n < asrc->n_points ; n++)
  {
    vsrc = &asrc->lv[n] ;
    vdst = &adst->lv[n+adst->n_points] ;
    memmove(vdst, vsrc, sizeof(LV)) ;
  }
  adst->n_points += asrc->n_points ;
  return(adst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
LABEL *
LabelCopy(LABEL *asrc, LABEL *adst)
{
  if (!adst)
    adst = LabelAlloc(asrc->n_points, asrc->subject_name, asrc->name) ;

  adst->n_points = asrc->n_points ;
  strcpy(adst->name, asrc->name) ;
  strcpy(adst->subject_name, asrc->subject_name) ;

  memmove(adst->lv, asrc->lv, asrc->n_points*sizeof(LABEL_VERTEX)) ;
  return(adst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
LabelRemoveDuplicates(LABEL *area)
{
  int    n1, n2, deleted = 0 ;
  LV     *lv1, *lv2 ;

  for (n1 = 0 ; n1 < area->n_points ; n1++)
  {
    lv1 = &area->lv[n1] ;
    if (lv1->deleted)
      continue ;
    for (n2 = n1+1 ; n2 < area->n_points ; n2++)
    {
      lv2 = &area->lv[n2] ;
      if (lv1->vno == lv2->vno)
      {
        deleted++ ;
        lv2->deleted = 1 ;
        lv1->stat++ ;
      }
    }
  }
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "%d duplicate vertices removed from label %s.\n",
            deleted, area->name) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
double
LabelArea(LABEL *area, MRI_SURFACE *mris)
{
  int      i ;
  double   total_area ;
  LV       *lv ;
  VERTEX   *v ;

  for (total_area = 0.0, i = 0 ; i < area->n_points ; i++)
  {
    lv = &area->lv[i] ;
    v = &mris->vertices[lv->vno] ;
    if (v->ripflag)
      continue ;
    total_area += v->area ;
  }

  return(total_area) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
LabelMean(LABEL *area, double *px, double *py, double *pz)
{
  int      i, n ;
  double   x, y, z ;
  LV       *lv ;

  for (x = y = z = 0.0, n = i = 0 ; i < area->n_points ; i++)
  {
    lv = &area->lv[i] ;
    
    x += lv->x ; y += lv->y ; z += lv->z ; n++ ;
  }

  *px = x / (double)n ; *py = y / (double)n ; *pz = z / (double)n ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
double
LabelVariance(LABEL *area, double ux, double uy, double uz)
{
  int      i, n ;
  double   xd, yd, zd, dsq ;
  LV       *lv ;

  for (dsq = 0.0, n = i = 0 ; i < area->n_points ; i++)
  {
    lv = &area->lv[i] ;

    xd = lv->x - ux ; yd = lv->y - uy ; zd = lv->z - uz ;
    dsq += xd*xd + yd*yd + zd*zd ;
    n++ ;
  }
  
  dsq /= (double)n ;
  return(dsq) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
LabelMark(LABEL *area, MRI_SURFACE *mris)
{
  int    n, vno ;
  VERTEX *v ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    vno = area->lv[n].vno ;
    v = &mris->vertices[vno] ;
    v->marked = 1 ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
LabelUnmark(LABEL *area, MRI_SURFACE *mris)
{
  int    n, vno ;
  VERTEX *v ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    vno = area->lv[n].vno ;
    v = &mris->vertices[vno] ;
    v->marked = 0 ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
LabelToOriginal(LABEL *area, MRI_SURFACE *mris)
{
  int     n, vno ;
  VERTEX  *v ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    vno = area->lv[n].vno ;
    v = &mris->vertices[vno] ;
    area->lv[n].x = v->origx ;
    area->lv[n].y = v->origy ;
    area->lv[n].z = v->origz ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
LABEL *
LabelFromMarkedSurface(MRI_SURFACE *mris)
{
  int    vno, npoints, n ;
  LABEL  *area ;
  VERTEX *v ;

  for (npoints = vno = 0 ; vno < mris->nvertices ; vno++)
    if (mris->vertices[vno].marked)
      npoints++ ;

  if (!npoints)
    return(NULL) ;
  area = LabelAlloc(npoints, NULL, NULL) ;
  for (n = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (!v->marked)
      continue ;
    area->lv[n].x = v->x ;
    area->lv[n].y = v->y ;
    area->lv[n].z = v->z ;
    area->lv[n].vno = vno ;
    area->lv[n].stat = v->stat ;
    n++ ;
  }
  area->n_points = npoints ;
  return(area) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
LabelMarkSurface(LABEL *area, MRI_SURFACE *mris)
{
  int     n, vno ;
  VERTEX  *v ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    vno = area->lv[n].vno ;
    if (vno < 0)
      continue ;
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->marked = 1 ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#include "mrishash.h"
int
LabelFillUnassignedVertices(MRI_SURFACE *mris, LABEL *area)
{
  int     n, i, vno, min_vno ;
  LV      *lv ;
  MHT     *mht ;
  MHBT    *bucket ;
  MHB     *bin ;
  VERTEX  *v ;
  float   dx, dy, dz, x, y, z, dist, min_dist ;

  for (i = n = 0 ; n < area->n_points ; n++)
  {
    lv = &area->lv[n] ;
    if (lv->vno >= 0 && lv->vno < mris->nvertices)
      continue ;
    i++ ;   /* count # of unassigned vertices */
  }
  if (i <= 0)
    return(NO_ERROR) ;  /* no work needed */

  fprintf(stderr,"%d unassigned vertices in label - building spatial LUT...\n",
          i) ;

  /* if we can't find a vertex within 10 mm of the point, something is wrong */
  mht = MHTfillVertexTableRes(mris, NULL, ORIGINAL_VERTICES, 10.0) ;
  fprintf(stderr, "assigning vertex numbers to label...\n") ;
  for (n = 0 ; n < area->n_points ; n++)
  {
    lv = &area->lv[n] ;
    if (lv->vno >= 0 && lv->vno <= mris->nvertices)
      continue ;
    bucket = MHTgetBucket(mht, lv->x, lv->y, lv->z) ;

    x = lv->x ; y = lv->y ; z = lv->z ;
    min_dist = 10000.0 ; min_vno = -1 ;
    if (Gdiag_no == 55)
      lv->vno = MRISfindClosestOriginalVertex(mris, lv->x, lv->y, lv->z) ;
    if (bucket)
    {
      for (bin = bucket->bins, i = 0 ; i < bucket->nused ; i++, bin++)
      {
        vno = bin->fno ; v = &mris->vertices[vno] ;
        if (vno == Gdiag_no)
          DiagBreak() ;
        dx = v->origx - x ; dy = v->origy - y ; dz = v->origz - z ;
        dist = sqrt(dx*dx + dy*dy + dz*dz) ;
        if (dist < min_dist)
        {
          min_dist = dist ;
          min_vno = vno ;
        }
      }
    }
    lv->vno = min_vno ;
  }
  LabelRemoveDuplicates(area) ;
  MHTfree(&mht) ;
  return(NO_ERROR) ;
}

LABEL *
LabelSphericalCombine(MRI_SURFACE *mris, LABEL *asrc, MRIS_HASH_TABLE *mht,
                      MRI_SURFACE *mris_dst, LABEL *adst)
{
  int              vno, n, nfilled, m ;
  VERTEX           *v, *vdst, *vn, *vsrc ;
  LABEL_VERTEX     *lv_dst ;
  MRIS_HASH_TABLE  *mht_src ;
  double           max_len ;

  if (!adst)
    adst = LabelClone(asrc) ;

  if (adst->max_points < asrc->n_points+adst->n_points)/* won't fit - expand */
  {
    LABEL   *atmp ;

    atmp = LabelAlloc(2*(asrc->n_points+adst->n_points),
                      asrc->subject_name, asrc->name) ;
    LabelCopy(adst, atmp) ;
    LabelFree(&adst) ;
    adst = atmp ;
  }

  for (n = 0 ; n < asrc->n_points ; n++)
  {
    vno = asrc->lv[n].vno ;
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    vdst = MHTfindClosestVertex(mht, mris_dst, v) ;
    if (!vdst)
    {
      ErrorPrintf(ERROR_BADPARM, "MRIScombine: cannot map vno %d", vno) ;
      continue ;
    }
    if (vdst-mris_dst->vertices == Gdiag_no)
      DiagBreak() ;
    lv_dst = &adst->lv[adst->n_points++] ;
    lv_dst->vno = vdst-mris_dst->vertices ;
    if (lv_dst->vno == 60008)
      DiagBreak() ;
    if (lv_dst->vno == 85592)
      DiagBreak() ;
    if (lv_dst->vno == Gdiag_no)
      DiagBreak() ;
    lv_dst->x = asrc->lv[n].x ;
    lv_dst->y = asrc->lv[n].y ;
    lv_dst->z = asrc->lv[n].z ;
    lv_dst->stat = 1 ;
  }

  MRIScomputeVertexSpacingStats(mris, NULL, NULL, &max_len, NULL,NULL);
  mht_src = 
    MHTfillVertexTableRes(mris, NULL, CURRENT_VERTICES, 2*max_len);

  MRISclearMarks(mris_dst) ;
  MRISclearMarks(mris) ;
  LabelMark(asrc, mris) ; LabelMark(adst, mris_dst) ;
  do
  {
    /* map the nbr of every point in the label back to src, and if in
       label add it
    */
    nfilled = 0 ;
    for (n = 0 ; n < adst->n_points ; n++)
    {
      v = &mris_dst->vertices[adst->lv[n].vno] ;
      for (m = 0 ; m < v->vnum ; m++)
      {
        vn = &mris_dst->vertices[v->v[m]] ;
        if (vn->marked)
          continue ;   /* already in label */
        vsrc = MHTfindClosestVertex(mht_src, mris, vn) ;
        if (vsrc == NULL)
          DiagBreak() ;
        if (vsrc - mris->vertices == 62644)
          DiagBreak() ;
        if (vsrc - mris->vertices == Gdiag_no)
          DiagBreak() ;
        if (vsrc->marked)
        {
          if (adst->n_points >= adst->max_points-1)
          {
            LABEL *atmp ;
            
            atmp = LabelAlloc(2*adst->n_points,adst->subject_name,adst->name);
            LabelCopy(adst, atmp) ;
            LabelFree(&adst) ;
            adst = atmp ;
          }

          vn->marked = 1 ;
          lv_dst = &adst->lv[adst->n_points++] ;
          lv_dst->x = adst->lv[n].x ;
          lv_dst->y = adst->lv[n].y ;
          lv_dst->z = adst->lv[n].z ;
          lv_dst->vno = v->v[m] ;
          lv_dst->stat = 1 ;
          if (lv_dst->vno == Gdiag_no)
            DiagBreak() ;
          if (lv_dst->vno == 55185)
            DiagBreak() ;
          if (lv_dst->vno == 85592)
            DiagBreak() ;
          nfilled++ ;
        }
      }
    }
  } while (nfilled != 0) ;
  return(adst) ;
}

int
LabelNormalizeStats(LABEL *area, float norm)
{
  LABEL_VERTEX   *lv ;
  int            n ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    lv = &area->lv[n] ;
    lv->stat /= norm ;
  }
  return(NO_ERROR) ;
}

/*-------------------------------------------------------------
  MaskSurfLabel() - removes vertices from a label based on the 
  value of those vertices in a surface mask.
  -------------------------------------------------------------*/
LABEL *MaskSurfLabel(LABEL *lbl, MRI *SurfMask, 
         float thresh, char *masksign, int frame)
{
  LABEL *msklbl;
  int n,vno, ok, noverlap;
  float maskval;

  if(frame >= SurfMask->nframes){
    printf("ERROR:  MaskSurfLabel: specified frame %d is greater "
     "than the number of frame %d\n",frame, SurfMask->nframes);
    return(NULL);
  }

  /* Check that the label is a surface-based label */
  ok = 0;
  for(n = 0; n < lbl->n_points; n++){
    if(lbl->lv[n].vno != 0) ok = 1;
  }
  if(!ok){
    printf("ERROR: MaskSurfLabel: all label vertices are 0, proably\n");
    printf("       not a label created from a surface\n");
    return(NULL);
  }

  /* -- Allocate the Target Label ---*/
  msklbl = LabelAlloc(lbl->n_points,lbl->subject_name,"labelfilenamehere");
  msklbl->n_points = lbl->n_points;

  /* Compute the number of overlapping vertices */
  noverlap = 0;
  for(n = 0; n < lbl->n_points; n++){
    vno = lbl->lv[n].vno;
    if(vno >= SurfMask->width){
      printf("ERROR: MaskSurfLabel: label vertex %d has vno %d which "
       "is larger than the number of vertices %d in mask\n",
       n,vno,SurfMask->width);
      return(NULL);
    }
    maskval = MRIgetVoxVal(SurfMask,vno,0,0,frame);
    if(!strcmp(masksign,"abs")) maskval = fabs(maskval);
    if(!strcmp(masksign,"neg")) maskval = -maskval;
    //printf("%4d  %6d  %g %g  %d\n",n,vno,maskval,thresh,maskval>thresh);
    if(maskval < thresh) continue;
    msklbl->lv[noverlap].vno = vno;
    msklbl->lv[noverlap].x = lbl->lv[noverlap].x;
    msklbl->lv[noverlap].y = lbl->lv[noverlap].y;
    msklbl->lv[noverlap].z = lbl->lv[noverlap].z;
    msklbl->lv[noverlap].stat = lbl->lv[noverlap].stat;
    noverlap ++;
  }
  if(noverlap==0){
    printf("WARNING: MaskSurfLabel: no overlap between label and mask\n");
    return(NULL);
  }
  msklbl->n_points = noverlap;

  return(msklbl);
}
