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



/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
LABEL *
LabelRead(char *subject_name, char *label_name)
{
  LABEL  *area ;
  char   fname[100], *cp, line[200], subjects_dir[100] ;
  FILE   *fp ;
  int    vno, nlines ;
  float  x, y, z ;

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
    sprintf(fname, "%s/%s/label/%s.label", subjects_dir,subject_name,
            label_name);
    strcpy(area->subject_name, subject_name) ;
  }
  else
    strcpy(fname, label_name) ;

  strcpy(area->name, label_name) ;

  /* read in the file */
  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorExit(ERROR_NOFILE, "%s: could not open label file %s",
                Progname, fname) ;

  cp = fgetl(line, 199, fp) ;
  if (!cp)
    ErrorExit(ERROR_BADFILE, "%s: empty label file %s", Progname, fname) ;
  if (!sscanf(cp, "%d", &area->n_points))
    ErrorExit(ERROR_BADFILE, "%s: could not scan # of lines from %s",
              Progname, fname) ;
  area->max_points = area->n_points ;
  area->lv = (LABEL_VERTEX *)calloc(area->n_points, sizeof(LABEL_VERTEX)) ;
  if (!area->lv)
    ErrorExit(ERROR_NOMEMORY, 
              "%s: LabelRead(%s) could not allocate %d-sized vector",
              Progname, label_name, sizeof(LV)*area->n_points) ;
  nlines = 0 ;
  while ((cp = fgetl(line, 199, fp)) != NULL)
  {
    if (sscanf(cp, "%d %f %f %f", &vno, &x, &y, &z) != 4)
      ErrorExit(ERROR_BADFILE, "%s: could not parse %dth line in %s",
                Progname, area->n_points, fname) ;
    area->lv[nlines].x = x ;
    area->lv[nlines].y = y ;
    area->lv[nlines].z = z ;
    area->lv[nlines].vno = vno ;
    nlines++ ;
  }

  fclose(fp) ;
  if (!nlines)
    ErrorExit(ERROR_BADFILE, "%s: no data in label file %s", Progname, fname);
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
  int     n, vno ;
  VERTEX  *v ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    vno =MRISfindClosestCanonicalVertex(mris,area->lv[n].x,area->lv[n].y,area->lv[n].z);
    v = &mris->vertices[vno] ;
    area->lv[n].vno = vno ;
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
  FILE  *fp ;
  int  n, num ;
  char   fname[100], *cp, subjects_dir[100] ;

  if (strlen(area->subject_name) > 0)
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
    strcpy(fname, label_name) ;

  for (num = n = 0 ; n < area->n_points ; n++)
    if (!area->lv[n].deleted)
      num++ ;

  fp = fopen(fname, "w") ;
  if (!fp)
    ErrorExit(ERROR_NOFILE, "%s: could not open label file %s",
                Progname, fname) ;

#if 1
  fprintf(fp, "#!ascii label %s, from subject %s\n", 
          area->name, area->subject_name);
#endif
  fprintf(fp, "%d\n", num) ;
  for (n = 0 ; n < area->n_points ; n++)
    if (!area->lv[n].deleted)
      fprintf(fp, "%d  %2.3f  %2.3f  %2.3f\n", area->lv[n].vno, area->lv[n].x, 
              area->lv[n].y, area->lv[n].z) ;
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
              Progname, label_name, sizeof(LV)*area->n_points) ;
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
    }
    fprintf(stderr, "%d vertices added.\n", nfilled) ;
    area->n_points += nfilled ;
  } while ((nfilled > 0) && (area->n_points < max_vertices)) ;
  MRISclearMarks(mris) ;
  fprintf(stderr, "%d vertices in label %s\n",area->n_points, area->name);
  return(NO_ERROR) ;
}

