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
  area->x = (float *)calloc(area->n_points, sizeof(float)) ;
  if (!area->x)
    ErrorExit(ERROR_NOMEMORY, 
              "%s: LabelRead(%s) could not allocate %d-sized vector",
              Progname, label_name, area->n_points) ;
  area->y = (float *)calloc(area->n_points, sizeof(float)) ;
  if (!area->y)
    ErrorExit(ERROR_NOMEMORY, 
              "%s: LabelRead(%s) could not allocate %d-sized vector",
              Progname, label_name, area->n_points) ;
  area->z = (float *)calloc(area->n_points, sizeof(float)) ;
  if (!area->z)
    ErrorExit(ERROR_NOMEMORY, 
              "%s: LabelRead(%s) could not allocate %d-sized vector",
              Progname, label_name, area->n_points) ;
  area->deleted = (unsigned char *)calloc(area->n_points, sizeof(float)) ;
  if (!area->deleted)
    ErrorExit(ERROR_NOMEMORY, 
              "%s: LabelRead(%s) could not allocate %d-sized vector",
              Progname, label_name, area->n_points) ;
  area->vno = (int *)calloc(area->n_points, sizeof(int)) ;
  if (!area->vno)
    ErrorExit(ERROR_NOMEMORY, 
              "%s: LabelRead(%s) could not allocate %d-sized vector",
              Progname, label_name, area->n_points) ;

  nlines = 0 ;
  while ((cp = fgetl(line, 199, fp)) != NULL)
  {
    if (sscanf(cp, "%d %f %f %f", &vno, &x, &y, &z) != 4)
      ErrorExit(ERROR_BADFILE, "%s: could not parse %dth line in %s",
                Progname, area->n_points, fname) ;
    area->x[nlines] = x ;
    area->y[nlines] = y ;
    area->z[nlines] = z ;
    area->vno[nlines] = vno ;
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
    fprintf(fp, "%d  %2.3f  %2.3f  %2.3f\n", area->vno[n], area->x[n], 
            area->y[n], area->z[n]) ;
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
  free(area->x) ; free(area->y) ; free(area->z) ; free(area->vno) ;
  free(area->deleted) ;
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
    vno = area->vno[n] ;
    v = &mris->vertices[vno] ;
    area->vno[n] = -1 ;      /* not associated with a vertex anymore */
    area->x[n] = v->cx ;
    area->y[n] = v->cy ;
    area->z[n] = v->cz ;
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
    vno =MRISfindClosestCanonicalVertex(mris,area->x[n],area->y[n],area->z[n]);
    v = &mris->vertices[vno] ;
    area->vno[n] = vno ;
    area->x[n] = v->cx ;
    area->y[n] = v->cy ;
    area->z[n] = v->cz ;
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
    vno = area->vno[n] ;
    if (vno >= 0)   /* already have associated vertex */
    {
      v = &mris->vertices[vno] ;
      area->x[n] = v->x ;
      area->y[n] = v->y ;
      area->z[n] = v->z ;
    }
    else    /* in canonical coordinate system - find closest vertex */
    {
      vno = MRISfindClosestVertex(mris, area->x[n], area->y[n], area->z[n]) ;
      v = &mris->vertices[vno] ;
      area->vno[n] = vno ;
      area->x[n] = v->x ;
      area->y[n] = v->y ;
      area->z[n] = v->z ;
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
LabelWrite(LABEL *area, char *fname)
{
  FILE  *fp ;
  int  n, num ;

  
  for (num = n = 0 ; n < area->n_points ; n++)
    if (!area->deleted[n])
      num++ ;

  fp = fopen(fname, "w") ;
  if (!fp)
    ErrorExit(ERROR_NOFILE, "%s: could not open label file %s",
                Progname, fname) ;

#if 0
  fprintf(fp, "# label %s, from subject %s\n", area->name, area->subject_name);
#endif
  for (n = 0 ; n < area->n_points ; n++)
    if (!area->deleted[n])
      fprintf(fp, "%d  %2.3f  %2.3f  %2.3f\n", area->vno[n], area->x[n], 
              area->y[n], area->z[n]) ;
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
    vno = area->vno[n] ;
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
    vno = area1->vno[n1] ;
    for (n2 = 0 ; n2 < area2->n_points ; n2++)
    {
      if (vno == area2->vno[n2])
      {
        area1->deleted[n1] = 1 ;
        break ;
      }
    }
  }
  return(NO_ERROR) ;
}
