/**
 * @file  label.c
 * @brief utilities for manipulating ROIs.
 *
 * utilities (surface and volume) for manipulating arbitrary lists of
 * vertices/voxels.
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2012/10/01 18:59:58 $
 *    $Revision: 1.110 $
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <errno.h>

#include "macros.h"
#include "minc_volume_io.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "mrisurf.h"
#include "label.h"
#include "mri.h"
#include "mrishash.h"
#include "fio.h"

extern const char* Progname;

static LABEL_VERTEX *labelFindVertexNumber(LABEL *area, int vno) ;
static Transform *labelLoadTransform(const char *subject_name,
                                     const char *sdir,
                                     General_transform *transform) ;

#define MAX_VERTICES 500000
/*-----------------------------------------------------
------------------------------------------------------*/
LABEL *LabelRead(const char *subject_name, const char *label_name)
{
  LABEL  *area ;
  char   fname[STRLEN], *cp, line[STRLEN], subjects_dir[STRLEN], lname[STRLEN];
  char   label_name0[STRLEN];
  FILE   *fp ;
  int    vno, nlines ;
  float  x, y, z, stat ;

  area = (LABEL *)calloc(1, sizeof(LABEL)) ;
  if (!area)
  {
    ErrorExit(ERROR_NOMEMORY,"%s: could not allocate LABEL struct.",Progname);
  }

  sprintf(label_name0,"%s",label_name); // keep a copy

  if (subject_name && !strchr(label_name, '/'))
  {
    // subject name passed and label name does not have /
    // so assume it is in subject/label
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: no subject's directory specified in environment "
                "(SUBJECTS_DIR)", Progname) ;
    strcpy(subjects_dir, cp) ;
    strcpy(lname, label_name) ;
    cp = strstr(lname, ".label") ;
    if(cp)
    {
      *cp = 0 ;
    }
    cp = strrchr(lname, '/') ;
    if(cp)
    {
      label_name = cp+1 ;
    }
    else
    {
      label_name = lname ;
    }
    sprintf(fname, "%s/%s/label/%s.label", subjects_dir,subject_name,
            label_name);
    strcpy(area->subject_name, subject_name) ;
  }
  else
  {
    strcpy(fname, label_name) ;
    cp = getenv("SUBJECTS_DIR") ;
    if(!cp)
    {
      strcpy(subjects_dir, ".") ;
    }
    else
    {
      strcpy(subjects_dir, cp) ;
    }
    strcpy(lname, label_name) ;
    cp = strstr(lname, ".label") ;
    if (cp == NULL)
    {
      sprintf(fname, "%s.label", lname);
    }
    else
    {
      strcpy(fname, label_name) ;
    }
  }

  // As a last resort, treat label_name0 as a full path name
  if(!fio_FileExistsReadable(fname) && fio_FileExistsReadable(label_name0))
  {
    sprintf(fname,"%s",label_name0);
  }

  //  printf("%s %s\n",label_name0,fname);

  strncpy(area->name, label_name, STRLEN-1) ;

  /* read in the file */
  errno = 0 ;
  fp = fopen(fname, "r") ;
  if (errno)
  {
    perror(NULL);
  }
  if (!fp)
    ErrorReturn(NULL, (ERROR_NOFILE, "%s: could not open label file %s",
                       Progname, fname)) ;


  cp = fgetl(line, STRLEN, fp) ;
  if (!cp)
    ErrorReturn(NULL,
                (ERROR_BADFILE, "%s: empty label file %s", Progname, fname)) ;
  if (!sscanf(cp, "%d", &area->n_points))
  {
    printf("\n%s\n",cp);
    ErrorReturn(NULL,
                (ERROR_BADFILE, "%s: could not scan # of lines from %s",
                 Progname, fname)) ;
  }
  area->max_points = area->n_points ;
  area->lv = (LABEL_VERTEX *)calloc(area->n_points, sizeof(LABEL_VERTEX)) ;
  if (!area->lv)
    ErrorExit(ERROR_NOMEMORY,
              "%s: LabelRead(%s) could not allocate %d-sized vector",
              Progname, label_name, sizeof(LV)*area->n_points) ;
  nlines = 0 ;
  while ((cp = fgetl(line, STRLEN, fp)) != NULL)
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

  fprintf(fp, "label %s, from subject %s %s\n",
          area->name, area->subject_name, area->space) ;
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
  strncpy (area->space, "coords=canonical", sizeof(area->space));
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
  LV      *lv ;

  mrisp = MRISPalloc(0.0f, 1) ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    MRISPcoordinate(mrisp, area->lv[n].x,area->lv[n].y,area->lv[n].z, &ui,&vi);
    MRISPvox(mrisp, ui, vi) = 1 ;
  }

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    ImageWrite(mrisp->Ip, "label.tif") ;
  }

  /* now go through the surface and build a list of the vertices and
     keep track of those which fall into volume regions marked as within
     the area.
     */
  for (nvertices = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    MRISPcoordinate(mrisp, v->cx, v->cy, v->cz, &ui, &vi) ;
    if (nvertices+1 >= MAX_VERTICES)
    {
      break ;
    }
    if (MRISPvox(mrisp, ui, vi) > 0.0f)
    {
      vlist[nvertices++] = vno ;
    }
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
    lv->x = v->cx ;
    lv->y = v->cy ;
    lv->z = v->cz ;
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
  double  xw, yw, zw ;
  LV      *lv ;

  /* first write it into a volume */
  mri = MRIalloc(VDIM,VDIM,VDIM, MRI_UCHAR) ;
  mri->inverse_linear_transform = mris->inverse_linear_transform ;
  mri->linear_transform = mris->linear_transform ;
  MRIsetResolution(mri, VRES, VRES, VRES) ;
  for (n = 0 ; n < area->n_points ; n++)
  {
    MRItalairachToVoxel(mri, area->lv[n].x,area->lv[n].y,area->lv[n].z,
                        &xw, &yw, &zw) ;
    xv = nint(xw) ;
    yv = nint(yw) ;
    zv = nint(zw) ;
    if (xv < 0)
    {
      xv = 0 ;
    }
    if (yv < 0)
    {
      yv = 0 ;
    }
    if (zv < 0)
    {
      zv = 0 ;
    }
    if (xv >= mri->width)
    {
      xv = mri->width-1 ;
    }
    if (yv >= mri->height)
    {
      yv = mri->height-1 ;
    }
    if (zv >= mri->depth)
    {
      zv = mri->depth-1 ;
    }
    MRIvox(mri, xv, yv, zv) = 1 ;
  }

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri, "label.mnc") ;
  }
  /* now go through the surface and build a list of the vertices and
     keep track of those which fall into volume regions marked as within
     the area.
     */
  for (nvertices = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    MRISorigVertexToVoxel(mris, v, mri, &xw, &yw, &zw) ;
    xv = nint(xw) ;
    yv = nint(yw) ;
    zv = nint(zw) ;
    if (xv < 0)
    {
      xv = 0 ;
    }
    if (yv < 0)
    {
      yv = 0 ;
    }
    if (zv < 0)
    {
      zv = 0 ;
    }
    if (xv >= mri->width)
    {
      xv = mri->width-1 ;
    }
    if (yv >= mri->height)
    {
      yv = mri->height-1 ;
    }
    if (zv >= mri->depth)
    {
      zv = mri->depth-1 ;
    }
    if (nvertices+1 >= MAX_VERTICES)
    {
      break ;
    }
    if (MRIvox(mri, xv, yv, zv) > 0)
    {
      vlist[nvertices++] = vno ;
    }
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
    lv->x = v->origx ;
    lv->y = v->origy ;
    lv->z = v->origz ;
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
    vno = MRIStalairachToVertex(mris,
                                area->lv[n].x,
                                area->lv[n].y,
                                area->lv[n].z);
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
  MRIS_HASH_TABLE *mht ;

  mht = MHTfillVertexTable(mris, NULL, CURRENT_VERTICES) ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    vno = area->lv[n].vno ;
    if (vno >= 0 && vno < mris->nvertices)/* already have associated vertex */
    {
      v = &mris->vertices[vno] ;
      area->lv[n].x = v->x ;
      area->lv[n].y = v->y ;
      area->lv[n].z = v->z ;
    }
    else    /* in canonical coordinate system - find closest vertex */
    {
#if 0
      vno =
        MRISfindClosestVertex(mris,
                              area->lv[n].x,
                              area->lv[n].y,
                              area->lv[n].z,
                              &dmin);
#endif
      v =
        MHTfindClosestVertexInTable(mht,
                                    mris,
                                    area->lv[n].x,
                                    area->lv[n].y,
                                    area->lv[n].z,
                                    0);
      if (v == NULL)
      {
        continue ;
      }
      vno = v - mris->vertices ;
      ;
      if (vno == Gdiag_no)
      {
        DiagBreak() ;
      }
      area->lv[n].vno = vno ;
      area->lv[n].x = v->x ;
      area->lv[n].y = v->y ;
      area->lv[n].z = v->z ;
    }
  }
  MHTfree(&mht) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
LabelWrite(LABEL *area, const char *label_name)
{
  FILE   *fp ;
  int    n, num, nbytes ;
  char   fname[STRLEN], *cp, subjects_dir[STRLEN], lname[STRLEN] ;

  strcpy(lname, label_name) ;
  cp = strrchr(lname, '.') ;
  if (cp && stricmp(cp, ".label") == 0)
  {
    *cp = 0 ;
  }
  label_name = lname ;
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
    {
      strcpy(fname, label_name) ;
    }
    else
    {
      sprintf(fname, "%s.label", label_name) ;
    }
  }

  for (num = n = 0 ; n < area->n_points ; n++)
    if (!area->lv[n].deleted)
    {
      num++ ;
    }

  printf("LabelWrite: saving to %s\n",fname);

  fp = fopen(fname, "w") ;
  if (!fp)
    ErrorReturn(ERROR_NOFILE, (ERROR_NO_FILE,
                               "%s: could not open label file %s",
                               Progname, fname)) ;

  nbytes = fprintf(fp, "#!ascii label %s , from subject %s vox2ras=TkReg %s\n",
                   area->name, area->subject_name, area->space);
  if(nbytes < 0)
  {
    printf("ERROR: writing to %s\n",fname);
    fclose(fp);
    return(1);
  }

  nbytes =   fprintf(fp, "%d\n", num) ;
  if(nbytes < 0)
  {
    printf("ERROR: writing to %s\n",fname);
    fclose(fp);
    return(1);
  }
  for (n = 0 ; n < area->n_points ; n++)
    if (!area->lv[n].deleted)
    {
      nbytes = fprintf(fp, "%d  %2.3f  %2.3f  %2.3f %10.10f\n",
                       area->lv[n].vno, area->lv[n].x,
                       area->lv[n].y, area->lv[n].z, area->lv[n].stat) ;
      if(nbytes < 0)
      {
        printf("ERROR: writing to %s\n",fname);
        fclose(fp);
        return(1);
      }
    }
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

  LabelToFlat(area, mris) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    v->ripflag = 1 ;
  }

  for (n = 0 ; n < area->n_points ; n++)
  {
    vno = area->lv[n].vno ;
    if (vno < 0 || vno >= mris->nvertices)
    {
      continue ;
    }
    if (vno == Gdiag_no)
    {
      DiagBreak() ;
    }
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
LabelRipRestOfSurfaceWithThreshold(LABEL *area,
                                   MRI_SURFACE *mris,
                                   float thresh)
{
  int    vno, n ;
  VERTEX *v ;

  LabelToFlat(area, mris) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    v->ripflag = 1 ;
  }

  for (n = 0 ; n < area->n_points ; n++)
  {
    if (area->lv[n].stat < thresh)
    {
      continue ;
    }
    vno = area->lv[n].vno ;
    if (vno < 0 || vno >= mris->nvertices)
    {
      continue ;
    }
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
        Parameters: two labels to be intersected

        Returns value: error code and intersection in area1

        Description: elements are not removed, only the deleted
                     flag is set
------------------------------------------------------*/
int
LabelIntersect(LABEL *area1, LABEL *area2)
{
  int n, vno ;
  int vmin, vmax,vnum;
  int * isec;
  if (area1->n_points == 0)
  {
    return(NO_ERROR) ;
  }

  vmin = area1->lv[0].vno;
  vmax = area1->lv[0].vno;
  for (n = 0 ; n < area1->n_points ; n++)
  {
    vno = area1->lv[n].vno ;
    if ( vno < vmin)
    {
      vmin = vno;
    }
    if ( vno > vmax)
    {
      vmax = vno;
    }
  }
  for (n = 0 ; n < area2->n_points ; n++)
  {
    vno = area2->lv[n].vno ;
    if ( vno < vmin)
    {
      vmin = vno;
    }
    if ( vno > vmax)
    {
      vmax = vno;
    }
  }
  vnum = vmax - vmin;
  isec = (int*) calloc(vnum,sizeof(int));
  if (!isec)
    ErrorExit(ERROR_NOMEMORY,
              "%s: could not allocate LabelIntersect struct.",Progname);
  for (n = 0 ; n<vnum; n++)
  {
    isec[n] = -1;
  }

  for (n = 0 ; n < area1->n_points ; n++)
  {
    vno = area1->lv[n].vno - vmin ;
    isec[vno] = n;
    area1->lv[n].deleted = 1;
  }
  for (n = 0 ; n < area2->n_points ; n++)
  {
    vno = area2->lv[n].vno -vmin;
    if (isec[vno] > -1)
    {
      area1->lv[isec[vno]].deleted = 0;
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
  char   *cp, subjects_dir[STRLEN] ;

  area = (LABEL *)calloc(1, sizeof(LABEL)) ;
  if (!area)
  {
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate LABEL struct.",Progname);
  }
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
  {
    if (label_name)
    {
      strncpy(area->name, label_name, STRLEN-1) ;
    }
  }
  strcpy(area->space, "");

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
  {
    max_vertices = area->max_points ;
  }

  MRISclearMarks(mris) ;
  max_curv = min_curv = 0 ;
  for (n = 0 ; n < nvertices ; n++)
  {
    v = &mris->vertices[vertex_list[n]] ;
    if (v->ripflag)
    {
      continue ;
    }
    v->marked = 1 ;
    lv = &area->lv[n] ;
    lv->vno = vertex_list[n] ;
    lv->x = v->x ;
    lv->y = v->y ;
    lv->z = v->z ;
    area->n_points++ ;
    if (v->curv > max_curv)
    {
      max_curv = v->curv ;
    }
    if (v->curv < min_curv)
    {
      min_curv = v->curv ;
    }
  }
  if (-min_curv > max_curv)
  {
    curv_thresh = min_curv ;
  }
  else
  {
    curv_thresh = max_curv ;
  }
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
      {
        continue ;
      }
      for (nv = 0 ; nv < v->vnum ; nv++)   /* go through neighbors */
      {
        vn = &mris->vertices[v->v[nv]] ;
        if (vn->ripflag || vn->marked)
        {
          continue ;
        }
        if (((curv_thresh > 0) && (vn->curv > curv_thresh)) ||
            ((curv_thresh < 0) && (vn->curv < curv_thresh)))
        {
          vn->marked = 1 ;
          lvn = &area->lv[area->n_points+nfilled] ;
          nfilled++ ;
          lvn->vno = v->v[nv] ;
          lvn->x = vn->x ;
          lvn->y = vn->y ;
          lvn->z = vn->z ;
        }
        if (area->n_points+nfilled >= area->max_points)
        {
          break ;
        }
      }
      if (area->n_points+nfilled >= area->max_points)
      {
        break ;
      }
    }
    fprintf(stderr, "%d vertices added.\n", nfilled) ;
    area->n_points += nfilled ;
  }
  while ((nfilled > 0) && (area->n_points < max_vertices)) ;
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
      v = &mris->vertices[lv->vno] ;
      v->marked = 2 ;
      if (v->ripflag)
      {
        continue ;
      }
      for (nv = 0 ; nv < v->vnum ; nv++)   /* go through neighbors */
      {
        vn = &mris->vertices[v->v[nv]] ;
        if (vn->ripflag)
        {
          continue ;
        }
        if (vn->marked == 1)  /* add it to the label */
        {
          vn->marked = 2 ;
          lvn = &area->lv[area->n_points+nfilled] ;
          nfilled++ ;
          lvn->vno = v->v[nv] ;
          lvn->x = vn->x ;
          lvn->y = vn->y ;
          lvn->z = vn->z ;
        }
        if (area->n_points+nfilled >= area->max_points)
        {
          break ;
        }
      }
      if (area->n_points+nfilled >= area->max_points)
      {
        break ;
      }
    }
    /*    fprintf(stderr, "%d vertices added.\n", nfilled) ;*/
    area->n_points += nfilled ;
  }
  while ((nfilled > 0) && (area->n_points < area->max_points)) ;
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
      v = &mris->vertices[lv->vno] ;
      v->marked = 1 ;
      if (v->ripflag)
      {
        continue ;
      }
      for (nv = 0 ; nv < v->vnum ; nv++)   /* go through neighbors */
      {
        vn = &mris->vertices[v->v[nv]] ;
        if (vn->ripflag || vn->marked)
        {
          continue ;
        }
        if (vn->annotation == annotation)  /* add it to the label */
        {
          vn->marked = 1 ;
          lvn = &area->lv[area->n_points+nfilled] ;
          nfilled++ ;
          lvn->vno = v->v[nv] ;
          lvn->x = vn->x ;
          lvn->y = vn->y ;
          lvn->z = vn->z ;
        }
        if (area->n_points+nfilled >= area->max_points)
        {
          break ;
        }
      }
      if (area->n_points+nfilled >= area->max_points)
      {
        break ;
      }
    }
    /*    fprintf(stderr, "%d vertices added.\n", nfilled) ;*/
    area->n_points += nfilled ;
  }
  while ((nfilled > 0) && (area->n_points < area->max_points)) ;
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
  {
    max_vertices = area->max_points ;
  }

  for (n = 0 ; n < nvertices ; n++)
  {
    v = &mris->vertices[vertex_list[n]] ;
    if (v->ripflag)
    {
      continue ;
    }
    v->marked = 1 ;
    lv = &area->lv[n] ;
    lv->vno = vertex_list[n] ;
    lv->x = v->x ;
    lv->y = v->y ;
    lv->z = v->z ;
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
      {
        continue ;
      }
      for (nv = 0 ; nv < v->vnum ; nv++)   /* go through neighbors */
      {
        vn = &mris->vertices[v->v[nv]] ;
        if (vn->ripflag || vn->marked)
        {
          continue ;
        }
        vn->marked = 1 ;
        lvn = &area->lv[area->n_points+nfilled] ;
        nfilled++ ;
        lvn->vno = v->v[nv] ;
        lvn->x = vn->x ;
        lvn->y = vn->y ;
        lvn->z = vn->z ;
        if (area->n_points+nfilled >= area->max_points)
        {
          break ;
        }
      }
      if (area->n_points+nfilled >= area->max_points)
      {
        break ;
      }
    }
    fprintf(stderr, "%d vertices added.\n", nfilled) ;
    area->n_points += nfilled ;
  }
  while ((nfilled > 0) && (area->n_points < max_vertices)) ;
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
labelLoadTransform(const char *subject_name,
                   const char *sdir,
                   General_transform *transform)
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
  double   x, y, z, xt, yt, zt ;
  VERTEX *v ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    lv = &area->lv[n] ;
    v = &mris->vertices[lv->vno] ;
    x = v->origx ;
    y = v->origy ;
    z = v->origz ;
    transform_point(area->linear_transform, x, y, z, &xt, &yt, &zt) ;
    lv->x = xt ;
    lv->y = yt ;
    lv->z = zt ;
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
    lv->x = v->cx ;
    lv->y = v->cy ;
    lv->z = v->cz ;
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
  {
    adst = LabelClone(asrc) ;
  }
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
  {
    adst = LabelAlloc(asrc->n_points, asrc->subject_name, asrc->name) ;
  }

  adst->n_points = asrc->n_points ;
  strcpy(adst->name, asrc->name) ;
  strcpy(adst->subject_name, asrc->subject_name) ;

  memmove(adst->lv, asrc->lv, asrc->n_points*sizeof(LABEL_VERTEX)) ;
  return(adst) ;
}
/*-----------------------------------------------------
  int LabelRemoveDuplicates(LABEL *area)
  Sets the 'deleted' flag of a label point if it is
  a duplicate. Does not actually remove duplicats!
  ------------------------------------------------------*/
int LabelRemoveDuplicates(LABEL *area)
{
  int    n1, n2, deleted = 0 ;
  LV     *lv1, *lv2 ;

  // loop thru each label point
  for (n1 = 0 ; n1 < area->n_points ; n1++)
  {
    lv1 = &area->lv[n1] ;
    if(lv1->deleted)
    {
      continue ;
    }
    // loop thru the remaining looking for duplicates
    for (n2 = n1+1 ; n2 < area->n_points ; n2++)
    {
      lv2 = &area->lv[n2] ;
      if(lv1->vno >= 0 && lv2->vno >= 0 && lv1->vno == lv2->vno)
      {
        deleted++ ;
        lv2->deleted = 1 ;
      }
    }
  }
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "%d duplicate vertices removed from label %s.\n",
            deleted, area->name) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------*/
double LabelArea(LABEL *area, MRI_SURFACE *mris)
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
    {
      continue ;
    }
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

    x += lv->x ;
    y += lv->y ;
    z += lv->z ;
    n++ ;
  }

  *px = x / (double)n ;
  *py = y / (double)n ;
  *pz = z / (double)n ;
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

    xd = lv->x - ux ;
    yd = lv->y - uy ;
    zd = lv->z - uz ;
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
    if (vno < 0 || vno >= mris->nvertices)
    {
      DiagBreak() ;
    }
    if (area->lv[n].deleted > 0)
    {
      continue ;
    }
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
LabelAddToMark(LABEL *area, MRI_SURFACE *mris, int val_to_add)
{
  int    n, vno ;
  VERTEX *v ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    vno = area->lv[n].vno ;
    if (vno < 0 || vno >= mris->nvertices)
    {
      DiagBreak() ;
    }
    if (area->lv[n].deleted > 0)
      continue ;

    v = &mris->vertices[vno] ;
    v->marked += val_to_add ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
LabelMarkWithThreshold(LABEL *area, MRI_SURFACE *mris, float thresh)
{
  int    n, vno ;
  VERTEX *v ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    if (area->lv[n].stat < thresh)
    {
      continue ;
    }
    vno = area->lv[n].vno ;
    v = &mris->vertices[vno] ;
    v->marked = 1 ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters: A label and a surface

        Returns value: An error if either parameter is NULL

        Description

        For every point in the label that is not deleted (doesn't have
        it's deleted flag on), the corresponding vertex in the surface
        has its marked flag turned on.

------------------------------------------------------*/
int
LabelMarkUndeleted(LABEL *area, MRI_SURFACE *mris)
{
  int    n, vno ;
  VERTEX *v ;

  if (NULL == area || NULL == mris)
    ErrorReturn( ERROR_BADPARM,
                 (ERROR_BADPARM,
                  "LabelMarkUndeleted: area or mris was NULL" ));

  /* For every point, if it's not deleted, get the vno. Check the
     vno. If it's good, set the vertex's marked flag.*/
  for (n = 0 ; n < area->n_points ; n++)
  {
    if (!area->lv[n].deleted)
    {
      vno = area->lv[n].vno ;
      if (vno < 0 || vno >= mris->nvertices)
      {
        printf ("LabelMarkUndeleted: label point %d had invalid vno %d\n",
                n, vno);
        continue;
      }
      v = &mris->vertices[vno] ;
      v->marked = 1 ;
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
LabelUnmark(LABEL *area, MRI_SURFACE *mris)
{
  int    n, vno ;
  VERTEX *v ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    vno = area->lv[n].vno ;
    if (vno < 0)
    {
      continue ;
    }
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
    if (vno < 0 || vno >= mris->nvertices)
    {
      continue;
    }
    v = &mris->vertices[vno] ;
    area->lv[n].x = v->origx ;
    area->lv[n].y = v->origy ;
    area->lv[n].z = v->origz ;
  }
  strncpy (area->space, "coords=orig", sizeof(area->space));
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters: allocated LABEL, valid MRI_SURFACE

        Returns value: error code

        Description: For every point in the label, gets the
                     corresponding white coords from the surface and
                     writes them in the label.
 ------------------------------------------------------*/
int
LabelToWhite(LABEL *area, MRI_SURFACE *mris)
{
  int     n, vno ;
  VERTEX  *v ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    vno = area->lv[n].vno ;
    if (vno < 0 || vno >= mris->nvertices)
    {
      continue;
    }
    v = &mris->vertices[vno] ;
    area->lv[n].x = v->whitex ;
    area->lv[n].y = v->whitey ;
    area->lv[n].z = v->whitez ;
  }
  strncpy (area->space, "coords=white", sizeof(area->space));
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
    {
      npoints++ ;
    }

  if (!npoints)
  {
    return(NULL) ;
  }
  area = LabelAlloc(npoints, NULL, NULL) ;
  for (n = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (!v->marked)
    {
      continue ;
    }
    if (vno == Gdiag_no)
    {
      DiagBreak() ;
    }
    area->lv[n].x = v->x ;
    area->lv[n].y = v->y ;
    area->lv[n].z = v->z ;
    area->lv[n].vno = vno ;
    area->lv[n].stat = v->stat ;
    if (FZERO(v->stat))
    {
      DiagBreak() ;
    }
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
LABEL *
LabelFromMarkValue(MRI_SURFACE *mris, int mark)
{
  int    vno, npoints, n ;
  LABEL  *area ;
  VERTEX *v ;

  for (npoints = vno = 0 ; vno < mris->nvertices ; vno++)
    if (mris->vertices[vno].marked == mark)
    {
      npoints++ ;
    }

  if (!npoints)
  {
    return(NULL) ;
  }
  area = LabelAlloc(npoints, NULL, NULL) ;
  for (n = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->marked != mark)
    {
      continue ;
    }
    if (vno == Gdiag_no)
    {
      DiagBreak() ;
    }
    area->lv[n].x = v->x ;
    area->lv[n].y = v->y ;
    area->lv[n].z = v->z ;
    area->lv[n].vno = vno ;
    area->lv[n].stat = v->stat ;
    if (FZERO(v->stat))
    {
      DiagBreak() ;
    }
    n++ ;
  }
  area->n_points = npoints ;
  return(area) ;
}
/*-----------------------------------------------------
------------------------------------------------------*/
int
LabelAddToSurfaceMark(LABEL *area, MRI_SURFACE *mris, int mark_to_add) 
{
  int     n, vno ;
  VERTEX  *v ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    vno = area->lv[n].vno ;
    if(vno < 0)
      continue ;

    if(vno >= mris->nvertices)
    {
      printf("ERROR: LabelMarkSurface: label point %d exceeds nvertices %d\n",
             vno,mris->nvertices);
      return(1);
    }
    v = &mris->vertices[vno] ;
    if(v->ripflag)
      continue ;

    v->marked += mark_to_add ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
------------------------------------------------------*/
int LabelMarkSurface(LABEL *area, MRI_SURFACE *mris)
{
  int     n, vno ;
  VERTEX  *v ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    vno = area->lv[n].vno ;
    if(vno < 0)
    {
      continue ;
    }
    if(vno >= mris->nvertices)
    {
      printf("ERROR: LabelMarkSurface: label point %d exceeds nvertices %d\n",
             vno,mris->nvertices);
      return(1);
    }
    v = &mris->vertices[vno] ;
    if(v->ripflag)
    {
      continue ;
    }
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
LabelIsCompletelyUnassigned(LABEL *area, int *unassigned)
{
  int i;

  if (NULL == area)
  {
    return ERROR_BADPARM;
  }
  if (NULL == unassigned)
  {
    return ERROR_BADPARM;
  }

  /* Go through the label. If a point vno is _not_ -1, return false,
     because we have an assigned point. If we get through the entire
     label without finding a vno that isn't -1, return true. */
  for (i = 0 ; i < area->n_points ; i++)
  {
    if (area->lv[i].vno != -1)
    {
      *unassigned = 0;
      return(NO_ERROR);
    }
  }

  *unassigned = 1;
  return(NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#include "mrishash.h"
int
LabelFillUnassignedVertices(MRI_SURFACE *mris, LABEL *area, int coords)
{
  int     n, i, vno, min_vno, nfilled = 0, max_vno ;
  LV      *lv ;
  MHT     *mht ;
  MHBT    *bucket ;
  MHB     *bin ;
  VERTEX  *v ;
  float   dx, dy, dz, x, y, z, dist, min_dist ;
  int     num_not_found;
  float   vx, vy, vz;
  double  max_spacing ;
  vx = vy = vz = -1;

  MRIScomputeVertexSpacingStats(mris, NULL, NULL, &max_spacing,
                                NULL,&max_vno, coords);

  for (i = n = 0 ; n < area->n_points ; n++)
  {
    lv = &area->lv[n] ;
    if (lv->vno >= 0 && lv->vno < mris->nvertices)
    {
      continue ;
    }
    i++ ;   /* count # of unassigned vertices */
  }
  if (i <= 0)
  {
    return(0) ;  /* no work needed */
  }

  fprintf(stderr,"%d unassigned vertices in label - building spatial LUT...\n",
          i) ;

  /* if we can't find a vertex within 10 mm of the point, something is wrong */
  mht = MHTfillVertexTableRes(mris, NULL, coords, 2*max_spacing) ;
  fprintf(stderr, "assigning vertex numbers to label...\n") ;
  num_not_found = 0;
  for (n = 0 ; n < area->n_points ; n++)
  {
    lv = &area->lv[n] ;
    if (lv->vno >= 0 && lv->vno <= mris->nvertices)
    {
      continue ;
    }
    nfilled++ ;
    bucket = MHTgetBucket(mht, lv->x, lv->y, lv->z) ;

    x = lv->x ;
    y = lv->y ;
    z = lv->z ;
    min_dist = 10000.0 ;
    min_vno = -1 ;
    if (bucket)
    {
      for (bin = bucket->bins, i = 0 ; i < bucket->nused ; i++, bin++)
      {
        vno = bin->fno ;
        v = &mris->vertices[vno] ;
        if (vno == Gdiag_no)
        {
          DiagBreak() ;
        }

        MRISgetCoords(v, coords, &vx, &vy, &vz);

        dx = vx - x;
        dy = vy - y;
        dz = vz - z;

        dist = sqrt(dx*dx + dy*dy + dz*dz) ;
        if (dist < min_dist)
        {
          min_dist = dist ;
          min_vno = vno ;
        }
      }
    }
    if (min_vno == -1)
    {
      switch (coords)
      {
      case ORIG_VERTICES:
        min_vno = MRISfindClosestOriginalVertex(mris, lv->x, lv->y, lv->z) ;
        break ;
      case WHITE_VERTICES:
        min_vno = MRISfindClosestWhiteVertex(mris, lv->x, lv->y, lv->z) ;
        break ;
      case CURRENT_VERTICES:
        min_vno = MRISfindClosestVertex(mris, lv->x, lv->y, lv->z, NULL) ;
        break ;
      case CANONICAL_VERTICES:
        min_vno = MRISfindClosestCanonicalVertex(mris, lv->x, lv->y, lv->z) ;
        break ;
      default:
        break ;
      }
    }
    if (min_vno == -1)
    {
      num_not_found++;
    }
    lv->vno = min_vno ;
  }
  LabelRemoveDuplicates(area) ;
  MHTfree(&mht) ;

  if (num_not_found > 0)
  {
    fprintf (stderr, "Couldn't assign %d vertices.\n", num_not_found);
  }


  return(nfilled) ;
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
  {
    adst = LabelClone(asrc) ;
  }

  if (adst->max_points < asrc->n_points+adst->n_points)/* won't fit - expand */
  {
    LABEL   *atmp ;

    atmp = LabelAlloc(2*(asrc->n_points+adst->n_points),
                      asrc->subject_name, asrc->name) ;
    LabelCopy(adst, atmp) ;
    LabelFree(&adst) ;
    adst = atmp ;
  }

  MRISclearMarks(mris_dst) ;
  for (n = 0 ; n < asrc->n_points ; n++)
  {
    vno = asrc->lv[n].vno ;
    if (vno == Gdiag_no)
    {
      DiagBreak() ;
    }
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    vdst = MHTfindClosestVertex(mht, mris_dst, v) ;
    if (!vdst)
    {
      ErrorPrintf(ERROR_BADPARM, "MRIScombine: cannot map vno %d", vno) ;
      continue ;
    }
    if (vdst->marked)  /* only add this vertex once */
    {
      continue ;
    }
    vdst->marked = 1 ;
    if (vdst-mris_dst->vertices == Gdiag_no)
    {
      DiagBreak() ;
    }
    lv_dst = labelFindVertexNumber(adst, vdst-mris_dst->vertices) ;
    if (lv_dst == NULL)
    {
      lv_dst = &adst->lv[adst->n_points++] ;
    }
    else if (lv_dst-adst->lv == Gdiag_no)
    {
      DiagBreak() ;
    }
    lv_dst->vno = vdst-mris_dst->vertices ;
    if (lv_dst->vno == 60008)
    {
      DiagBreak() ;
    }
    if (lv_dst->vno == 85592)
    {
      DiagBreak() ;
    }
    if (lv_dst->vno == Gdiag_no)
    {
      DiagBreak() ;
    }
    lv_dst->x = asrc->lv[n].x ;
    lv_dst->y = asrc->lv[n].y ;
    lv_dst->z = asrc->lv[n].z ;
    lv_dst->stat += asrc->lv[n].stat ;
  }

  MRIScomputeVertexSpacingStats(mris, NULL, NULL, &max_len, NULL,
                                NULL,CURRENT_VERTICES);
  mht_src =
    MHTfillVertexTableRes(mris, NULL, CURRENT_VERTICES, 2*max_len);

  MRISclearMarks(mris) ;
  LabelMarkStats(asrc, mris) ;

  do
  {
    /* map the nbr of every point in the label back to src, and if in
       label add it
    */
    nfilled = 0 ;
    for (n = 0 ; n < adst->n_points ; n++)
    {
      v = &mris_dst->vertices[adst->lv[n].vno] ;
      if (adst->lv[n].vno == Gdiag_no)
      {
        DiagBreak() ;
      }
      if (n == Gdiag_no)
      {
        DiagBreak() ;
      }
      if (v->marked == 0)   /* hasn't been processed for this surface yet */
      {
        vsrc = MHTfindClosestVertex(mht_src, mris, v) ;
        if (vsrc && vsrc->marked)   /* in label */
        {
          adst->lv[n].stat += vsrc->stat ;
          v->marked = 1 ;
        }
      }
      for (m = 0 ; m < v->vnum ; m++)
      {
        vn = &mris_dst->vertices[v->v[m]] ;
        if (vn->marked)
        {
          continue ;  /* already in label */
        }
        vsrc = MHTfindClosestVertex(mht_src, mris, vn) ;
        if (vsrc == NULL)
        {
          DiagBreak() ;
        }
        if (vsrc - mris->vertices == 62644)
        {
          DiagBreak() ;
        }
        if (vsrc - mris->vertices == Gdiag_no)
        {
          DiagBreak() ;
        }
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
          lv_dst = labelFindVertexNumber(adst, v->v[m]) ;
          if (lv_dst == NULL)
          {
            lv_dst = &adst->lv[adst->n_points++] ;
          }
          else if (lv_dst-adst->lv == Gdiag_no)
          {
            DiagBreak() ;
          }
          lv_dst->x = adst->lv[n].x ;
          lv_dst->y = adst->lv[n].y ;
          lv_dst->z = adst->lv[n].z ;
          lv_dst->vno = v->v[m] ;
          lv_dst->stat += vsrc->stat ;
          if (lv_dst->vno == Gdiag_no)
          {
            DiagBreak() ;
          }
          if (lv_dst->vno == 55185)
          {
            DiagBreak() ;
          }
          if (lv_dst->vno == 85592)
          {
            DiagBreak() ;
          }
          nfilled++ ;
        }
      }
    }
  }
  while (nfilled != 0) ;

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

  if (frame >= SurfMask->nframes)
  {
    printf("ERROR:  MaskSurfLabel: specified frame %d is greater "
           "than the number of frame %d\n",frame, SurfMask->nframes);
    return(NULL);
  }

  /* Check that the label is a surface-based label */
  ok = 0;
  for (n = 0; n < lbl->n_points; n++)
  {
    if (lbl->lv[n].vno != 0)
    {
      ok = 1;
    }
  }
  if (!ok)
  {
    printf("ERROR: MaskSurfLabel: all label vertices are 0, proably\n");
    printf("       not a label created from a surface\n");
    return(NULL);
  }

  /* -- Allocate the Target Label ---*/
  msklbl = LabelAlloc(lbl->n_points,lbl->subject_name,"labelfilenamehere");
  msklbl->n_points = lbl->n_points;

  /* Compute the number of overlapping vertices */
  noverlap = 0;
  for (n = 0; n < lbl->n_points; n++)
  {
    vno = lbl->lv[n].vno;
    if (vno >= SurfMask->width)
    {
      printf("ERROR: MaskSurfLabel: label vertex %d has vno %d which "
             "is larger than the number of vertices %d in mask\n",
             n,vno,SurfMask->width);
      return(NULL);
    }
    maskval = MRIgetVoxVal(SurfMask,vno,0,0,frame);
    if (!strcmp(masksign,"abs"))
    {
      maskval = fabs(maskval);
    }
    if (!strcmp(masksign,"neg"))
    {
      maskval = -maskval;
    }
    //printf("%4d  %6d  %g %g  %d\n",n,vno,maskval,thresh,maskval>thresh);
    if (maskval < thresh)
    {
      continue;
    }
    msklbl->lv[noverlap].vno = vno;
    msklbl->lv[noverlap].x = lbl->lv[noverlap].x;
    msklbl->lv[noverlap].y = lbl->lv[noverlap].y;
    msklbl->lv[noverlap].z = lbl->lv[noverlap].z;
    msklbl->lv[noverlap].stat = lbl->lv[noverlap].stat;
    noverlap ++;
  }
  if (noverlap==0)
  {
    printf("WARNING: MaskSurfLabel: no overlap between label and mask\n");
    return(NULL);
  }
  msklbl->n_points = noverlap;

  return(msklbl);
}

int
LabelErode(LABEL *area, MRI_SURFACE *mris, int num_times)
{
  int n, num_new_lvs, label_vno, vno, neighbor_index, neighbor_vno, found;
  LV* new_lv;
  VERTEX *vn ;

  if (NULL == area)
  {
    ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM,"LabelErode: NULL label"));
  }
  if (NULL == mris)
  {
    ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM,"LabelErode: NULL mris"));
  }
  if (num_times < 1)
  {
    ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM,"LabelErode: num_times < 1"));
  }

  for (n = 0; n < num_times; n++)
  {
    /* Create an array of LVs the same size as the current label. */
    new_lv = (LV*) calloc( area->n_points, sizeof(LV) );
    if (NULL == new_lv)
      ErrorReturn(ERROR_NOMEMORY,(ERROR_NOMEMORY,
                                  "LabelErode: couldn't allocate new_lv"));
    num_new_lvs = 0;

    MRISclearMarks(mris) ;
    LabelMark(area, mris) ; // all vertices in label now have v->marked==1

    /* For each vertex in the label... */
    for (label_vno = 0; label_vno < area->n_points; label_vno++)
    {
      vno = area->lv[label_vno].vno;
      if (vno == Gdiag_no)
      {
        DiagBreak() ;
      }

      // check to see if we should not add this label
      // (if one of it's nbrs is not in label)
      found = 0 ;
      for (neighbor_index = 0;
           neighbor_index < mris->vertices[vno].vnum; neighbor_index++)
      {
        neighbor_vno = mris->vertices[vno].v[neighbor_index];
        if (neighbor_vno == Gdiag_no)
        {
          DiagBreak() ;
        }

        /* Look for neighbor_vno in the label. */
        vn = &mris->vertices[neighbor_vno] ;
        if (vn->marked == 0) // found a nbr not in the label
        {
          if (neighbor_vno == Gdiag_no)
          {
            DiagBreak() ;
          }
          found = 1;
          break;
        }
      }

      if (found == 0) // all nbrs on - add it to the new label
      {
        memmove (&new_lv[num_new_lvs], &area->lv[label_vno],sizeof(LV));
        num_new_lvs++;
      }
    }

    /* Point the label's lv to the new one and update the number of
    points. */
    free (area->lv);
    area->lv = (LV*) realloc (new_lv, num_new_lvs * sizeof(LV));
    area->n_points = num_new_lvs;
    area->max_points = num_new_lvs;
  }

  MRISclearMarks(mris) ;
  return (NO_ERROR);
}

int
LabelDilate(LABEL *area, MRI_SURFACE *mris, int num_times)
{
  int    n, num_new_lvs, vno, neighbor_index, neighbor_vno,found;
  LV     *new_lv;
  VERTEX *vn, *v ;

  if (NULL == area)
  {
    ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM,"LabelDilate: NULL label"));
  }
  if (NULL == mris)
  {
    ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM,"LabelDilate: NULL mris"));
  }
  if (num_times < 1)
  {
    ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM,"LabelDilate: num_times < 1"));
  }

  for (n = 0; n< num_times; n++ )
  {
    MRISclearMarks(mris) ;
    LabelMarkStats(area, mris) ; // all vertices in label now have v->marked==1

    /* Allocate an LV array the size of the surface. */
    new_lv = (LV*) calloc( mris->nvertices, sizeof(LV) );
    if (NULL == new_lv)
      ErrorReturn(ERROR_NOMEMORY,(ERROR_NOMEMORY,
                                  "LabelDilate: couldn't allocate new_lv"));
    num_new_lvs = 0;

    /* Copy the existing lvs over first and increment our count. */
    memmove (new_lv, area->lv, area->n_points * sizeof(LV));
    num_new_lvs = area->n_points;


    /* For each vertex in the label... */
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (vno == Gdiag_no)
      {
        DiagBreak() ;
      }
      if (v->marked == 1) // already in label
      {
        continue ;
      }

      // Check its neighbors. If any are in the label, add it
      found = 0 ;
      for (neighbor_index = 0; neighbor_index < v->vnum; neighbor_index++)
      {
        /* Look for neighbor_vno in the label. */
        neighbor_vno = mris->vertices[vno].v[neighbor_index];
        vn = &mris->vertices[neighbor_vno] ;
        if (vn->marked > 0)
        {
          if (neighbor_vno == Gdiag_no)
          {
            DiagBreak() ;
          }
          found = 1;
          break;
        }
      }

      // add it if at least one nbr was found that was in label
      if (found)
      {
        if (vno == Gdiag_no)
        {
          DiagBreak() ;
        }
        new_lv[num_new_lvs].vno = vno ;
        new_lv[num_new_lvs].x = v->x;
        new_lv[num_new_lvs].y = v->y;
        new_lv[num_new_lvs].z = v->z;
        new_lv[num_new_lvs].stat = v->stat;
        num_new_lvs++;
      }
    }

    /* Point the label's lv to the new one and update the number of
      points. */
    free (area->lv);
    area->lv = (LV*) realloc (new_lv, num_new_lvs * sizeof(LV));
    area->n_points = num_new_lvs;
    area->max_points = num_new_lvs;
  }

  return (NO_ERROR);
}


static LABEL_VERTEX *
labelFindVertexNumber(LABEL *area, int vno)
{
  int          n ;
  LABEL_VERTEX *lv ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    lv = &area->lv[n] ;
    if (lv->vno == vno)
    {
      return(lv) ;
    }
  }
  return(NULL) ;
}

int
LabelSetStat(LABEL *area, float stat)
{
  int  n ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    area->lv[n].stat = stat ;
  }

  return(NO_ERROR) ;
}

int
LabelMarkStats(LABEL *area, MRI_SURFACE *mris)
{
  int    n, vno ;
  VERTEX *v ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    if (area->lv[n].deleted > 0)
    {
      continue ;
    }
    vno = area->lv[n].vno ;
    v = &mris->vertices[vno] ;
    v->marked = 1 ;
    v->stat = area->lv[n].stat ;
  }
  return(NO_ERROR) ;
}

LABEL *
LabelFillHoles(LABEL *area_src, MRI_SURFACE *mris, int coords)
{
  MRI    *mri ;
  int    i, dst_index, vno, xi, yi, zi, xk, yk, zk, found, x2, y2, z2,nchanged;
  double xw, yw, zw, xv, yv, zv ;
  VERTEX *v ;
  LABEL  *area_dst ;
  float  vx, vy, vz, dist, dx, dy, dz;

  vx = vy = vz = -1;
  mri = MRIalloc(256,256,256,MRI_UCHAR) ;
  
  useVolGeomToMRI(&mris->vg, mri);
  area_dst = LabelAlloc(mris->nvertices, mris->subject_name, area_src->name) ;
  LabelCopy(area_src, area_dst) ;
  LabelFillUnassignedVertices(mris, area_dst, coords) ;
  LabelMarkSurface(area_dst, mris) ;

  for (i = 0 ; i < area_src->n_points ; i++)
  {
    xw = area_src->lv[i].x  ;
    yw = area_src->lv[i].y  ;
    zw = area_src->lv[i].z ;
    // MRIworldToVoxel(mri, xw, yw, zw, &xv, &yv, &zv) ;
#if 0
    if (mris->useRealRAS)
    {
      MRIworldToVoxel(mri, xw, yw, zw, &xv, &yv, &zv) ;
    }
    else
    {
      MRIsurfaceRASToVoxel(mri, xw, yw, zw, &xv, &yv, &zv) ;
    }
#else
    MRISsurfaceRASToVoxel(mris, mri, xw, yw, zw, &xv, &yv, &zv) ;
#endif
    MRIvox(mri, nint(xv), nint(yv), nint(zv)) = 1 ;
  }

  dst_index = area_dst->n_points ;
  do
  {
    nchanged = 0 ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (vno == Gdiag_no)
      {
        DiagBreak() ;
      }
      if (v->ripflag || v->marked)   /* already in label */
      {
        continue ;
      }
      MRISvertexCoordToVoxel(mris, v, mri, coords, &xv, &yv, &zv) ;
      xi = nint(xv) ;
      yi = nint(yv) ;
      zi = nint(zv) ;
      if (xi >= mri->width)
      {
        xi = mri->width-1;
      }
      if (yi >= mri->height)
      {
        yi = mri->height-1;
      }
      if (zi >= mri->depth)
      {
        zi = mri->depth-1;
      }
      if (MRIvox(mri, xi, yi, zi) == 1)
      {
        MRISgetCoords(v, coords, &vx, &vy, &vz);
        area_dst->lv[dst_index].vno = vno ;
        area_dst->lv[dst_index].x = vx ;
        area_dst->lv[dst_index].y = vy ;
        area_dst->lv[dst_index].z = vz ;
        v->marked = 1 ;
        dst_index++ ;
        nchanged++ ;
        area_dst->n_points++ ;
      }
      if (MRIneighborsOn(mri, xi, yi, zi, 1) > 0)
      {
        found = 0 ;
        for (xk = -1 ; !found && xk <= 1 ; xk++)
          for (yk = -1 ; !found && yk <= 1 ; yk++)
            for (zk = -1 ; !found && zk <= 1 ; zk++)
            {
              x2 = mri->xi[xi+xk] ;
              y2 = mri->xi[yi+yk] ;
              z2 = mri->xi[zi+zk] ;
              if (MRIvox(mri, x2, y2, z2) == 1)
              {
                dx = xv-(xi+xk) ;
                dy = yv-(yi+yk) ;
                dz = zv-(zi+zk) ;
                dist = sqrt(SQR(dx)+SQR(dy)+SQR(dz)) ;
                if (dist < sqrt(.5*.5+.5*.5))
                {
                  found = 1 ;
                  MRISgetCoords(v, coords, &vx, &vy, &vz);
                  area_dst->lv[dst_index].vno = vno ;
                  area_dst->lv[dst_index].x = vx ;
                  area_dst->lv[dst_index].y = vy ;
                  area_dst->lv[dst_index].z = vz ;
                  dst_index++ ;
                  v->marked = 1 ;
                  nchanged++ ;
                  area_dst->n_points++ ;
                  break ;
                }
              }
            }
      }
    }
  }
  while (nchanged > 0) ;

  MRIfree(&mri) ;
  return(area_dst) ;
}
LABEL *
LabelFillHolesWithOrig(LABEL *area_src, MRI_SURFACE *mris)
{
  MRI    *mri ;
  int    i, dst_index, vno ;
  double xw, yw, zw, xv, yv, zv ;
  VERTEX *v ;
  LABEL  *area_dst ;

  mri = MRIalloc(256,256,256,MRI_UCHAR) ;
  area_dst = LabelAlloc(mris->nvertices, mris->subject_name, area_src->name) ;
  LabelCopy(area_src, area_dst) ;
  LabelFillUnassignedVertices(mris, area_dst, ORIG_VERTICES) ;
  LabelMarkSurface(area_dst, mris) ;

  for (i = 0 ; i < area_src->n_points ; i++)
  {
    xw = area_src->lv[i].x  ;
    yw = area_src->lv[i].y  ;
    zw = area_src->lv[i].z  ;
    // MRIworldToVoxel(mri, xw, yw, zw, &xv, &yv, &zv) ;
    if (mris->useRealRAS)
    {
      MRIworldToVoxel(mri, xw, yw, zw, &xv, &yv, &zv) ;
    }
    else
    {
      MRIsurfaceRASToVoxel(mri, xw, yw, zw, &xv, &yv, &zv) ;
    }
    MRIvox(mri, nint(xv), nint(yv), nint(zv)) = 1 ;
  }

  dst_index = area_dst->n_points ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->marked)   /* already in label */
    {
      continue ;
    }
    MRISorigVertexToVoxel(mris, v, mri, &xv, &yv, &zv) ;
    if (MRIvox(mri, nint(xv), nint(yv), nint(zv)) == 1)
    {
      area_dst->lv[dst_index].vno = vno ;
      area_dst->lv[dst_index].x = v->origx ;
      area_dst->lv[dst_index].y = v->origy ;
      area_dst->lv[dst_index].z = v->origz ;
      dst_index++ ;
      area_dst->n_points++ ;
    }
  }

  MRIfree(&mri) ;
  return(area_dst) ;
}

/*---------------------------------------------------------------
  LabelHasVertex() - returns -1 if the vertex is not in the label,
  otherwise returns the number of the label point that corresponds
  to the vertex number.
  ---------------------------------------------------------------*/
int LabelHasVertex(int vtxno, LABEL *lb)
{
  int n;
  for (n = 0; n < lb->n_points; n++)
    if (lb->lv[n].vno == vtxno)
    {
      return(n);
    }
  return(-1);
}
/*---------------------------------------------------------------
  LabelRealloc() - reallocates the number of label vertices
  to have max_points. If something goes wrong, returns 1 without
  changing anything. Otherwise returns 0.
  ---------------------------------------------------------------*/
int LabelRealloc(LABEL *lb, int max_points)
{
  LV *lvtmp;

  if (max_points <= lb->max_points)
  {
    return(0);
  }

  lvtmp = realloc(lb->lv,sizeof(LV)*max_points);
  if (lvtmp == NULL)
  {
    return(1);
  }
  lb->max_points = max_points;
  lb->lv = lvtmp;

  return(0);
}
/*------------------------------------------------------------
  LabelfromASeg() - creates a label from a segmentation given
  the numeric segmentation code. If there are no voxels that
  match the code, then it silently returns NULL. Note: this
  generates the same result as the code in mri_cor2label.c.
  ------------------------------------------------------------*/
LABEL *LabelfromASeg(MRI *aseg, int segcode)
{
  int c,r,s,nlabel,v;
  LABEL *lb;
  MATRIX *Vox2RAS, *crs, *ras;
  //char labelfile[100];

  // Count number of label points first
  nlabel = 0;
  for (c=0; c < aseg->width; c++)
  {
    for (r=0; r < aseg->height; r++)
    {
      for (s=0; s < aseg->depth; s++)
      {
        v = (int)MRIgetVoxVal(aseg,c,r,s,0);
        if (v == segcode)
        {
          nlabel++;
        }
      }
    }
  }
  //printf("Found %d voxels in label\n",nlabel);
  if (nlabel==0)
  {
    return(NULL);
  }

  lb = LabelAlloc(nlabel,"dontknow","dontcare");

  Vox2RAS = MRIxfmCRS2XYZtkreg(aseg);
  //printf("ASeg Vox2RAS---------------------------\n");
  //MatrixPrint(stdout,Vox2RAS);
  //printf("---------------------------\n");
  crs = MatrixAlloc(4,1,MATRIX_REAL);
  crs->rptr[4][1] = 1;
  ras = MatrixAlloc(4,1,MATRIX_REAL);

  nlabel = 0;
  for (c=0; c < aseg->width; c++)
  {
    for (r=0; r < aseg->height; r++)
    {
      for (s=0; s < aseg->depth; s++)
      {
        v = (int)MRIgetVoxVal(aseg,c,r,s,0);
        if (v == segcode)
        {
          crs->rptr[1][1] = c;
          crs->rptr[2][1] = r;
          crs->rptr[3][1] = s;
          ras = MatrixMultiply(Vox2RAS,crs,ras);
          lb->lv[nlabel].x = ras->rptr[1][1];
          lb->lv[nlabel].y = ras->rptr[2][1];
          lb->lv[nlabel].z = ras->rptr[3][1];
	  lb->lv[nlabel].vno = -1 ;   // not assigned yet
          nlabel++;
        }
      }
    }
  }
  lb->n_points = nlabel;
  MatrixFree(&Vox2RAS);
  MatrixFree(&crs);
  MatrixFree(&ras);

  //sprintf(labelfile,"./tmp-%d.label",segcode);
  //LabelWrite(lb,labelfile);

  return(lb);
}

/*---------------------------------------------------------------
  VertexIsInLabel() - returns a 1 if the given vertex number is
  in the label. Label must be surface-based.
  ---------------------------------------------------------------*/
int VertexIsInLabel(int vtxno, LABEL *label)
{
  int n;

  for (n=0; n<label->n_points; n++)
    if (label->lv[n].vno == vtxno)
    {
      return(1);
    }

  return(0);
}
/*---------------------------------------------------------------
  LabelBoundary() - returns a label of all the points in the input
  label on the edge/boundary. Label must be surface-based.
  ---------------------------------------------------------------*/
LABEL *LabelBoundary(LABEL *label, MRIS *surf)
{
  int n,nnbrs,nthnbr,vtxno,nbrvtxno;
  LABEL *boundary;

  boundary = LabelAlloc(label->n_points,"","");
  boundary->n_points=0;

  for (n=0; n<label->n_points; n++)
  {
    vtxno = label->lv[n].vno;
    nnbrs = surf->vertices[vtxno].vnum;
    for (nthnbr = 0; nthnbr < nnbrs; nthnbr++)
    {
      nbrvtxno = surf->vertices[vtxno].v[nthnbr];
      if (! VertexIsInLabel(nbrvtxno, label))
      {
        boundary->lv[boundary->n_points].vno = vtxno;
        boundary->lv[boundary->n_points].x = label->lv[n].x;
        boundary->lv[boundary->n_points].y = label->lv[n].y;
        boundary->lv[boundary->n_points].z = label->lv[n].z;
        boundary->lv[boundary->n_points].stat = label->lv[n].stat;
        boundary->n_points ++;
        break;
      }
    }
  }

  return(boundary);
}

/*---------------------------------------------------------------
  LabelFitXYZ() - fits a plane to the xyz of the label. If order=1,
  then the shape is a plane. If order=2, then it is a quad, 3 for
  cubic. Returns the regression coefficients.
  ---------------------------------------------------------------*/
MATRIX *LabelFitXYZ(LABEL *label, int order)
{
  int n;
  MATRIX *X, *Xt, *XtX, *iXtX, *iXtXXt, *y, *beta;

  if (order > 3)
  {
    printf("ERROR: LabelFitXYZ(): order = %d, must be <= 3\n",order);
    return(NULL);
  }

  X = MatrixAlloc(label->n_points,1+2*order,MATRIX_REAL);
  y = MatrixAlloc(label->n_points,1,MATRIX_REAL);

  // Load the design matrix
  for (n=0; n<label->n_points; n++)
  {
    X->rptr[n+1][1] = 1; // Intercept
    X->rptr[n+1][2] = label->lv[n].x; // Slope with x
    X->rptr[n+1][3] = label->lv[n].y; // Slope with y
    if (order > 1)
    {
      X->rptr[n+1][4] = pow(label->lv[n].x,2); // Slope with x^2
      X->rptr[n+1][5] = pow(label->lv[n].y,2); // Slope with y^2
    }
    if (order > 2)
    {
      X->rptr[n+1][6] = pow(label->lv[n].x,3); // Slope with x^3
      X->rptr[n+1][7] = pow(label->lv[n].y,3); // Slope with y^3
    }
    y->rptr[n+1][1] = label->lv[n].z; // Fit z
  }

  Xt     = MatrixTranspose(X,NULL);
  XtX    = MatrixMultiply(Xt,X,NULL);
  iXtX   = MatrixInverse(XtX,NULL);
  iXtXXt = MatrixMultiply(iXtX,Xt,NULL);
  beta   = MatrixMultiply(iXtXXt,y,NULL);

  MatrixFree(&X);
  MatrixFree(&Xt);
  MatrixFree(&XtX);
  MatrixFree(&iXtX);
  MatrixFree(&iXtXXt);
  MatrixFree(&y);

  //printf("beta -----------------\n");
  //MatrixPrint(stdout,beta);
  //printf("----------------------\n");

  return(beta);
}
LABEL *
LabelInFOV(MRI_SURFACE *mris, MRI *mri, float pad)
{
  LABEL        *area ;
  int          vno, npoints ;
  VERTEX       *v ;
  LABEL_VERTEX *lv ;
  double       xv, yv, zv, nvox ;

  // convert to voxels
  nvox = pad / MAX(MAX(mri->xsize, mri->ysize), mri->zsize) ;

  for (npoints = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    if (vno == Gdiag_no)
    {
      DiagBreak() ;
    }
#if 0
    if (mris->useRealRAS)
    {
      MRIworldToVoxel(mri, v->x, v->y, v->z, &xv, &yv, &zv);
    }
    else
    {
      MRIsurfaceRASToVoxel(mri, v->x, v->y, v->z, &xv, &yv, &zv);
    }
#else
    MRISsurfaceRASToVoxel(mris, mri, v->x, v->y, v->z, &xv, &yv, &zv) ;
#endif
    if (xv < nvox || yv < nvox || zv < nvox ||
        xv > (mri->width-1)-nvox ||
        yv > (mri->height-1)-nvox ||
        zv > (mri->depth-1)-nvox)
    {
      continue ;
    }
    npoints++ ;
  }

  area = LabelAlloc(npoints, NULL, NULL) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    if (vno == Gdiag_no)
    {
      DiagBreak() ;
    }
#if 0
    if (mris->useRealRAS)
    {
      MRIworldToVoxel(mri, v->x, v->y, v->z, &xv, &yv, &zv);
    }
    else
    {
      MRIsurfaceRASToVoxel(mri, v->x, v->y, v->z, &xv, &yv, &zv);
    }
#else
    MRISsurfaceRASToVoxel(mris, mri, v->x, v->y, v->z, &xv, &yv, &zv) ;
#endif
    if (xv < nvox || yv < nvox || zv < nvox ||
        xv > (mri->width-1)-nvox ||
        yv > (mri->height-1)-nvox ||
        zv > (mri->depth-1)-nvox)
    {
      continue ;
    }
    lv = &area->lv[area->n_points] ;
    area->n_points++ ;
    lv->vno = vno ;
    lv->x = v->x ;
    lv->y = v->y ;
    lv->z = v->z ;
  }

  // LabelWrite(area, "./lh.all.label");
  return(area) ;
}
LABEL *
LabelTranslate(LABEL *area, LABEL *area_offset, float dx, float dy, float dz)
{
  int i ;

  if (area_offset == NULL)
  {
    area_offset = LabelCopy(area, NULL) ;
  }


  for (i = 0 ; i < area->n_points ; i++)
  {
    area_offset->lv[i].x = area->lv[i].x  + dx ;
    area_offset->lv[i].y = area->lv[i].y  + dy ;
    area_offset->lv[i].z = area->lv[i].z  + dz ;
  }
  return(area_offset) ;
}

int
LabelUnassign(LABEL *area)
{
  int i ;

  for (i = 0 ; i < area->n_points ; i++)
  {
    area->lv[i].vno = -1 ;
  }
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------
  MRISlabelInvert(MRIS *surf, LABEL *label)
  Creates a label with all the vertices NOT in the input label.
  ---------------------------------------------------------------*/
LABEL *MRISlabelInvert(MRIS *surf, LABEL *label)
{
  MRI *tmpmri;
  LABEL *invlabel;
  int n,vtxno,ninv,nlabel;

  // Create binary mask of label
  tmpmri = MRISlabel2Mask(surf,label,NULL);

  // Count number of points. Need to do this in case duplicates
  nlabel = 0;
  for(vtxno=0; vtxno < surf->nvertices; vtxno++)
    if(MRIgetVoxVal(tmpmri,vtxno,0,0,0) > 0.5)
    {
      nlabel++;
    }

  // Alloc inverse label
  ninv = surf->nvertices-nlabel;
  invlabel = LabelAlloc(ninv,label->subject_name, NULL);
  invlabel->n_points = ninv;

  // Assign label points to vtxs NOT in the mask
  n = 0;
  for(vtxno=0; vtxno < surf->nvertices; vtxno++)
  {
    if(MRIgetVoxVal(tmpmri,vtxno,0,0,0) < 0.5)
    {
      invlabel->lv[n].vno = vtxno;
      invlabel->lv[n].x = surf->vertices[vtxno].x;
      invlabel->lv[n].y = surf->vertices[vtxno].y;
      invlabel->lv[n].z = surf->vertices[vtxno].z;
      invlabel->lv[n].stat = 0.0;
      n++;
    }
  }
  MRIfree(&tmpmri);
  return(invlabel);
}
int
LabelCopyStatsToSurface(LABEL *area, MRI_SURFACE *mris, int which)
{
  int    n, vno ;
  VERTEX *v ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    vno = area->lv[n].vno ;
    v = &mris->vertices[vno] ;
    v->marked = 1 ;
    switch (which)
    {
    case VERTEX_VALS:
      v->val = area->lv[n].stat ;
      break ;
    case VERTEX_CURV:
      v->curv = area->lv[n].stat ;
      break ;
    case VERTEX_STATS:
      v->stat = area->lv[n].stat ;
      break ;
    default:
      ErrorExit(ERROR_BADPARM,
                "LabelCopyStatsToSurface: unsupported which = %d", which);
    }
  }
  return(NO_ERROR) ;
}

int
LabelThreshold(LABEL *area, float thresh)
{
  int          n, ndel ;
  LABEL_VERTEX *lv ;

  for (ndel = n = 0 ; n < area->n_points ; n++)
  {
    lv = &area->lv[n] ;
    if (lv->stat < thresh)
    {
      lv->deleted = 1 ;
      ndel++ ;
    }
  }
  return(ndel) ;
}

int
LabelCropAnterior(LABEL *area, float anterior_dist)
{
  int    n ;
  float  amax ;

  amax = -100000;
  for (n = 0 ; n < area->n_points ; n++)
  {
    if (area->lv[n].y > amax)
    {
      amax = area->lv[n].y ;
    }
  }

  amax -= anterior_dist ;
  printf("cropping all vertices with Y > %2.0f\n", amax) ;
  for (n = 0 ; n < area->n_points ; n++)
  {
    if (area->lv[n].y < amax)
    {
      area->lv[n].deleted = 1 ;
    }
  }
  return(NO_ERROR) ;
}

int
LabelCropPosterior(LABEL *area, float anterior_dist)
{
  int    n ;
  float  amax ;

  amax = -100000;
  for (n = 0 ; n < area->n_points ; n++)
  {
    if (area->lv[n].y > amax)
    {
      amax = area->lv[n].y ;
    }
  }

  amax += anterior_dist ;
  printf("cropping all vertices with Y > %2.0f\n", amax) ;
  for (n = 0 ; n < area->n_points ; n++)
  {
    if (area->lv[n].y > amax)
    {
      area->lv[n].deleted = 1 ;
    }
  }
  return(NO_ERROR) ;
}

int
LabelMaskSurface(LABEL *area, MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *v ;

  LabelMarkSurface(area, mris) ;  /* mark all points in label */

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->marked || v->ripflag)
    {
      continue ;
    }
    v->stat = v->val = v->imag_val = v->val2 = v->valbak = v->val2bak = 0.0 ;
  }
  MRISclearMarks(mris) ;
  return(NO_ERROR) ;
}
int
LabelMaskSurfaceVolume(LABEL *area, MRI *mri, float nonmask_val)
{
  int    vno, n, x, y, z, f ;
  uchar  *marked ;
  long   nvox ;

  nvox = mri->width*mri->height*mri->depth*mri->nframes ;
  marked = (uchar *)calloc(nvox, sizeof(uchar)) ;
  for (n = 0 ; n < area->n_points ; n++)
  {
    vno = area->lv[n].vno ;
    if (vno >= 0)
      marked[vno] = 1 ;
  }
  for (n = x = 0 ; x < mri->width ; x++)
    for (y = 0 ; y < mri->height ; y++)
      for (z = 0 ; z < mri->depth ; z++)
	for (f = 0 ; f < mri->nframes ; f++, n++)
	{
	  if (marked[n])
	    continue ;
	  MRIsetVoxVal(mri, x, y, z, f, nonmask_val) ;
	}

  free(marked) ;
  return(NO_ERROR) ;
}
int
LabelCentroid(LABEL *area, MRI_SURFACE *mris, double *px, double *py, double *pz)
{
  int    vno, num, n  ;
  VERTEX *v ;
  double xc, yc, zc ;

  xc = yc = zc = 0.0 ;
  for (num = n = 0 ; n < area->n_points ; n++)
  {
    vno = area->lv[n].vno ;
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
    {
      DiagBreak() ;
    }
    num++ ;
    xc += v->x ;
    yc += v->y ;
    zc += v->z ;
  }
  *px = xc / num ;
  *py = yc / num ;
  *pz = zc / num ;
  return(NO_ERROR) ;
}

int
LabelFillVolume(MRI *mri, LABEL *area, int fillval)
{
  int           n ;
  double        xv, yv, zv ;
  LABEL_VERTEX  *lv ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    lv = &area->lv[n] ;
    MRIsurfaceRASToVoxel(mri, lv->x, lv->y, lv->z, &xv, &yv, &zv) ;
    MRIsetVoxVal(mri, nint(xv), nint(yv), nint(zv), 0, fillval) ;
  }
  return(NO_ERROR) ;
}

int
LabelSetVals(MRI_SURFACE *mris, LABEL *area, float fillval)
{
  int           n ;
  LABEL_VERTEX  *lv ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    lv = &area->lv[n] ;
    if (lv->deleted || lv->vno < 0)
      continue ;
    mris->vertices[lv->vno].val = fillval ;
  }
  return(NO_ERROR) ;
}
