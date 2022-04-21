/**
 * @brief utilities for manipulating ROIs.
 *
 * utilities (surface and volume) for manipulating arbitrary lists of
 * vertices/voxels.
 */
/*
 * Original Author: Bruce Fischl
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mri.h"
#include "mrisurf.h"
#include "mrishash_internals.h"

#include "diag.h"
#include "error.h"
#include "fio.h"
#include "macros.h"
#include "minc.h"
#include "proto.h"
#include "utils.h"
#include "surfcluster.h"

#include "label.h"

extern const char *Progname;

static int labelGetVoxelCoords(LABEL *area, LABEL_VERTEX *lv, float *px, float *py, float *pz);
static int labelGetSurfaceRasCoords(LABEL *area, LABEL_VERTEX *lv, float *px, float *py, float *pz);
static int update_vertex_indices(LABEL *area);
;
static LABEL_VERTEX *labelFindVertexNumber(LABEL *area, int vno);
static Transform *labelLoadTransform(const char *subject_name, const char *sdir, General_transform *transform);
#define MAX_VERTICES 500000
/*-----------------------------------------------------
------------------------------------------------------*/
LABEL *LabelReadFrom(const char *subject_name, FILE *fp)
{
  LABEL *area;
  char line[STRLEN], subjects_dir[STRLEN], *cp, *str;
  int vno, nlines;
  float x, y, z, stat;

  area = (LABEL *)calloc(1, sizeof(LABEL));
  if (!area) {
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate LABEL struct.", Progname);
  }
  cp = fgets(line, STRLEN, fp);  // read comment line
  if (cp == NULL) return (NULL);
  str = strstr(cp, "vox2ras=");
  if (str) {
    if (*(cp + strlen(cp) - 1) == '\n') *(cp + strlen(cp) - 1) = 0;
    sprintf(area->space, "%s", str + strlen("vox2ras="));
  }

  if (strstr(area->space, "voxel"))
    area->coords = LABEL_COORDS_VOXEL;
  else if (strstr(cp, "scanner"))
    area->coords = LABEL_COORDS_SCANNER_RAS;
  else
    area->coords = LABEL_COORDS_TKREG_RAS;

  cp = fgetl(line, STRLEN, fp);
  if (!cp) ErrorReturn(NULL, (ERROR_BADFILE, "%s: empty label", Progname));
  if (!sscanf(cp, "%d", &area->n_points)) {
    printf("\n%s\n", cp);
    ErrorReturn(NULL, (ERROR_BADFILE, "%s: could not scan # of lines from label file", Progname));
  }
  area->max_points = area->n_points;
  area->lv = (LABEL_VERTEX *)calloc(area->n_points, sizeof(LABEL_VERTEX));
  if (!area->lv)
    ErrorExit(
          ERROR_NOMEMORY, "%s: LabelReadFrom could not allocate %d-sized vector", Progname, sizeof(LV) * area->n_points);
  nlines = 0;
  while ((cp = fgetl(line, STRLEN, fp)) != NULL) {
    if (sscanf(cp, "%d %f %f %f %f", &vno, &x, &y, &z, &stat) != 5)
      ErrorReturn(NULL, (ERROR_BADFILE, "%s: could not parse %dth line '%s' in label file", Progname, nlines + 1, cp));
    area->lv[nlines].x = x;
    area->lv[nlines].y = y;
    area->lv[nlines].z = z;
    area->lv[nlines].stat = stat;
    area->lv[nlines].vno = vno;
    nlines++;
    if (nlines == area->n_points) break;
  }

  if (!nlines) ErrorReturn(NULL, (ERROR_BADFILE, "%s: no data in label file", Progname));
  if (subject_name) {
    cp = getenv("SUBJECTS_DIR");
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: no subject's directory specified in environment "
                "(SUBJECTS_DIR)",
                Progname);
    strncpy(subjects_dir, cp, STRLEN - 1);
    strncpy(area->subject_name, subject_name, STRLEN - 1);
    area->linear_transform = labelLoadTransform(subject_name, subjects_dir, &area->transform);
    area->inverse_linear_transform = get_inverse_linear_transform_ptr(&area->transform);
  }
  return (area);
}

LABEL *LabelRead(const char *subject_name, const char *label_name)
{
  LABEL *area;
  char *cp, subjects_dir[STRLEN], lname[STRLEN];
  char label_name0[STRLEN];
  FILE *fp;
  std::string fname;

  sprintf(label_name0, "%s", label_name);  // keep a copy

  if (subject_name && !strchr(label_name, '/')) {
    // subject name passed and label name does not have /
    // so assume it is in subject/label
    cp = getenv("SUBJECTS_DIR");
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: no subject's directory specified in environment "
                "(SUBJECTS_DIR)",
                Progname);
    strcpy(subjects_dir, cp);
    strcpy(lname, label_name);
    cp = strstr(lname, ".label");
    if (cp) *cp = 0;

    cp = strrchr(lname, '/');
    if (cp) {
      label_name = cp + 1;
    } else {
      label_name = lname;
    }

    fname = std::string(subjects_dir) + "/" + std::string(subject_name)
      + "/label/" + std::string(label_name) + ".label";
  }
  else {
    fname = label_name;
    cp = getenv("SUBJECTS_DIR");
    if (!cp) {
      strcpy(subjects_dir, ".");
    } else {
      strcpy(subjects_dir, cp);
    }
    strcpy(lname, label_name);
    cp = strstr(lname, ".label");
    if (cp == NULL) {
      fname = std::string(lname) + ".label";
    } else {
      fname = label_name;
    }
  }

  // As a last resort, treat label_name0 as a full path name
  if (!fio_FileExistsReadable(fname.c_str()) && fio_FileExistsReadable(label_name0)) {
    fname = label_name0;
  }

  //  printf("%s %s\n",label_name0,fname);

  /* read in the file */
  errno = 0;
  fp = fopen(fname.c_str(), "r");
  if (errno) perror(NULL);

  if (!fp) ErrorReturn(NULL, (ERROR_NOFILE, "%s: could not open label file %s", Progname, fname.c_str()));

  area = LabelReadFrom(subject_name, fp);
  if (area) {
    strncpy(area->name, fname.c_str(), STRLEN-1);
  }
  fclose(fp);
  return (area);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int LabelDump(FILE *fp, LABEL *area)
{
  int n;

  fprintf(fp, "label %s, from subject %s %s\n", area->name, area->subject_name, area->space);
  for (n = 0; n < area->n_points; n++)
    fprintf(fp, "%d  %2.3f  %2.3f  %2.3f\n", area->lv[n].vno, area->lv[n].x, area->lv[n].y, area->lv[n].z);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int LabelFree(LABEL **parea)
{
  LABEL *area;

  area = *parea;
  *parea = NULL;
  if (area->vertex_label_ind) free(area->vertex_label_ind);

  free(area->lv);
  free(area);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int LabelToCanonical(LABEL *area, MRI_SURFACE *mris)
{
  int n, vno;
  VERTEX *v;

  for (n = 0; n < area->n_points; n++) {
    vno = area->lv[n].vno;
    v = &mris->vertices[vno];
    area->lv[n].vno = -1; /* not associated with a vertex anymore */
    area->lv[n].x = v->cx;
    area->lv[n].y = v->cy;
    area->lv[n].z = v->cz;
  }
  strncpy(area->space, "TkReg coords=canonical", sizeof(area->space));
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int LabelToCurrent(LABEL *area, MRI_SURFACE *mris)
{
  int n, vno;
  VERTEX *v;

  for (n = 0; n < area->n_points; n++) {
    vno = area->lv[n].vno;
    v = &mris->vertices[vno];
    area->lv[n].vno = -1; /* not associated with a vertex anymore */
    area->lv[n].x = v->x;
    area->lv[n].y = v->y;
    area->lv[n].z = v->z;
  }
  strncpy(area->space, "TkReg coords=canonical", sizeof(area->space));
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int LabelFromCanonical(LABEL *area, MRI_SURFACE *mris)
{
  int n, vno, ui, vi, vlist[500000], nvertices;
  VERTEX *v;
  MRI_SP *mrisp;
  LV *lv;

  mrisp = MRISPalloc(0.0f, 1);

  for (n = 0; n < area->n_points; n++) {
    MRISPcoordinate(mrisp, area->lv[n].x, area->lv[n].y, area->lv[n].z, &ui, &vi);
    MRISPvox(mrisp, ui, vi) = 1;
  }

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    ImageWrite(mrisp->Ip, "label.tif");
  }

  /* now go through the surface and build a list of the vertices and
     keep track of those which fall into volume regions marked as within
     the area.
     */
  for (nvertices = vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    MRISPcoordinate(mrisp, v->cx, v->cy, v->cz, &ui, &vi);
    if (nvertices + 1 >= MAX_VERTICES) {
      break;
    }
    if (MRISPvox(mrisp, ui, vi) > 0.0f) {
      vlist[nvertices++] = vno;
    }
  }

  if (area->max_points < nvertices) {
    free(area->lv);
    area->max_points = nvertices;
    area->lv = (LABEL_VERTEX *)calloc(nvertices, sizeof(LABEL_VERTEX));
    if (!area->lv)
      ErrorExit(ERROR_NOMEMORY,
                "%s: LabelFromCanonical could not allocate %d-sized vector",
                Progname,
                sizeof(LV) * nvertices);
  }
  area->n_points = nvertices;
  for (n = 0; n < nvertices; n++) {
    lv = &area->lv[n];
    lv->vno = vlist[n];
    v = &mris->vertices[vlist[n]];
    lv->x = v->cx;
    lv->y = v->cy;
    lv->z = v->cz;
  }
  MRISPfree(&mrisp);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#if 0
#define VRES 2.0f
#define VSIZE 256
#define VDIM (VSIZE / VRES)
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
int LabelFromTalairach(LABEL *area, MRI_SURFACE *mris)
{
  int n, vno;
  VERTEX *v;

  for (n = 0; n < area->n_points; n++) {
    DiagHeartbeat((float)(n + 1) / (float)area->n_points);
    vno = MRIStalairachToVertex(mris, area->lv[n].x, area->lv[n].y, area->lv[n].z);
    v = &mris->vertices[vno];
    area->lv[n].vno = vno;
    area->lv[n].x = v->origx;
    area->lv[n].y = v->origy;
    area->lv[n].z = v->origz;
  }
  return (NO_ERROR);
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int LabelToFlat(LABEL *area, MRI_SURFACE *mris)
{
  int n, vno;
  MRIS_HASH_TABLE* mht = MHTcreateVertexTable(mris, CURRENT_VERTICES);

  for (n = 0; n < area->n_points; n++) {
    vno = area->lv[n].vno;
    if (vno >= 0 && vno < mris->nvertices) /* already have associated vertex */
    {
      VERTEX *v = &mris->vertices[vno];
      area->lv[n].x = v->x;
      area->lv[n].y = v->y;
      area->lv[n].z = v->z;
    }
    else /* in canonical coordinate system - find closest vertex */
    {
      vno = MHTfindVnoOfClosestVertexInTable(mht, mris, area->lv[n].x, area->lv[n].y, area->lv[n].z, 0);
      if (vno < 0) {
        continue;
      }
      VERTEX *v = &mris->vertices[vno];

      if (vno == Gdiag_no) {
        DiagBreak();
      }
      area->lv[n].vno = vno;
      area->lv[n].x = v->x;
      area->lv[n].y = v->y;
      area->lv[n].z = v->z;
    }
  }
  MHTfree(&mht);
  return (NO_ERROR);
}
int LabelWriteInto(LABEL *area, FILE *fp)
{
  int n, num, nbytes;

  for (num = n = 0; n < area->n_points; n++)
    if (!area->lv[n].deleted) {
      num++;
    }

  nbytes = fprintf(fp, "#!ascii label %s , from subject %s vox2ras=%s\n", area->name, area->subject_name, area->space);
  if (nbytes < 0) {
    printf("ERROR: writing to label file 1\n");
    fclose(fp);
    return (1);
  }

  nbytes = fprintf(fp, "%d\n", num);
  if (nbytes < 0) {
    printf("ERROR: writing to label file 2\n");
    fclose(fp);
    return (1);
  }
  for (n = 0; n < area->n_points; n++)
    if (!area->lv[n].deleted) {
      nbytes = fprintf(fp,
                       "%d  %2.3f  %2.3f  %2.3f %10.10f\n",
                       area->lv[n].vno,
                       area->lv[n].x,
                       area->lv[n].y,
                       area->lv[n].z,
                       area->lv[n].stat);
      if (nbytes < 0) {
        printf("ERROR: writing to label file 3\n");
        fclose(fp);
        return (1);
      }
    }
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int LabelWrite(LABEL *area, const char *label_name)
{
  char *cp, subjects_dir[STRLEN], lname[STRLEN];
  FILE *fp;
  int ret;
  std::string fname;

  strcpy(lname, label_name);
  cp = strrchr(lname, '.');
  if (cp && stricmp(cp, ".label") == 0) {
    *cp = 0;
  }
  label_name = lname;
  cp = strrchr(lname, '/');
  if ((cp == NULL) && strlen(area->subject_name) > 0) {
    cp = getenv("SUBJECTS_DIR");
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: no subject's directory specified in environment "
                "(SUBJECTS_DIR)",
                Progname);
    strcpy(subjects_dir, cp);
    fname = std::string(subjects_dir) + '/' + std::string(area->subject_name) +
      "/label/" + std::string(label_name) + ".label";
  }
  else {
    cp = strrchr(lname, '.');
    if (cp && stricmp(cp, ".label") == 0) {
      fname = label_name;
    }
    else {
      fname = std::string(label_name) + ".label";
    }
  }

  fp = fopen(fname.c_str(), "w");
  if (!fp) ErrorReturn(ERROR_NOFILE, (ERROR_NO_FILE, "%s: could not open label file %s", Progname, fname.c_str()));

  ret = LabelWriteInto(area, fp);
  fclose(fp);
  return (ret);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int LabelRipRestOfSurface(LABEL *area, MRI_SURFACE *mris)
{
  int vno, n;
  VERTEX *v;

  LabelToFlat(area, mris);
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    v->ripflag = 1;
  }

  for (n = 0; n < area->n_points; n++) {
    vno = area->lv[n].vno;
    if (vno < 0 || vno >= mris->nvertices) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    v = &mris->vertices[vno];
    v->ripflag = 0;
  }
  MRISsetRipInFacesWithRippedVertices(mris);
  MRISremoveRipped(mris);
  return (NO_ERROR);
}

/*!
  \fn int LabelRip(MRI_SURFACE *mris, const LABEL *area, const int Outside)
  Use label to rip vertices and faces. Outside should be 0 or 1.
  If Outside=1, then vertices outside the label are ripped.
  If Outside=0, then vertices inside the label are ripped.
 */
int LabelRip(MRI_SURFACE *mris, const LABEL *area, const int Outside)
{
  int vno, n;
  VERTEX *v;

  // Make sure ripflag is set to a known value for all vertices
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    v->ripflag = Outside;
  }
  // Rip or unrip vertices in the label
  for (n = 0; n < area->n_points; n++) {
    vno = area->lv[n].vno;
    if(vno < 0 || vno >= mris->nvertices) continue;
    v = &mris->vertices[vno];
    v->ripflag = !Outside;
  }
  // Now ripped the faces
  MRISsetRipInFacesWithRippedVertices(mris);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int LabelRipRestOfSurfaceWithThreshold(LABEL *area, MRI_SURFACE *mris, float thresh)
{
  int vno, n;
  VERTEX *v;

  LabelToFlat(area, mris);
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    v->ripflag = 1;
  }

  for (n = 0; n < area->n_points; n++) {
    if (area->lv[n].stat < thresh) {
      continue;
    }
    vno = area->lv[n].vno;
    if (vno < 0 || vno >= mris->nvertices) {
      continue;
    }
    v = &mris->vertices[vno];
    v->ripflag = 0;
  }
  MRISsetRipInFacesWithRippedVertices(mris);
  MRISremoveRipped(mris);
  return (NO_ERROR);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int LabelRemoveOverlap(LABEL *area1, LABEL *area2)
{
  int n1, n2, vno;

  for (n1 = 0; n1 < area1->n_points; n1++) {
    vno = area1->lv[n1].vno;
    for (n2 = 0; n2 < area2->n_points; n2++) {
      if (vno == area2->lv[n2].vno) {
        area1->lv[n1].deleted = 1;
        break;
      }
    }
  }
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters: two labels to be intersected

        Returns value: error code and intersection in area1

        Description: elements are not removed, only the deleted
                     flag is set
------------------------------------------------------*/
int LabelIntersect(LABEL *area1, LABEL *area2)
{
  int n, vno;
  int vmin, vmax, vnum;
  int *isec;
  if (area1->n_points == 0) {
    return (NO_ERROR);
  }

  vmin = area1->lv[0].vno;
  vmax = area1->lv[0].vno;
  for (n = 0; n < area1->n_points; n++) {
    vno = area1->lv[n].vno;
    if (vno < vmin) {
      vmin = vno;
    }
    if (vno > vmax) {
      vmax = vno;
    }
  }
  for (n = 0; n < area2->n_points; n++) {
    vno = area2->lv[n].vno;
    if (vno < vmin) {
      vmin = vno;
    }
    if (vno > vmax) {
      vmax = vno;
    }
  }
  vnum = vmax - vmin;
  isec = (int *)calloc(vnum, sizeof(int));
  if (!isec) ErrorExit(ERROR_NOMEMORY, "%s: could not allocate LabelIntersect struct.", Progname);
  for (n = 0; n < vnum; n++) {
    isec[n] = -1;
  }

  for (n = 0; n < area1->n_points; n++) {
    vno = area1->lv[n].vno - vmin;
    isec[vno] = n;
    area1->lv[n].deleted = 1;
  }
  for (n = 0; n < area2->n_points; n++) {
    vno = area2->lv[n].vno - vmin;
    if (isec[vno] > -1) {
      area1->lv[isec[vno]].deleted = 0;
    }
  }

  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
LABEL *LabelAlloc(int max_points, const char *subject_name, const char *label_name)
{
  LABEL *area;
  char *cp, subjects_dir[STRLEN];
  int  n ;

  area = (LABEL *)calloc(1, sizeof(LABEL));
  if (!area) {
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate LABEL struct.", Progname);
  }
  if (subject_name) {
    cp = getenv("SUBJECTS_DIR");
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: no subject's directory specified in environment "
                "(SUBJECTS_DIR)",
                Progname);
    strcpy(subjects_dir, cp);
    strcpy(area->subject_name, subject_name);
  }
  else {
    if (label_name) {
      strncpy(area->name, label_name, STRLEN - 1);
    }
  }

  area->n_points = 0;
  area->max_points = max_points;
  area->lv = (LABEL_VERTEX *)calloc(area->max_points, sizeof(LABEL_VERTEX));
  if (!area->lv)
    ErrorExit(ERROR_NOMEMORY,
              "%s: LabelAlloc(%s) could not allocate %d-sized vector",
              Progname,
              label_name ? label_name : "",
              sizeof(LV) * area->n_points);
  for (n = 0 ; n < area->max_points ; n++)
    area->lv[n].vno = -1 ;   // mark them as unassigned (0 is a valid vertex index)
  strcpy(area->space, "TkReg");
  area->coords = LABEL_COORDS_TKREG_RAS;  // would like to use scanner RAS, but need an volume or header
  return (area);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int LabelCurvFill(LABEL *area, int *vertex_list, int nvertices, int max_vertices, MRI_SURFACE *mris)
{
  if (!max_vertices) {
    max_vertices = area->max_points;
  }

  MRISclearMarks(mris);
  
  float max_curv = 0, min_curv = 0;
  int n;
  for (n = 0; n < nvertices; n++) {
    VERTEX * const v = &mris->vertices[vertex_list[n]];
    if (v->ripflag) {
      continue;
    }
    v->marked = 1;
    LV* const lv = &area->lv[n];
    lv->vno = vertex_list[n];
    lv->x = v->x;
    lv->y = v->y;
    lv->z = v->z;
    area->n_points++;
    if (v->curv > max_curv) {
      max_curv = v->curv;
    }
    if (v->curv < min_curv) {
      min_curv = v->curv;
    }
  }
  
  float curv_thresh;
  if (-min_curv > max_curv) {
    curv_thresh = min_curv;
  }
  else {
    curv_thresh = max_curv;
  }
  curv_thresh *= 0.1; /* 10% of max */
  
  fprintf(stderr, "starting fill with curvature threshold = %2.3f\n", curv_thresh);
  int nfilled;
  do {
    nfilled = 0;

    for (n = 0; n < area->n_points && area->n_points + nfilled < area->max_points; n++) {
      LV const * const lv = &area->lv[n];
      
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[lv->vno];
      VERTEX                * const v  = &mris->vertices[lv->vno];
      if (v->ripflag) {
        continue;
      }
      
      int nv;
      for (nv = 0; nv < vt->vnum; nv++) /* go through neighbors */
      {
        VERTEX * const vn = &mris->vertices[vt->v[nv]];
        if (vn->ripflag || vn->marked) {
          continue;
        }
        if (((curv_thresh > 0) && (vn->curv > curv_thresh)) || ((curv_thresh < 0) && (vn->curv < curv_thresh))) {
          vn->marked = 1;
          LV* const lvn = &area->lv[area->n_points + nfilled];
          nfilled++;
          lvn->vno = vt->v[nv];
          lvn->x = vn->x;
          lvn->y = vn->y;
          lvn->z = vn->z;
        }
        if (area->n_points + nfilled >= area->max_points) {
          break;
        }
      }
      if (area->n_points + nfilled >= area->max_points) {
        break;
      }
    }
    fprintf(stderr, "%d vertices added.\n", nfilled);
    area->n_points += nfilled;
  } while ((nfilled > 0) && (area->n_points < max_vertices));
  MRISclearMarks(mris);
  fprintf(stderr, "%d vertices in label %s\n", area->n_points, area->name);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int LabelFillMarked(LABEL *area, MRI_SURFACE *mris)
{
  int nfilled;
  do {
    nfilled = 0;
    int n;
    for (n = 0; n < area->n_points && area->n_points + nfilled < area->max_points; n++) {
      LV const * const lv = &area->lv[n];
      
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[lv->vno];
      VERTEX                * const v  = &mris->vertices         [lv->vno];
      v->marked = 2;
      if (v->ripflag) {
        continue;
      }
      int nv;
      for (nv = 0; nv < vt->vnum; nv++) /* go through neighbors */
      {
        VERTEX * const vn = &mris->vertices[vt->v[nv]];
        if (vn->ripflag) {
          continue;
        }
        if (vn->marked == 1) /* add it to the label */
        {
          vn->marked = 2;
          LV * const lvn = &area->lv[area->n_points + nfilled];
          nfilled++;
          lvn->vno = vt->v[nv];
          lvn->x = vn->x;
          lvn->y = vn->y;
          lvn->z = vn->z;
        }
        if (area->n_points + nfilled >= area->max_points) {
          break;
        }
      }
      if (area->n_points + nfilled >= area->max_points) {
        break;
      }
    }
    /*    fprintf(stderr, "%d vertices added.\n", nfilled) ;*/
    area->n_points += nfilled;
  } while ((nfilled > 0) && (area->n_points < area->max_points));
#if 0
  fprintf(stderr, "%d vertices in label %s\n",area->n_points, area->name);
#endif
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int LabelFillAnnotated(LABEL *area, MRI_SURFACE *mris)
{
  int const annotation = mris->vertices[area->lv[0].vno].annotation;
  int nfilled;
  do {
    nfilled = 0;
    int n;
    for (n = 0; n < area->n_points && area->n_points + nfilled < area->max_points; n++) {
      LV * const lv = &area->lv[n];
      VERTEX_TOPOLOGY * const vt = &mris->vertices_topology[lv->vno];
      VERTEX          * const v  = &mris->vertices         [lv->vno];
      v->marked = 1;
      if (v->ripflag) {
        continue;
      }
      int nv;
      for (nv = 0; nv < vt->vnum; nv++) /* go through neighbors */
      {
        VERTEX * const vn = &mris->vertices[vt->v[nv]];
        if (vn->ripflag || vn->marked) {
          continue;
        }
        if (vn->annotation == annotation) /* add it to the label */
        {
          vn->marked = 1;
          LV * const lvn = &area->lv[area->n_points + nfilled];
          nfilled++;
          lvn->vno = vt->v[nv];
          lvn->x = vn->x;
          lvn->y = vn->y;
          lvn->z = vn->z;
        }
        if (area->n_points + nfilled >= area->max_points) {
          break;
        }
      }
      if (area->n_points + nfilled >= area->max_points) {
        break;
      }
    }
    /*    fprintf(stderr, "%d vertices added.\n", nfilled) ;*/
    area->n_points += nfilled;
  } while ((nfilled > 0) && (area->n_points < area->max_points));
#if 0
  fprintf(stderr, "%d vertices in label %s\n",area->n_points, area->name);
#endif
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int LabelFillAll(LABEL *area, int *vertex_list, int nvertices, int max_vertices, MRI_SURFACE *mris)
{
  int n, nfilled, nv;
  LV *lv, *lvn;

  if (!max_vertices) {
    max_vertices = area->max_points;
  }

  for (n = 0; n < nvertices; n++) {
    VERTEX * const v = &mris->vertices[vertex_list[n]];
    if (v->ripflag) {
      continue;
    }
    v->marked = 1;
    lv = &area->lv[n];
    lv->vno = vertex_list[n];
    lv->x = v->x;
    lv->y = v->y;
    lv->z = v->z;
    area->n_points++;
  }
  MRISclearMarks(mris);
  do {
    nfilled = 0;

    for (n = 0; n < area->n_points && area->n_points + nfilled < area->max_points; n++) {
      lv = &area->lv[n];
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[lv->vno];
      VERTEX          const * const v  = &mris->vertices         [lv->vno];
      if (v->ripflag) {
        continue;
      }
      for (nv = 0; nv < vt->vnum; nv++) /* go through neighbors */
      {
        VERTEX * const vn = &mris->vertices[vt->v[nv]];
        if (vn->ripflag || vn->marked) {
          continue;
        }
        vn->marked = 1;
        lvn = &area->lv[area->n_points + nfilled];
        nfilled++;
        lvn->vno = vt->v[nv];
        lvn->x = vn->x;
        lvn->y = vn->y;
        lvn->z = vn->z;
        if (area->n_points + nfilled >= area->max_points) {
          break;
        }
      }
      if (area->n_points + nfilled >= area->max_points) {
        break;
      }
    }
    fprintf(stderr, "%d vertices added.\n", nfilled);
    area->n_points += nfilled;
  } while ((nfilled > 0) && (area->n_points < max_vertices));
  MRISclearMarks(mris);
  fprintf(stderr, "%d vertices in label %s\n", area->n_points, area->name);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static Transform *labelLoadTransform(const char *subject_name, const char *sdir, General_transform *transform)
{
  std::string xform_fname;

  xform_fname = std::string(sdir) + '/' + std::string(subject_name) + "/mri/transforms/talairach.xfm";
  if (FileExists(xform_fname.c_str()) == 0) return (NULL);
  if (input_transform_file(xform_fname.c_str(), transform) != OK) {
    ErrorReturn(NULL, (ERROR_NOFILE, "%s: could not load transform file '%s'", Progname, xform_fname.c_str()));
  }

  return (get_linear_transform_ptr(transform));
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int LabelTalairachTransform(LABEL *area, MRI_SURFACE *mris)
{
  int n;
  LV *lv;
  double x, y, z, xt, yt, zt;
  VERTEX *v;

  for (n = 0; n < area->n_points; n++) {
    lv = &area->lv[n];
    v = &mris->vertices[lv->vno];
    x = v->origx;
    y = v->origy;
    z = v->origz;
    transform_point(area->linear_transform, x, y, z, &xt, &yt, &zt);
    lv->x = xt;
    lv->y = yt;
    lv->z = zt;
  }
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int LabelSphericalTransform(LABEL *area, MRI_SURFACE *mris)
{
  int n;
  LV *lv;
  VERTEX *v;

  for (n = 0; n < area->n_points; n++) {
    lv = &area->lv[n];
    v = &mris->vertices[lv->vno];
    lv->x = v->cx;
    lv->y = v->cy;
    lv->z = v->cz;
  }
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MATRIX *LabelCovarianceMatrix(LABEL *area, MATRIX *mCov)
{
  MATRIX *mInputs;
  int i;
  LV *lv;

  mInputs = MatrixAlloc(area->n_points, 3, MATRIX_REAL);
  for (i = 0; i < area->n_points; i++) {
    lv = &area->lv[i];
    *MATRIX_RELT(mInputs, i + 1, 1) = lv->x;
    *MATRIX_RELT(mInputs, i + 1, 2) = lv->y;
    *MATRIX_RELT(mInputs, i + 1, 3) = lv->z;
  }

  mCov = MatrixCovariance(mInputs, mCov, NULL);

  return (mCov);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
LABEL *LabelCombine(LABEL *asrc, LABEL *adst)
{
  LABEL *atmp;
  int n;
  LV *vsrc, *vdst;

  if (!adst) {
    adst = LabelClone(asrc);
  }
  if (adst->max_points < asrc->n_points + adst->n_points) /* won't fit - expand */
  {
    atmp = LabelAlloc(2 * (asrc->n_points + adst->n_points), asrc->subject_name, asrc->name);
    LabelCopy(adst, atmp);
    LabelFree(&adst);
    adst = atmp;
  }

  for (n = 0; n < asrc->n_points; n++) {
    vsrc = &asrc->lv[n];
    vdst = &adst->lv[n + adst->n_points];
    memmove(vdst, vsrc, sizeof(LV));
  }
  adst->n_points += asrc->n_points;
  return (adst);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
LABEL *LabelCopy(LABEL *asrc, LABEL *adst)
{
  if (!adst) {
    adst = LabelAlloc(asrc->max_points, asrc->subject_name, asrc->name);
  }

  adst->mris = asrc->mris;
  adst->mri_template = asrc->mri_template;
  adst->mht = asrc->mht;
  adst->coords = asrc->coords;
  strcpy(adst->space, asrc->space);
  if (asrc->vertex_label_ind) {
    adst->vertex_label_ind = (int *)calloc(asrc->max_points, sizeof(int));
    if (adst->vertex_label_ind == NULL)
      ErrorExit(ERROR_NOMEMORY, "LabelCompact: could not allocate %d long label index array", asrc->max_points);
    memmove(adst->vertex_label_ind, asrc->vertex_label_ind, asrc->max_points * sizeof(int));
  }

  adst->n_points = asrc->n_points;
  strcpy(adst->name, asrc->name);
  strcpy(adst->subject_name, asrc->subject_name);

  memmove(adst->lv, asrc->lv, asrc->n_points * sizeof(LABEL_VERTEX));
  return (adst);
}
/*-----------------------------------------------------
  int LabelRemoveDuplicates(LABEL *area)
  Sets the 'deleted' flag of a label point if it is
  a duplicate. Does not actually remove duplicats!
  ------------------------------------------------------*/
int LabelRemoveDuplicates(LABEL *area)
{
  int n1, n2, deleted = 0;
  LV *lv1, *lv2;

  // loop thru each label point
  for (n1 = 0; n1 < area->n_points; n1++) {
    lv1 = &area->lv[n1];
    if (lv1->deleted) {
      continue;
    }
    // loop thru the remaining looking for duplicates
    for (n2 = n1 + 1; n2 < area->n_points; n2++) {
      lv2 = &area->lv[n2];
      if (lv1->vno >= 0 && lv2->vno >= 0 && lv1->vno == lv2->vno) {
        deleted++;
        lv2->deleted = 1;
      }
      else if (lv1->vno < 0 && lv2->vno < 0) {
        if (FEQUAL(lv1->x, lv2->x) && FEQUAL(lv1->y, lv2->y) && FEQUAL(lv1->z, lv2->z)) {
          deleted++;
          lv2->deleted = 1;
        }
      }
    }
  }

  if (Gdiag & DIAG_SHOW) fprintf(stderr, "%d duplicate vertices removed from label %s.\n", deleted, area->name);
  return (NO_ERROR);
}

/*-----------------------------------------------------
  int LabelRemoveDuplicates(LABEL *area)
  Sets the 'deleted' flag of a label point if it is
  a duplicate. Does not actually remove duplicats!
  ------------------------------------------------------*/
LABEL *LabelRemoveAlmostDuplicates(LABEL *area, double dist, LABEL *ldst)
{
  int n1, n2, deleted = 0;
  LV *lv1, *lv2;

  // loop thru each label point
  for (n1 = 0; n1 < area->n_points; n1++) {
    lv1 = &area->lv[n1];
    if (lv1->deleted) {
      continue;
    }
    // loop thru the remaining looking for duplicates
    for (n2 = n1 + 1; n2 < area->n_points; n2++) {
      lv2 = &area->lv[n2];
      if (lv1->vno >= 0 && lv2->vno >= 0 && lv1->vno == lv2->vno) {
        deleted++;
        lv2->deleted = 1;
      }
      else if (lv1->vno < 0 && lv2->vno < 0) {
        if (fabs(lv1->x - lv2->x) < dist && fabs(lv1->y - lv2->y) < dist && fabs(lv1->z - lv2->z) < dist) {
          deleted++;
          lv2->deleted = 1;
        }
      }
    }
  }
  if (Gdiag & DIAG_SHOW) fprintf(stderr, "%d duplicate vertices removed from label %s.\n", deleted, area->name);
  ldst = LabelCompact(area, ldst);
  return (ldst);
}
LABEL *LabelCompact(LABEL *lsrc, LABEL *ldst)
{
  int i, n;
  for (i = n = 0; i < lsrc->n_points; i++)
    if (lsrc->lv[i].deleted == 0) n++;

  ldst = LabelRealloc(ldst, n);
  if (ldst != lsrc && lsrc->vertex_label_ind) {
    ldst->vertex_label_ind = (int *)calloc(lsrc->max_points, sizeof(int));
    if (ldst->vertex_label_ind == NULL)
      ErrorExit(ERROR_NOMEMORY, "LabelCompact: could not allocate %d long label index array", lsrc->max_points);
  }
  if (ldst->vertex_label_ind != NULL)
    for (i = 0; i < ldst->max_points; i++) ldst->vertex_label_ind[i] = -1;

  for (i = n = 0; i < lsrc->n_points; i++)
    if (lsrc->lv[i].deleted == 0) {
      ldst->lv[n].deleted = 0;
      ldst->lv[n].x = lsrc->lv[i].x;
      ldst->lv[n].y = lsrc->lv[i].y;
      ldst->lv[n].z = lsrc->lv[i].z;
      ldst->lv[n].vno = lsrc->lv[i].vno;
      ldst->lv[n].stat = lsrc->lv[i].stat;
      if (lsrc->lv[i].vno >= 0 && lsrc->lv[i].vno < ldst->max_points) ldst->vertex_label_ind[lsrc->lv[i].vno] = n;
      n++;
    }
  ldst->n_points = n;
  return (ldst);
}

/*-----------------------------------------------------*/
double LabelArea(LABEL *area, MRI_SURFACE *mris)
{
  int i;
  double total_area;
  LV *lv;
  VERTEX *v;

  for (total_area = 0.0, i = 0; i < area->n_points; i++) {
    lv = &area->lv[i];
    v = &mris->vertices[lv->vno];
    if (v->ripflag || lv->deleted) continue;

    total_area += v->area;
  }

  return (total_area);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int LabelMean(LABEL *area, double *px, double *py, double *pz)
{
  int i, n;
  double x, y, z;
  LV *lv;

  for (x = y = z = 0.0, n = i = 0; i < area->n_points; i++) {
    lv = &area->lv[i];

    x += lv->x;
    y += lv->y;
    z += lv->z;
    n++;
  }

  *px = x / (double)n;
  *py = y / (double)n;
  *pz = z / (double)n;
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
double LabelMeanIntensity(LABEL *area, MRI *mri)
{
  int i;
  double x, y, z, mean, val;
  LV *lv;
  LABEL *area2;

  area2 = LabelToVoxel(area, mri, NULL);

  for (mean = x = y = z = 0.0, i = 0; i < area->n_points; i++) {
    lv = &area2->lv[i];

    x = lv->x;
    y = lv->y;
    z = lv->z;
    MRIsampleVolume(mri, x, y, z, &val);
    mean += val;
  }

  if (i > 0) mean /= (double)i;
  LabelFree(&area2);
  return (mean);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
double LabelVariance(LABEL *area, double ux, double uy, double uz)
{
  int i, n;
  double xd, yd, zd, dsq;
  LV *lv;

  for (dsq = 0.0, n = i = 0; i < area->n_points; i++) {
    lv = &area->lv[i];

    xd = lv->x - ux;
    yd = lv->y - uy;
    zd = lv->z - uz;
    dsq += xd * xd + yd * yd + zd * zd;
    n++;
  }

  dsq /= (double)n;
  return (dsq);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int LabelMark(LABEL *area, MRI_SURFACE *mris)
{
  int n, vno;
  VERTEX *v;

  for (n = 0; n < area->n_points; n++) {
    vno = area->lv[n].vno;
    if (vno < 0 || vno >= mris->nvertices) {
      DiagBreak();
    }
    if (area->lv[n].deleted > 0) {
      continue;
    }
    v = &mris->vertices[vno];
    v->marked = 1;
  }
  return (NO_ERROR);
}
int LabelMark2(LABEL *area, MRI_SURFACE *mris)
{
  int n, vno;
  VERTEX *v;

  for (n = 0; n < area->n_points; n++) {
    vno = area->lv[n].vno;
    if (vno < 0 || vno >= mris->nvertices) {
      DiagBreak();
    }
    if (area->lv[n].deleted > 0) {
      continue;
    }
    v = &mris->vertices[vno];
    v->marked2 = 1;
  }
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int LabelAddToMark(LABEL *area, MRI_SURFACE *mris, int val_to_add)
{
  int n, vno;
  VERTEX *v;

  for (n = 0; n < area->n_points; n++) {
    vno = area->lv[n].vno;
    if (vno < 0 || vno >= mris->nvertices) {
      DiagBreak();
    }
    if (area->lv[n].deleted > 0) continue;

    v = &mris->vertices[vno];
    v->marked += val_to_add;
  }
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int LabelMarkWithThreshold(LABEL *area, MRI_SURFACE *mris, float thresh)
{
  int n, vno;
  VERTEX *v;

  for (n = 0; n < area->n_points; n++) {
    if (area->lv[n].stat < thresh) {
      continue;
    }
    vno = area->lv[n].vno;
    v = &mris->vertices[vno];
    v->marked = 1;
  }
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters: A label and a surface

        Returns value: An error if either parameter is NULL

        Description

        For every point in the label that is not deleted (doesn't have
        it's deleted flag on), the corresponding vertex in the surface
        has its marked flag turned on.

------------------------------------------------------*/
int LabelMarkUndeleted(LABEL *area, MRI_SURFACE *mris)
{
  int n, vno;
  VERTEX *v;

  if (NULL == area || NULL == mris)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "LabelMarkUndeleted: area or mris was NULL"));

  /* For every point, if it's not deleted, get the vno. Check the
     vno. If it's good, set the vertex's marked flag.*/
  for (n = 0; n < area->n_points; n++) {
    if (!area->lv[n].deleted) {
      vno = area->lv[n].vno;
      if (vno < 0 || vno >= mris->nvertices) {
        printf("LabelMarkUndeleted: label point %d had invalid vno %d\n", n, vno);
        continue;
      }
      v = &mris->vertices[vno];
      v->marked = 1;
    }
  }
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int LabelUnmark(LABEL *area, MRI_SURFACE *mris)
{
  int n, vno;
  VERTEX *v;

  for (n = 0; n < area->n_points; n++) {
    vno = area->lv[n].vno;
    if (vno < 0) {
      continue;
    }
    v = &mris->vertices[vno];
    v->marked = 0;
  }
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int LabelToOriginal(LABEL *area, MRI_SURFACE *mris)
{
  int n, vno;
  VERTEX *v;

  for (n = 0; n < area->n_points; n++) {
    vno = area->lv[n].vno;
    if (vno < 0 || vno >= mris->nvertices) {
      continue;
    }
    v = &mris->vertices[vno];
    area->lv[n].x = v->origx;
    area->lv[n].y = v->origy;
    area->lv[n].z = v->origz;
  }
  strncpy(area->space, "TkReg coords=orig", sizeof(area->space));
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters: allocated LABEL, valid MRI_SURFACE

        Returns value: error code

        Description: For every point in the label, gets the
                     corresponding white coords from the surface and
                     writes them in the label.
 ------------------------------------------------------*/
int LabelToWhite(LABEL *area, MRI_SURFACE *mris)
{
  int n, vno;
  VERTEX *v;

  for (n = 0; n < area->n_points; n++) {
    vno = area->lv[n].vno;
    if (vno < 0 || vno >= mris->nvertices) {
      continue;
    }
    v = &mris->vertices[vno];
    area->lv[n].x = v->whitex;
    area->lv[n].y = v->whitey;
    area->lv[n].z = v->whitez;
  }
  strncpy(area->space, "TkReg coords=white", sizeof(area->space));
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
LABEL *LabelFromMarkedSurface(MRI_SURFACE *mris)
{
  int vno, npoints, n;
  LABEL *area;
  VERTEX *v;

  for (npoints = vno = 0; vno < mris->nvertices; vno++)
    if (mris->vertices[vno].marked) {
      npoints++;
    }

  if (!npoints) {
    return (NULL);
  }
  area = LabelAlloc(npoints, NULL, NULL);
  for (n = vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (!v->marked) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    area->lv[n].x = v->x;
    area->lv[n].y = v->y;
    area->lv[n].z = v->z;
    area->lv[n].vno = vno;
    area->lv[n].stat = v->stat;
    if (FZERO(v->stat)) {
      DiagBreak();
    }
    n++;
  }
  area->n_points = npoints;
  return (area);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
LABEL *LabelFromMarkValue(MRI_SURFACE *mris, int mark)
{
  int vno, npoints, n;
  LABEL *area;
  VERTEX *v;

  for (npoints = vno = 0; vno < mris->nvertices; vno++)
    if (mris->vertices[vno].marked == mark) {
      npoints++;
    }

  if (!npoints) {
    return (NULL);
  }
  area = LabelAlloc(npoints, NULL, NULL);
  for (n = vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->marked != mark) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    area->lv[n].x = v->x;
    area->lv[n].y = v->y;
    area->lv[n].z = v->z;
    area->lv[n].vno = vno;
    area->lv[n].stat = v->stat;
    if (FZERO(v->stat)) {
      DiagBreak();
    }
    n++;
  }
  area->n_points = npoints;
  return (area);
}
/*-----------------------------------------------------
------------------------------------------------------*/
int LabelAddToSurfaceMark(LABEL *area, MRI_SURFACE *mris, int mark_to_add)
{
  int n, vno;
  VERTEX *v;

  for (n = 0; n < area->n_points; n++) {
    vno = area->lv[n].vno;
    if (vno < 0) continue;

    if (vno >= mris->nvertices) {
      printf("ERROR: LabelMarkSurface: label point %d exceeds nvertices %d\n", vno, mris->nvertices);
      return (1);
    }
    v = &mris->vertices[vno];
    if (v->ripflag) continue;

    v->marked += mark_to_add;
  }
  return (NO_ERROR);
}
/*-----------------------------------------------------
------------------------------------------------------*/
int LabelMarkSurface(LABEL *area, MRI_SURFACE *mris)
{
  int n, vno;
  VERTEX *v;

  for (n = 0; n < area->n_points; n++) {
    vno = area->lv[n].vno;
    if (vno < 0 || area->lv[n].deleted) {
      continue;
    }
    if (vno >= mris->nvertices) {
      printf("ERROR: LabelMarkSurface: label point %d exceeds nvertices %d\n", vno, mris->nvertices);
      return (1);
    }
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->marked = 1;
  }
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int LabelIsCompletelyUnassigned(LABEL *area, int *unassigned)
{
  int i;

  if (NULL == area) {
    return ERROR_BADPARM;
  }
  if (NULL == unassigned) {
    return ERROR_BADPARM;
  }

  /* Go through the label. If a point vno is _not_ -1, return false,
     because we have an assigned point. If we get through the entire
     label without finding a vno that isn't -1, return true. */
  for (i = 0; i < area->n_points; i++) {
    if (area->lv[i].vno != -1) {
      *unassigned = 0;
      return (NO_ERROR);
    }
  }

  *unassigned = 1;
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#include "mrishash.h"
int LabelFillUnassignedVertices(MRI_SURFACE *mris, LABEL *area, int coords)
{
  int n, i, vno, min_vno, nfilled = 0, max_vno;
  LV *lv;
  MHT *mht;
  VERTEX *v;
  float dx, dy, dz, x, y, z, dist, min_dist;
  int num_not_found;
  float vx, vy, vz;
  double max_spacing;
  vx = vy = vz = x = y = z = -1;

  MRIScomputeVertexSpacingStats(mris, NULL, NULL, &max_spacing, NULL, &max_vno, coords);

  for (i = n = 0; n < area->n_points; n++) {
    lv = &area->lv[n];
    if (lv->vno >= 0 && lv->vno < mris->nvertices) {
      continue;
    }
    i++; /* count # of unassigned vertices */
  }
  if (i <= 0) {
    return (0); /* no work needed */
  }

  fprintf(stderr, "%d unassigned vertices in label - building spatial LUT...\n", i);

  /* if we can't find a vertex within 10 mm of the point, something is wrong */
  mht = MHTcreateVertexTable_Resolution(mris, coords, 2 * max_spacing);
  fprintf(stderr, "assigning vertex numbers to label...\n");
  num_not_found = 0;
  if (area->mri_template == NULL)
  {
    area->mri_template = MRIalloc(mris->vg.width, mris->vg.height, mris->vg.depth, MRI_UCHAR) ;
    useVolGeomToMRI(&mris->vg, area->mri_template);
  }

  if (area->mris == NULL)
    area->mris = mris ;
  for (n = 0; n < area->n_points; n++) {
    lv = &area->lv[n];
    if (lv->vno >= 0 && lv->vno <= mris->nvertices) {
      continue;
    }
    nfilled++;
    labelGetSurfaceRasCoords(area, lv, &x, &y, &z);

    {
      min_dist = 10000.0;
      min_vno = -1;

      MHBT *bucket = MHTacqBucket(mht, x, y, z);
      if (bucket) {
        MHB *bin;
        for (bin = bucket->bins, i = 0; i < bucket->nused; i++, bin++) {
          vno = bin->fno;
          v = &mris->vertices[vno];
          if (vno == Gdiag_no) {
            DiagBreak();
          }

          MRISgetCoords(v, coords, &vx, &vy, &vz);

          dx = vx - x;
          dy = vy - y;
          dz = vz - z;

          dist = sqrt(dx * dx + dy * dy + dz * dz);
          if (dist < min_dist) {
            min_dist = dist;
            min_vno = vno;
          }
        }
      }
      MHTrelBucket(&bucket);
    }

    if (min_vno == -1) {
      switch (coords) {
      case ORIG_VERTICES:
        min_vno = MRISfindClosestOriginalVertex(mris, lv->x, lv->y, lv->z);
        break;
      case WHITE_VERTICES:
        min_vno = MRISfindClosestWhiteVertex(mris, lv->x, lv->y, lv->z);
        break;
      case CURRENT_VERTICES:
        min_vno = MRISfindClosestVertex(mris, lv->x, lv->y, lv->z, NULL, CURRENT_VERTICES);
        break;
      case CANONICAL_VERTICES:
        min_vno = MRISfindClosestCanonicalVertex(mris, lv->x, lv->y, lv->z);
        break;
      default:
        break;
      }
    }
    if (min_vno == -1) {
      num_not_found++;
    }
    lv->vno = min_vno;
  }
  LabelRemoveDuplicates(area);
  MHTfree(&mht);

  if (num_not_found > 0) {
    fprintf(stderr, "Couldn't assign %d vertices.\n", num_not_found);
  }

  return (nfilled);
}

LABEL *LabelSphericalCombine(MRI_SURFACE *mris, LABEL *asrc, MRIS_HASH_TABLE *mht, MRI_SURFACE *mris_dst, LABEL *adst)
{
  int vno, n, nfilled, m;
  LABEL_VERTEX *lv_dst;
  MRIS_HASH_TABLE *mht_src;
  double max_len;

  if (!adst) {
    adst = LabelClone(asrc);
  }

  if (adst->max_points < asrc->n_points + adst->n_points) /* won't fit - expand */
  {
    LABEL *atmp;

    atmp = LabelAlloc(2 * (asrc->n_points + adst->n_points), asrc->subject_name, asrc->name);
    LabelCopy(adst, atmp);
    LabelFree(&adst);
    adst = atmp;
  }

  MRISclearMarks(mris_dst);
  for (n = 0; n < asrc->n_points; n++) {
    vno = asrc->lv[n].vno;
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    VERTEX const * const v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    float min_dist;
    VERTEX * const vdst = MHTfindClosestVertex2(mht, mris_dst, mris, v, &min_dist);
    if (!vdst) {
      ErrorPrintf(ERROR_BADPARM, "MRIScombine: cannot map vno %d", vno);
      continue;
    }
    if (vdst->marked) /* only add this vertex once */
    {
      continue;
    }
    vdst->marked = 1;
    if (vdst - mris_dst->vertices == Gdiag_no) {
      DiagBreak();
    }
    lv_dst = labelFindVertexNumber(adst, vdst - mris_dst->vertices);
    if (lv_dst == NULL) {
      lv_dst = &adst->lv[adst->n_points++];
    }
    else if (lv_dst - adst->lv == Gdiag_no) {
      DiagBreak();
    }
    lv_dst->vno = vdst - mris_dst->vertices;
    if (lv_dst->vno == 60008) {
      DiagBreak();
    }
    if (lv_dst->vno == 85592) {
      DiagBreak();
    }
    if (lv_dst->vno == Gdiag_no) {
      DiagBreak();
    }
    lv_dst->x = asrc->lv[n].x;
    lv_dst->y = asrc->lv[n].y;
    lv_dst->z = asrc->lv[n].z;
    lv_dst->stat += asrc->lv[n].stat;
  }

  MRIScomputeVertexSpacingStats(mris, NULL, NULL, &max_len, NULL, NULL, CURRENT_VERTICES);
  mht_src = MHTcreateVertexTable_Resolution(mris, CURRENT_VERTICES, 2 * max_len);

  MRISclearMarks(mris);
  LabelMarkStats(asrc, mris);

  do {
    /* map the nbr of every point in the label back to src, and if in
       label add it
    */
    nfilled = 0;
    for (n = 0; n < adst->n_points; n++) {
      VERTEX_TOPOLOGY const * const vt = &mris_dst->vertices_topology[adst->lv[n].vno];
      VERTEX                * const v  = &mris_dst->vertices         [adst->lv[n].vno];
      if (adst->lv[n].vno == Gdiag_no) {
        DiagBreak();
      }
      if (n == Gdiag_no) {
        DiagBreak();
      }
      if (v->marked == 0) /* hasn't been processed for this surface yet */
      {
        float min_dist;
        VERTEX const * const vsrc = MHTfindClosestVertex2(mht_src, mris, mris_dst, v, &min_dist);
        if (vsrc && vsrc->marked) /* in label */
        {
          adst->lv[n].stat += vsrc->stat;
          v->marked = 1;
        }
      }
      for (m = 0; m < vt->vnum; m++) {
        VERTEX * const vn = &mris_dst->vertices[vt->v[m]];
        if (vn->marked) {
          continue; /* already in label */
        }
        float min_dist;
        VERTEX const * const vsrc = MHTfindClosestVertex2(mht_src, mris, mris_dst, vn, &min_dist);
        if (vsrc == NULL) {
          DiagBreak();
        }
        if (vsrc - mris->vertices == 62644) {
          DiagBreak();
        }
        if (vsrc - mris->vertices == Gdiag_no) {
          DiagBreak();
        }
        if (vsrc->marked) {
          if (adst->n_points >= adst->max_points - 1) {
            LABEL *atmp;

            atmp = LabelAlloc(2 * adst->n_points, adst->subject_name, adst->name);
            LabelCopy(adst, atmp);
            LabelFree(&adst);
            adst = atmp;
          }

          vn->marked = 1;
          lv_dst = labelFindVertexNumber(adst, vt->v[m]);
          if (lv_dst == NULL) {
            lv_dst = &adst->lv[adst->n_points++];
          }
          else if (lv_dst - adst->lv == Gdiag_no) {
            DiagBreak();
          }
          lv_dst->x = adst->lv[n].x;
          lv_dst->y = adst->lv[n].y;
          lv_dst->z = adst->lv[n].z;
          lv_dst->vno = vt->v[m];
          lv_dst->stat += vsrc->stat;
          if (lv_dst->vno == Gdiag_no) {
            DiagBreak();
          }
          if (lv_dst->vno == 55185) {
            DiagBreak();
          }
          if (lv_dst->vno == 85592) {
            DiagBreak();
          }
          nfilled++;
        }
      }
    }
  } while (nfilled != 0);

  return (adst);
}

int LabelNormalizeStats(LABEL *area, float norm)
{
  LABEL_VERTEX *lv;
  int n;

  for (n = 0; n < area->n_points; n++) {
    lv = &area->lv[n];
    lv->stat /= norm;
  }
  return (NO_ERROR);
}

/*-------------------------------------------------------------
  MaskSurfLabel() - removes vertices from a label based on the
  value of those vertices in a surface mask.
  -------------------------------------------------------------*/
LABEL *MaskSurfLabel(LABEL *lbl, MRI *SurfMask, float thresh, const char *masksign, int frame)
{
  LABEL *msklbl;
  int n, vno, ok, noverlap;
  float maskval;

  if (frame >= SurfMask->nframes) {
    printf(
          "ERROR:  MaskSurfLabel: specified frame %d is greater "
          "than the number of frame %d\n",
          frame,
          SurfMask->nframes);
    return (NULL);
  }

  /* Check that the label is a surface-based label */
  ok = 0;
  for (n = 0; n < lbl->n_points; n++) {
    if (lbl->lv[n].vno != 0) {
      ok = 1;
    }
  }
  if (!ok) {
    printf("ERROR: MaskSurfLabel: all label vertices are 0, proably\n");
    printf("       not a label created from a surface\n");
    return (NULL);
  }

  /* -- Allocate the Target Label ---*/
  msklbl = LabelAlloc(lbl->n_points, lbl->subject_name, "labelfilenamehere");
  msklbl->n_points = lbl->n_points;

  /* Compute the number of overlapping vertices */
  noverlap = 0;
  for (n = 0; n < lbl->n_points; n++) {
    vno = lbl->lv[n].vno;
    if (vno >= SurfMask->width) {
      printf(
            "ERROR: MaskSurfLabel: label vertex %d has vno %d which "
            "is larger than the number of vertices %d in mask\n",
            n,
            vno,
            SurfMask->width);
      return (NULL);
    }
    maskval = MRIgetVoxVal(SurfMask, vno, 0, 0, frame);
    if (!strcmp(masksign, "abs")) {
      maskval = fabs(maskval);
    }
    if (!strcmp(masksign, "neg")) {
      maskval = -maskval;
    }
    // printf("%4d  %6d  %g %g  %d\n",n,vno,maskval,thresh,maskval>thresh);
    if (maskval < thresh) {
      continue;
    }
    msklbl->lv[noverlap].vno = vno;
    msklbl->lv[noverlap].x = lbl->lv[noverlap].x;
    msklbl->lv[noverlap].y = lbl->lv[noverlap].y;
    msklbl->lv[noverlap].z = lbl->lv[noverlap].z;
    msklbl->lv[noverlap].stat = lbl->lv[noverlap].stat;
    noverlap++;
  }
  if (noverlap == 0) {
    printf("WARNING: MaskSurfLabel: no overlap between label and mask\n");
    return (NULL);
  }
  msklbl->n_points = noverlap;

  return (msklbl);
}

int LabelErode(LABEL *area, MRI_SURFACE *mris, int num_times)
{
  int n, label_vno, vno, neighbor_index, neighbor_vno, found_nbr_off;
  VERTEX *vn;

  if (NULL == area) {
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "LabelErode: NULL label"));
  }
  if (NULL == mris) {
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "LabelErode: NULL mris"));
  }
  if (num_times < 1) {
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "LabelErode: num_times < 1"));
  }

  MRISclearMarks(mris);
  for (n = 0; n < num_times; n++) {
    LabelMark(area, mris);  // all vertices in label now have v->marked==1

    /* For each vertex in the label... */
    for (label_vno = 0; label_vno < area->n_points; label_vno++) {
      if (area->lv[label_vno].deleted) continue;
      vno = area->lv[label_vno].vno;
      if (vno < 0) continue;
      if (vno >= mris->nvertices)
        ErrorExit(ERROR_BADPARM, "LabelErode: label vertex %d too big for surface (%d)", vno, mris->nvertices);
      if (vno == Gdiag_no) DiagBreak();

      // check to see if we should not add this label
      // (if one of it's nbrs is not in label)
      found_nbr_off = 0;
      for (neighbor_index = 0; neighbor_index < mris->vertices_topology[vno].vnum; neighbor_index++) {
        neighbor_vno = mris->vertices_topology[vno].v[neighbor_index];
        if (neighbor_vno == Gdiag_no) DiagBreak();
        if (neighbor_vno < 0) continue;
        if (neighbor_vno >= mris->nvertices)
          ErrorExit(ERROR_BADPARM,
                    "LabelErode: neighbor label vertex %d too big for surface (%d)",
                    neighbor_vno,
                    mris->nvertices);

        /* Look for neighbor_vno in the label. */
        vn = &mris->vertices[neighbor_vno];
        if (vn->marked == 0)  // found a nbr not in the label
        {
          if (neighbor_vno == Gdiag_no) DiagBreak();

          found_nbr_off = 1;
          break;
        }
      }

      if (found_nbr_off)  // all nbrs on - add it to the new label
        area->lv[label_vno].deleted = 1;
    }
    LabelUnmark(area, mris);
  }

  update_vertex_indices(area);
  //  printf("area->max_points = %d\n",area->max_points) ;
  MRISclearMarks(mris);
  return (NO_ERROR);
}

int LabelDilate(LABEL *area, MRI_SURFACE *mris, int num_times, int coords)
{
  int n, neighbor_index, neighbor_vno, found, vno;

  //  printf("LabelDilate(%d, %d)\n", num_times, coords) ;
  if (NULL == area) {
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "LabelDilate: NULL label"));
  }
  if (NULL == mris) {
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "LabelDilate: NULL mris"));
  }
  if (num_times < 1) {
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "LabelDilate: num_times < 1"));
  }

  MRISclearMarks(mris);
  for (n = 0; n < num_times; n++) {
    LabelMark(area, mris);  // all vertices in label now have v->marked==1

    /* For each vertex in the label... */
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX          const * const v  = &mris->vertices         [vno];
      if (vno == Gdiag_no) DiagBreak();

      if (v->marked == 1)  // already in label
        continue;

      // Check its neighbors. If any are in the label, add it
      found = 0;
      for (neighbor_index = 0; neighbor_index < vt->vnum; neighbor_index++) {
        /* Look for neighbor_vno in the label. */
        neighbor_vno = mris->vertices_topology[vno].v[neighbor_index];
        VERTEX const * const vn = &mris->vertices[neighbor_vno];
        if (vn->marked > 0) {
          if (neighbor_vno == Gdiag_no) DiagBreak();

          found = 1;
          break;
        }
      }

      // add it if at least one nbr was found that was in label
      if (found) {
        LV *lv;
        int n;

        if (area->n_points >= area->max_points) LabelRealloc(area, nint(area->max_points * 1.5));

        if (vno == Gdiag_no) DiagBreak();

        n = area->n_points++;
        lv = &area->lv[n];
        lv->vno = vno;
        MRISgetCoords(v, coords, &lv->x, &lv->y, &lv->z);
        if (area->vertex_label_ind) area->vertex_label_ind[vno] = n;
        if (area->mris && area->mri_template) {
          double xv, yv, zv;
          MRISsurfaceRASToVoxel((MRIS *)area->mris, area->mri_template, lv->x, lv->y, lv->z, &xv, &yv, &zv);
          lv->xv = nint(xv);
          lv->yv = nint(yv);
          lv->zv = nint(zv);
        }
        //	printf("LabelDilate: added vertex %d (%d)\n", vno, n) ;
      }
    }

    /* Point the label's lv to the new one and update the number of
      points. */
    LabelUnmark(area, mris);
  }

  update_vertex_indices(area);

  //  printf("area->max_points = %d\n",area->max_points) ;
  return (NO_ERROR);
}

static LABEL_VERTEX *labelFindVertexNumber(LABEL *area, int vno)
{
  int n;
  LABEL_VERTEX *lv;

  for (n = 0; n < area->n_points; n++) {
    lv = &area->lv[n];
    if (lv->vno == vno) {
      return (lv);
    }
  }
  return (NULL);
}

int LabelSetStat(LABEL *area, float stat)
{
  int n;

  for (n = 0; n < area->n_points; n++) {
    area->lv[n].stat = stat;
  }

  return (NO_ERROR);
}

int LabelMarkStats(LABEL *area, MRI_SURFACE *mris)
{
  int n, vno;
  VERTEX *v;

  for (n = 0; n < area->n_points; n++) {
    if (area->lv[n].deleted > 0) {
      continue;
    }
    vno = area->lv[n].vno;
    v = &mris->vertices[vno];
    v->marked = 1;
    v->stat = area->lv[n].stat;
  }
  return (NO_ERROR);
}

LABEL *LabelFillHoles(LABEL *area_src, MRI_SURFACE *mris, int coords)
{
  MRI *mri;
  int i, dst_index, vno, xi, yi, zi, xk, yk, zk, found, x2, y2, z2, nchanged;
  double xw, yw, zw, xv, yv, zv;
  VERTEX *v;
  LABEL *area_dst;
  float vx, vy, vz, dist, dx, dy, dz;

  vx = vy = vz = -1;
  mri = MRIalloc(256, 256, 256, MRI_UCHAR);

  useVolGeomToMRI(&mris->vg, mri);
  area_dst = LabelAlloc(mris->nvertices, mris->subject_name, area_src->name);
  LabelCopy(area_src, area_dst);
  LabelFillUnassignedVertices(mris, area_dst, coords);
  LabelMarkSurface(area_dst, mris);

  for (i = 0; i < area_src->n_points; i++) {
    xw = area_src->lv[i].x;
    yw = area_src->lv[i].y;
    zw = area_src->lv[i].z;
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
    MRISsurfaceRASToVoxel(mris, mri, xw, yw, zw, &xv, &yv, &zv);
#endif
    MRIvox(mri, nint(xv), nint(yv), nint(zv)) = 1;
  }

  dst_index = area_dst->n_points;
  do {
    nchanged = 0;
    for (vno = 0; vno < mris->nvertices; vno++) {
      v = &mris->vertices[vno];
      if (vno == Gdiag_no) {
        DiagBreak();
      }
      if (v->ripflag || v->marked) /* already in label */
      {
        continue;
      }
      MRISvertexCoordToVoxel(mris, v, mri, coords, &xv, &yv, &zv);
      xi = nint(xv);
      yi = nint(yv);
      zi = nint(zv);
      if (xi >= mri->width) {
        xi = mri->width - 1;
      }
      if (yi >= mri->height) {
        yi = mri->height - 1;
      }
      if (zi >= mri->depth) {
        zi = mri->depth - 1;
      }
      if (MRIvox(mri, xi, yi, zi) == 1) {
        MRISgetCoords(v, coords, &vx, &vy, &vz);
        area_dst->lv[dst_index].vno = vno;
        area_dst->lv[dst_index].x = vx;
        area_dst->lv[dst_index].y = vy;
        area_dst->lv[dst_index].z = vz;
        v->marked = 1;
        dst_index++;
        nchanged++;
        area_dst->n_points++;
      }
      if (MRIneighborsOn(mri, xi, yi, zi, 1) > 0) {
        found = 0;
        for (xk = -1; !found && xk <= 1; xk++)
          for (yk = -1; !found && yk <= 1; yk++)
            for (zk = -1; !found && zk <= 1; zk++) {
              x2 = mri->xi[xi + xk];
              y2 = mri->xi[yi + yk];
              z2 = mri->xi[zi + zk];
              if (MRIvox(mri, x2, y2, z2) == 1) {
                dx = xv - (xi + xk);
                dy = yv - (yi + yk);
                dz = zv - (zi + zk);
                dist = sqrt(SQR(dx) + SQR(dy) + SQR(dz));
                if (dist < sqrt(.5 * .5 + .5 * .5)) {
                  found = 1;
                  MRISgetCoords(v, coords, &vx, &vy, &vz);
                  area_dst->lv[dst_index].vno = vno;
                  area_dst->lv[dst_index].x = vx;
                  area_dst->lv[dst_index].y = vy;
                  area_dst->lv[dst_index].z = vz;
                  dst_index++;
                  v->marked = 1;
                  nchanged++;
                  area_dst->n_points++;
                  break;
                }
              }
            }
      }
    }
  } while (nchanged > 0);

  MRIfree(&mri);
  return (area_dst);
}
LABEL *LabelFillHolesWithOrig(LABEL *area_src, MRI_SURFACE *mris)
{
  MRI *mri;
  int i, dst_index, vno;
  double xw, yw, zw, xv, yv, zv;
  VERTEX *v;
  LABEL *area_dst;

  mri = MRIalloc(256, 256, 256, MRI_UCHAR);
  area_dst = LabelAlloc(mris->nvertices, mris->subject_name, area_src->name);
  LabelCopy(area_src, area_dst);
  LabelFillUnassignedVertices(mris, area_dst, ORIG_VERTICES);
  LabelMarkSurface(area_dst, mris);

  for (i = 0; i < area_src->n_points; i++) {
    xw = area_src->lv[i].x;
    yw = area_src->lv[i].y;
    zw = area_src->lv[i].z;
    // MRIworldToVoxel(mri, xw, yw, zw, &xv, &yv, &zv) ;
    if (mris->useRealRAS) {
      MRIworldToVoxel(mri, xw, yw, zw, &xv, &yv, &zv);
    }
    else {
      MRIsurfaceRASToVoxel(mri, xw, yw, zw, &xv, &yv, &zv);
    }
    MRIvox(mri, nint(xv), nint(yv), nint(zv)) = 1;
  }

  dst_index = area_dst->n_points;
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag || v->marked) /* already in label */
    {
      continue;
    }
    MRISorigVertexToVoxel(mris, v, mri, &xv, &yv, &zv);
    if (MRIvox(mri, nint(xv), nint(yv), nint(zv)) == 1) {
      area_dst->lv[dst_index].vno = vno;
      area_dst->lv[dst_index].x = v->origx;
      area_dst->lv[dst_index].y = v->origy;
      area_dst->lv[dst_index].z = v->origz;
      dst_index++;
      area_dst->n_points++;
    }
  }

  MRIfree(&mri);
  return (area_dst);
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
    if (lb->lv[n].vno == vtxno) {
      return (n);
    }
  return (-1);
}
/*---------------------------------------------------------------
  LabelRealloc() - reallocates the number of label vertices
  to have max_points. If something goes wrong, returns 1 without
  changing anything. Otherwise returns 0.
  ---------------------------------------------------------------*/
LABEL *LabelRealloc(LABEL *lb, int max_points)
{
  LV *lvtmp;

  if (lb == NULL) return (LabelAlloc(max_points, NULL, NULL));
  if (max_points <= lb->max_points) {
    return (lb);
  }

  lvtmp = (LV *)realloc(lb->lv, sizeof(LV) * max_points);
  if (lvtmp == NULL) {
    return (lb);
  }
  memset(lvtmp + lb->max_points, 0, (max_points - lb->max_points) * sizeof(LABEL_VERTEX));
  lb->max_points = max_points;
  lb->lv = lvtmp;

  return (lb);
}
/*------------------------------------------------------------
  LabelfromASeg() - creates a label from a segmentation given
  the numeric segmentation code. If there are no voxels that
  match the code, then it silently returns NULL. Note: this
  generates the same result as the code in mri_cor2label.c.
  ------------------------------------------------------------*/
LABEL *LabelfromASeg(MRI *aseg, int segcode)
{
  int c, r, s, nlabel, v;
  LABEL *lb;
  MATRIX *Vox2RAS, *crs, *ras;
  // char labelfile[100];

  // Count number of label points first
  nlabel = 0;
  for (c = 0; c < aseg->width; c++) {
    for (r = 0; r < aseg->height; r++) {
      for (s = 0; s < aseg->depth; s++) {
        v = (int)MRIgetVoxVal(aseg, c, r, s, 0);
        if (v == segcode) {
          nlabel++;
        }
      }
    }
  }
  // printf("Found %d voxels in label\n",nlabel);
  if (nlabel == 0) {
    return (NULL);
  }

  lb = LabelAlloc(nlabel, "dontknow", "dontcare");

  Vox2RAS = MRIxfmCRS2XYZtkreg(aseg);
  // printf("ASeg Vox2RAS---------------------------\n");
  // MatrixPrint(stdout,Vox2RAS);
  // printf("---------------------------\n");
  crs = MatrixAlloc(4, 1, MATRIX_REAL);
  crs->rptr[4][1] = 1;
  ras = MatrixAlloc(4, 1, MATRIX_REAL);

  nlabel = 0;
  for (c = 0; c < aseg->width; c++) {
    for (r = 0; r < aseg->height; r++) {
      for (s = 0; s < aseg->depth; s++) {
        v = (int)MRIgetVoxVal(aseg, c, r, s, 0);
        if (v == segcode) {
          crs->rptr[1][1] = c;
          crs->rptr[2][1] = r;
          crs->rptr[3][1] = s;
          ras = MatrixMultiply(Vox2RAS, crs, ras);
          lb->lv[nlabel].x = ras->rptr[1][1];
          lb->lv[nlabel].y = ras->rptr[2][1];
          lb->lv[nlabel].z = ras->rptr[3][1];
          lb->lv[nlabel].vno = -1;  // not assigned yet
          nlabel++;
        }
      }
    }
  }
  lb->n_points = nlabel;
  MatrixFree(&Vox2RAS);
  MatrixFree(&crs);
  MatrixFree(&ras);

  // sprintf(labelfile,"./tmp-%d.label",segcode);
  // LabelWrite(lb,labelfile);

  return (lb);
}

/*---------------------------------------------------------------
  VertexIsInLabel() - returns a 1 if the given vertex number is
  in the label. Label must be surface-based.
  ---------------------------------------------------------------*/
int VertexIsInLabel(int vtxno, LABEL *label)
{
  int n;

  for (n = 0; n < label->n_points; n++)
    if (label->lv[n].vno == vtxno) {
      return (1);
    }

  return (0);
}
/*---------------------------------------------------------------
  LabelBoundary() - returns a label of all the points in the input
  label on the edge/boundary. Label must be surface-based.
  ---------------------------------------------------------------*/
LABEL *LabelBoundary(LABEL *label, MRIS *surf)
{
  int n, nnbrs, nthnbr, vtxno, nbrvtxno;
  LABEL *boundary;

  boundary = LabelAlloc(label->n_points, "", "");
  boundary->n_points = 0;

  for (n = 0; n < label->n_points; n++) {
    vtxno = label->lv[n].vno;
    nnbrs = surf->vertices_topology[vtxno].vnum;
    for (nthnbr = 0; nthnbr < nnbrs; nthnbr++) {
      nbrvtxno = surf->vertices_topology[vtxno].v[nthnbr];
      if (!VertexIsInLabel(nbrvtxno, label)) {
        boundary->lv[boundary->n_points].vno = vtxno;
        boundary->lv[boundary->n_points].x = label->lv[n].x;
        boundary->lv[boundary->n_points].y = label->lv[n].y;
        boundary->lv[boundary->n_points].z = label->lv[n].z;
        boundary->lv[boundary->n_points].stat = label->lv[n].stat;
        boundary->n_points++;
        break;
      }
    }
  }

  return (boundary);
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

  if (order > 3) {
    printf("ERROR: LabelFitXYZ(): order = %d, must be <= 3\n", order);
    return (NULL);
  }

  X = MatrixAlloc(label->n_points, 1 + 2 * order, MATRIX_REAL);
  y = MatrixAlloc(label->n_points, 1, MATRIX_REAL);

  // Load the design matrix
  for (n = 0; n < label->n_points; n++) {
    X->rptr[n + 1][1] = 1;               // Intercept
    X->rptr[n + 1][2] = label->lv[n].x;  // Slope with x
    X->rptr[n + 1][3] = label->lv[n].y;  // Slope with y
    if (order > 1) {
      X->rptr[n + 1][4] = pow(label->lv[n].x, 2);  // Slope with x^2
      X->rptr[n + 1][5] = pow(label->lv[n].y, 2);  // Slope with y^2
    }
    if (order > 2) {
      X->rptr[n + 1][6] = pow(label->lv[n].x, 3);  // Slope with x^3
      X->rptr[n + 1][7] = pow(label->lv[n].y, 3);  // Slope with y^3
    }
    y->rptr[n + 1][1] = label->lv[n].z;  // Fit z
  }

  Xt = MatrixTranspose(X, NULL);
  XtX = MatrixMultiply(Xt, X, NULL);
  iXtX = MatrixInverse(XtX, NULL);
  iXtXXt = MatrixMultiply(iXtX, Xt, NULL);
  beta = MatrixMultiply(iXtXXt, y, NULL);

  MatrixFree(&X);
  MatrixFree(&Xt);
  MatrixFree(&XtX);
  MatrixFree(&iXtX);
  MatrixFree(&iXtXXt);
  MatrixFree(&y);

  // printf("beta -----------------\n");
  // MatrixPrint(stdout,beta);
  // printf("----------------------\n");

  return (beta);
}
LABEL *LabelInFOV(MRI_SURFACE *mris, MRI *mri, float pad)
{
  LABEL *area;
  int vno, npoints;
  VERTEX *v;
  LABEL_VERTEX *lv;
  double xv, yv, zv, nvox;

  // convert to voxels
  nvox = pad / MAX(MAX(mri->xsize, mri->ysize), mri->zsize);

  for (npoints = vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
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
    MRISsurfaceRASToVoxel(mris, mri, v->x, v->y, v->z, &xv, &yv, &zv);
#endif
    if (xv < nvox || yv < nvox || zv < nvox || xv > (mri->width - 1) - nvox || yv > (mri->height - 1) - nvox ||
        zv > (mri->depth - 1) - nvox) {
      continue;
    }
    npoints++;
  }

  area = LabelAlloc(npoints, NULL, NULL);
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
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
    MRISsurfaceRASToVoxel(mris, mri, v->x, v->y, v->z, &xv, &yv, &zv);
#endif
    if (xv < nvox || yv < nvox || zv < nvox || xv > (mri->width - 1) - nvox || yv > (mri->height - 1) - nvox ||
        zv > (mri->depth - 1) - nvox) {
      continue;
    }
    lv = &area->lv[area->n_points];
    area->n_points++;
    lv->vno = vno;
    lv->x = v->x;
    lv->y = v->y;
    lv->z = v->z;
  }

  // LabelWrite(area, "./lh.all.label");
  return (area);
}
LABEL *LabelTranslate(LABEL *area, LABEL *area_offset, float dx, float dy, float dz)
{
  int i;

  if (area_offset == NULL) {
    area_offset = LabelCopy(area, NULL);
  }

  for (i = 0; i < area->n_points; i++) {
    area_offset->lv[i].x = area->lv[i].x + dx;
    area_offset->lv[i].y = area->lv[i].y + dy;
    area_offset->lv[i].z = area->lv[i].z + dz;
  }
  return (area_offset);
}

int LabelUnassign(LABEL *area)
{
  int i;

  for (i = 0; i < area->n_points; i++) {
    area->lv[i].vno = -1;
  }
  return (NO_ERROR);
}

/*---------------------------------------------------------------
  MRISlabelInvert(MRIS *surf, LABEL *label)
  Creates a label with all the vertices NOT in the input label.
  ---------------------------------------------------------------*/
LABEL *MRISlabelInvert(MRIS *surf, LABEL *label)
{
  MRI *tmpmri;
  LABEL *invlabel;
  int n, vtxno, ninv, nlabel;

  // Create binary mask of label
  tmpmri = MRISlabel2Mask(surf, label, NULL);

  // Count number of points. Need to do this in case duplicates
  nlabel = 0;
  for (vtxno = 0; vtxno < surf->nvertices; vtxno++)
    if (MRIgetVoxVal(tmpmri, vtxno, 0, 0, 0) > 0.5) {
      nlabel++;
    }

  // Alloc inverse label
  ninv = surf->nvertices - nlabel;
  invlabel = LabelAlloc(ninv, label->subject_name, NULL);
  invlabel->n_points = ninv;

  // Assign label points to vtxs NOT in the mask
  n = 0;
  for (vtxno = 0; vtxno < surf->nvertices; vtxno++) {
    if (MRIgetVoxVal(tmpmri, vtxno, 0, 0, 0) < 0.5) {
      invlabel->lv[n].vno = vtxno;
      invlabel->lv[n].x = surf->vertices[vtxno].x;
      invlabel->lv[n].y = surf->vertices[vtxno].y;
      invlabel->lv[n].z = surf->vertices[vtxno].z;
      invlabel->lv[n].stat = 0.0;
      n++;
    }
  }
  MRIfree(&tmpmri);
  return (invlabel);
}
int LabelCopyStatsToSurface(LABEL *area, MRI_SURFACE *mris, int which)
{
  int n, vno;
  VERTEX *v;

  for (n = 0; n < area->n_points; n++) {
    vno = area->lv[n].vno;
    v = &mris->vertices[vno];
    v->marked = 1;
    switch (which) {
    case VERTEX_VALS:
      v->val = area->lv[n].stat;
      break;
    case VERTEX_CURV:
      v->curv = area->lv[n].stat;
      break;
    case VERTEX_STATS:
      v->stat = area->lv[n].stat;
      break;
    default:
      ErrorExit(ERROR_BADPARM, "LabelCopyStatsToSurface: unsupported which = %d", which);
    }
  }
  return (NO_ERROR);
}

int LabelThreshold(LABEL *area, float thresh)
{
  int n, ndel;
  LABEL_VERTEX *lv;

  for (ndel = n = 0; n < area->n_points; n++) {
    lv = &area->lv[n];
    if (lv->stat < thresh) {
      lv->deleted = 1;
      ndel++;
    }
  }
  return (ndel);
}

int LabelCropAnterior(LABEL *area, float anterior_dist)
{
  int n;
  float amax;

  amax = -100000;
  for (n = 0; n < area->n_points; n++) {
    if (area->lv[n].y > amax) {
      amax = area->lv[n].y;
    }
  }

  amax -= anterior_dist;
  printf("cropping all vertices with Y > %2.0f\n", amax);
  for (n = 0; n < area->n_points; n++) {
    if (area->lv[n].y < amax) {
      area->lv[n].deleted = 1;
    }
  }
  return (NO_ERROR);
}

int LabelCropPosterior(LABEL *area, float anterior_dist)
{
  int n;
  float amax;

  amax = -100000;
  for (n = 0; n < area->n_points; n++) {
    if (area->lv[n].y > amax) {
      amax = area->lv[n].y;
    }
  }

  amax += anterior_dist;
  printf("cropping all vertices with Y > %2.0f\n", amax);
  for (n = 0; n < area->n_points; n++) {
    if (area->lv[n].y > amax) {
      area->lv[n].deleted = 1;
    }
  }
  return (NO_ERROR);
}

int
LabelMaskSurfaceValues(LABEL *label, MRI_SURFACE *mris)
{
  return(LabelMaskSurface(label, mris)) ;
}


int
LabelMaskSurfaceCurvature(LABEL *area, MRI_SURFACE *mris) 
{
  int vno;
  VERTEX *v;

  LabelMarkSurface(area, mris); /* mark all points in label */

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->marked || v->ripflag) {
      continue;
    }
    v->curv = 0.0;
  }
  MRISclearMarks(mris);
  return (NO_ERROR);
}

int LabelMaskSurface(LABEL *area, MRI_SURFACE *mris)
{
  int vno;
  VERTEX *v;

  LabelMarkSurface(area, mris); /* mark all points in label */

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->marked || v->ripflag) {
      continue;
    }
    v->stat = v->annotation = v->val = v->imag_val = v->val2 = v->valbak = v->val2bak = 0.0;
  }
  MRISclearMarks(mris);
  return (NO_ERROR);
}
int LabelMaskSurfaceVolume(LABEL *area, MRI *mri, float nonmask_val)
{
  int vno, n, x, y, z, f;
  uchar *marked;
  long nvox;

  nvox = mri->width * mri->height * mri->depth * mri->nframes;
  marked = (uchar *)calloc(nvox, sizeof(uchar));
  for (n = 0; n < area->n_points; n++) {
    vno = area->lv[n].vno;
    if (vno >= 0) marked[vno] = 1;
  }
  for (n = x = 0; x < mri->width; x++)
    for (y = 0; y < mri->height; y++)
      for (z = 0; z < mri->depth; z++)
        for (f = 0; f < mri->nframes; f++, n++) {
          if (marked[n]) continue;
          MRIsetVoxVal(mri, x, y, z, f, nonmask_val);
        }

  free(marked);
  return (NO_ERROR);
}
int LabelCentroid(LABEL *area, MRI_SURFACE *mris, double *px, double *py, double *pz, int *pvno)
{
  int vno, num, n;
  VERTEX *v;
  double xc, yc, zc;

  xc = yc = zc = 0.0;
  for (num = n = 0; n < area->n_points; n++) {
    vno = area->lv[n].vno;
    v = &mris->vertices[vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    num++;
    xc += v->x;
    yc += v->y;
    zc += v->z;
  }
  *px = xc / num;
  *py = yc / num;
  *pz = zc / num;
  if (pvno != NULL) {
    int vno_closest;
    double min_dist, dist;

    xc /= num;
    yc /= num;
    zc /= num;
    vno_closest = -1;
    min_dist = 1e10;
    for (n = 0; n < area->n_points; n++) {
      vno = area->lv[n].vno;
      v = &mris->vertices[vno];
      dist = SQR(v->x - xc) + SQR(v->y - yc) + SQR(v->z - zc);
      if (dist < min_dist) {
        min_dist = dist;
        vno_closest = vno;
      }
    }
    *pvno = vno_closest;
  }
  return (NO_ERROR);
}

int LabelFillVolume(MRI *mri, LABEL *area, int fillval)
{
  int n;
  double xv, yv, zv;
  LABEL_VERTEX *lv;

  for (n = 0; n < area->n_points; n++) {
    lv = &area->lv[n];
    MRIsurfaceRASToVoxel(mri, lv->x, lv->y, lv->z, &xv, &yv, &zv);
    MRIsetVoxVal(mri, nint(xv), nint(yv), nint(zv), 0, fillval);
  }
  return (NO_ERROR);
}

int LabelSetVals(MRI_SURFACE *mris, LABEL *area, float fillval)
{
  int n;
  LABEL_VERTEX *lv;

  for (n = 0; n < area->n_points; n++) {
    lv = &area->lv[n];
    if (lv->deleted || lv->vno < 0) continue;
    mris->vertices[lv->vno].val = fillval;
  }
  return (NO_ERROR);
}
/*
  convert the label coords from tkreg (surface) RAS to scanner RAS. Note that this assumes that the
  label coords are in tkreg space
*/
LABEL *LabelToScannerRAS(LABEL *lsrc, MRI *mri, LABEL *ldst)
{
  int i;
  MATRIX *M_surface_to_RAS = RASFromSurfaceRAS_(mri);
  VECTOR *v1, *v2;

  if (ldst != lsrc)
    ldst = LabelCopy(lsrc, ldst) ;

  if (lsrc->coords == LABEL_COORDS_SCANNER_RAS) // already in the right space
    return(ldst) ;
  if (lsrc->coords == LABEL_COORDS_VOXEL)
    LabelToSurfaceRAS(lsrc, mri, ldst) ;

  M_surface_to_RAS = RASFromSurfaceRAS_(mri);
  v1 = VectorAlloc(4, MATRIX_REAL);
  v2 = VectorAlloc(4, MATRIX_REAL);
  VECTOR_ELT(v1, 4) = 1.0;
  VECTOR_ELT(v2, 4) = 1.0;
  for (i = 0; i < lsrc->n_points; i++) {
    V3_X(v1) = ldst->lv[i].x;
    V3_Y(v1) = ldst->lv[i].y;
    V3_Z(v1) = ldst->lv[i].z;
    MatrixMultiply(M_surface_to_RAS, v1, v2);
    ldst->lv[i].x = V3_X(v2);
    ldst->lv[i].y = V3_Y(v2);
    ldst->lv[i].z = V3_Z(v2);
  }
  strncpy(ldst->space, "scanner", sizeof(ldst->space));
  ldst->coords = LABEL_COORDS_SCANNER_RAS;
  VectorFree(&v1);
  VectorFree(&v2);
  MatrixFree(&M_surface_to_RAS);
  return (ldst);
}
/*
  convert the label coords to tkreg (surface) RAS

*/
LABEL *LabelToSurfaceRAS(LABEL *lsrc, MRI *mri, LABEL *ldst)
{
  int i;
  MATRIX *M_surface_to_RAS = RASFromSurfaceRAS_(mri), *M_surface_from_RAS;
  VECTOR *v1, *v2;

  if (ldst != lsrc)
    ldst = LabelCopy(lsrc, ldst) ;
  switch (lsrc->coords)
  {
  case LABEL_COORDS_TKREG_RAS:
    return(ldst) ;  // already done
  case LABEL_COORDS_VOXEL:
    return(LabelVoxelToSurfaceRAS(lsrc, mri, ldst)) ;
  case LABEL_COORDS_SCANNER_RAS:
    break ; // will be done below
  default:
    ErrorExit(ERROR_UNSUPPORTED, "LabelToSurfaceRAS: unsupported coords (was %d)",lsrc->coords) ;
    break ;
  }
  M_surface_from_RAS = MatrixInverse(M_surface_to_RAS, NULL);

  v1 = VectorAlloc(4, MATRIX_REAL);
  v2 = VectorAlloc(4, MATRIX_REAL);
  VECTOR_ELT(v1, 4) = 1.0;
  VECTOR_ELT(v2, 4) = 1.0;
  for (i = 0; i < lsrc->n_points; i++) {
    V3_X(v1) = ldst->lv[i].x;
    V3_Y(v1) = ldst->lv[i].y;
    V3_Z(v1) = ldst->lv[i].z;
    MatrixMultiply(M_surface_from_RAS, v1, v2);
    ldst->lv[i].x = V3_X(v2);
    ldst->lv[i].y = V3_Y(v2);
    ldst->lv[i].z = V3_Z(v2);
  }
  strcpy(ldst->space, "TkReg");
  ldst->coords = LABEL_COORDS_TKREG_RAS;
  VectorFree(&v1);
  VectorFree(&v2);
  MatrixFree(&M_surface_to_RAS);
  MatrixFree(&M_surface_from_RAS);
  return (ldst);
}

/*
  convert the label coords from tkreg (surface) RAS to voxels. Note that this assumes that the
  label coords are in tkreg space
*/
LABEL *LabelToVoxel(LABEL *lsrc, MRI *mri, LABEL *ldst)
{
  int i;
  MATRIX *M_surface_to_vox ;
  VECTOR *v1, *v2;

  if (ldst != lsrc)
    ldst = LabelCopy(lsrc, ldst) ;
  if (lsrc->coords == LABEL_COORDS_VOXEL)  // alread done
    return(ldst) ;

  if (lsrc->coords == LABEL_COORDS_SCANNER_RAS) // already in the right space
    LabelFromScannerRAS(lsrc, mri, ldst);  // make sure it starts as tkreg ras

  M_surface_to_vox = voxelFromSurfaceRAS_(mri);

  v1 = VectorAlloc(4, MATRIX_REAL);
  v2 = VectorAlloc(4, MATRIX_REAL);
  VECTOR_ELT(v1, 4) = 1.0;
  VECTOR_ELT(v2, 4) = 1.0;
  for (i = 0; i < lsrc->n_points; i++) {
    V3_X(v1) = ldst->lv[i].x;
    V3_Y(v1) = ldst->lv[i].y;
    V3_Z(v1) = ldst->lv[i].z;
    MatrixMultiply(M_surface_to_vox, v1, v2);
    ldst->lv[i].x = V3_X(v2);
    ldst->lv[i].y = V3_Y(v2);
    ldst->lv[i].z = V3_Z(v2);
    ldst->lv[i].xv = nint(V3_X(v2));
    ldst->lv[i].yv = nint(V3_Y(v2));
    ldst->lv[i].zv = nint(V3_Z(v2));
  }
  strncpy(ldst->space, "voxel", sizeof(ldst->space));
  ldst->coords = LABEL_COORDS_VOXEL;
  VectorFree(&v1);
  VectorFree(&v2);
  MatrixFree(&M_surface_to_vox);
  return (ldst);
}
LABEL *LabelClone(LABEL *a)
{
  LABEL *l;
  l = LabelAlloc(a->max_points, a->subject_name, a->name);
  l->mri_template = a->mri_template;
  l->coords = a->coords;
  l->mht = a->mht;
  l->mris = a->mris;
  l->avg_stat = a->avg_stat;
  strcpy(l->space, a->space);
  return (l);
}

LABEL *LabelTransform(LABEL *lsrc, TRANSFORM *xform, MRI *mri, LABEL *ldst)
{
  int i;
  MATRIX *M;
  VECTOR *v1, *v2;

  if (ldst != lsrc)
    ldst = LabelCopy(lsrc, ldst) ;

  if (xform->type != LINEAR_RAS_TO_RAS)
    ErrorExit(ERROR_NOFILE, "LabelTransform: unsupported type %d. Must be RAS->RAS", xform->type);
  if (strstr(ldst->space, "scanner") == NULL)
    ErrorExit(ERROR_NOFILE, "LabelTransform: label must be in scanner RAS not %s", ldst->space);

  M = ((LTA *)(xform->xform))->xforms[0].m_L;

  v1 = VectorAlloc(4, MATRIX_REAL);
  v2 = VectorAlloc(4, MATRIX_REAL);
  VECTOR_ELT(v1, 4) = 1.0;
  VECTOR_ELT(v2, 4) = 1.0;
  for (i = 0; i < lsrc->n_points; i++) {
    V3_X(v1) = ldst->lv[i].x;
    V3_Y(v1) = ldst->lv[i].y;
    V3_Z(v1) = ldst->lv[i].z;
    MatrixMultiply(M, v1, v2);
    ldst->lv[i].x = V3_X(v2);
    ldst->lv[i].y = V3_Y(v2);
    ldst->lv[i].z = V3_Z(v2);
  }
  strncpy(ldst->space, "scanner", sizeof(ldst->space));
  ldst->coords = LABEL_COORDS_SCANNER_RAS ;

  return (ldst);
}

LABEL *LabelVoxelToSurfaceRAS(LABEL *lsrc, MRI *mri, LABEL *ldst)
{
  int i;
  double xs, ys, zs;

  if (ldst != lsrc)
    ldst = LabelCopy(lsrc, ldst) ;

  for (i = 0; i < lsrc->n_points; i++) {
    if (lsrc->mris)
      MRISsurfaceRASToVoxel((MRIS *)lsrc->mris, mri, lsrc->lv[i].x, lsrc->lv[i].y, lsrc->lv[i].z, &xs, &ys, &zs);
    else
      MRIvoxelToSurfaceRAS(mri, lsrc->lv[i].x, lsrc->lv[i].y, lsrc->lv[i].z, &xs, &ys, &zs);
    ldst->lv[i].x = xs;
    ldst->lv[i].y = ys;
    ldst->lv[i].z = zs;
  }

  ldst->coords = LABEL_COORDS_TKREG_RAS;
  strcpy(ldst->space, "TkReg");
  return (ldst);
}

/*!
  \fn LABEL *LabelBaryFill(MRIS *mris, LABEL *srclabel, double delta)
  \brief Fills a surface-based label to points between the vertices
  using barycentric coordinates.  The resulting label is not
  surface-based since there will be points off the mesh. delta
  controls how finely the face is filled. delta should be between 0
  and 1.
 */
LABEL *LabelBaryFill(MRIS *mris, LABEL *srclabel, double delta)
{
  MRI *mri;
  LABEL *outlabel;
  int k, ok, n, vtxno, nPerTriangle, nTriangles, nthlab;
  VERTEX *v1, *v2, *v3;
  double x, y, z, l1, l2, l3;

  mri = MRISlabel2Mask(mris, srclabel, NULL);

  nPerTriangle = 0;
  for (l1 = 0; l1 <= 1; l1 += .1) {
    for (l2 = 0; l2 <= 1; l2 += .1) {
      if (l1 + l2 > 1) continue;
      nPerTriangle++;
    }
  }

  nTriangles = 0;
  for (k = 0; k < mris->nfaces; k++) {
    ok = 0;
    for (n = 0; n < 3; n++) {
      vtxno = mris->faces[k].v[n];
      if (MRIgetVoxVal(mri, vtxno, 0, 0, 0) < 0.5) continue;
      ok = 1;
      break;
    }
    if (!ok) continue;
    nTriangles++;
  }
  // printf("nPerTriangle = %d, nTriangles = %d, %d\n",nPerTriangle,nTriangles,nPerTriangle*nTriangles);

  outlabel = LabelAlloc(nPerTriangle * nTriangles, srclabel->subject_name, NULL);
  outlabel->n_points = nPerTriangle * nTriangles;
  nthlab = 0;
  for (k = 0; k < mris->nfaces; k++) {
    ok = 0;
    for (n = 0; n < 3; n++) {
      vtxno = mris->faces[k].v[n];
      if (MRIgetVoxVal(mri, vtxno, 0, 0, 0) < 0.5) continue;
      ok = 1;
      break;
    }
    if (!ok) continue;
    v1 = &(mris->vertices[mris->faces[k].v[0]]);
    v2 = &(mris->vertices[mris->faces[k].v[1]]);
    v3 = &(mris->vertices[mris->faces[k].v[2]]);
    for (l1 = 0; l1 <= 1; l1 += .1) {
      for (l2 = 0; l2 <= 1; l2 += .1) {
        if (l1 + l2 > 1) continue;
        l3 = 1 - (l1 + l2);
        x = l1 * v1->x + l2 * v2->x + l3 * v3->x;
        y = l1 * v1->y + l2 * v2->y + l3 * v3->y;
        z = l1 * v1->z + l2 * v2->z + l3 * v3->z;
        outlabel->lv[nthlab].x = x;
        outlabel->lv[nthlab].y = y;
        outlabel->lv[nthlab].z = z;
        nthlab++;
      }
    }
  }

  MRIfree(&mri);

  return (outlabel);
}
LABEL *LabelSampleToSurface(MRI_SURFACE *mris, LABEL *area, MRI *mri_template, int coords)
{
  int n, i, vno, min_vno, nfilled = 0, max_vno;
  LV *lv;
  static MHT *mht = NULL;
  static MRI_SURFACE *mris_cached = NULL;
  float dx, dy, dz, x, y, z, dist, min_dist;
  int num_not_found, nchanged, num_brute_force;
  float vx, vy, vz;
  double max_spacing;
  static MRI *mri = NULL, *mri_cached = NULL;
  LABEL *area_dst;
  double xw, yw, zw, xv, yv, zv;

  LabelInit(area, mri_template, mris, coords) ;
  printf("LabelSampleToSurface(%d vertices)\n", area->n_points);
  xw = yw = zw = vx = vy = vz = x = y = z = -1;  // remove compiler warnings

  for (i = n = 0; n < area->n_points; n++) {
    lv = &area->lv[n];
    if (lv->vno >= 0 && lv->vno < mris->nvertices)  // already assigned
      continue;

    i++; /* count # of unassigned vertices */
  }
  area_dst = LabelAlloc(mris->nvertices, area->subject_name, area->name);
  LabelCopy(area, area_dst);

  if (i <= 0) {
    printf("LabelSampleToSurface: no hole filling needed, returning (%d vertices)\n", area_dst->n_points);
    return (area_dst); /* no work needed */
  }

  fprintf(stderr, "%d unassigned vertices in label...\n", i);

  /* if we can't find a vertex within 10 mm of the point, something is wrong */
  if (mris != mris_cached) {
    fprintf(stderr, "building spatial LUT...\n");
    MRIScomputeVertexSpacingStats(mris, NULL, NULL, &max_spacing, NULL, &max_vno, coords);
    if (mht) MHTfree(&mht);
    mris_cached = mris;
    mht = MHTcreateVertexTable_Resolution(mris, coords, 2 * max_spacing);
  }
  fprintf(stderr, "assigning vertex numbers to label...\n");
  num_not_found = num_brute_force = 0;
  for (n = 0; n < area_dst->n_points; n++) {
    lv = &area_dst->lv[n];
    if (lv->vno >= 0 && lv->vno <= mris->nvertices)  // vertex already assigned
      continue;

    nfilled++;
    labelGetSurfaceRasCoords(area, lv, &x, &y, &z);

    {
      min_dist = 10000.0;
      min_vno  = -1;

      MHBT *bucket = MHTacqBucket(mht, x, y, z);  // x,y,z now in surface ras
      if (bucket) {
        MHB *bin;
        for (bin = bucket->bins, i = 0; i < bucket->nused; i++, bin++)  // find min dist vertex
        {
          vno = bin->fno;
          VERTEX const * const v = &mris->vertices[vno];
          if (vno == Gdiag_no) DiagBreak();

          MRISgetCoords(v, coords, &vx, &vy, &vz);

          dx = vx - x;
          dy = vy - y;
          dz = vz - z;

          dist = sqrt(dx * dx + dy * dy + dz * dz);
          if (dist < min_dist) {
            min_dist = dist;
            min_vno  = vno;
          }
        }
        MHTrelBucket(&bucket);
      }
    }

    if (min_vno == 0) DiagBreak();
    if (min_vno == -1) {
      num_brute_force++;
      switch (coords) {
      case ORIG_VERTICES:
        min_vno = MRISfindClosestOriginalVertex(mris, lv->x, lv->y, lv->z);
        break;
      case WHITE_VERTICES:
        min_vno = MRISfindClosestWhiteVertex(mris, lv->x, lv->y, lv->z);
        break;
      case CURRENT_VERTICES:
        min_vno = MRISfindClosestVertex(mris, lv->x, lv->y, lv->z, NULL, CURRENT_VERTICES);
        break;
      case CANONICAL_VERTICES:
        min_vno = MRISfindClosestCanonicalVertex(mris, lv->x, lv->y, lv->z);
        break;
      default:
        break;
      }
    }
    if (min_vno == -1) {
      num_not_found++;
    }
    if (min_vno == Gdiag_no) DiagBreak();
    lv->vno = min_vno;
    area_dst->vertex_label_ind[min_vno] = n;
  }
  //  LabelRemoveDuplicates(area) ;
  //  MHTfree(&mht) ;

  if (num_not_found > 0) fprintf(stderr, "Couldn't assign %d vertices.\n", num_not_found);
  if (num_brute_force > 0) fprintf(stderr, "%d vertices required brute force search.\n", num_brute_force);

  // now sample from surface back into volume to see if there are any holes
  if (mri_template != mri_cached) {
    if (mri) MRIfree(&mri);
    mri = MRIclone(mri_template, NULL);
    mri_cached = mri_template;
  }
  for (i = 0; i < area_dst->n_points; i++) {
    labelGetVoxelCoords(area_dst, &area_dst->lv[i], &vx, &vy, &vz);
    if (nint(vx) == Gx && nint(vy) == Gy && nint(vz) == Gz) DiagBreak();
    MRIsetVoxVal(mri, nint(vx), nint(vy), nint(vz), 0, 1);
  }

  //  MRISclearMarks(mris) ;

  do {
    int nbr, vno2;
    VERTEX *vn;

    nchanged = 0;
    LabelMarkSurface(area_dst, mris);
    for (i = 0; i < area_dst->n_points; i++) {
      vno = area_dst->lv[i].vno;
      if (vno == Gdiag_no) DiagBreak();
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX          const * const v  = &mris->vertices         [vno];
      for (nbr = 0; nbr < vt->vnum; nbr++) {
        vno2 = vt->v[nbr];
        if (vno2 == Gdiag_no) DiagBreak();
        vn = &mris->vertices[vno2];
        if (vn->marked != 0)  // already in the label
          continue;

        MRISgetCoords(v, coords, &xw, &yw, &zw);
        MRISsurfaceRASToVoxel(mris, mri, xw, yw, zw, &xv, &yv, &zv);
        if (MRIgetVoxVal(mri, nint(xv), nint(yv), nint(zv), 0) > 0) {
          LABEL_VERTEX *lv;

          printf("adding vertex %d at index %d to fill label hole\n", vno2, area_dst->n_points);
          lv = &area_dst->lv[area_dst->n_points];
          lv->vno = vno2;
          lv->x = xw;
          lv->y = yw;
          lv->z = zw;
          lv->stat = vn->val;
          vn->marked = 1;
          area_dst->vertex_label_ind[vno2] = area_dst->n_points;
          nchanged++;
          area_dst->n_points++;
        }
      }
    }
    LabelUnmark(area_dst, mris);
  } while (nchanged > 0);
  for (i = 0; i < area_dst->n_points; i++) {
    labelGetVoxelCoords(area_dst, area_dst->lv + i, &vx, &vy, &vz);
    if (nint(vx) == Gx && nint(vy) == Gy && nint(vz) == Gz) DiagBreak();
    MRIsetVoxVal(mri, nint(vx), nint(vy), nint(vz), 0, 0);  // turn them all off
  }
  //  MRIfree(&mri) ;
  printf("LabelSampleToSurface: returning (%d vertices)\n", area_dst->n_points);
  return (area_dst);
}

int LabelInit(LABEL *area, MRI *mri_template, MRI_SURFACE *mris, int coords)
{
  double xv, yv, zv, x, y, z, min_dist, vx, vy, vz;
  int n, min_vno, i, vno;
  LV *lv;
  VERTEX *v;

  area->mri_template = mri_template;
  if (mris == NULL)
  {
    if (area->coords != LABEL_COORDS_SCANNER_RAS)
      LabelToScannerRAS(area, mri_template, area);

    for (n = 0; n < area->n_points; n++) {
      lv = &area->lv[n];

      if (lv->deleted) continue;

      // always in scanner ras now
      MRIscannerRASToVoxel(mri_template, lv->x, lv->y, lv->z, &xv, &yv, &zv);
      lv->xv = nint(xv);  lv->yv = nint(yv); lv->zv = nint(zv);
    }  // for loop

    return (NO_ERROR);   // end case where no mris is passed
  }  // mris == NULL

  x = y = z = -1;
  LabelRealloc(area, mris->nvertices);  // allocate enough room in the label for the whole surface
  area->mris = mris;

  // create a volume of indices into the label. Voxels < 0 are not mapped
  area->vertex_label_ind = (int *)calloc(mris->nvertices, sizeof(int));
  if (area->vertex_label_ind == NULL)
    ErrorExit(ERROR_NOMEMORY, "LabelInit: could not allocate %d-long vertex index array", mris->nvertices);
  for (n = 0; n < mris->nvertices; n++) area->vertex_label_ind[n] = -1;  // means that this vertex is not in th elabel

  double vxl_spacing = (mri_template->xsize < mri_template->ysize ? mri_template->xsize : mri_template->ysize);
  vxl_spacing = (vxl_spacing < mri_template->zsize ? vxl_spacing : mri_template->zsize);
  area->mht = MHTcreateVertexTable_Resolution(mris, coords, vxl_spacing);

  // map unassigned vertices to surface locations
  for (n = 0; n < area->n_points; n++) {
    lv = &area->lv[n];

    if (lv->deleted) continue;
    if (lv->vno >= 0 && lv->vno <= mris->nvertices)  // vertex already assigned
    {
      area->vertex_label_ind[lv->vno] = n;
      continue;
    }

    // vertex not assigned to this label vertex - assign all surface vertices at this location to the label
    // note that this includes the closest vertex, even if it does not sit in the same voxel

    vx = lv->x;
    vy = lv->y;
    vz = lv->z;

    switch (area->coords) {
    case LABEL_COORDS_SCANNER_RAS:
      MRIworldToVoxel(mri_template, vx, vy, vz, &xv, &yv, &zv);
      break;
    case LABEL_COORDS_TKREG_RAS:
      MRISsurfaceRASToVoxel(mris, mri_template, vx, vy, vz, &xv, &yv, &zv);
      break;
    case LABEL_COORDS_VOXEL:
      xv = lv->x;
      yv = lv->y;
      zv = lv->z;
      break;
    default:
      ErrorExit(ERROR_UNSUPPORTED, "LabelInit: area->coords %d unsupported", area->coords);
      break;
    }
    lv->xv = nint(xv);
    lv->yv = nint(yv);
    lv->zv = nint(zv);

    if (area->mris)
      MRISsurfaceRASFromVoxel((MRIS *)area->mris, area->mri_template, xv, yv, zv, &vx, &vy, &vz);  // surfaceRASToVoxel BRF!!
    else
      MRIvoxelToSurfaceRAS(area->mri_template, xv, yv, zv, &vx, &vy, &vz);

    // must use surface coords for finding vertex
    if (area->mris)
      MRISsurfaceRASFromVoxel((MRIS *)area->mris, area->mri_template, xv, yv, zv, &vx, &vy, &vz);  // surfaceRASToVoxel BRF!!
    else
      MRIvoxelToSurfaceRAS(area->mri_template, xv, yv, zv, &vx, &vy, &vz);
    x = y = z = -1;

    //  printf("LabelAddVoxel(%d, %d, %d): coords = %2.1f, %2.1f, %2.1f\n", xv, yv, zv, vx, vy, vz) ;

    min_dist = 10000.0;
    min_vno = -1;

    float distance;
    min_vno = MHTfindClosestVertexNoXYZ( (MHT *)area->mht, (MRIS *)area->mris, vx,vy,vz, &distance );

    lv->vno = min_vno;
    if (min_vno >= 0 && area->vertex_label_ind[min_vno] < 0)  // found one that isn't in label
    {
      area->vertex_label_ind[min_vno] = n;
      //      printf("LabelAddVoxel(%d, %d, %d): added min_dist vno %d at %d\n", xv, yv, zv, min_vno, n) ;

      // now add other vertices that also map to this voxel
      VERTEX_TOPOLOGY const * const min_vt = &((MRI_SURFACE *)((MRIS *)area->mris))->vertices_topology[min_vno];
      for (i = 0; i < min_vt->vnum; i++)  // find min dist vertex
      {
        vno = min_vt->v[i];
        if (area->vertex_label_ind[vno] >= 0) continue;  // already in the label
        v = &((MRI_SURFACE *)((MRIS *)area->mris))->vertices[vno];
        if (vno == Gdiag_no) DiagBreak();

        MRISgetCoords(v, coords, &x, &y, &z);
        MRISsurfaceRASToVoxel((MRIS *)area->mris, area->mri_template, x, y, z, &vx, &vy, &vz);
        if ((xv == nint(vx)) && (yv == nint(vy)) && (zv == nint(vz))) {
          LabelAddVertex(area, vno, coords);
        }
      }
    }
  }
  return (NO_ERROR);
}

int LabelAddVoxel(LABEL *area, int xv, int yv, int zv, int coords, int *vertices, int *pnvertices)
{
  int n, min_vno, i, vno;
  LV *lv;
  double min_dist, x, y, z, vx = 0, vy = 0, vz = 0;
  VERTEX *v;

  x = y = z = 0;  // quick bug fix. should review whether code is correct

  if (pnvertices) *pnvertices = 0;

  if (area->n_points >= area->max_points) LabelRealloc(area, nint(1.5 * area->n_points));
  n = area->n_points++;
  lv = &area->lv[n];
  lv->xv = xv;
  lv->yv = yv;
  lv->zv = zv;
  lv->vno = -1;
  if (area->mri_template) {
    switch (area->coords) {
    case LABEL_COORDS_TKREG_RAS:
      if (area->mris)
        MRISsurfaceRASFromVoxel((MRIS *)area->mris, area->mri_template, xv, yv, zv, &vx, &vy, &vz);  // surfaceRASToVoxel BRF!!
      else
        MRIvoxelToSurfaceRAS(area->mri_template, xv, yv, zv, &vx, &vy, &vz);
      break;
    case LABEL_COORDS_SCANNER_RAS:
      MRIvoxelToWorld(area->mri_template, xv, yv, zv, &vx, &vy, &vz);
      break;
    default:
      ErrorExit(ERROR_UNSUPPORTED, "LabelAddVoxel: coords %d not supported yet", area->coords);
    }
    lv->x = vx;
    lv->y = vy;
    lv->z = vz;
  }

  if (area->mht == NULL) return (NO_ERROR);

  // must use surface coords for finding vertex
  if (area->mris)
    MRISsurfaceRASFromVoxel((MRIS *)area->mris, area->mri_template, xv, yv, zv, &vx, &vy, &vz);  // surfaceRASToVoxel BRF!!
  else
    MRIvoxelToSurfaceRAS(area->mri_template, xv, yv, zv, &vx, &vy, &vz);
  x = y = z = -1;

  //  printf("LabelAddVoxel(%d, %d, %d): coords = %2.1f, %2.1f, %2.1f\n", xv, yv, zv, vx, vy, vz) ;

  min_dist = 10000.0;
  min_vno = -1;

  float distance;
  min_vno = MHTfindClosestVertexNoXYZ( (MHT *)area->mht, (MRIS *)area->mris, vx,vy,vz, &distance );

  lv->vno = min_vno;
  if (min_vno >= 0)
  {
    if (area->vertex_label_ind[min_vno] < 0)  // found one that isn't in label
    {
      area->vertex_label_ind[min_vno] = n;
      if (pnvertices) {
        int n = *pnvertices;
        vertices[n] = min_vno;
        (*pnvertices)++;
      }
      //      printf("LabelAddVoxel(%d, %d, %d): added min_dist vno %d at %d\n", xv, yv, zv, min_vno, n) ;
    }
    else  //
    {
      n = area->vertex_label_ind[min_vno];
      LV *lv2 = &area->lv[n];
      if (FEQUAL(lv->x, lv2->x) && FEQUAL(lv->y, lv2->y) && FEQUAL(lv->z, lv2->z))
      {
        lv->deleted = 1;
      }
      else
      {
        lv->vno = -1;
        for ( i = 0; i < area->n_points-1; i++)
        {
          LV *lv2 = &area->lv[i];
          if ( lv2->vno < 0 && !lv2->deleted)
          {
            if (FEQUAL(lv->x, lv2->x) && FEQUAL(lv->y, lv2->y) && FEQUAL(lv->z, lv2->z))
            {
              lv->deleted = 1;
              break;
            }
          }
        }
      }
    }
  }
  else
  {
    for ( i = 0; i < area->n_points-1; i++)
    {
      LV *lv2 = &area->lv[i];
      if ( lv2->vno < 0 && !lv2->deleted)
      {
        if (FEQUAL(lv->x, lv2->x) && FEQUAL(lv->y, lv2->y) && FEQUAL(lv->z, lv2->z))
        {
          lv->deleted = 1;
          break;
        }
      }
    }
    return (NO_ERROR);
  }
  //  else
  //    printf("min dist vertex %d already in label\n", min_vno) ;

  // now add other vertices that also map to this voxel
  VERTEX_TOPOLOGY const * min_vt = &((MRI_SURFACE *)(area->mris))->vertices_topology[min_vno];
  for (i = 0; i < min_vt->vnum; i++)  // find min dist vertex
  {
    vno = min_vt->v[i];
    if (area->vertex_label_ind[vno] >= 0) continue;  // already in the label
    v = &((MRI_SURFACE *)(area->mris))->vertices[vno];
    if (vno == Gdiag_no) DiagBreak();

    MRISgetCoords(v, coords, &x, &y, &z);
    MRISsurfaceRASToVoxel((MRIS *)area->mris, area->mri_template, x, y, z, &vx, &vy, &vz);
    if ((xv == nint(vx)) && (yv == nint(vy)) && (zv == nint(vz))) {
      if (pnvertices) {
        int n = *pnvertices;
        vertices[n] = vno;
        (*pnvertices)++;
      }
      LabelAddVertex(area, vno, coords);
    }
  }

  return (NO_ERROR);
}

int LabelDeleteVoxel(LABEL *area, int xv, int yv, int zv, int *vertices, int *pnvertices)
{
  int n, ndeleted;
  LV *lv;
#if 0
  MATRIX *m_vox2ras ;
  VECTOR *v1, *v2 ;

  if (area->coords == LABEL_COORDS_SCANNER_RAS)
  {
    m_vox2ras = MRIxfmCRS2XYZ( area->mri_template, 0 );  // Native Vox2RAS Matrix
    v1 = VectorAlloc(4, MATRIX_REAL);
    v1->rptr[4][1] = 1.0f ;
    V3_X(v1) = xv ; V3_Y(v1) = yv ; V3_Z(v1) = zv ;
    v2 = MatrixMultiply(m_vox2ras, v1, NULL) ;
    xv = nint(V3_X(v2)) ; yv = nint(V3_Y(v2)) ; zv = nint(V3_Z(v2)) ;
  }
  else if (area->coords == LABEL_COORDS_TKREG_RAS)
    ErrorExit(ERROR_UNSUPPORTED, "LabelDeleteVoxel: label coords tkreg unsupported\n") ;
#endif

  for (ndeleted = n = 0; n < area->n_points; n++) {
    lv = &area->lv[n];
    if (lv->deleted) continue;
    if (lv->xv == xv && lv->yv == yv && lv->zv == zv) {
      lv->deleted = 1;
      if (lv->vno >= 0) {
        if (vertices) vertices[ndeleted] = lv->vno;
        ndeleted++;
        if (area->vertex_label_ind) area->vertex_label_ind[lv->vno] = -1;
      }
    }
  }

  //  MatrixFree(&m_vox2ras) ; MatrixFree(&v1) ; MatrixFree(&v2) ;
  if (pnvertices) *pnvertices = ndeleted;
  //  printf("LabelDeleteVoxel(%d, %d, %d): %d deleted\n", xv, yv, zv, ndeleted) ;
  return (ndeleted);
}

/*!
  \fn LABEL *LabelAddPoint(LABEL *label, LV *srclv)
  \brief Adds a label point (ie, label vertex, not the same as a
  surface vertex) to a label (creates the label if label==NULL). This
  will reallocate the label if necessary. Does not check whether the
  label point is there or not. It does not distinguish between surface
  and volume label points. All it does is to add the passed label
  point (srclv) to the end of the list of label points, copying over
  the entire LV structure using memcpy(), and increments
  n_points. This differs from LabelAddVertex(). LabelAddVertex()
  specifically adds a surface vertex to the label computing the xyz
  from the surface coords and updates vertex_label_ind; this function
  does neither of those.
*/
LABEL *LabelAddPoint(LABEL *label, LV *lv)
{
  // alloc label if nec
  if (label == NULL) label = LabelAlloc(100, NULL, NULL);
  // realloc label if nec
  if (label->n_points >= label->max_points) LabelRealloc(label, nint(label->max_points * 1.5));
  // copy structure
  memcpy(&(label->lv[label->n_points]), lv, sizeof(LV));
  // increment the number of points
  label->n_points++;
  return (label);
}
/*!
  \fn int LabelAddVertex(LABEL *area, int vno, int coords)
  Adds a surface vertex to a label (note: a surface vertex is
  different from a label vertex). Reallocs label if needed.  Sets the
  coordinates based on the coordinates of the given vertex and the
  coords argument. The voxel coordinates are computed. If the label
  already exists as indicated by area->vertex_label_ind[vno], then
  nothing is done and -1 is returned. See also LabelDeleteVertex() and
  LabelAddPoint().
 */
int LabelAddVertex(LABEL *area, int vno, int coords)
{
  LV *lv;
  int n;
  double xv, yv, zv, x, y, z;
  VERTEX *v;

  x = y = z = -1;

  if (area->vertex_label_ind[vno] >= 0) return (-1);  // already in the label
  if (area->n_points >= area->max_points) LabelRealloc(area, nint(1.5 * area->n_points));

  n = area->n_points++;  // n is the number of points before incr
  lv = &area->lv[n];
  v = &((MRI_SURFACE *)(area->mris))->vertices[vno];
  MRISgetCoords(v, coords, &x, &y, &z);
  MRISsurfaceRASToVoxel((MRIS *)area->mris, area->mri_template, v->x, v->y, v->z, &xv, &yv, &zv);
  if (area->mri_template && area->coords == LABEL_COORDS_SCANNER_RAS)
    MRIvoxelToWorld(area->mri_template, xv, yv, zv, &x, &y, &z);

  lv->vno = vno;
  lv->x = x;
  lv->y = y;
  lv->z = z;
  lv->xv = nint(xv);
  lv->yv = nint(yv);
  lv->zv = nint(zv);

  area->vertex_label_ind[vno] = n;
  return (NO_ERROR);
}

int LabelDeleteVertex(LABEL *area, int vno, int coords)
{
  int n;

  if (area->vertex_label_ind[vno] < 0) return (-1);
  n = area->vertex_label_ind[vno];
  area->vertex_label_ind[vno] = -1;
  area->lv[n].deleted = 1;
  return (NO_ERROR);
}
static int update_vertex_indices(LABEL *area)
{
  int vno, n;
  MRI_SURFACE *mris = (MRI_SURFACE *)(area->mris);

  if (area->vertex_label_ind == NULL || mris == NULL) return (NO_ERROR);

  for (vno = 0; vno < mris->nvertices; vno++) area->vertex_label_ind[vno] = -1;

  for (n = 0; n < area->n_points; n++)
    if (area->lv[n].deleted == 0) area->vertex_label_ind[area->lv[n].vno] = n;
  return (NO_ERROR);
}
LABEL *LabelApplyMatrix(LABEL *lsrc, MATRIX *m, LABEL *ldst)
{
  int n;
  VECTOR *v1, *v2;

  v1 = VectorAlloc(4, MATRIX_REAL);
  v2 = VectorAlloc(4, MATRIX_REAL);
  v1->rptr[4][1] = 1.0f;
  v2->rptr[4][1] = 1.0f;

  if (lsrc != ldst)
    ldst = LabelCopy(lsrc, ldst) ;

  for (n = 0; n < lsrc->n_points; n++) {
    V3_X(v1) = ldst->lv[n].x;
    V3_Y(v1) = ldst->lv[n].y;
    V3_Z(v1) = ldst->lv[n].z;
    MatrixMultiply(m, v1, v2);
    ldst->lv[n].x = V3_X(v2);
    ldst->lv[n].y = V3_Y(v2);
    ldst->lv[n].z = V3_Z(v2);
  }

  VectorFree(&v1);
  VectorFree(&v2);
  return (ldst);
}
static int labelGetSurfaceRasCoords(LABEL *area, LABEL_VERTEX *lv, float *px, float *py, float *pz)
{
  double xv, yv, zv, xw, yw, zw;

  switch (area->coords) {
  case LABEL_COORDS_TKREG_RAS:
    *px = lv->x;
    *py = lv->y;
    *pz = lv->z;
    break;
  case LABEL_COORDS_SCANNER_RAS:
    MRIworldToVoxel(area->mri_template, lv->x, lv->y, lv->z, &xv, &yv, &zv);
    MRISsurfaceRASFromVoxel((MRIS *)area->mris, area->mri_template, xv, yv, zv, &xw, &yw, &zw);
    *px = (float)xw;
    *py = (float)yw;
    *pz = (float)zw;  // double->float (uggh)
    break;
  default:
    ErrorExit(ERROR_UNSUPPORTED, "labelGetSurfaceRAScoords: coords %d not supported yet", area->coords);
    break;
  }
  return (NO_ERROR);
}

static int labelGetVoxelCoords(LABEL *area, LABEL_VERTEX *lv, float *px, float *py, float *pz)
{
  double xv = 0, yv = 0, zv = 0;

  switch (area->coords) {
  case LABEL_COORDS_TKREG_RAS:
    if (area->mris)
      MRISsurfaceRASToVoxel((MRIS *)area->mris, area->mri_template, lv->x, lv->y, lv->z, &xv, &yv, &zv);
    else
      MRIvoxelToSurfaceRAS(area->mri_template, lv->x, lv->y, lv->z, &xv, &yv, &zv);
    break;
  case LABEL_COORDS_SCANNER_RAS:
    MRIworldToVoxel(area->mri_template, lv->x, lv->y, lv->z, &xv, &yv, &zv);
    break;
  default:
    ErrorExit(ERROR_UNSUPPORTED, "labelGetSurfaceRAScoords: coords %d not supported yet", area->coords);
    break;
  }
  *px = (float)xv;
  *py = (float)yv;
  *pz = (float)zv;
  return (NO_ERROR);
}

double
LabelAverageVal(LABEL *area, MRI_SURFACE *mris)
{
  int vno, n, num ;
  VERTEX *v ;
  double avg ;

  for (avg = 0.0, num = n = 0; n < area->n_points; n++) {
    if (area->lv[n].deleted)
      continue ;
    vno = area->lv[n].vno ;
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    num++ ;
    avg += v->val ;
  }
  if (num > 0)
    avg /= num ;
  return(avg) ;
}
float
LabelMaxStat(LABEL *area) 
{
  int     n ;
  float   max_stat ;

  for (max_stat = -1e10, n = 0; n < area->n_points; n++) 
  {
    if (area->lv[n].deleted)
      continue ;
    if (area->lv[n].stat > max_stat)
      max_stat = area->lv[n].stat;
  }
  return(max_stat) ;
}
LABEL  *
LabelFromSurface(MRI_SURFACE *mris, int which, double thresh) 
{
  int    vno, npoints ;
  VERTEX *v ;
  LABEL  *area ;
  double val ;


  for (npoints = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;

    switch (which)
    {
    case VERTEX_VALS:
      val = v->val ;
      break;
    case VERTEX_CURV:
      val = v->curv ;
      break;
    case VERTEX_STATS:
      val = v->stat ;
      break;
    default:
      ErrorExit(ERROR_UNSUPPORTED, "LabelFromSurface: unsupported which %d", which);
      val = 0 ;
      break ;
    }
    if (val >= thresh)
      npoints++ ;
  }

  area = LabelAlloc(npoints, NULL, NULL) ;
  for (npoints = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;

    switch (which)
    {
    case VERTEX_VALS:
      val = v->val ;
      break;
    case VERTEX_CURV:
      val = v->curv ;
      break;
    case VERTEX_STATS:
      val = v->stat ;
      break;
    default:
      ErrorExit(ERROR_UNSUPPORTED, "LabelFromSurface: unsupported which %d", which);
      val = 0 ;
      break ;
    }
    if (val >= thresh)
    {
      LABEL_VERTEX *lv ;

      lv = &area->lv[npoints] ;
      area->n_points++ ;
      npoints++ ;
      lv->stat = val ;
      lv->vno = vno ;
      lv->x = v->x ;
      lv->y = v->y ;
      lv->z = v->z ;
    }
  }

  return(area) ;
}

LABEL *LabelRemoveHolesSurf(MRIS *surf, LABEL *lb)
{
  LABEL *lbinv, *lbinvNoIslands, *lbNoHoles;
  int n;

  // Invert the label
  lbinv = MRISlabelInvert(surf, lb);

  // Remove the island from the inverted label
  lbinvNoIslands = LabelRemoveIslandsSurf(surf, lbinv);

  // Invert it back
  lbNoHoles = MRISlabelInvert(surf, lbinvNoIslands);

  LabelFree(&lbinv);
  LabelFree(&lbinvNoIslands);

  // Copy the original stat back (not possible in a whol)
  double *stat = (double *)calloc(sizeof(double),surf->nvertices);
  for(n=0; n < lb->n_points; n++)
    stat[lb->lv[n].vno] = lb->lv[n].stat;
  for(n=0; n < lbNoHoles->n_points; n++)
    lbNoHoles->lv[n].stat = stat[lbNoHoles->lv[n].vno];
  free(stat);
  
  return(lbNoHoles);
}


LABEL *LabelRemoveIslandsSurf(MRIS *surf, LABEL *lb)
{
  int n,vtxno,NClusters,MaxSize,nMax;
  SURFCLUSTERSUM *scs;
  LABEL *lbcluster=NULL;

  // Set val to 0 in al vertices
  for (vtxno = 0; vtxno < surf->nvertices; vtxno++)
    surf->vertices[vtxno].val = 0;

  for(n=0; n < lb->n_points; n++)
    surf->vertices[lb->lv[n].vno].val = 1;

  scs = sclustMapSurfClusters(surf,0.5,2,+1,0,&NClusters,NULL,NULL);
  printf("-------------------------------------\n");
  DumpSurfClusterSum(stdout, scs, NClusters);
  printf("-------------------------------------\n");

  MaxSize = 0;
  nMax = 0;
  for(n=0; n < NClusters; n++){
    printf("%d \n",scs[n].nmembers);
    if(MaxSize < scs[n].nmembers){
      MaxSize = scs[n].nmembers;
      nMax = n;
    }
  }

  printf("n = %d %d\n",nMax,MaxSize);
  lbcluster = LabelAlloc(scs[nMax].nmembers+1, lb->subject_name, NULL);
  lbcluster->coords = lb->coords;
  lbcluster->n_points = scs[nMax].nmembers;
  strcpy(lbcluster->space,lb->space);
  n = 0;
  for (vtxno = 0; vtxno < surf->nvertices; vtxno++){
    if(surf->vertices[vtxno].undefval != nMax+1) continue;
    //printf("%d %d  %d\n",vtxno,surf->vertices[vtxno].undefval,n);
    lbcluster->lv[n].vno = vtxno;
    lbcluster->lv[n].x = surf->vertices[vtxno].x;
    lbcluster->lv[n].y = surf->vertices[vtxno].y;
    lbcluster->lv[n].z = surf->vertices[vtxno].z;
    lbcluster->lv[n].stat = surf->vertices[vtxno].stat;
    n++;
  }
  free(scs);

  return(lbcluster);
}

LABEL *LabelRemoveHolesAndIslandsSurf(MRIS *surf, LABEL *lb)
{
  int n;
  printf("Removing label holes\n");
  LABEL *tmplabel = LabelRemoveHolesSurf(surf, lb);
  printf("Removing label islands\n");
  LABEL *tmplabel2 = LabelRemoveIslandsSurf(surf, tmplabel);
  LabelFree(&tmplabel);

  // Copy the original stat back (not possible in a whol)
  double *stat = (double *)calloc(sizeof(double),surf->nvertices);
  for(n=0; n < lb->n_points; n++)
    stat[lb->lv[n].vno] = lb->lv[n].stat;
  for(n=0; n < tmplabel2->n_points; n++)
    tmplabel2->lv[n].stat = stat[tmplabel2->lv[n].vno];

  free(stat);
  return(tmplabel2);
}
