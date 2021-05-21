/**
 * @brief utilities for surface-based parcellations
 *
 * utilities for surface-based parcellations (see Fischl et al.,
 * Cerebral Cortex)
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <map>
#include <vector>

#include "colortab.h"
#include "const.h"
#include "diag.h"
#include "error.h"
#include "label.h"
#include "mri2.h"
#include "log.h"
#include "fio.h"
#include "mrisurf.h"
#include "utils.h"
#include "tags.h"

#define ANNOTATION_SRC
#include "annotation.h"
#undef ANNOTATION_SRC

typedef struct
{
  int index;
  int r, g, b;
  int annotation;
  char name[100];
} ATABLE_ELT;

static ATABLE_ELT *atable;
static int num_entries = 0;

/*-----------------------------------------------*/
int print_annotation_table(FILE *fp)
{
  int n;
  if (num_entries <= 0) {
    read_annotation_table();
  }

  for (n = 0; n < num_entries; n++) {
    fprintf(fp, "%3d   %s\n", atable[n].index, atable[n].name);
  }
  return (0);
}

/*-----------------------------------------------*/
int print_annotation_colortable(FILE *fp)
{
  int n;
  if (num_entries <= 0) {
    read_annotation_table();
  }

  for (n = 0; n < num_entries; n++)
    fprintf(
        fp, "%3d   %-40s  %3d %3d %3d  0\n", atable[n].index, atable[n].name, atable[n].r, atable[n].g, atable[n].b);
  return (0);
}


std::vector<int> readAnnotationIntoVector(const std::string& filename)
{
  FILE *file = fopen(filename.c_str(), "r");
  if (!file) fs::fatal() << "could not open " << filename;

  std::map<int,int> annotmap;
  int num = freadInt(file);

  for (int i = 0; i < num; i++) {
    int vno = freadInt(file);
    int v = freadInt(file);
    annotmap[vno] = v;
  }

  fclose(file);

  int maxv = annotmap.rbegin()->first;
  std::vector<int> annotlist(maxv, 0);
  for (auto const& elt : annotmap) annotlist[elt.first] = elt.second;
  return annotlist;
}


int read_named_annotation_table(const char *name)
{
  FILE *fp;
  char fname[STRLEN], line[STRLEN];
  const char *cp; 
  int i;

  if (num_entries) {
    return (NO_ERROR); /* already read */
  }

  cp = strchr(name, '/');
  if (!cp) /* no path - use same one as mris was read from */
  {
    cp = getenv("FREESURFER_HOME");
    if (!cp) {
      cp = ".";
    }
    sprintf(fname, "%s/%s", cp, name);
  }
  else {
    cp = ""; /* use path in name */
    sprintf(fname, "%s", name);
  }

  fp = fopen(fname, "r");
  if (!fp) {
    fprintf(stderr, "could not open translation file %s\n", fname);
    return (ERROR_NO_FILE);
  }

  num_entries = 0;
  do {
    cp = fgetl(line, 199, fp);
    if (!cp) {
      break;
    }
    num_entries++;
  } while (cp && !feof(fp));

  rewind(fp);

  atable = (ATABLE_ELT *)calloc(num_entries, sizeof(ATABLE_ELT));
  for (i = 0; i < num_entries; i++) {
    cp = fgetl(line, 199, fp);
    if (!cp) {
      break;
    }
    sscanf(cp, "%d %s %d %d %d %*d", &atable[i].index, atable[i].name, &atable[i].r, &atable[i].g, &atable[i].b);
    atable[i].annotation = atable[i].r + (atable[i].g << 8) + (atable[i].b << 16);
  }
  return (NO_ERROR);
}

/*-----------------------------------------------*/
int read_annotation_table(void)
{
  FILE *fp;
  char fname[STRLEN], line[STRLEN];
  const char *cp;
  int i;
  extern char *annotation_table_file;

  if (num_entries) {
    return (NO_ERROR); /* already read */
  }

  cp = getenv("FREESURFER_HOME");
  if (!cp) {
    cp = ".";
  }

  if (annotation_table_file == NULL) {
    sprintf(fname, "%s/Simple_surface_labels2009.txt", cp);
  }
  else {
    sprintf(fname, "%s", annotation_table_file);
  }

  fp = fopen(fname, "r");
  if (!fp) {
    fprintf(stderr, "could not open translation file %s\n", fname);
    return (ERROR_NO_FILE);
  }

  num_entries = 0;
  do {
    cp = fgetl(line, 199, fp);
    if (!cp) {
      break;
    }
    num_entries++;
  } while (cp && !feof(fp));

  rewind(fp);

  atable = (ATABLE_ELT *)calloc(num_entries, sizeof(ATABLE_ELT));
  for (i = 0; i < num_entries; i++) {
    cp = fgetl(line, 199, fp);
    if (!cp) {
      break;
    }
    sscanf(cp, "%d %s %d %d %d %*d", &atable[i].index, atable[i].name, &atable[i].r, &atable[i].g, &atable[i].b);
    atable[i].annotation = atable[i].r + (atable[i].g << 8) + (atable[i].b << 16);
  }
  return (NO_ERROR);
}
int annotation_to_index(int annotation)
{
  int i;

  if (num_entries <= 0) {
    read_annotation_table();
  }

  for (i = 0; i < num_entries; i++) {
    if (atable[i].annotation == annotation) {
      return (atable[i].index);
    }
  }

  return (-1);
}

const char *index_to_name(int index)
{
  int i;

  if (num_entries <= 0) {
    read_annotation_table();
  }

  if (num_entries < 0) {
    static char name[100];

    sprintf(name, "%d", index);
    return (name);
  }

  for (i = 0; i < num_entries; i++) {
    if (atable[i].index == index) {
      return (atable[i].name);
    }
  }

  return ("NOT_FOUND");
}

int index_to_annotation(int index)
{
  int i;

  if (num_entries <= 0) {
    read_annotation_table();
  }

  for (i = 0; i < num_entries; i++) {
    if (atable[i].index == index) {
      return (atable[i].annotation);
    }
  }

  return (-1);
}

const char *annotation_to_name(int annotation, int *pindex)
{
  int i;

  if (num_entries <= 0) {
    read_annotation_table();
  }

  if (num_entries < 0) {
    static char name[100];

    if (pindex) {
      *pindex = -1;
    }
    sprintf(name, "%d", annotation);
    return (name);
  }

  for (i = 0; i < num_entries; i++) {
    if (atable[i].annotation == annotation) {
      if (pindex) {
        *pindex = atable[i].index;
      }
      return (atable[i].name);
    }
  }
  if (pindex) {
    *pindex = -1;
  }
  return ("NOT_FOUND");
}
/*------------------------------------------------------------
  annotation2label() - converts an annotation into a label
  given the index of the annotation in the color table.
  If no vertices with the index can be found, returns NULL.
------------------------------------------------------------*/
LABEL *annotation2label(int annotid, MRIS *Surf)
{
  int npoints, vtxno, annot, vtxannotid;
  VERTEX *vtx;
  LABEL *label;

  // Count number of points in the label
  npoints = 0;
  for (vtxno = 0; vtxno < Surf->nvertices; vtxno++) {
    vtx = &(Surf->vertices[vtxno]);
    annot = Surf->vertices[vtxno].annotation;
    // Given this annotation, find its index in the ctab
    if (Surf->ct) {
      CTABfindAnnotation(Surf->ct, annot, &vtxannotid);
    }
    else {
      vtxannotid = annotation_to_index(annot);
    }
    if (vtxannotid == annotid) {
      npoints++;
    }
  }
  if (npoints == 0) {
    return (NULL);
  }

  // Allocate the label
  label = LabelAlloc(npoints, "", "");
  label->n_points = npoints;

  // Fill the label
  npoints = 0;
  for (vtxno = 0; vtxno < Surf->nvertices; vtxno++) {
    vtx = &(Surf->vertices[vtxno]);
    annot = Surf->vertices[vtxno].annotation;
    if (Surf->ct) {
      CTABfindAnnotation(Surf->ct, annot, &vtxannotid);
    }
    else {
      vtxannotid = annotation_to_index(annot);
    }
    if (vtxannotid == annotid) {
      label->lv[npoints].vno = vtxno;
      label->lv[npoints].x = vtx->x;
      label->lv[npoints].y = vtx->y;
      label->lv[npoints].z = vtx->z;
      npoints++;
    }
  }
  return (label);
}

int set_atable_from_ctable(COLOR_TABLE *pct)
{
  CTE *cte;
  int i;

  if (pct == NULL) {
    return (ERROR_BAD_PARM);
  }

  if (num_entries > 0)  // atable already set
  {
    return (NO_ERROR);
  }

  num_entries = pct->nentries;

  if (num_entries <= 0) {
    return (ERROR_BAD_PARM);
  }

  atable = (ATABLE_ELT *)calloc(num_entries, sizeof(ATABLE_ELT));
  for (i = 0; i < num_entries; i++) {
    cte = pct->entries[i];
    if (NULL != cte) {
      atable[i].index = i;
      CTABcopyName(pct, i, atable[i].name, sizeof(atable[i].name));
      CTABrgbAtIndexi(pct, i, &atable[i].r, &atable[i].g, &atable[i].b);
      CTABannotationAtIndex(pct, i, &atable[i].annotation);
    }
  }

  return (NO_ERROR);
}
int MRISdivideAnnotation(MRI_SURFACE *mris, int *nunits)
{
  int *done, vno, index, nadded, i, num, j, annot, new_annot;
  VERTEX *v;
  COLOR_TABLE *ct;
  int rgb_scale = 30;  // need to pass as arg

  MRIScomputeMetricProperties(mris);
  MRIScomputeSecondFundamentalForm(mris);
  done = (int *)calloc(mris->ct->nentries, sizeof(int));
  if (done == NULL)
    ErrorExit(ERROR_NOMEMORY, "ERROR: MRISdivideAnnotation: could not allocate %d index table\n", mris->ct->nentries);

  MRISclearMarks(mris);
  MRISsetNeighborhoodSizeAndDist(mris, 2);
  for (nadded = vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag || v->annotation <= 0) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    
    CTABfindAnnotation(mris->ct, v->annotation, &index);
    if (index <= 0 || done[index])  // don't do unknown (index = 0)
    {
      continue;
    }
    if (index == Gdiag_no) {
      DiagBreak();
    }
    if (Gdiag_no >= 0 && v->annotation == mris->vertices[Gdiag_no].annotation)
      DiagBreak();
    num = MRISdivideAnnotationUnit(mris, v->annotation, nunits[index]);
    nadded += (num + 1);
    done[index] = 1 + num;
  }

  if (DIAG_VERBOSE_ON)
    printf("allocating new colortable with %d additional units...\n", nadded);
  ct = CTABalloc(mris->ct->nentries + nadded);
  index = mris->ct->nentries;
  for (i = 0; i < mris->ct->nentries; i++) {
    if (i == Gdiag_no) {
      DiagBreak();
    }
    if (mris->ct->entries[i] == NULL) {
      continue;
    }
    *(ct->entries[i]) = *(mris->ct->entries[i]);
    for (j = 0; done[i] > 1 && j < done[i]; j++) {
      int offset, new_index, ri, gi, bi, found, old_index;

      *(ct->entries[index]) = *(mris->ct->entries[i]);
      auto cx = snprintf(ct->entries[index]->name, STRLEN, "%s_div%d", ct->entries[i]->name, j + 1);
      if( (cx<0) || (cx>STRLEN) ) {
	std::cerr << __FUNCTION__ << ": snprintf returned error value" << std::endl;
      }
      offset = j;
      found = 0;
      do {
        ri = (ct->entries[i]->ri + rgb_scale * offset) % 256;
        gi = (ct->entries[i]->gi + rgb_scale * offset) % 256;
        bi = (ct->entries[i]->bi + rgb_scale * offset) % 256;
        CTABfindRGBi(ct, ri, gi, bi, &new_index);
        if (new_index < 0)  // couldn't find this r,g,b set - can use it for new entry
        {
          ct->entries[index]->ri = ri;
          ct->entries[index]->gi = gi;
          ct->entries[index]->bi = bi;

          ct->entries[index]->rf = (float)ri / 255.0f;
          ct->entries[index]->gf = (float)gi / 255.0f;
          ct->entries[index]->bf = (float)bi / 255.0f;

          CTABannotationAtIndex(ct, i, &annot);
          CTABannotationAtIndex(ct, index, &new_annot);
	  CTABfindAnnotation(mris->ct, new_annot, &old_index);
	  if (old_index >= 0)
	    continue ;  // not unique
          found = 1;

          // translate old annotations to new ones
          for (vno = 0; vno < mris->nvertices; vno++) {
	    if (vno == Gdiag_no)
	      DiagBreak() ;
            v = &mris->vertices[vno];
            if (v->ripflag || v->marked != j || v->annotation != annot) {
              continue;
            }
	    if (vno == Gdiag_no)
	      DiagBreak() ;
            v->annotation = new_annot;
	    v->marked = -1 ;   // don't process it again
          }
        }
        else {
          offset++;
          found = 0;
        }
      } while (!found && offset < 256);
      index++;
    }
  }
  CTABfree(&mris->ct);
  mris->ct = ct;
  return (NO_ERROR);
}

/*
  will split up parcellation units based on area and write the index of the
  new unit into the marked field (marked=0->old annot, marked=1-> 1st new
  subdivision, etc..)
  and will also put the continuous measure of how far along the eigendirection
  the
  vertex is into v->curv (in [0 1]).

  return the # of additional parcellation units that have been added
*/
int MRISdivideAnnotationUnit(MRI_SURFACE *mris, int annot, int nunits)
{
  int vno, num, min_vno;
  VERTEX *v, *vc;
  float cx, cy, cz, dist, min_dist, evalues[3], u, w, dx, dy, dz, min_dot, max_dot, e1x, e1y, e1z, dot;
  MATRIX *m_obs, *m_obs_T, *m_cov, *m_eig;

  if (nunits < 2) {
    return (0);
  }

  // compute centroid of annotation
  cx = cy = cz = 0;
  for (num = vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag || v->annotation != annot) {
      continue;
    }
    cx += v->x;
    cy += v->y;
    cz += v->z;
    num++;
  }
  if (num == 0)  // unused parcellation
  {
    return (0);
  }

  cx /= num;
  cy /= num;
  cz /= num;

  // find vertex in annotation closest to centroid
  min_dist = 100000;
  min_vno = -1;
  vc = NULL;
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag || v->annotation != annot) {
      continue;
    }
    dist = sqrt(SQR(v->x - cx) + SQR(v->y - cy) + SQR(v->z - cz));
    if (dist < min_dist) {
      min_vno = vno;
      min_dist = dist;
      vc = v;
    }
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    printf("using v %d as closest (%2.3f mm) from centroid (%2.1f, %2.1f, %2.1f)\n", min_vno, min_dist, cx, cy, cz);

  // now compute eigensystem around this vertex
  m_obs = MatrixAlloc(2, num, MATRIX_REAL);
  for (num = vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag || v->annotation != annot) {
      continue;
    }
    dx = v->x - cx;
    dy = v->y - cy;
    dz = v->z - cz;
    u = vc->e1x * dx + vc->e1y * dy + vc->e1z * dz;
    w = vc->e2x * dx + vc->e2y * dy + vc->e2z * dz;
    *MATRIX_RELT(m_obs, 1, num + 1) = u;
    *MATRIX_RELT(m_obs, 2, num + 1) = w;
    num++;
  }

  m_obs_T = MatrixTranspose(m_obs, NULL);
  m_cov = MatrixMultiply(m_obs, m_obs_T, NULL);
  m_eig = MatrixEigenSystem(m_cov, evalues, NULL);
  e1x = *MATRIX_RELT(m_eig, 1, 1) * vc->e1x + *MATRIX_RELT(m_eig, 2, 1) * vc->e2x;
  e1y = *MATRIX_RELT(m_eig, 1, 1) * vc->e1y + *MATRIX_RELT(m_eig, 2, 1) * vc->e2y;
  e1z = *MATRIX_RELT(m_eig, 1, 1) * vc->e1z + *MATRIX_RELT(m_eig, 2, 1) * vc->e2z;
  if (e1y < 0)  // orient them from posterior to anterior
  {
    e1x *= -1;
    e1y *= -1;
    e1z *= -1;
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    MatrixPrint(stdout, m_eig);
  }
  dist = sqrt(e1x * e1x + e1y * e1y + e1z * e1z);
  e1x /= dist;
  e1y /= dist;
  e1z /= dist;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    printf("principle eigendirection = (%2.3f, %2.3f, %2.3f)\n", e1x, e1y, e1z);
  }
  min_dot = 10000;
  max_dot = -min_dot;
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag || v->annotation != annot) {
      continue;
    }
    dx = v->x - cx;
    dy = v->y - cy;
    dz = v->z - cz;
    dot = dx * e1x + dy * e1y + dz * e1z;
    if (dot > max_dot) {
      max_dot = dot;
    }
    if (dot < min_dot) {
      min_dot = dot;
    }
  }

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag || v->annotation != annot) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    dx = v->x - cx;
    dy = v->y - cy;
    dz = v->z - cz;
    dot = dx * e1x + dy * e1y + dz * e1z;
    v->curv = (dot - min_dot) / (max_dot - min_dot);
    v->marked = (int)(nunits * v->curv);
    if (v->marked >= nunits) {
      v->marked = nunits - 1;  // one with max_dot will be just too big
    }
  }

  MatrixFree(&m_obs);
  MatrixFree(&m_cov);
  MatrixFree(&m_obs_T);
  return (nunits - 1);
}

/*!
  \fn int MRISmergeAnnotations(MRIS *mris, int nparcs, char **parcnames, char
  *newparcname)
  \param mris - surface structure
  \param nparcs - number of parcs to merge
  \param parcnames - names of parcs to merge
  \param newparcname - name of new parcellation
  \brief Merges parcellations into a single parcellation with the new name. The
  color
  will be that of the first parcellation in the list of names.
  Example:
    parcnames[0] = "caudalmiddlefrontal";
    parcnames[1] = "rostralmiddlefrontal";
    MRISmergeAnnotations(surf, 2, parcnames, "frontal");
*/
int MRISmergeAnnotations(MRIS *mris, int nparcs, std::vector<std::string> parcnames, const char *newparcname)
{
  int err, nthparc, parcid, nnewparcs, nthnewparc, m, match;
  int vtxno, *annotlist;
  COLOR_TABLE *ct;

  if (nparcs == 1) {
    printf("ERROR: nparcs must be > 1\n");
    return (1);
  }

  printf("MRISmergeAnnotations: parcCount=%d, newparcname=%s\n", nparcs, newparcname);

  // Make sure each parc name is in the parcellation
  // Get the list of annotation numbers too
  annotlist = (int *)calloc(nparcs, sizeof(int));
  for (m = 0; m < nparcs; m++) {
    err = CTABfindName(mris->ct, parcnames[m].c_str(), &parcid);
    if (err) {
      printf("ERROR: cannot find %s in annotation\n", parcnames[m].c_str());
      return (1);
    }
    CTABannotationAtIndex(mris->ct, parcid, &(annotlist[m]));
  }

  // Create a new color table
  // The merged parc gets the same color as the first listed parc
  nnewparcs = mris->ct->nentries - nparcs + 1;
  ct = CTABalloc(nnewparcs);
  nthnewparc = 0;
  for (nthparc = 0; nthparc < mris->ct->nentries; nthparc++) {
    // This checks whether the nth parc is in the list to merge
    match = 0;
    for (m = 0; m < nparcs; m++) {
      if (!strcmp(mris->ct->entries[nthparc]->name, parcnames[m].c_str())) {
        match = 1;
        break;
      }
    }
    if (match && m != 0) {
      continue;
    }
    // Gets here if it is not in the list, or if it is in the list
    // and it is the first in the list

    // Copy colors to the new color table
    *(ct->entries[nthnewparc]) = *(mris->ct->entries[nthparc]);

    // If it is the first in the list, change its name
    if (m == 0) {
      sprintf(ct->entries[nthnewparc]->name, "%s", newparcname);
    }

    nthnewparc++;
  }

  // Now change the vertex annotation values of the parcs
  // in the list to that of the 1st list member
  for (vtxno = 0; vtxno < mris->nvertices; vtxno++) {
    for (m = 0; m < nparcs; m++) {
      if (mris->vertices[vtxno].annotation == annotlist[m]) {
        mris->vertices[vtxno].annotation = annotlist[0];
        break;
      }
    }
  }

  CTABfree(&mris->ct);
  mris->ct = ct;
  free(annotlist);

  return (NO_ERROR);
}

/*---------------------------------------------------------------*/
/*!
  \fn MRI *MRISannot2seg(MRIS *surf, int base)
  \brief Constructs a 'segmentation' MRI from an annotation. The
  segmentation index is the annotation ID + base. See mri_annnotation2label
  for more details on setting the base. The returned
  MRI is a volume-encoded surface. If an annotation is not found in
  the color table, it is given a value equal to the base (assumes
  base is "unknown").
*/
MRI *MRISannot2seg(MRIS *surf, int base)
{
  int k, annot, annotid;
  MRI *seg;

  seg = MRIalloc(surf->nvertices, 1, 1, MRI_INT);

  for (k = 0; k < surf->nvertices; k++) {
    annot = surf->vertices[k].annotation;
    if (surf->ct) {
      CTABfindAnnotation(surf->ct, annot, &annotid);
    }
    else {
      annotid = annotation_to_index(annot);
    }
    if (annotid == -1) {
      annotid = 0;  // Assume base is unknown
    }
    MRIsetVoxVal(seg, k, 0, 0, 0, annotid + base);
  }
  return (seg);
}
/*---------------------------------------------------------------*/
/*!
  \fn MRI *MRISannot2border(MRIS *surf)
  Creates a binary overlay that is 1 for vertices as the border
  of parcellations and 0 everywhere else.
*/
MRI *MRISannot2border(MRIS *surf)
{
  MRI* const border = MRIalloc(surf->nvertices, 1, 1, MRI_INT);
  int k;
  for (k = 0; k < surf->nvertices; k++) {
    int const nnbrs = surf->vertices_topology[k].vnum;
    int const annot = surf->vertices[k].annotation;
    int isborder = 0;
    int nthnbr;
    for (nthnbr = 0; nthnbr < nnbrs; nthnbr++) {
      int const knbr     = surf->vertices_topology[k].v[nthnbr];
      int const nbrannot = surf->vertices[knbr].annotation;
      if (nbrannot != annot) {
        isborder = 1;
        break;
      }
    }
    MRIsetVoxVal(border, k, 0, 0, 0, isborder);
  }
  return (border);
}

/*---------------------------------------------------------------*/
/*!
  \fn int MRISaparc2lobes(MRIS *surf)
  Merges aparc labels into lobes.

      Lobar division type can be specified by the 'a_lobeDivisionType' --
      see 'mris_annotation2Label.c' for the enumerated types.
*/
int MRISaparc2lobes(MRIS *surf, int a_lobeDivisionType)
{
  int parcCount = 0;
  std::vector<std::string> parcnames(64);

  parcnames[0] = "caudalmiddlefrontal";
  parcnames[1] = "superiorfrontal";
  parcnames[2] = "rostralmiddlefrontal";
  parcnames[3] = "parsopercularis";
  parcnames[4] = "parstriangularis";
  parcnames[5] = "parsorbitalis";
  parcnames[6] = "lateralorbitofrontal";
  parcnames[7] = "medialorbitofrontal";
  parcnames[8] = "paracentral";
  parcnames[9] = "frontalpole";
  switch (a_lobeDivisionType) {
    case 0:
      parcCount = 10;
      break;
    case 1:
    case 2:
      parcnames[10] = "precentral";
      parcCount = 11;
      break;
  }
  MRISmergeAnnotations(surf, parcCount, parcnames, "frontal");

  switch (a_lobeDivisionType) {
    case 0:
    case 1:
      parcnames[0] = "superiortemporal";
      parcnames[1] = "entorhinal";
      parcnames[2] = "temporalpole";
      parcnames[3] = "fusiform";
      parcnames[4] = "inferiortemporal";
      parcnames[5] = "middletemporal";
      parcnames[6] = "parahippocampal";
      parcnames[7] = "bankssts";
      parcnames[8] = "transversetemporal";
      parcCount = 9;
      break;
    case 2:
      parcnames[0] = "superiortemporal";
      parcnames[1] = "inferiortemporal";
      parcnames[2] = "middletemporal";
      parcnames[3] = "bankssts";
      parcnames[4] = "transversetemporal";
      parcCount = 5;
      break;
  }
  MRISmergeAnnotations(surf, parcCount, parcnames, "temporal");

  switch (a_lobeDivisionType) {
    case 2:
      parcnames[0] = "entorhinal";
      parcnames[1] = "temporalpole";
      parcnames[2] = "fusiform";
      parcnames[3] = "parahippocampal";
      MRISmergeAnnotations(surf, 4, parcnames, "parahippocampalgyrus");
      break;
    default:
      break;
  }

  parcnames[0] = "supramarginal";
  parcnames[1] = "inferiorparietal";
  parcnames[2] = "superiorparietal";
  parcnames[3] = "precuneus";
  switch (a_lobeDivisionType) {
    case 0:
      parcCount = 4;
      break;
    case 1:
    case 2:
      parcnames[4] = "postcentral";
      parcCount = 5;
      break;
  }
  MRISmergeAnnotations(surf, parcCount, parcnames, "parietal");

  parcnames[0] = "pericalcarine";
  parcnames[1] = "cuneus";
  parcnames[2] = "lingual";
  parcnames[3] = "lateraloccipital";
  MRISmergeAnnotations(surf, 4, parcnames, "occipital");

  parcnames[0] = "isthmuscingulate";
  parcnames[1] = "posteriorcingulate";
  parcnames[2] = "caudalanteriorcingulate";
  parcnames[3] = "rostralanteriorcingulate";
  MRISmergeAnnotations(surf, 4, parcnames, "cingulate");

  return (0);
}

/* -------------------------------------------------------------------
   int MRISfbirnAnnot(MRIS *surf) - creates an annotation for fBIRN
   by dividing some aparcs or other geometric manipulations. Based on
   conversations with Jim Fallon.
   -----------------------------------------------------------------*/
int MRISfbirnAnnot(MRIS *surf)
{
  int *nunits;
  int area32p, area32v, superiorfrontal, medialorbitofrontal;
  int rostralanteriorcingulate, rostralmiddlefrontal;
  int index;
  COLOR_TABLE *ct;
  MRI *mri;

  // Create new entries in the CTAB for area32p and area32v
  ct = CTABaddEntry(surf->ct, "area32p");
  CTABfree(&surf->ct);
  surf->ct = ct;
  ct = CTABaddEntry(surf->ct, "area32v");
  CTABfree(&surf->ct);
  surf->ct = ct;

  // Get indices into CTAB for labels of interest
  CTABfindName(surf->ct, "area32p", &area32p);
  CTABfindName(surf->ct, "area32v", &area32v);
  CTABfindName(surf->ct, "superiorfrontal", &superiorfrontal);
  CTABfindName(surf->ct, "medialorbitofrontal", &medialorbitofrontal);
  CTABfindName(surf->ct, "rostralanteriorcingulate", &rostralanteriorcingulate);
  CTABfindName(surf->ct, "rostralmiddlefrontal", &rostralmiddlefrontal);

  /* Create area32v by creating a mask of the region bordering
     MOF and RA Cingulate, then dilating 12 times and constraining
     to be in MOF. The new area is carved out of MOF.*/
  mri = MRISfbirnMask_MOF_RACing(surf);
  MRISdilateConfined(surf, mri, medialorbitofrontal, 12, area32v);
  MRIfree(&mri);

  /* Create area32p by creating a mask of the region bordering
     SFG and several Cingulates, then dilating 12 times and constraining
     to be in SFG. The new area is carved out of MOF.*/
  mri = MRISfbirnMask_SFG_Cing(surf);
  MRISdilateConfined(surf, mri, superiorfrontal, 12, area32p);
  MRIfree(&mri);

  // Now, divide up some units
  nunits = (int *)calloc(surf->ct->nentries, sizeof(int));
  nunits[rostralanteriorcingulate] = 3;
  nunits[rostralmiddlefrontal] = 3;
  nunits[superiorfrontal] = 5;
  nunits[area32p] = 2;
  // nunits[area32v] = 2;

  MRISdivideAnnotation(surf, nunits);

  CTABfindName(surf->ct, "area32p", &index);
  sprintf(surf->ct->entries[index]->name, "%s", "area32p_pseudo");
  CTABfindName(surf->ct, "area32p_div2", &index);
  sprintf(surf->ct->entries[index]->name, "%s", "area32a_pseudo");
  CTABfindName(surf->ct, "area32v", &index);
  sprintf(surf->ct->entries[index]->name, "%s", "area32v_pseudo");
  // CTABfindName(surf->ct, "area32v_div2", &index);
  // sprintf(surf->ct->entries[index]->name, "%s","area32v_pseudo");

  CTABfindName(surf->ct, "rostralanteriorcingulate", &index);
  sprintf(surf->ct->entries[index]->name, "%s", "area24d_pseudo");
  CTABfindName(surf->ct, "rostralanteriorcingulate_div2", &index);
  sprintf(surf->ct->entries[index]->name, "%s", "area24pg_pseudo");
  CTABfindName(surf->ct, "rostralanteriorcingulate_div3", &index);
  sprintf(surf->ct->entries[index]->name, "%s", "area24v_pseudo");

  free(nunits);

  return (0);
}

/*---------------------------------------------------------------*/
/*!
  \fn double *MRISannotDice(MRIS *surf1, MRIS *surf2, int *nsegs, int
  **segidlist)
  \brief Computes dice coefficient for each parcellation unit. surf1
  and surf2 should be the same surface with different parcellations
  loaded.  *nsegs returns the number of parcellations. **segidlist is
  a list of the parcellation id numbers (usually just 0 to nsegs-1.
*/
double *MRISannotDice(MRIS *surf1, MRIS *surf2, int *nsegs, int **segidlist)
{
  MRI *seg1, *seg2;
  int k, id1, id2, k1 = 0, k2 = 0, vtxno;
  int nsegid1, *segidlist1;
  int nsegid2, *segidlist2;
  double *area1, *area2, *area12, *dice;
  *nsegs = -1;

  // Create a seg from the 1st annot
  seg1 = MRISannot2seg(surf1, 0);
  // Extract a unique, sorted list of the ids
  segidlist1 = MRIsegIdList(seg1, &nsegid1, 0);

  // Create a seg from the 2nd annot
  seg2 = MRISannot2seg(surf2, 0);
  // Extract a unique, sorted list of the ids
  segidlist2 = MRIsegIdList(seg1, &nsegid2, 0);

  if (nsegid1 != nsegid2) {
    printf("ERROR: MRISannotDice(): nsegs do not match %d %d\n", nsegid1, nsegid2);
    return (NULL);
  }
  // Note: segidlist1 and 2 should be the same too
  printf("MRISannotDice(): found %d segs\n", nsegid1);
  *nsegs = nsegid1;

  area1 = (double *)calloc(nsegid1, sizeof(double));
  area2 = (double *)calloc(nsegid1, sizeof(double));
  area12 = (double *)calloc(nsegid1, sizeof(double));

  for (vtxno = 0; vtxno < surf1->nvertices; vtxno++) {
    // id at vtxno for 1st annot
    id1 = MRIgetVoxVal(seg1, vtxno, 0, 0, 0);
    // determine its index in the segidlist
    for (k = 0; k < nsegid1; k++) {
      if (id1 == segidlist1[k]) {
        k1 = k;
        break;
      }
    }
    // id at vtxno for 2nd annot
    id2 = MRIgetVoxVal(seg2, vtxno, 0, 0, 0);
    // determine its index in the segidlist
    for (k = 0; k < nsegid1; k++) {
      if (id2 == segidlist2[k]) {
        k2 = k;
        break;
      }
    }
    // accum areas
    area1[k1] += surf1->vertices[vtxno].area;
    area2[k2] += surf1->vertices[vtxno].area;
    if (id1 == id2) {
      area12[k1] += surf1->vertices[vtxno].area;
    }
  }

  // Compute dice for each area
  dice = (double *)calloc(nsegid1, sizeof(double));
  for (k = 0; k < nsegid1; k++) {
    dice[k] = area12[k] / ((area1[k] + area2[k]) / 2.0);
  }

  MRIfree(&seg1);
  MRIfree(&seg2);
  free(area1);
  free(area2);
  free(area12);
  free(segidlist2);
  *segidlist = segidlist1;

  return (dice);
}


/*!
  \brief Reads an annotation file into a 1D seg overlay. Expects that the
  file has an embedded lookup table so that it can convert annotation values
  to segmentation indices.
*/
MRI *readAnnotationIntoSeg(const std::string& filename)
{
  FILE *f = fopen(filename.c_str(), "r");
  if (!f) fs::fatal() << "could not open annot file " << filename;

  // first int is nvertices
  int nvertices = freadInt(f);
  MRI *overlay = new MRI({nvertices, 1, 1}, MRI_INT);

  // read annotations into seg
  for (int n = 0; n < nvertices; n++) {
    int vno = freadInt(f);
    int val = freadInt(f);
    MRIsetVoxVal(overlay, vno, 0, 0, 0, val);
  }

  // check for embedded colortable
  int tag = freadInt(f);
  if (!feof(f) && TAG_OLD_COLORTABLE == tag) overlay->ct = CTABreadFromBinary(f);

  fclose(f);

  // since we're converting from annot to seg index values, we
  // should expect that the file contains a colortable
  if (overlay->ct) {
    for (int vno = 0; vno < nvertices; vno++) {
      int annot = MRIgetVoxVal(overlay, vno, 0, 0, 0);
      CTABfindAnnotation(overlay->ct, annot, &annot);
      MRIsetVoxVal(overlay, vno, 0, 0, 0, annot);
    }
  } else {
    fs::warning() << "reading annotation without embedded lookup table";
  }

  return overlay;
}


/*!
  \brief Writes an annotation file from a 1D seg overlay. Expects that the
  seg MRI has a lookup table so that it can convert segmentation indices to
  annotation values.
*/
void writeAnnotationFromSeg(const MRI *overlay, const std::string& filename)
{
  FILE *f = fopen(filename.c_str(), "wb");
  if (!f) fs::fatal() << "could not write annot file " << filename;

  // first int is nvertices
  int nvertices = overlay->width;
  fwriteInt(nvertices, f);

  // convert seg values back into annotation values (if lookup table exists)
  for (int vno = 0; vno < nvertices; vno++) {
    // round in case it's a floating point overlay
    int annot = std::round(MRIgetVoxVal(overlay, vno, 0, 0, 0));
    if (overlay->ct) {
      if (annot < 0) {
         // set all negative seg values to 0 and don't search for annotation
        annot = 0;
      } else {
        // set annotation from seg value
        CTABannotationAtIndex(overlay->ct, annot, &annot);
      }
    }
    fwriteInt(vno, f);
    fwriteInt(annot, f);
  }

  // embed lookup table
  if (overlay->ct) {
    fwriteInt(TAG_OLD_COLORTABLE, f);
    CTABwriteIntoBinary(overlay->ct, f);
  } else {
    fs::warning() << "writing annotation without embedded lookup table";
  }

  fclose(f);
}
