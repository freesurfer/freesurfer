/**
 * @file  annotation.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2007/05/08 03:47:42 $
 *    $Revision: 1.18 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "utils.h"
#include "const.h"
#include "error.h"
#include "mrisurf.h"
#include "label.h"
#include "colortab.h"
#include "diag.h"

#define ANNOTATION_SRC
#include "annotation.h"
#undef ANNOTATION_SRC

typedef struct
{
  int    index ;
  int    r, g, b ;
  int    annotation ;
  char   name[100] ;
}
ATABLE_ELT ;

static ATABLE_ELT *atable ;
static int num_entries = 0 ;

/*-----------------------------------------------*/
int print_annotation_table(FILE *fp)
{
  int n;
  if (num_entries <= 0) read_annotation_table() ;

  for (n = 0; n < num_entries; n++)
    fprintf(fp,"%3d   %s\n",atable[n].index,atable[n].name);
  return(0);
}

/*-----------------------------------------------*/
int print_annotation_colortable(FILE *fp)
{
  int n;
  if (num_entries <= 0) read_annotation_table() ;

  for (n = 0; n < num_entries; n++)
    fprintf(fp,"%3d   %-40s  %3d %3d %3d  0\n",atable[n].index,atable[n].name,
            atable[n].r,atable[n].g,atable[n].b);
  return(0);
}

int
read_named_annotation_table(char *name)
{

  FILE  *fp ;
  char  *cp, fname[STRLEN], line[STRLEN] ;
  int   i ;

  if (num_entries)
    return(NO_ERROR) ;   /* already read */

  cp = strchr(name, '/') ;
  if (!cp)                 /* no path - use same one as mris was read from */
  {
    cp = getenv("FREESURFER_HOME") ;
    if (!cp)
      cp = "." ;
    sprintf(fname, "%s/%s", cp, name) ;
  }
  else
  {
    cp = "" ;  /* use path in name */
    sprintf(fname, "%s", name) ;
  }

  fp = fopen(fname, "r") ;
  if (!fp)
  {
    fprintf(stderr, "could not open translation file %s\n", fname) ;
    return(ERROR_NO_FILE) ;
  }

  num_entries = 0 ;
  do
  {
    cp = fgetl(line, 199, fp) ;
    if (!cp)
      break ;
    num_entries++ ;
  }
  while (cp && !feof(fp)) ;

  rewind(fp) ;

  atable = (ATABLE_ELT *)calloc(num_entries, sizeof(ATABLE_ELT)) ;
  for (i = 0 ; i < num_entries ; i++)
  {
    cp = fgetl(line, 199, fp) ;
    if (!cp)
      break ;
    sscanf(cp, "%d %s %d %d %d %*d",
           &atable[i].index,
           atable[i].name,
           &atable[i].r,
           &atable[i].g,
           &atable[i].b) ;
    atable[i].annotation = atable[i].r+(atable[i].g << 8)+(atable[i].b << 16);
  }
  return(NO_ERROR) ;
}

/*-----------------------------------------------*/
int
read_annotation_table(void)
{
  FILE  *fp ;
  char  *cp, fname[STRLEN], line[STRLEN] ;
  int   i ;
  extern char *annotation_table_file;

  if (num_entries)
    return(NO_ERROR) ;   /* already read */

  cp = getenv("FREESURFER_HOME") ;
  if (!cp)
    cp = "." ;

  if (annotation_table_file == NULL)
    sprintf(fname, "%s/Simple_surface_labels2002.txt", cp) ;
  else
    sprintf(fname, "%s",annotation_table_file);

  fp = fopen(fname, "r") ;
  if (!fp)
  {
    fprintf(stderr, "could not open translation file %s\n", fname) ;
    return(ERROR_NO_FILE) ;
  }

  num_entries = 0 ;
  do
  {
    cp = fgetl(line, 199, fp) ;
    if (!cp)
      break ;
    num_entries++ ;
  }
  while (cp && !feof(fp)) ;

  rewind(fp) ;

  atable = (ATABLE_ELT *)calloc(num_entries, sizeof(ATABLE_ELT)) ;
  for (i = 0 ; i < num_entries ; i++)
  {
    cp = fgetl(line, 199, fp) ;
    if (!cp)
      break ;
    sscanf(cp, "%d %s %d %d %d %*d",
           &atable[i].index,
           atable[i].name,
           &atable[i].r,
           &atable[i].g,
           &atable[i].b) ;
    atable[i].annotation = atable[i].r+(atable[i].g << 8)+(atable[i].b << 16);
  }
  return(NO_ERROR) ;
}
int
annotation_to_index(int annotation)
{
  int   i ;

  if (num_entries <= 0)
    read_annotation_table() ;

  for (i = 0 ; i < num_entries ; i++)
  {
    if (atable[i].annotation == annotation)
      return(atable[i].index) ;
  }

  return(-1) ;
}

char *
index_to_name(int index)
{
  int   i ;

  if (num_entries <= 0)
    read_annotation_table() ;

  if (num_entries < 0)
  {
    static char name[100] ;

    sprintf(name, "%d", index) ;
    return(name) ;
  }

  for (i = 0 ; i < num_entries ; i++)
  {
    if (atable[i].index == index)
    {
      return(atable[i].name) ;
    }
  }

  return("NOT_FOUND") ;
}

int
index_to_annotation(int index)
{
  int   i ;

  if (num_entries <= 0)
    read_annotation_table() ;

  for (i = 0 ; i < num_entries ; i++)
  {
    if (atable[i].index == index)
      return(atable[i].annotation) ;
  }

  return(-1) ;
}

char *
annotation_to_name(int annotation, int *pindex)
{
  int   i ;

  if (num_entries <= 0)
    read_annotation_table() ;

  if (num_entries < 0)
  {
    static char name[100] ;

    if (pindex)
      *pindex = -1 ;
    sprintf(name, "%d", annotation) ;
    return(name) ;
  }

  for (i = 0 ; i < num_entries ; i++)
  {
    if (atable[i].annotation == annotation)
    {
      if (pindex)
        *pindex = atable[i].index ;
      return(atable[i].name) ;
    }
  }
  if (pindex)
    *pindex = -1 ;
  return("NOT_FOUND") ;
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
  for (vtxno = 0; vtxno < Surf->nvertices; vtxno++)
  {
    vtx = &(Surf->vertices[vtxno]);
    annot = Surf->vertices[vtxno].annotation;
    // Given this annotation, find its index in the ctab
    if (Surf->ct)
      CTABfindAnnotation(Surf->ct, annot, &vtxannotid);
    else
      vtxannotid = annotation_to_index(annot);
    if (vtxannotid == annotid) npoints++;
  }
  if (npoints==0) return(NULL);

  // Allocate the label
  label = LabelAlloc(npoints,"","");
  label->n_points = npoints;

  // Fill the label
  npoints = 0;
  for (vtxno = 0; vtxno < Surf->nvertices; vtxno++)
  {
    vtx = &(Surf->vertices[vtxno]);
    annot = Surf->vertices[vtxno].annotation;
    if (Surf->ct)
      CTABfindAnnotation(Surf->ct, annot, &vtxannotid);
    else
      vtxannotid = annotation_to_index(annot);
    if (vtxannotid == annotid)
    {
      label->lv[npoints].vno = vtxno;
      label->lv[npoints].x = vtx->x;
      label->lv[npoints].y = vtx->y;
      label->lv[npoints].z = vtx->z;
      npoints++;
    }
  }
  return(label);
}


int set_atable_from_ctable(COLOR_TABLE *pct)
{
  CTE *cte;
  int i;

  if (pct == NULL)
    return(ERROR_BAD_PARM);

  if (num_entries > 0) // atable already set
    return(NO_ERROR);

  num_entries = pct->nentries;

  if (num_entries <= 0)
    return(ERROR_BAD_PARM);

  atable = (ATABLE_ELT *)calloc(num_entries, sizeof(ATABLE_ELT)) ;
  for (i = 0 ; i < num_entries ; i++)
  {
    cte = pct->entries[i];
    if (NULL != cte)
    {
      atable[i].index = i;
      CTABcopyName(pct, i, atable[i].name, sizeof(atable[i].name));
      CTABrgbAtIndexi(pct, i, &atable[i].r, &atable[i].g, &atable[i].b );
      CTABannotationAtIndex(pct, i, &atable[i].annotation);
    }
  }

  return(NO_ERROR) ;
}
int MRISdivideAnnotation(MRI_SURFACE *mris, int *nunits) {
  int   *done, vno, index, nadded, i, num, j, annot, new_annot ;
  VERTEX *v ;
  COLOR_TABLE *ct ;
  int rgb_scale = 30 ; // need to pass as arg


  MRIScomputeMetricProperties(mris) ;
  MRIScomputeSecondFundamentalForm(mris) ;
  done = (int *)calloc(mris->ct->nentries, sizeof(int)) ;
  if (done == NULL)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d index table",
              Progname, mris->ct->nentries) ;

  MRISclearMarks(mris) ;
  MRISsetNeighborhoodSize(mris, 2) ;
  for (nadded = vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->annotation <= 0)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    CTABfindAnnotation(mris->ct, v->annotation, &index) ;
    if (index <= 0 || done[index])  // don't do unknown (index = 0)
      continue ;
    if (index == Gdiag_no)
      DiagBreak() ;
#if 0
    if (stricmp("postcentral", mris->ct->entries[index]->name))
      continue ;
#endif
    num = MRISdivideAnnotationUnit(mris, v->annotation, nunits[index]) ;
    nadded += num ;
    done[index] = 1+num ;
  }

  printf("allocating new colortable with %d additional units...\n", nadded) ;
  ct = CTABalloc(mris->ct->nentries+nadded) ;
  index = mris->ct->nentries ;
  for (i = 0 ; i < mris->ct->nentries ; i++) {
    if (i == Gdiag_no)
      DiagBreak() ;
    *(ct->entries[i]) = *(mris->ct->entries[i]) ;
    for (j = 1 ; j < done[i] ; j++) {
      int offset, new_index, ri, gi, bi, found ;

      *(ct->entries[index]) = *(mris->ct->entries[i]) ;
      sprintf(ct->entries[index]->name, "%s_div%d", ct->entries[i]->name, j+1) ;
      offset = j ;
      found = 0 ;
      do {
        ri = (ct->entries[i]->ri+rgb_scale*offset) % 256 ;
        gi = (ct->entries[i]->gi+rgb_scale*offset) % 256 ;
        bi = (ct->entries[i]->bi+rgb_scale*offset) % 256 ;
        CTABfindRGBi(ct, ri, gi, bi, &new_index) ;
        if (new_index < 0)  // couldn't find this r,g,b set - can use it for new entry
        {
          ct->entries[index]->ri = ri ;
          ct->entries[index]->gi = gi ;
          ct->entries[index]->bi = bi ;

          ct->entries[index]->rf = (float)ri/255.0f;
          ct->entries[index]->gf = (float)gi/255.0f;
          ct->entries[index]->bf = (float)bi/255.0f;
          found = 1 ;

          CTABannotationAtIndex(ct, i, &annot) ;
          CTABannotationAtIndex(ct, index, &new_annot) ;
          // translate old annotations to new ones
          for (vno = 0 ; vno < mris->nvertices ; vno++) {
            v = &mris->vertices[vno] ;
            if (v->ripflag || v->marked != j || v->annotation != annot)
              continue ;
            v->annotation = new_annot ;
          }
        }
        else {
          offset++ ;
          found = 0 ;
        }
      } while (!found) ;
      index++ ;
    }
  }
  CTABfree(&mris->ct) ;
  mris->ct = ct ;
  return(NO_ERROR) ;
}

/*
  will split up parcellation units based on area and write the index of the
  new unit into the marked field (marked=0->old annot, marked=1-> 1st new subdivision, etc..)
  and will also put the continuous measure of how far along the eigendirection the
  vertex is into v->curv (in [0 1]).

  return the # of additional parcellation units that have been added
*/
int
MRISdivideAnnotationUnit(MRI_SURFACE *mris, int annot, int nunits) {
  int    vno, num, min_vno ;
  VERTEX *v, *vc ;
  float  cx, cy, cz, dist, min_dist, evalues[3], u, w, dx, dy, dz,
  min_dot, max_dot, e1x, e1y, e1z, dot, mx ;
  MATRIX *m_obs, *m_obs_T, *m_cov, *m_eig ;


  if (nunits < 2)
    return(0);

  // compute centroid of annotation
  cx = cy = cz = 0 ;
  for (num = vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->annotation != annot)
      continue ;
    cx += v->x ;
    cy += v->y ;
    cz += v->z ;
    num++ ;
  }
  if (num == 0) // unused parcellation
    return(0) ;

  cx /= num ;
  cy /= num ;
  cz /= num ;

  // find vertex in annotation closest to centroid
  min_dist = 100000 ;
  min_vno = -1 ;
  vc = NULL ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->annotation != annot)
      continue ;
    dist = sqrt(SQR(v->x-cx)+SQR(v->y-cy)+SQR(v->z-cz));
    if (dist < min_dist) {
      min_vno = vno;
      min_dist = dist ;
      vc = v ;
    }
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    printf("using v %d as closest (%2.3f mm) from centroid (%2.1f, %2.1f, %2.1f)\n",
           min_vno, min_dist, cx, cy, cz) ;

  // now compute eigensystem around this vertex
  m_obs = MatrixAlloc(2, num, MATRIX_REAL) ;
  for (num = vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->annotation != annot)
      continue ;
    dx = v->x - cx ;
    dy = v->y - cy ;
    dz = v->z - cz ;
    u = vc->e1x * dx + vc->e1y * dy + vc->e1z*dz ;
    w = vc->e2x * dx + vc->e2y * dy + vc->e2z*dz ;
    *MATRIX_RELT(m_obs, 1, num+1) = u ;
    *MATRIX_RELT(m_obs, 2, num+1) = w ;
    num++ ;
  }

  m_obs_T = MatrixTranspose(m_obs, NULL) ;
  m_cov = MatrixMultiply(m_obs,m_obs_T, NULL) ;
  m_eig = MatrixEigenSystem(m_cov, evalues, NULL) ;
  e1x = *MATRIX_RELT(m_eig, 1,1) * vc->e1x + *MATRIX_RELT(m_eig, 2,1) * vc->e2x;
  e1y = *MATRIX_RELT(m_eig, 1,1) * vc->e1y + *MATRIX_RELT(m_eig, 2,1) * vc->e2y;
  e1z = *MATRIX_RELT(m_eig, 1,1) * vc->e1z + *MATRIX_RELT(m_eig, 2,1) * vc->e2z;
  if (fabs(e1x) > fabs(e1y) &&  fabs(e1x) > fabs(e1z))
    mx = e1x ;
  else if (fabs(e1y) > fabs(e1z))
    mx = e1y ;
  else
    mx = e1z ;
  //  if (mx < 0)
  if (e1y < 0)  // orient them from posterior to anterior
  {
    e1x *= -1 ;
    e1y *= -1 ;
    e1z *= -1 ;
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    MatrixPrint(stdout, m_eig) ;
  dist = sqrt(e1x*e1x + e1y*e1y + e1z*e1z) ;
  e1x /= dist ;
  e1y /= dist ;
  e1z /= dist ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    printf("principle eigendirection = (%2.3f, %2.3f, %2.3f)\n",e1x,e1y,e1z) ;
  min_dot = 10000 ;
  max_dot = -min_dot ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->annotation != annot)
      continue ;
    dx = v->x - cx ;
    dy = v->y - cy ;
    dz = v->z - cz ;
    dot = dx*e1x + dy*e1y + dz*e1z ;
    if (dot > max_dot)
      max_dot = dot ;
    if (dot < min_dot)
      min_dot = dot ;
  }

  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->annotation != annot)
      continue ;
    if  (vno == Gdiag_no)
      DiagBreak() ;
    dx = v->x - cx ;
    dy = v->y - cy ;
    dz = v->z - cz ;
    dot = dx*e1x + dy*e1y + dz*e1z ;
    v->curv = (dot - min_dot) / (max_dot-min_dot) ;
    v->marked = (int)(nunits * v->curv) ;
    if (v->marked >= nunits)
      v->marked = nunits-1 ;  // one with max_dot will be just too big
  }

  MatrixFree(&m_obs) ;
  MatrixFree(&m_cov) ;
  MatrixFree(&m_obs_T) ;
  return(nunits-1) ;
}

/*!
  \fn int MRISmergeAnnotations(MRIS *mris, int nparcs, char **parcnames, char *newparcname)
  \param mris - surface structure
  \param nparcs - number of parcs to merge
  \param parcnames - names of parcs to merge
  \param newparcname - name of new parcellation
  \brief Merges parcellations into a single parcellation with the new name. The color
  will be that of the first parcellation in the list of names.
  Example:
    parcnames[0] = "caudalmiddlefrontal";
    parcnames[1] = "rostralmiddlefrontal";
    MRISmergeAnnotations(surf, 2, parcnames, "frontal");
*/
int MRISmergeAnnotations(MRIS *mris, int nparcs, char **parcnames, char *newparcname)
{
  int err, nthparc, parcid, nnewparcs, nthnewparc, m, match;
  int vtxno, *annotlist;
  COLOR_TABLE *ct ;
  VERTEX *vtx;
  
  if(nparcs == 1){
    printf("ERROR: nparcs must be > 1\n");
    return(1);
  }

  // Make sure each parc name is in the parcellation
  // Get the list of annotation numbers too
  annotlist = (int *) calloc(nparcs,sizeof(int));
  for(m = 0; m < nparcs; m++){
    err = CTABfindName(mris->ct, parcnames[m], &parcid);
    if(err){
      printf("ERROR: cannot find %s in annotation\n",parcnames[m]);
      return(1);
    }
    CTABannotationAtIndex(mris->ct, parcid, &(annotlist[m])) ;
  }

  // Create a new color table
  // The merged parc gets the same color as the first listed parc
  nnewparcs = mris->ct->nentries - nparcs + 1;
  ct = CTABalloc(nnewparcs);
  nthnewparc = 0;
  for(nthparc = 0; nthparc < mris->ct->nentries; nthparc++){

    // This checks whether the nth parc is in the list to merge
    match = 0;
    for(m = 0; m < nparcs; m++){
      if(!strcmp(mris->ct->entries[nthparc]->name,parcnames[m])){
	match = 1;
	break;
      }
    }
    if(match && m != 0) continue;
    // Gets here if it is not in the list, or if it is in the list
    // and it is the first in the list

    // Copy colors to the new color table
    *(ct->entries[nthnewparc]) = *(mris->ct->entries[nthparc]);

    // If it is the first in the list, change its name
    if(m == 0)  sprintf(ct->entries[nthnewparc]->name,"%s",newparcname);

    nthnewparc ++;
  }

  // Now change the vertex annotation values of the parcs 
  // in the list to that of the 1st list member
  for(vtxno = 0; vtxno < mris->nvertices; vtxno++){
    vtx = &(mris->vertices[vtxno]);
    for(m = 0; m < nparcs; m++){
      if(mris->vertices[vtxno].annotation == annotlist[m]){
	mris->vertices[vtxno].annotation = annotlist[0];
	break;
      }
    }
  }

  CTABfree(&mris->ct) ;
  mris->ct = ct ;
  free(annotlist);

  return(NO_ERROR);
}
