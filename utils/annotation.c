#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "utils.h"
#include "const.h"
#include "error.h"
#include "mrisurf.h"
#include "label.h"
#include "colortab.h"

#define ANNOTATION_SRC
#include "annotation.h"
#undef ANNOTATION_SRC

typedef struct
{
  int    index ;
  int    r, g, b ;
  int    annotation ;
  char   name[100] ;
} ATABLE_ELT ;

static ATABLE_ELT *atable ;
static int num_entries = 0 ;

/*-----------------------------------------------*/
int print_annotation_table(FILE *fp)
{
  int n;
  if (num_entries <= 0) read_annotation_table() ;

  for(n = 0; n < num_entries; n++)
    fprintf(fp,"%3d   %s\n",atable[n].index,atable[n].name);
  return(0);
}

/*-----------------------------------------------*/
int print_annotation_colortable(FILE *fp)
{
  int n;
  if (num_entries <= 0) read_annotation_table() ;

  for(n = 0; n < num_entries; n++)
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
  else{
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
  } while (cp && !feof(fp)) ;

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

  if(annotation_table_file == NULL)
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
  } while (cp && !feof(fp)) ;

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
  given the annotation identifier. If no vertices with the
  identifier can be found, returns NULL.
------------------------------------------------------------*/
LABEL *annotation2label(int annotid, MRIS *Surf)
{
  int npoints, vtxno, annot, vtxannotid;
  VERTEX *vtx;
  LABEL *label;

  // Count number of points in the label 
  npoints = 0;
  for(vtxno = 0; vtxno < Surf->nvertices; vtxno++){
    vtx = &(Surf->vertices[vtxno]);
    annot = Surf->vertices[vtxno].annotation;
    if (Surf->ct)
      CTABfindAnnotation(Surf->ct, annot, &vtxannotid);
    else
      vtxannotid = annotation_to_index(annot);
    if(vtxannotid == annotid) npoints++;
  }
  if(npoints==0) return(NULL);

  // Allocate the label
  label = LabelAlloc(npoints,"","");
  label->n_points = npoints;

  // Fill the label
  npoints = 0;
  for(vtxno = 0; vtxno < Surf->nvertices; vtxno++){
    vtx = &(Surf->vertices[vtxno]);
    annot = Surf->vertices[vtxno].annotation;
    if (Surf->ct)
      CTABfindAnnotation(Surf->ct, annot, &vtxannotid);
    else
      vtxannotid = annotation_to_index(annot);
    if(vtxannotid == annotid){
      label->lv[npoints].vno = vtxno;
      label->lv[npoints].x = vtx->x;
      label->lv[npoints].y = vtx->y;
      label->lv[npoints].z = vtx->z;
      npoints++;
    }
  }
  return(label);
}


int set_atable_from_ctable(COLOR_TABLE *pct){
  CTE *cte;
  int i;

  if(pct == NULL)
    return(ERROR_BAD_PARM);
  
  if (num_entries > 0) // atable already set
    return(NO_ERROR);

  num_entries = pct->nentries;

  if(num_entries <= 0)
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
