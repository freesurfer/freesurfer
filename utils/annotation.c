#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "utils.h"
#include "const.h"
#include "error.h"
#include "annotation.h"

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
		cp = getenv("MRI_DIR") ;
		if (!cp)
			cp = "." ;
	}
	else
		cp = "" ;  /* use path in name */

  sprintf(fname, "%s/%s", cp, name) ;
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

  if (num_entries)
    return(NO_ERROR) ;   /* already read */

  cp = getenv("MRI_DIR") ;
  if (!cp)
    cp = "." ;

  sprintf(fname, "%s/christophe_parc.txt", cp) ;
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
