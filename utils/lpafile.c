#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "lpafile.h"
#include "utils.h"
#include "error.h"
#include "proto.h"

static int lpafFillEntries(LP_ANSWER_FILE *lpaf, char *fname, int entryno) ;
static void lpafDump(FILE *fp, LP_ANSWER_FILE *lpaf) ;

LP_ANSWER_FILE *
LPAFcreate(char *out_fname, int argc, char *argv[])
{
  LP_ANSWER_FILE *lpaf ;
  int            i, nentries, nfiles, entryno ;

  nfiles = 0 ;
  for (i = 0 ; i < argc ; i++)
  {
    nentries = FileNumberOfEntries(argv[i]) ;
    fprintf(stderr, "argv[%d] = '%s', # %d\n", i, argv[i], nentries) ;
    nfiles += nentries ;
  }

  fprintf(stderr, "total of %d files\n", nfiles) ;
  if (nfiles <= 0)
    ErrorReturn(NULL, (ERROR_NO_FILE, "LPAFcreate: no valid files specified"));
  
  lpaf = (LP_ANSWER_FILE *)calloc(1, sizeof(*lpaf)) ;
  if (!lpaf)
    ErrorExit(ERROR_NO_MEMORY, "LPAFcreate: allocation failed") ;
  strcpy(lpaf->fname, out_fname) ;
  lpaf->fp = fopen(lpaf->fname, "a+b") ;
  if (!lpaf->fp)
    ErrorReturn(NULL,
                (ERROR_NO_FILE,"LPAFcreate: could not open %s",out_fname));
  lpaf->nfiles = nfiles ;
  lpaf->last_written = lpaf->current = 0 ;
  lpaf->filelist = (char **)calloc(nfiles, sizeof(char *)) ;
  if (!lpaf->filelist)
    ErrorExit(ERROR_NO_MEMORY, "LPAFcreate: allocation failed") ;
  lpaf->coords = (LP_BOX *)calloc(nfiles, sizeof(LP_BOX)) ;
  if (!lpaf->coords)
    ErrorExit(ERROR_NO_MEMORY, "LPAFcreate: allocation failed") ;
  lpaf->last_written = 0 ;

  /* now fill out filelist array */
  entryno = 0 ;
  for (i = 0 ; i < argc ; i++)
  {
    nentries = lpafFillEntries(lpaf, argv[i], entryno) ;
    entryno += nentries ;
  }

  for (i = 0 ; i < lpaf->nfiles ; i++)
    lpaf->coords[i].fpos = -1L ;

  lpafDump(stderr, lpaf) ;

  return(lpaf) ;
}
static int
lpafFillEntries(LP_ANSWER_FILE *lpaf, char *fname, int entryno)
{
  int  nentries, type, i, num ;
  char buf[100], *base_name, line[200], *cp ;
  FILE *fp ;

  nentries = FileNumberOfEntries(fname) ;
  type = FileType(fname) ;
  base_name = FileFullName(fname) ;

  for (i = 0 ; i < nentries ; i++)
  {
    switch (type)
    {
    case LIST_FILE:
      fp = fopen(base_name, "rb") ;
      if (!fp)
        ErrorReturn(0, (ERROR_NO_FILE, "lpafFillEntries: could not open %s\n",
                        base_name)) ;
      cp = fgetl(line, 199, fp) ;
      nentries = 0 ;
      while (cp)
      {
        sscanf(cp, "%s", buf) ;
        num = lpafFillEntries(lpaf, buf, entryno+nentries) ;
        nentries += num ;
        cp = fgetl(line, 199, fp) ;
      }
      fclose(fp) ;
      break ;
    default:
      sprintf(buf, "%s:%d", base_name, i) ;
      lpaf->filelist[entryno+i] = (char *)calloc(strlen(buf)+1, sizeof(char));
      strcpy(lpaf->filelist[entryno+i], buf) ;
      break ;
    }
  }
  return(nentries) ;
}



static void
lpafDump(FILE *fp, LP_ANSWER_FILE *lpaf)
{
  int i ;

  if (strcmp(lpaf->fname, "test"))
    return ;
  fprintf(fp, "lpaf has %d entries:\n", lpaf->nfiles) ;
  for (i = 0 ; i < lpaf->nfiles ; i++)
    fprintf(fp, "%s\n", lpaf->filelist[i]) ;
}
int
LPAFwrite(LPAF *lpaf, int current)
{
  LP_BOX *lpb ;
  int    i ;

  lpb = &lpaf->coords[current] ;
  lpb->fpos = ftell(lpaf->fp) ;
  fprintf(lpaf->fp, "%s (%d, %d) ", lpaf->filelist[current], lpb->xc, lpb->yc);
  for (i = 0 ; i < NPOINTS ; i++)
    fprintf(lpaf->fp, "(%d %d) ", lpb->xp[i], lpb->yp[i]) ;
  fflush(lpaf->fp) ;
  return(0) ;
}

int
LPAFread(LPAF *lpaf, int current)
{
  LP_BOX *lpb ;
  int    i ;
  char   fname[100] ;

  lpb = &lpaf->coords[current] ;
  if (lpb->fpos < 0)
    return(0) ;     /* hasn't been written yet */

  if (fseek(lpaf->fp, lpb->fpos, SEEK_SET) < 0)
    ErrorReturn(-1, (ERROR_BADFILE, "LPAFread could not seek to %ld",
                     lpb->fpos)) ;
  fscanf(lpaf->fp, "%s (%d, %d) ", fname, &lpb->xc, &lpb->yc) ;
  for (i = 0 ; i < NPOINTS ; i++)
    fscanf(lpaf->fp, "(%d %d) ", &lpb->xp[i], &lpb->yp[i]) ;

  fprintf(stderr, "read %s: (%d, %d)\n", fname, lpb->xc, lpb->yc);
  for (i = 0 ; i < NPOINTS ; i++)
    fprintf(stderr, "(%d, %d) ", lpb->xp[i], lpb->yp[i]) ;
  fprintf(stderr, "\n") ;

  return(1) ;
}

int
LPAFset(LPAF *lpaf, int current, int *xp, int *yp, int xc, int yc)
{
  LP_BOX *lpb ;
  int    i ;

  lpb = &lpaf->coords[current] ;
  lpb->xc = xc ;
  lpb->yc = yc ;
  for (i = 0 ; i < NPOINTS ; i++)
  {
    lpb->xp[i] = xp[i] ;
    lpb->yp[i] = yp[i] ;
  }
  return(1) ;
}

