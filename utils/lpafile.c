/**
 * @file  lpafile.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:45 $
 *    $Revision: 1.11 $
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
#include <unistd.h>
#include <hipl_format.h>
#include <memory.h>

#include "lpafile.h"
#include "utils.h"
#include "machine.h"
#include "error.h"
#include "proto.h"
#include "hipsu.h"
#include "image.h"
#include "macros.h"

static int lpafFillEntries(LP_ANSWER_FILE *lpaf, char *fname, int entryno) ;
static void lpafDump(FILE *fp, LP_ANSWER_FILE *lpaf) ;
static void lpafAllocParms(IMAGE *I) ;

#define NFLUSH    10     /* flush after every 10th write */
#define INIT_VAL  10000  /* initial value of parameters */

int
LPAFreadImageAnswer(LPAF *lpaf, int current)
{
  char   *fullname, fname[100] ;
  FILE   *fp ;
  IMAGE  Iheader ;
  int    i, ecode, frame, current_frame ;
  LP_BOX *lpb ;
  struct extpar *xp ;
#ifdef _MSDOS
  long   *parms ;
#else
  int    *parms ;
#endif

  fullname = lpaf->filelist[current] ;

  ImageUnpackFileName(fullname, &current_frame, &i, fname) ;

  fp = fopen(fname, "rb") ;
  if (!fp)
    ErrorReturn(-1,(ERROR_NO_FILE,"LPAFreadImageAnswer(%d): could not open %s",
                    current, fname)) ;

  ecode = fread_header(fp, &Iheader, fname) ;
  fclose(fp) ;
  if (ecode)
    ErrorReturn(-2, (ERROR_BADFILE,
                     "LPAFreadImageAnswer(%s): could not read header",fname));

  if (Iheader.numparam < Iheader.num_frame)
    return(0) ;

  /* read answer from header */
#if 0
  fprintf(stderr, "reading lp values from %dth entry in image file\n",
          current_frame);
#endif
  lpb = &lpaf->coords[current] ;
  for (frame = 0, xp = Iheader.params ; xp ; xp = xp->nextp)
    if (frame++ == current_frame)
      break ;

  /*
   if hips file created on Sun, then the parameters are actually longs.
  */
#ifndef _MSDOS
  parms = xp->val.v_pi ;
#else
  parms = (long *)xp->val.v_pi ;
#endif

#ifndef _MSDOS
  if (parms[0] < 0 || parms[0] >= Iheader.cols)
  {
    parms[0] = swapInt(parms[0]) ;
    parms[1] = swapInt(parms[1]) ;
    for (i = 0 ; i < NPOINTS ; i++)
    {
      parms[2+2*i] = swapInt(parms[2*i]) ;
      parms[2+2*i+1] = swapInt(parms[2*i+1]) ;
    }
  }
#else
  if (parms[0] < 0 || parms[0] >= (long)Iheader.cols)
  {
    parms[0] = swapLong(parms[0]) ;
    parms[1] = swapLong(parms[1]) ;
    for (i = 0 ; i < NPOINTS ; i++)
    {
      parms[2+2*i] = swapLong(parms[2*i]) ;
      parms[2+2*i+1] = swapLong(parms[2*i+1]) ;
    }
  }
#endif

  if ((int)parms[0] == INIT_VAL)  /* not yet written with real value */
    return(0) ;

  lpb->xc = (int)parms[0] ;
  lpb->yc  = (int)parms[1] ;
  for (i = 0 ; i < NPOINTS ; i++)
  {
    lpb->xp[i] = (int)parms[2+2*i] ;
    lpb->yp[i] = (int)parms[2+2*i+1] ;
  }

  if (lpb->xc < 0 || lpb->xc >= Iheader.cols ||
      lpb->yc < 0 || lpb->xc >= Iheader.rows )
    return(0) ;

  return(1) ;
}

int
LPAFresetImageAnswer(LPAF *lpaf, int current)
{
  char   *fullname, tmpname[100], fname[100] ;
  IMAGE  Iheader, *I ;
  FILE   *infp, *outfp ;
  int    ecode, frame, nframes, *parms, i, current_frame, type ;
  LP_BOX *lpb ;
  struct extpar *xp ;

  fullname = lpaf->filelist[current] ;

  ImageUnpackFileName(fullname, &current_frame, &type, fname) ;
  if (type != HIPS_IMAGE)
    return(0) ;

  infp = fopen(fname, "rb") ;
  if (!infp)
    ErrorReturn(-1,(ERROR_NO_FILE,
                    "LPAFwriteImageAnswer(%d): could not open %s",
                    current, fname)) ;
  ecode = fread_header(infp, &Iheader, fname) ;
  if (ecode)
  {
    fclose(infp) ;
    ErrorReturn(-2, (ERROR_BADFILE,
                     "LPAFwriteImageAnswer(%d): could not read header", current));
  }
  fclose(infp) ;
  if (Iheader.numparam == 0)  /* must make room for header in image file */
  {
    lpafAllocParms(&Iheader) ;

    /* now copy the old image file to a new one which has room for parms */
    nframes = Iheader.num_frame ;
    strcpy(tmpname, FileTmpName(NULL)) ;
    outfp = fopen(tmpname, "wb") ;
    Iheader.num_frame = 0 ;  /* ImageAppend will bump num_frame on each call */
    ecode = fwrite_header(outfp, &Iheader, tmpname) ;
    fclose(outfp) ;   /* flush file */

    fprintf(stderr, "rewriting image file to make room for parms...\n") ;
    for (frame = 0 ; frame < nframes ; frame++)
    {
      I = ImageReadFrames(fname, frame, 1) ;
      if (!I)
        ErrorExit(ERROR_BADFILE, "LPwriteImageAnswer: could not read frame");
      ImageAppend(I, tmpname) ;
      ImageFree(&I) ;
    }
    FileRename(tmpname, fname) ;
    Iheader.num_frame = nframes ;  /* reset correct # of frames */
  }

  /* now put answer into header */
  lpb = &lpaf->coords[current] ;

  for (frame = 0, xp = Iheader.params ; xp ; xp = xp->nextp)
    if (frame++ == current_frame)
      break ;

  parms = xp->val.v_pi ;
  parms[0] = INIT_VAL ;
  parms[1] = INIT_VAL ;
  for (i = 0 ; i < NPOINTS ; i++)
  {
    parms[2+2*i] = INIT_VAL ;
    parms[2+2*i+1] = INIT_VAL ;
  }
  ImageUpdateHeader(&Iheader, fname) ;
  free_hdrcon(&Iheader) ;
  return(1) ;
}

int
LPAFwriteImageAnswer(LPAF *lpaf, int current)
{
  char   *fullname, tmpname[100], fname[100] ;
  IMAGE  Iheader, *I ;
  FILE   *infp, *outfp ;
  int    ecode, frame, nframes, *parms, i, current_frame, type ;
  LP_BOX *lpb ;
  struct extpar *xp ;

  fullname = lpaf->filelist[current] ;

  ImageUnpackFileName(fullname, &current_frame, &type, fname) ;
  if (type != HIPS_IMAGE)
    return(0) ;

  infp = fopen(fname, "rb") ;
  if (!infp)
    ErrorReturn(-1,(ERROR_NO_FILE,
                    "LPAFwriteImageAnswer(%d): could not open %s",
                    current, fname)) ;
  ecode = fread_header(infp, &Iheader, fname) ;
  if (ecode)
  {
    fclose(infp) ;
    ErrorReturn(-2, (ERROR_BADFILE,
                     "LPAFwriteImageAnswer(%d): could not read header", current));
  }
  fclose(infp) ;
  if (Iheader.numparam == 0)  /* must make room for header in image file */
  {
    lpafAllocParms(&Iheader) ;

    /* now copy the old image file to a new one which has room for parms */
    nframes = Iheader.num_frame ;
    strcpy(tmpname, FileTmpName(NULL)) ;
    outfp = fopen(tmpname, "wb") ;
    Iheader.num_frame = 0 ;  /* ImageAppend will bump num_frame on each call */
    ecode = fwrite_header(outfp, &Iheader, tmpname) ;
    fclose(outfp) ;   /* flush file */

    fprintf(stderr, "rewriting image file to make room for parms...\n") ;
    for (frame = 0 ; frame < nframes ; frame++)
    {
      I = ImageReadFrames(fname, frame, 1) ;
      if (!I)
        ErrorExit(ERROR_BADFILE, "LPwriteImageAnswer: could not read frame");
      ImageAppend(I, tmpname) ;
      ImageFree(&I) ;
    }
    FileRename(tmpname, fname) ;
    Iheader.num_frame = nframes ;  /* reset correct # of frames */
  }

  /* now put answer into header */
  lpb = &lpaf->coords[current] ;

  for (frame = 0, xp = Iheader.params ; xp ; xp = xp->nextp)
    if (frame++ == current_frame)
      break ;

  parms = xp->val.v_pi ;
  parms[0] = lpb->xc ;
  parms[1] = lpb->yc ;
  for (i = 0 ; i < NPOINTS ; i++)
  {
    parms[2+2*i] = lpb->xp[i] ;
    parms[2+2*i+1] = lpb->yp[i] ;
  }
  ImageUpdateHeader(&Iheader, fname) ;
  free_hdrcon(&Iheader) ;
  return(1) ;
}

LP_ANSWER_FILE *
LPAFcreate(char *out_fname, int argc, char *argv[])
{
  LP_ANSWER_FILE *lpaf ;
  int            i, nentries, nfiles, entryno ;

  nfiles = 0 ;
  for (i = 0 ; i < argc ; i++)
  {
    nentries = FileNumberOfEntries(argv[i]) ;
    nfiles += nentries ;
  }

  if (nfiles <= 0)
    ErrorReturn(NULL, (ERROR_NO_FILE, "LPAFcreate: no valid files specified"));

  lpaf = (LP_ANSWER_FILE *)calloc(1, sizeof(*lpaf)) ;
  if (!lpaf)
    ErrorExit(ERROR_NO_MEMORY, "LPAFcreate: allocation failed") ;
  if (out_fname)   /* write output to file as well as into image file */
  {
    strcpy(lpaf->fname, out_fname) ;
    unlink(lpaf->fname) ;
    lpaf->fp = fopen(lpaf->fname, "a+b") ;
    if (!lpaf->fp)
      ErrorReturn(NULL,
                  (ERROR_NO_FILE,"LPAFcreate: could not open %s",out_fname));
  }
  else
    lpaf->fp = NULL ;

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
      sprintf(buf, "%s#%d", base_name, i) ;
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


/* filename (centroid) (x, y) ... */
#define FILE_FMT  "%s (%3d, %3d) (%3d, %3d) (%3d, %3d) (%3d, %3d) (%3d, %3d)\n"

int
LPAFwrite(LPAF *lpaf, int current)
{
  LP_BOX *lpb ;

  LPAFwriteImageAnswer(lpaf, current) ;

  if (lpaf->fp)
  {
    lpb = &lpaf->coords[current] ;
    if (lpb->fpos >= 0L)   /* overwrite previous entry */
    {
      if (fseek(lpaf->fp, lpb->fpos, SEEK_SET) < 0)
        ErrorReturn(-1, (ERROR_BADFILE, "LPAFwrite could not seek to %ld",
                         lpb->fpos)) ;
    }
    else
      lpb->fpos = ftell(lpaf->fp) ;

    if (current > lpaf->last_written)
      lpaf->last_written = current ;
    else   /* write out rest of file */
    {}

    fprintf(lpaf->fp, FILE_FMT,
            lpaf->filelist[current], lpb->xc, lpb->yc,
            lpb->xp[0], lpb->yp[0], lpb->xp[1], lpb->yp[1],
            lpb->xp[2], lpb->yp[2], lpb->xp[3], lpb->yp[3]) ;

    if (lpaf->flush++ >= NFLUSH)
    {
      fflush(lpaf->fp) ;
      lpaf->flush = 0 ;
    }
  }
  return(0) ;
}

int
LPAFread(LPAF *lpaf, int current)
{
  LP_BOX *lpb ;
  char   line[300], *cp ;

  lpb = &lpaf->coords[current] ;

  if (LPAFreadImageAnswer(lpaf, current) <= 0)  /* not in image file */
  {
    if (!lpaf->fp || (lpb->fpos < 0))
      return(0) ;     /* hasn't been written yet */

    if (fseek(lpaf->fp, lpb->fpos, SEEK_SET) < 0)
      ErrorReturn(-1, (ERROR_BADFILE, "LPAFread could not seek to %ld",
                       lpb->fpos)) ;

    cp = fgetl(line, 299, lpaf->fp) ;
    if (!cp)
      ErrorReturn(-1, (ERROR_BADFILE, "LPAFread: could not read line")) ;

    if (sscanf(cp, FILE_FMT,
               lpaf->filelist[current], &lpb->xc, &lpb->yc,
               &lpb->xp[0], &lpb->yp[0], &lpb->xp[1], &lpb->yp[1],
               &lpb->xp[2], &lpb->yp[2], &lpb->xp[3], &lpb->yp[3]) != 11)
      ErrorReturn(-1,
                  (ERROR_BADFILE, "LPAFread: could not scan all parms from %s",
                   cp)) ;

  }

#if 0
  {
    int i ;

    for (i = 0 ; i < NPOINTS ; i++)
      fprintf(stderr, "(%d, %d) ", lpb->xp[i], lpb->yp[i]) ;
    fprintf(stderr, "\n") ;
  }
#endif

  return(abs(lpb->xc) < INIT_VAL) ;  /* handles garbages as well as unwritten */
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

#define NPARAMS   10

static void
lpafAllocParms(IMAGE *I)
{
  int           frame, nframes, i ;
  struct extpar *params ;

  nframes = I->num_frame ;
  I->numparam = nframes ;
  I->params = params = (struct extpar *)calloc(nframes, sizeof(*params)) ;
  if (!params)
    ErrorExit(ERROR_NO_MEMORY, "lpafAllocParms: could not allocate (%d, %d)",
              nframes, sizeof(*params)) ;
  I->paramdealloc = TRUE ;

  /* set up link-list and allocate enough space */
  for (frame = 0 ; frame < nframes ; frame++)
  {
    if (frame < (nframes - 1))
      params[frame].nextp = &params[frame+1] ;
    else
      params[frame].nextp = NULL ;

    params[frame].name = STRCPALLOC("LP_ANSWER") ;
    params[frame].format = PFINT ;
    params[frame].count = NPARAMS ; /* 4 corners and centroid */
    params[frame].dealloc = TRUE ;
    params[frame].val.v_pi = (int *)calloc(NPARAMS, sizeof(int)) ;
    for (i = 0 ; i < NPARAMS ; i++)
      params[frame].val.v_pi[i] = INIT_VAL ;
  }
}

