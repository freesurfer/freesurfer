/**
 * @file  tags.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2007/01/09 08:03:46 $
 *    $Revision: 1.7 $
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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fio.h"
#include "error.h"
#include "tags.h"

/* $Id : $Exp */

int
TAGskip(FILE *fp, int tag, long long len)
{
#if 1
  unsigned char *buf ;
  int ret ;

  buf = (unsigned char *)calloc(len, sizeof(unsigned char)) ;
  if (buf==NULL)
    ErrorExit(ERROR_NOMEMORY,
              "TAGskip: failed to calloc %u bytes!\n",len);
  ret = fread(buf, sizeof(unsigned char), len, fp) ;
  free(buf) ;
  return(ret) ;
#else
  return(fseek(fp, len, SEEK_CUR)) ;  // doesn't work for gzipped files
#endif
}

int TAGreadStart(FILE *fp, long long *plen)
{
  int  tag ;

  tag = freadInt(fp) ;
  if (feof(fp))
    return(0) ;
  switch (tag)
  {
  case TAG_OLD_MGH_XFORM:
    *plen = (long long)freadInt(fp) ;  /* sorry - backwards compatibility
                                                            with Tosa's stuff */
    *plen = *plen -1 ; // doesn't include null
    break ;
  case TAG_OLD_SURF_GEOM:    // these don't take lengths at all
  case TAG_OLD_USEREALRAS:
  case TAG_OLD_COLORTABLE:
    *plen = 0 ;
    break ;
  default:
    *plen = freadLong(fp) ;
  }

  return(tag) ;
}

int
TAGwriteStart(FILE *fp, int tag, long long *phere, long long len)
{
  long here ;

  fwriteInt(tag, fp) ;

  here = ftell(fp) ;
  *phere = (long long)here ;
  fwriteLong(len, fp) ;

  return(NO_ERROR) ;
}

int TAGwrite(FILE *fp, int tag, void *buf, long long len)
{
  long long here ;

  TAGwriteStart(fp, tag, &here, len) ;
  fwrite(buf, sizeof(char), len, fp) ;
  TAGwriteEnd(fp, here) ;
  return(NO_ERROR) ;
}

int
TAGwriteEnd(FILE *fp, long long there)
{
  long long here ;

  here = ftell(fp) ;
#if 0
  fseek(fp, there, SEEK_SET) ;
  fwriteLong((long long)(here-(there+sizeof(long long))), fp) ;
  fseek(fp, here, SEEK_SET) ;
#endif

  return(NO_ERROR) ;
}

int
TAGmakeCommandLineString(int argc, char **argv, char *cmd_line)
{
  int i ;

  cmd_line[0] = 0 ;
  for (i = 0 ; i < argc ; i++)
  {
    strcat(cmd_line, argv[i]) ;
    strcat(cmd_line, " ") ;
  }
  return(NO_ERROR) ;
}

int
TAGwriteCommandLine(FILE *fp, char *cmd_line)
{
  long long here ;

  TAGwriteStart(fp, TAG_CMDLINE, &here, strlen(cmd_line)+1) ;
  fprintf(fp, "%s", cmd_line) ;
  TAGwriteEnd(fp, here) ;
  return(NO_ERROR) ;
}

int TAGwriteAutoAlign(FILE *fp, MATRIX *M)
{
  long long here ;
  char buf[16*100];
  long long len;

  bzero(buf,16*100);
  sprintf(buf,"AutoAlign %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf",
	  M->rptr[1][1],M->rptr[1][2],M->rptr[1][3],M->rptr[1][4],
	  M->rptr[2][1],M->rptr[2][2],M->rptr[2][3],M->rptr[2][4],
	  M->rptr[3][1],M->rptr[3][2],M->rptr[3][3],M->rptr[3][4],
	  M->rptr[4][1],M->rptr[4][2],M->rptr[4][3],M->rptr[4][4]);
  len = strlen(buf);
  TAGwriteStart(fp, TAG_AUTO_ALIGN, &here, len) ;
  fwrite(buf, sizeof(char), len, fp) ;

  TAGwriteEnd(fp, here) ;
  return(NO_ERROR) ;
}

MATRIX *TAGreadAutoAlign(FILE *fp)
{
  int c,r;
  char buf[1000];
  MATRIX *M;

  M = MatrixAlloc(4,4,MATRIX_REAL);
  fscanf(fp,"%s",buf); // get past "AutoAlign" string
  for(r=1; r<=4; r++){
    for(c=1; c<=4; c++){
      fscanf(fp,"%f",&(M->rptr[r][c]));
    }
  }
  return(M);
}

