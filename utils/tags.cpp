/**
 * @brief utils for adding tags (meta info) to mgz/h files
 *
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

#include "error.h"
#include "fio.h"
#include "machine.h"
#include "tags.h"

#define MATRIX_STRLEN (4 * 4 * 100)

/* For this to skip the tagid and data correctly, data-length is needed after tagid.
 * The format is 'tagid data-length data'.
 */
int TAGskip(FILE *fp, int tag, long long len)
{
#if 1
  unsigned char *buf;
  int ret;

  buf = (unsigned char *)calloc(len, sizeof(unsigned char));
  if (buf == NULL) ErrorExit(ERROR_NOMEMORY, "TAGskip: tag=%d, failed to calloc %u bytes!\n", tag, len);
  ret = fread(buf, sizeof(unsigned char), len, fp);
  free(buf);
  return (ret);
#else
  return (fseek(fp, len, SEEK_CUR));  // doesn't work for gzipped files
#endif
}

int TAGreadStart(FILE *fp, long long *plen)
{
  int tag;

  tag = freadInt(fp);
  if (feof(fp)) return (0);
  switch (tag) {
    case TAG_OLD_MGH_XFORM:
      *plen = (long long)freadInt(fp); /* sorry - backwards compatibility
                                          with Tosa's stuff */
      *plen = *plen - 1;               // doesn't include null
      break;
    case TAG_OLD_SURF_GEOM:  // these don't take lengths at all
    case TAG_OLD_USEREALRAS:
    case TAG_OLD_COLORTABLE:
      *plen = 0;  // these tags have no data-length output after tagid
      break;
    default:
      *plen = freadLong(fp);  // read data-length for the tagid 
  }

  return (tag);
}

/* This function is used to output 'tagid data-length'.
 *   1. write tag
 *   2. remember file position after tag write
 *   3. write data length (long)
 */
int TAGwriteStart(FILE *fp, int tag, long long *phere, long long len)
{
  long here;

  fwriteInt(tag, fp);

  here = ftell(fp);
  *phere = (long long)here;
  fwriteLong(len, fp);

  return (NO_ERROR);
}

int TAGwrite(FILE *fp, int tag, void *buf, long long len)
{
  long long here;

  TAGwriteStart(fp, tag, &here, len);
  fwrite(buf, sizeof(char), len, fp);
  TAGwriteEnd(fp, here);
  return (NO_ERROR);
}

int TAGwriteEnd(FILE *fp, long long there)
{
  long long here;

  // ??? why get current file position ???
  here = ftell(fp);
#if 0
  fseek(fp, there, SEEK_SET) ;
  fwriteLong((long long)(here-(there+sizeof(long long))), fp) ;
  fseek(fp, here, SEEK_SET) ;
#endif

  return (NO_ERROR);
}

int TAGmakeCommandLineString(int argc, char **argv, char *cmd_line)
{
  int i;

  cmd_line[0] = 0;
  for (i = 0; i < argc; i++) {
    strcat(cmd_line, argv[i]);
    strcat(cmd_line, " ");
  }
  return (NO_ERROR);
}

int TAGwriteCommandLine(FILE *fp, char *cmd_line)
{
  long long here;

  TAGwriteStart(fp, TAG_CMDLINE, &here, strlen(cmd_line) + 1);
  fprintf(fp, "%s", cmd_line);
  TAGwriteEnd(fp, here);
  return (NO_ERROR);
}

int TAGwriteAutoAlign(FILE *fp, MATRIX *M)
{
  long long here;
  char buf[MATRIX_STRLEN];
  long long len;

  bzero(buf, MATRIX_STRLEN);
  sprintf(buf,
          "AutoAlign %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf",
          M->rptr[1][1],
          M->rptr[1][2],
          M->rptr[1][3],
          M->rptr[1][4],
          M->rptr[2][1],
          M->rptr[2][2],
          M->rptr[2][3],
          M->rptr[2][4],
          M->rptr[3][1],
          M->rptr[3][2],
          M->rptr[3][3],
          M->rptr[3][4],
          M->rptr[4][1],
          M->rptr[4][2],
          M->rptr[4][3],
          M->rptr[4][4]);
  len = strlen(buf);
  TAGwriteStart(fp, TAG_AUTO_ALIGN, &here, len);
  fwrite(buf, sizeof(char), len, fp);

  TAGwriteEnd(fp, here);
  return (NO_ERROR);
}

MATRIX *TAGreadAutoAlign(FILE *fp)
{
  int c, r;
  char buf[MATRIX_STRLEN];
  MATRIX *M;

  M = MatrixAlloc(4, 4, MATRIX_REAL);
  fscanf(fp, "%s", buf);  // get past "AutoAlign" string
  for (r = 1; r <= 4; r++) {
    for (c = 1; c <= 4; c++) {
      fscanf(fp, "%f", &(M->rptr[r][c]));
    }
  }
  return (M);
}

/* zlib support */

/* This is the zlib version of TAGskip().
 * For this to skip the tagid and data correctly, data-length is needed after tagid.
 * The format is 'tagid data-length data'.
 */
int znzTAGskip(znzFile fp, int tag, long long len)
{
#if 1
  unsigned char *buf;
  int ret;

  buf = (unsigned char *)calloc(len, sizeof(unsigned char));
  if (buf == NULL) ErrorExit(ERROR_NOMEMORY, "znzTAGskip: tag=%d, failed to calloc %u bytes!\n", tag, len);
  ret = znzread(buf, sizeof(unsigned char), len, fp);
  free(buf);
  return (ret);
#else
  return (znzseek(fp, len, SEEK_CUR));  // doesn't work for gzipped files
#endif
}

int znzTAGreadStart(znzFile fp, long long *plen, int tagwithzerolen)
{
  int tag;

  if (znzeof(fp)) return (0);
  tag = znzreadInt(fp);
  if (znzeof(fp)) return (0);

  if (tagwithzerolen && tagwithzerolen == tag)
  {
    /* This is to handle following situation: 
     * TAG_MGH_XFORM is used in both mgz and m3z, but in different format.
     * in mgz, data-length is output after TAG_MGH_XFORM
     * in m3z, no data-length is output after TAG_MGH_XFORM
     *
     * in __m3zRead(), the function is called as znzTAGreadStart(file, &len, TAG_MGH_XFORM);
     */
    *plen = 0;
    return tag;
  }
  
  switch (tag) {
    case TAG_OLD_MGH_XFORM:  // the tag doesn't look like being used
      *plen = (long long)znzreadInt(fp); /* sorry - backwards compatibility
                                            with Tosa's stuff */
      *plen = *plen - 1;                 // doesn't include null
      break;
    case TAG_OLD_SURF_GEOM:  // these don't take lengths at all
    case TAG_OLD_USEREALRAS:
    case TAG_OLD_COLORTABLE:
    case TAG_GCAMORPH_GEOM:
    case TAG_GCAMORPH_TYPE:
    case TAG_GCAMORPH_LABELS:
      *plen = 0;  // these tags have no data-length output after tagid
      break;
    default:
      *plen = znzreadLong(fp);  // read data-length for the tagid
  }

  return (tag);
}

/* This is the zlib version of TAGwriteStart().
 * It is used to output 'tagid data-length'.
 *   1. write tag
 *   2. remember file position after tag write
 *   3. write data length (long)
 */
int znzTAGwriteStart(znzFile fp, int tag, long long *phere, long long len)
{
  long here;

  znzwriteInt(tag, fp);

  here = znztell(fp);
  *phere = (long long)here;
  znzwriteLong(len, fp);

  return (NO_ERROR);
}

/*
  Note: TAGwrite()  is not compatible with znzreadFloatEx()
  because TAGwrite() does not do a byte order swap. Use
  znzTAGreadFloat() instead of znzreadFloatEx().
 */
int znzTAGwrite(znzFile fp, int tag, void *buf, long long len)
{
  long long here;
  znzTAGwriteStart(fp, tag, &here, len);
  znzwrite(buf, sizeof(char), len, fp);
  znzTAGwriteEnd(fp, here);
  return (NO_ERROR);
}

int znzTAGwriteEnd(znzFile fp, long long there)
{
  long long here;

  here = znztell(fp);
#if 0
  fseek(fp, there, SEEK_SET) ;
  fwriteLong((long long)(here-(there+sizeof(long long))), fp) ;
  fseek(fp, here, SEEK_SET) ;
#endif

  return (NO_ERROR);
}

int znzTAGwriteCommandLine(znzFile fp, char *cmd_line)
{
  long long here;

  znzTAGwriteStart(fp, TAG_CMDLINE, &here, strlen(cmd_line) + 1);
  znzprintf(fp, "%s", cmd_line);
  znzTAGwriteEnd(fp, here);
  return (NO_ERROR);
}

int znzWriteMatrix(znzFile fp, MATRIX *M, int tag)
{
  long long here;
  char buf[MATRIX_STRLEN];

  bzero(buf, MATRIX_STRLEN);
  sprintf(buf,
          "Matrix %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf",
          M->rptr[1][1],
          M->rptr[1][2],
          M->rptr[1][3],
          M->rptr[1][4],
          M->rptr[2][1],
          M->rptr[2][2],
          M->rptr[2][3],
          M->rptr[2][4],
          M->rptr[3][1],
          M->rptr[3][2],
          M->rptr[3][3],
          M->rptr[3][4],
          M->rptr[4][1],
          M->rptr[4][2],
          M->rptr[4][3],
          M->rptr[4][4]);
  znzTAGwriteStart(fp, tag, &here, MATRIX_STRLEN);
  here = znztell(fp);
  znzwrite(buf, sizeof(char), MATRIX_STRLEN, fp);
  here = znztell(fp);

  znzTAGwriteEnd(fp, here);
  return (NO_ERROR);
}

/*
  Note: znzTAGreadFloat() is compatible with
  TAGwrite(). znzreadFloatEx() is not compatible with TAGwrite()
  because it performs a byte order swap. Use znzTAGreadFloat() instead
  of znzreadFloatEx().
 */
int znzTAGreadFloat(float *pf, znzFile fp)
{
  int ret;
  ret = znzreadFloatEx(pf, fp);
#if (BYTE_ORDER == LITTLE_ENDIAN)
  // Now unswap if needed
  *pf = swapFloat(*pf);
#endif
  return ret;
}

/*
  \fn MATRIX *znzReadMatrix(znzFile fp)
  Different than znzReadAutoAlignMatrix(znzFile fp) in that
  this function includes a znzTAGreadStart() which is not
  included in znzReadAutoAlignMatrix()
 */
MATRIX *znzReadMatrix(znzFile fp)
{
  char buf[MATRIX_STRLEN];
  MATRIX *M;
  char ch[100];
  int ret;
  long long here, len;

  /* no fscanf equivalent in zlib!! have to hack it */
  znzTAGreadStart(fp, &len);
  M = MatrixAlloc(4, 4, MATRIX_REAL);
  here = znztell(fp);  // not used?
  ret = znzread(buf, sizeof(unsigned char), MATRIX_STRLEN, fp);
  here = znztell(fp);  // not used?
  sscanf(buf,
         "%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
         ch,
         &(M->rptr[1][1]),
         &(M->rptr[1][2]),
         &(M->rptr[1][3]),
         &(M->rptr[1][4]),
         &(M->rptr[2][1]),
         &(M->rptr[2][2]),
         &(M->rptr[2][3]),
         &(M->rptr[2][4]),
         &(M->rptr[3][1]),
         &(M->rptr[3][2]),
         &(M->rptr[3][3]),
         &(M->rptr[3][4]),
         &(M->rptr[4][1]),
         &(M->rptr[4][2]),
         &(M->rptr[4][3]),
         &(M->rptr[4][4]));

  return (M);
}

int znzWriteAutoAlignMatrix(znzFile fp, MATRIX *M)
{
  // This does not appear to be used
  long long here;
  char buf[MATRIX_STRLEN];
  long long len;

  bzero(buf, MATRIX_STRLEN);
  sprintf(buf,
          "AutoAlign %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf",
          M->rptr[1][1],
          M->rptr[1][2],
          M->rptr[1][3],
          M->rptr[1][4],
          M->rptr[2][1],
          M->rptr[2][2],
          M->rptr[2][3],
          M->rptr[2][4],
          M->rptr[3][1],
          M->rptr[3][2],
          M->rptr[3][3],
          M->rptr[3][4],
          M->rptr[4][1],
          M->rptr[4][2],
          M->rptr[4][3],
          M->rptr[4][4]);
  len = strlen(buf);
  znzTAGwriteStart(fp, TAG_AUTO_ALIGN, &here, len);
  znzwrite(buf, sizeof(char), len, fp);

  znzTAGwriteEnd(fp, here);
  return (NO_ERROR);
}

/*
 \fn MATRIX *znzReadAutoAlignMatrix(znzFile fp)
 Compatible with znzWriteMatrix(znzFile fp, MATRIX *M)
*/
MATRIX *znzReadAutoAlignMatrix(znzFile fp)
{
  MATRIX *M;
  char buf[MATRIX_STRLEN];

  /* no fscanf equivalent in zlib!! have to hack it */
  znzread(buf, sizeof(unsigned char), MATRIX_STRLEN, fp);
  M = MatrixAlloc(4, 4, MATRIX_REAL);
  char ch[100];
  sscanf(buf,
         "%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
         ch,
         &(M->rptr[1][1]),
         &(M->rptr[1][2]),
         &(M->rptr[1][3]),
         &(M->rptr[1][4]),
         &(M->rptr[2][1]),
         &(M->rptr[2][2]),
         &(M->rptr[2][3]),
         &(M->rptr[2][4]),
         &(M->rptr[3][1]),
         &(M->rptr[3][2]),
         &(M->rptr[3][3]),
         &(M->rptr[3][4]),
         &(M->rptr[4][1]),
         &(M->rptr[4][2]),
         &(M->rptr[4][3]),
         &(M->rptr[4][4]));

  return (M);
}
