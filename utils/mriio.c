/*
 *       FILE NAME:   mriio.c
 *
 *       DESCRIPTION: utilities for reading/writing MRI data structure
 *
 *       AUTHOR:      Bruce Fischl
 *       DATE:        4/12/97
 *
*/

/*-----------------------------------------------------
                    INCLUDE FILES
-------------------------------------------------------*/
#define USE_ELECTRIC_FENCE 1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <memory.h>
#include <sys/stat.h>
#include <ctype.h>

#include "error.h"
#include "proto.h"
#include "mri.h"
#include "macros.h"
#include "diag.h"
#include "volume_io.h"
#include "region.h"
#include "machine.h"
#include "analyze.h"
#include "fio.h"
#include "utils.h"
#include "mri_identify.h"
#include "fio.h"
#include "matfile.h"
#include "math.h"

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/

#define MM_PER_METER  1000.0f
#define INFO_FNAME    "COR-.info"
#define MAX_DIM       4

#define IMAX3(a, b, c)  (a > b ? (a > c ? 0 : 2) : (b > c ? 1 : 2))

/*-----------------------------------------------------
                    STATIC DATA
-------------------------------------------------------*/

static char MIfspace[] = "frame space" ;
static int data_size[] = { 1, 4, 4, 4, 2 };

/*-----------------------------------------------------
                    STATIC PROTOTYPES
-------------------------------------------------------*/

#if 0
static void image_to_int_buffer(int *buf, MRI *mri, int slice);
static void image_to_long_buffer(long *buf, MRI *mri, int slice);
static void image_to_float_buffer(float *buf, MRI *mri, int slice);
#endif
static void short_buffer_to_image(short *buf, MRI *mri, int slice, int frame) ;
static void image_to_short_buffer(short *buf, MRI *mri, int slice);
static void int_buffer_to_image(int *buf, MRI *mri, int slice, int frame) ;
static void long_buffer_to_image(long *buf, MRI *mri, int slice, int frame) ;
static void float_buffer_to_image(float *buf, MRI *mri, int slice, int frame) ;
static void buffer_to_image(BUFTYPE *buf, MRI *mri, int slice, int frame) ;
static void image_to_buffer(BUFTYPE *buf, MRI *mri, int slice) ;
static MRI *sdtRead(char *fname, int read_volume, int fream) ;
static MRI *mncRead(char *fname, int read_volume, int frame) ;
static MRI *mghRead(char *fname, int read_volume, int frame) ;
static MRI *genesisRead(char *fname, int read_volume, int frame) ;
static MRI *siemensRead(char *fname, int read_volume, int frame) ;
static MRI *gelxRead(char *fname, int read_volume, int frame) ;
static MRI *brikRead(char *fname, int read_volume, int frame) ;
static int brikWrite(MRI *mri, char *fname, int frame) ;
static MRI *bshortRead(char *fname, int read_volume, int frame) ;
static int bshortWrite(MRI *mri, char *fname, int frame) ;
static int mghWrite(MRI *mri, char *fname, int frame) ;
static int mghAppend(MRI *mri, char *fname, int frame) ;
static int mncWrite(MRI *mri, char *fname, int frame) ;
static MRI *analyzeRead(char *fname, int read_volume, int frame) ;
static int analyzeWrite(MRI *mri, char *fname, int frame) ;
static void flipAnalyzeHeader(dsr *header);
static int write_analyze_header(char *fpref, MRI *mri) ;
static int write_analyze_mat(char *fpref, MRI *mri) ;
static int write_analyze_image(char *fpref, MRI *mri) ;
static void get_file_limits(char *format_string, int initial, int *min, int *max, int ge_check_flag);

/*-----------------------------------------------------
                    GLOBAL FUNCTIONS
-------------------------------------------------------*/
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIunpackFileName(char *inFname, int *pframe, int *ptype, char *outFname)
{
  char *number = NULL, *at = NULL, buf[STRLEN] ;
  struct stat stat_buf;

  strcpy(outFname, inFname) ;
  number = strrchr(outFname, '#') ;
  at = strrchr(outFname, '@');

  if(at)
    *at = '\0';

  if (number)   /* '#' in filename indicates frame # */
  {
    if (sscanf(number+1, "%d", pframe) < 1)
      *pframe = -1 ;
    *number = 0 ;
  }
  else
    *pframe = -1 ;

  if (at)
  {
    at = StrUpper(strcpy(buf, at+1)) ;
    if (!strcmp(at, "MNC"))
      *ptype = MRI_MINC_FILE ;
    else if (!strcmp(at, "MINC"))
      *ptype = MRI_MINC_FILE ;
    else if (!strcmp(at, "BRIK"))
      *ptype = BRIK_FILE ;
    else if (!strcmp(at, "SIEMENS"))
      *ptype = SIEMENS_FILE ;
    else if (!strcmp(at, "MGH"))
      *ptype = MRI_MGH_FILE ;
    else if (!strcmp(at, "MR"))
      *ptype = GENESIS_FILE ;
    else if (!strcmp(at, "GE"))
      *ptype = GE_LX_FILE ;
    else if (!strcmp(at, "IMG"))
      *ptype = MRI_ANALYZE_FILE ;
    else if (!strcmp(at, "COR"))
      *ptype = MRI_CORONAL_SLICE_DIRECTORY ;
    else if (!strcmp(at, "BSHORT"))
      *ptype = BSHORT_FILE;
    else if (!strcmp(at, "SDT"))
      *ptype = SDT_FILE;
    else
      ErrorExit(ERROR_UNSUPPORTED, "unknown file type %s", at);
  }
  else  /* no '@' found */
  {

    *ptype = -1;


    if(is_genesis(outFname))
      *ptype = GENESIS_FILE;
    else if(is_ge_lx(outFname))
      *ptype = GE_LX_FILE;
    else if(is_brik(outFname))
      *ptype = BRIK_FILE;
    else if(is_siemens(outFname))
      *ptype = SIEMENS_FILE;
    else if(is_analyze(outFname))
      *ptype = MRI_ANALYZE_FILE;
    else if(is_sdt(outFname))
      *ptype = SDT_FILE;
    else if(is_mgh(outFname))
      *ptype = MRI_MGH_FILE;
    else if(is_mnc(outFname))
      *ptype = MRI_MINC_FILE;
    else if(is_bshort(outFname))
      *ptype = BSHORT_FILE;
    else 
    {
      if(stat(outFname, &stat_buf) < 0)
      {
        ErrorReturn(ERROR_BADFILE, (ERROR_BAD_FILE, "can't stat file %s", outFname));
      }
      if(S_ISDIR(stat_buf.st_mode))
        *ptype = MRI_CORONAL_SLICE_DIRECTORY;
    }

    if(*ptype == -1)
      ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "unrecognized file type for file %s", outFname));

  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           read an input stream of raw bytes and organize it into
           an MRI structure
------------------------------------------------------*/
MRI *
MRIreadRaw(FILE *fp, int width, int height, int depth, int type)
{
  MRI     *mri ;
  BUFTYPE *buf ;
  int     slice, pixels ;
  int     i;

  mri = MRIalloc(width, height, depth, type) ;
  if (!mri)
    return(NULL) ;

  pixels = width*height ;
  buf = (BUFTYPE *)calloc(pixels, data_size[type]) ;

  /* every width x height pixels should be another slice */
  for (slice = 0 ; slice < depth ; slice++)
  {
    if (fread(buf, data_size[type], pixels, fp) != pixels)
      ErrorReturn(NULL,
                  (ERROR_BADFILE, "%s: could not read %dth slice (%d)",
                   Progname, slice, pixels)) ;
    if(type == 0)
      buffer_to_image(buf, mri, slice, 0) ;
    if(type == 1)
    {
      for(i = 0;i < pixels;i++)
        ((int *)buf)[i] = orderIntBytes(((int *)buf)[i]);
      int_buffer_to_image((int *)buf, mri, slice, 0);
    }
    if(type == 2)
    {
      for(i = 0;i < pixels;i++)
        ((long *)buf)[i] = orderLongBytes(((long *)buf)[i]);
      long_buffer_to_image((long *)buf, mri, slice, 0);
    }
    if(type == 3)
    {
      for(i = 0;i < pixels;i++)
        ((float *)buf)[i] = orderFloatBytes(((float *)buf)[i]);
      float_buffer_to_image((float *)buf, mri, slice, 0);
    }
    if(type == 4)
    {
      for(i = 0;i < pixels;i++)
        ((short *)buf)[i] = orderShortBytes(((short *)buf)[i]);
      short_buffer_to_image((short *)buf, mri, slice, 0);
    }
  }

  MRIinitHeader(mri) ;
  free(buf) ;
  return(mri) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Allocate an MRI data structure and read it in from data
           contained in a set of files in the directory 'fpref'
------------------------------------------------------*/
MRI *
MRIread(char *fpref)
{
  MRI     *mri ;
  int     slice, row, type, frame ;
  char    fname[2*STR_LEN] ;
  long    bytes ;
  BUFTYPE *buf ;   /* tmp space to read in whole buffer */
  FILE    *fp ;
#if USE_ELECTRIC_FENCE
  BUFTYPE  *slice_ptr ;  
#endif

  if(MRIunpackFileName(fpref, &frame, &type, fname) != NO_ERROR)
    return(NULL);

  switch (type)
  {
  case BSHORT_FILE:
    mri = bshortRead(fname, 1, frame);
    break;
  case BRIK_FILE:
    mri = brikRead(fname, 1, frame) ;
    break ;
  case SIEMENS_FILE:
    mri = siemensRead(fname, 1, frame) ;
    break ;
  case GENESIS_FILE:
    mri = genesisRead(fname, 1, frame) ;
    break ;
  case GE_LX_FILE:
    mri = gelxRead(fname, 1, frame) ;
    break ;
  case MRI_MGH_FILE:
    mri = mghRead(fname, 1, frame) ;
    if (!mri)
      return(NULL) ;
    break ;
  case MRI_MINC_FILE:
    mri = mncRead(fname, 1, frame) ;
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      printf("MNC dimensions: %d,%d,%d, slice direction: %d\n", 
             mri->xdir, mri->ydir, mri->zdir, mri->slice_direction);
    if (!mri)
      return(NULL) ;
   break ;
  case MRI_ANALYZE_FILE:
    mri = analyzeRead(fname, 1, frame) ;
    if (!mri)
      return(NULL) ;
    break ;
  case SDT_FILE:
    mri = sdtRead(fname, 1, frame) ;
    if (!mri)
      return(NULL) ;
    break;
  default:   /* coronal slice data */
    mri = MRIreadInfo(fname) ;
    if (!mri)
      return(NULL) ;
    MRIallocIndices(mri) ;
    mri->slices = (BUFTYPE ***)calloc(mri->depth, sizeof(BUFTYPE **)) ;
    if (!mri->slices)
      ErrorExit(ERROR_NO_MEMORY, 
                "MRIread(%s): could not allocate %d slices\n", mri->depth) ;
    
    bytes = (long)mri->width * (long)mri->height ;
    buf = (BUFTYPE *)calloc(bytes, sizeof(BUFTYPE)) ;
    if (!buf)
      ErrorExit(ERROR_NO_MEMORY, 
                "MRIread(%s): could not allocate %d bytes for buf\n",
                fpref,bytes);
    
    for (slice = 0 ; slice < mri->depth ; slice++)
    {
      /* allocate pointer to array of rows */
      mri->slices[slice] = (BUFTYPE **)calloc(mri->height, sizeof(BUFTYPE *)) ;
      if (!mri->slices[slice])
        ErrorExit(ERROR_NO_MEMORY, 
                  "MRIread(%s): could not allocate %d bytes for %dth slice\n",
                  fpref, mri->height*sizeof(BUFTYPE *), slice) ;
      
#if USE_ELECTRIC_FENCE

      if (slice == 255)
        DiagBreak() ;
      slice_ptr = (BUFTYPE *)calloc(mri->height * mri->width, sizeof(BUFTYPE));
      if (!slice_ptr)
        ErrorExit(ERROR_NO_MEMORY, 
                  "MRIread(%s): could not allocate %dth slice\n",fpref, slice);
      for (row = 0 ; row < mri->height ; row++)
      {
        mri->slices[slice][row] = slice_ptr+row*mri->width ;
      }
#else
      /* allocate each row */
      for (row = 0 ; row < mri->height ; row++)
      {
        mri->slices[slice][row] = 
          (BUFTYPE *)calloc(mri->width, sizeof(BUFTYPE));
        if (!mri->slices[slice][row])
          ErrorExit(ERROR_NO_MEMORY, 
                    "MRIread(%s): could not allocate %dth row in %dth slice\n",
                    fpref, slice, row) ;
      }
#endif

      sprintf(fname,"%s/COR-%03d", fpref, slice+mri->imnr0) ;

      fp = fopen(fname, "rb") ;
      if (!fp)
      {
        MRIfree(&mri) ;
        ErrorReturn(NULL, 
                    (ERROR_NO_FILE,
                     "MRIread(%s): could not open slice file '%s'",
                     fpref, fname)) ;       
      }
      if (fread(buf, sizeof(BUFTYPE), bytes, fp) != bytes)
      {
        fclose(fp) ;
        MRIfree(&mri) ;
        ErrorReturn(NULL, 
                    (ERROR_NO_FILE,
                     "MRIread(%s): could not read slice file '%s'",
                     fpref, fname)) ;
      }
      fclose(fp) ;
      buffer_to_image(buf, mri, slice, 0) ;
    }
    free(buf) ;
    break ;
  }

  if(mri == NULL)
    return(NULL);

  if (!MRIisValid(mri))
    MRIflipByteOrder(mri, mri) ;

  return(mri) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Allocate an MRI data structure and read in its
           header information from the file COR-.info in the
           directory specified by 'fpref'
------------------------------------------------------*/
MRI *
MRIreadInfo(char *fpref)
{
  MRI     *mri ;
  FILE    *fp;
  char    cmd[STRLEN], fname[STRLEN];
  int     imnr0, imnr1, width, height, depth, ptype, type, frame ;
  char    info_line[STRLEN];

  MRIunpackFileName(fpref, &frame, &type, fname) ;
  switch (type)
  {
  case BSHORT_FILE:
    mri = bshortRead(fname, 0, frame) ;
    if(!mri)
      return(NULL);
    break;
  case BRIK_FILE:
    mri = brikRead(fname, 0, frame) ;
    if(!mri)
      return(NULL);
    break ;
  case SIEMENS_FILE:
    mri = siemensRead(fname, 0, frame) ;
    if(!mri)
      return(NULL);
    break ;
  case GENESIS_FILE:
    mri = genesisRead(fname, 0, frame) ;
    if (!mri)
      return(NULL) ;
    break ;
  case GE_LX_FILE:
    mri = gelxRead(fname, 0, frame) ;
    if (!mri)
      return(NULL) ;
    break ;
  case MRI_MGH_FILE:
    mri = mghRead(fname, 0, frame) ;
    if (!mri)
      return(NULL) ;
    break ;
  case MRI_MINC_FILE:
    mri = mncRead(fname, 0, frame) ;
    if (!mri)
      return(NULL) ;
    break ;
  case MRI_ANALYZE_FILE:
    mri = analyzeRead(fname, 0, frame) ;
    if (!mri)
      return(NULL) ;
    break ;
  case SDT_FILE:
    mri = sdtRead(fname, 0, frame) ;
    if (!mri)
      return(NULL) ;
    break ;
  default:   /* coronal slice data */
    strcpy(fpref, fname);
    sprintf(fname,"%s/%s",fpref, INFO_FNAME);
    fp = fopen(fname,"r");
    if (fp == NULL) 
      ErrorReturn(NULL,
                  (ERROR_NO_FILE, "File %s not found.\n",fname));
    
    fscanf(fp,"%*s %d",&imnr0);
    fscanf(fp,"%*s %d",&imnr1);
    fscanf(fp,"%*s %d",&ptype);
    fscanf(fp,"%*s %d",&width);
    fscanf(fp,"%*s %d",&height);
    depth = imnr1-imnr0+1;
    mri = MRIallocHeader(width, height, depth, MRI_UCHAR) ;
    if (!mri)
      ErrorExit(ERROR_NO_MEMORY, "MRIreadInfo: could not allocate MRI\n") ;
    mri->imnr0 = imnr0 ;
    mri->imnr1 = imnr1 ;
    mri->ptype = ptype ;
    strncpy(mri->fname, fpref, STR_LEN) ;

    fscanf(fp,"%*s %f",&mri->fov);
    fscanf(fp,"%*s %f",&mri->thick);
    fscanf(fp,"%*s %f",&mri->ps);
    fscanf(fp,"%*s %f",&mri->location); /* locatn */
    fscanf(fp,"%*s %f",&mri->xstart); /* strtx */
    fscanf(fp,"%*s %f",&mri->xend); /* endx */
    fscanf(fp,"%*s %f",&mri->ystart); /* strty */
    fscanf(fp,"%*s %f",&mri->yend); /* endy */
    fscanf(fp,"%*s %f",&mri->zstart); /* strtz */
    fscanf(fp,"%*s %f",&mri->zend); /* endz */
    fscanf(fp,"%*s %f",&mri->tr); 
    fscanf(fp,"%*s %f",&mri->te); 
    fscanf(fp,"%*s %f",&mri->ti);

    mri->linear_transform = NULL;
    mri->ras_good_flag = 0;

    while(fgets(info_line, STRLEN, fp) != NULL)
    {
      if(strncmp(info_line, "ras_good_flag", 13) == 0)
        sscanf(info_line, "%*s %d", &mri->ras_good_flag);
      else if(strncmp(info_line, "x_ras", 5) == 0)
        sscanf(info_line, "%*s %f %f %f", &mri->x_r, &mri->x_a, &mri->x_s);
      else if(strncmp(info_line, "y_ras", 5) == 0)
        sscanf(info_line, "%*s %f %f %f", &mri->y_r, &mri->y_a, &mri->y_s);
      else if(strncmp(info_line, "z_ras", 5) == 0)
        sscanf(info_line, "%*s %f %f %f", &mri->z_r, &mri->z_a, &mri->z_s);
      else if(strncmp(info_line, "c_ras", 5) == 0)
        sscanf(info_line, "%*s %f %f %f", &mri->c_r, &mri->c_a, &mri->c_s);
      else if(strncmp(info_line, "xform", 5) == 0 || strncmp(info_line, "transform", 9) == 0)
      {

        sscanf(info_line, "%s %s", cmd, mri->transform_fname);

        if (*mri->transform_fname != '/') /* relative path, add prefix */
          sprintf(fname, "%s/%s", fpref, mri->transform_fname) ;
        else
          strcpy(fname, mri->transform_fname) ; /* absolute path */
        FileNameAbsolute(fname, mri->transform_fname) ;
        if (!FileExists(mri->transform_fname))  /* try typical location */
          sprintf(mri->transform_fname,"%s/../transforms/talairach.xfm",fpref);

        if (FileExists(mri->transform_fname))
        {
          if (input_transform_file(mri->transform_fname, &mri->transform)!=OK)
            ErrorPrintf(ERROR_NO_MEMORY, 
                        "MRIreadInfo: could not read xform file '%s'\n", 
                        mri->transform_fname) ;
          else
          {
            mri->linear_transform = get_linear_transform_ptr(&mri->transform) ;
            mri->inverse_linear_transform = 
              get_inverse_linear_transform_ptr(&mri->transform) ;
            mri->free_transform = 1 ;
          }
        }
        else
        {
          mri->linear_transform = NULL ;
#if 0
          ErrorPrintf(ERROR_NO_MEMORY, 
                      "MRIreadInfo: could not read xform file '%s'\n", 
                      mri->transform_fname) ;
#endif
        }
    
      }
      else
      {
      /* do nothing -- toro has been known to return a blank line in the first pass */
      }
    }

    mri->fov *= MM_PER_METER ;
    mri->ps *= MM_PER_METER ;
    mri->thick *= MM_PER_METER ;
    mri->xstart *= MM_PER_METER ;
    mri->xend *= MM_PER_METER ;
    mri->ystart *= MM_PER_METER ;
    mri->yend *= MM_PER_METER ;
    mri->zstart *= MM_PER_METER ;
    mri->zend *= MM_PER_METER ;
    mri->ysize = mri->xsize = mri->ps ;
    mri->zsize = mri->thick ;
    fclose(fp);
    mri->depth = mri->imnr1-mri->imnr0+1;
    break ;
  }

  mri->slice_direction = MRI_CORONAL ;  /* should read from data file */
  return(mri) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Write an MRI header and a set of data files to
          the directory specified by 'fpref'
------------------------------------------------------*/
int
MRIwrite(MRI *mri, char *fpref)
{
  int      slice, type, frame ;
  BUFTYPE  *buf ;
  long     bytes ;
  char     fname[2*STR_LEN] ;
  FILE     *fp ;

  MRIunpackFileName(fpref, &frame, &type, fname) ;
  if (type == MRI_MINC_FILE)
    return(mncWrite(mri, fname, frame)) ;
  else if (type == MRI_ANALYZE_FILE)
    return(analyzeWrite(mri, fname, frame)) ;
  else if (type == MRI_MGH_FILE)
    return(mghWrite(mri, fname, frame)) ;
  else if (type == BRIK_FILE)
    return(brikWrite(mri, fname, frame)) ;
  else if (type == BSHORT_FILE)
    return(bshortWrite(mri, fname, frame)) ;
  else if (type == GENESIS_FILE)
  {
    ErrorReturn(0, (ERROR_UNSUPPORTED, "Genesis write unsupported"));
  }
  else if(type == GE_LX_FILE)
  {
    ErrorReturn(0, (ERROR_UNSUPPORTED, "GE LX write unsupported"));
  }
  else if(type == SIEMENS_FILE)
  {
    ErrorReturn(0, (ERROR_UNSUPPORTED, "Siemens write unsupported"));
  }

  MRIwriteInfo(mri, fpref) ;
  bytes = (long)mri->width * (long)mri->height ;
  buf = (BUFTYPE *)calloc(bytes, sizeof(BUFTYPE)) ;
  if (!buf)
    ErrorExit(ERROR_NO_MEMORY, 
              "MRIwrite(%s): could not allocate %d bytes for buf",fpref,bytes);

  for (slice = 0 ; slice < mri->depth ; slice++)
  {
    image_to_buffer(buf, mri, slice) ;
    sprintf(fname,"%s/COR-%03d", fpref, slice+mri->imnr0) ;
    fp = fopen(fname, "wb") ;
    if (!fp)
      ErrorReturn(ERROR_NO_FILE, 
                (ERROR_NO_FILE, "MRIwrite(%s): could not open slice file '%s'",
                  fpref, fname)) ;       

    if (fwrite(buf, sizeof(BUFTYPE), bytes, fp) != bytes)
    {
      fclose(fp) ;
      ErrorReturn(ERROR_BAD_FILE, 
              (ERROR_BAD_FILE, "MRIwrite(%s): could not write slice file '%s'",
                  fpref, fname)) ;
    }
    fclose(fp) ;
  }

  free(buf) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Write an MRI header and a set of data files to
          the directory specified by 'fpref'
------------------------------------------------------*/
int
MRIappend(MRI *mri, char *fpref)
{
  int      type, frame ;
  char     fname[STRLEN] ;

  MRIunpackFileName(fpref, &frame, &type, fname) ;
  if (type == MRI_MGH_FILE)
    return(mghAppend(mri, fname, frame)) ;
  else
    ErrorReturn(ERROR_UNSUPPORTED, 
                (ERROR_UNSUPPORTED, "MRIappend(%s): file type not supported",
                fname)) ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Write the MRI header information to the file
           COR-.info in the directory specified by 'fpref'
------------------------------------------------------*/
int
MRIwriteInfo(MRI *mri, char *fpref)
{
  FILE    *fp;
  char    fname[STRLEN];

  sprintf(fname,"%s/%s",fpref, INFO_FNAME);
  fp = fopen(fname,"w");
  if (fp == NULL) 
    ErrorReturn(ERROR_NO_FILE, 
                (ERROR_NO_FILE,
                 "MRIwriteInfo(%s): could not open %s.\n", fpref, fname)) ;

  fprintf(fp, "%s %d\n", "imnr0", mri->imnr0);
  fprintf(fp, "%s %d\n", "imnr1", mri->imnr1);
  fprintf(fp, "%s %d\n", "ptype", 
          mri->slice_direction == MRI_CORONAL ? 2 :
          mri->slice_direction == MRI_HORIZONTAL ? 0 : 1) ;
  fprintf(fp, "%s %d\n", "x", mri->width);
  fprintf(fp, "%s %d\n", "y", mri->height);
  fprintf(fp, "%s %f\n", "fov", mri->fov/MM_PER_METER);
  fprintf(fp, "%s %f\n", "thick", mri->ps/MM_PER_METER);
  fprintf(fp, "%s %f\n", "psiz", mri->ps/MM_PER_METER);
  fprintf(fp, "%s %f\n", "locatn", mri->location); /* locatn */
  fprintf(fp, "%s %f\n", "strtx", mri->xstart/MM_PER_METER); /* strtx */
  fprintf(fp, "%s %f\n", "endx", mri->xend/MM_PER_METER); /* endx */
  fprintf(fp, "%s %f\n", "strty", mri->ystart/MM_PER_METER); /* strty */
  fprintf(fp, "%s %f\n", "endy", mri->yend/MM_PER_METER); /* endy */
  fprintf(fp, "%s %f\n", "strtz", mri->zstart/MM_PER_METER); /* strtz */
  fprintf(fp, "%s %f\n", "endz", mri->zend/MM_PER_METER); /* endz */
  fprintf(fp, "%s %f\n", "tr", mri->tr) ;
  fprintf(fp, "%s %f\n", "te", mri->te) ;
  fprintf(fp, "%s %f\n", "ti", mri->ti) ;
  if (mri->linear_transform)
  {
    char fname[STRLEN] ;

/* 
   this won't work for relative paths which are not the same for the
   destination directory as for the the source directory.
   */
    sprintf(fname,"%s", mri->transform_fname);
    fprintf(fp, "xform %s\n", fname) ;

#if 0
    /* doesn't work - I don't know why */
    if (output_transform_file(fname, "talairach xfm", &mri->transform) != OK)
      ErrorPrintf(ERROR_BADFILE, "MRIwriteInfo(%s): xform write failed",fpref);
#endif

  }

  fprintf(fp, "%s %d\n", "ras_good_flag", mri->ras_good_flag);
  fprintf(fp, "%s %f %f %f\n", "x_ras", mri->x_r, mri->x_a, mri->x_s);
  fprintf(fp, "%s %f %f %f\n", "y_ras", mri->y_r, mri->y_a, mri->y_s);
  fprintf(fp, "%s %f %f %f\n", "z_ras", mri->z_r, mri->z_a, mri->z_s);
  fprintf(fp, "%s %f %f %f\n", "c_ras", mri->c_r, mri->c_a, mri->c_s);

  fclose(fp);

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Copy the data from a flat buffer to an array
           of rows in the MRI data structure (after reading)
------------------------------------------------------*/
int
MRIfileType(char *fname)
{
  int file_type = MRI_CORONAL_SLICE_DIRECTORY ;
  char *dot, buf[STRLEN], *number ;

  if (*fname == '@')
    return(LIST_FILE) ;

  strcpy(buf, fname) ;
  dot = strrchr(buf, '@') ;
  number = strchr(buf, '#') ;
  if (number)
    *number = 0 ;  /* don't consider : part of extension */
  if (!dot)
    dot = strrchr(buf, '.') ;

  if (dot)
  {
    dot++ ;
    StrUpper(buf) ;
    if (!strcmp(dot, "LST"))
      return(LIST_FILE) ;
    else if (!strcmp(dot, "MNC"))
      return(MRI_MINC_FILE) ;
    else if (!strcmp(dot, "MGH"))
      return(MRI_MGH_FILE) ;
  }

  return(file_type) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIfromVolume(Volume volume, int start_frame, int end_frame)
{
  MRI   *mri ;
  int   type, width, height, depth, x, y, z, ystep, y1, ndim, nframes,
        sizes[MAX_DIM], frame ;
  Real  separations[MAX_DIM] ;

/*
   the MNC coordinate system is related to coronal slices in the following way:

   MNC   CORONAL
   z  -->  y
   y  -->  z
   x  -->  x
*/
  ndim = get_volume_n_dimensions(volume) ;
  get_volume_sizes(volume, sizes) ;
  if (ndim > 3)
    nframes = sizes[3] ;
  else
    nframes = 1 ;
  if (start_frame < 0)
    start_frame = 0 ;
  else if (start_frame >= nframes)
    start_frame = nframes-1 ;
  if (end_frame >= nframes)
    end_frame = nframes - 1 ;

  width = sizes[0] ;
  height = sizes[2] ;
  depth = sizes[1] ;
  get_volume_separations(volume, separations) ;
  ystep = nint(separations[2]) ;
  switch (volume->nc_data_type)
  {
  case NC_BYTE:  type = MRI_UCHAR ; break ;
  case NC_SHORT: type = MRI_SHORT ; break ;
  case NC_FLOAT: type = MRI_FLOAT ; break ;
  default:  
    ErrorReturn(NULL, 
                (ERROR_UNSUPPORTED, "MRIfromVolume: unsupported MNC type %d",
                 volume->nc_data_type)) ;
    break ;
  }
  if (volume_is_alloced(volume))
    mri = MRIallocSequence(width, height, depth, type,end_frame-start_frame+1);
  else
    mri = MRIallocHeader(width, height, depth, type) ;

  if (!mri)
    ErrorExit(ERROR_NO_MEMORY, "MRIfromVolume: could not allocate MRI\n") ;

  if (volume_is_alloced(volume)) switch (type)
  {
  case MRI_UCHAR:
    for (frame = start_frame ; frame <= end_frame ; frame++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        for (y = 0 ; y < height ; y++)
        {
          if (ystep < 0)
            y1 = y ;
          else
            y1 = height - y - 1 ;
          for (x = 0 ; x < width ; x++)
          {
            if (ndim <= 3)
              MRIvox(mri, x, y1, z) = 
                GET_MULTIDIM_TYPE_3D(volume->array, BUFTYPE, x, z, y) ;
            else
              MRIseq_vox(mri, x, y1, z, frame-start_frame) = 
                GET_MULTIDIM_TYPE_4D(volume->array, BUFTYPE, x, z, y, frame) ;
          }
        }
      }
    }
    break ;
  case MRI_SHORT:
    for (frame = start_frame ; frame <= end_frame ; frame++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        for (y = 0 ; y < height ; y++)
        {
          if (ystep < 0)
            y1 = y ;
          else
            y1 = height - y - 1 ;
          for (x = 0 ; x < width ; x++)
          {
            if (ndim <= 3)
              MRISvox(mri, x, y1, z) = 
                GET_MULTIDIM_TYPE_3D(volume->array, short, x, z, y) ;
            else
              MRISseq_vox(mri, x, y1, z, frame-start_frame) = 
                GET_MULTIDIM_TYPE_4D(volume->array, short, x, z, y, frame) ;
          }
        }
      }
    }
    break ;
  case MRI_FLOAT:
    for (frame = start_frame ; frame <= end_frame ; frame++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        for (y = 0 ; y < height ; y++)
        {
          if (ystep < 0)
            y1 = y ;
          else
            y1 = height - y - 1 ;
          for (x = 0 ; x < width ; x++)
          {
            if (ndim <= 3)
              MRIFvox(mri, x, y1, z) = 
                GET_MULTIDIM_TYPE_3D(volume->array, float, x, z, y) ;
            else
              MRIFseq_vox(mri, x, y1, z, frame-start_frame) = 
                GET_MULTIDIM_TYPE_4D(volume->array, float, x, z, y, frame) ;
          }
        }
      }
    }
    break ;
  default:
    break ;
  }
  mri->slice_direction = MRI_CORONAL ;
  mri->imnr0 = 0 ;
  mri->imnr1 = depth - 1 ;
  mri->xsize = separations[0] ;
  mri->ysize = separations[2] ;
  mri->zsize = separations[1] ;
  mri->xstart = (float)volume->world_space_for_translation_voxel[0] ;
  mri->zstart = (float)volume->world_space_for_translation_voxel[1] ;
  mri->ystart = (float)volume->world_space_for_translation_voxel[2] ;

  if (ystep < 0)
    mri->yinvert = -1 ;
  else
    mri->yinvert = 1 ;

  mri->xend = mri->xstart + width * mri->xsize ;
  mri->yend = mri->ystart + height * mri->ysize ;
  mri->zend = mri->zstart + depth * mri->zsize ;
  mri->ps = mri->thick = mri->zsize ; 
  mri->fov = (mri->xend - mri->xstart) / mri->xsize ;
  return(mri) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
Volume
MRItoVolume(MRI *mri)
{
  Volume        volume ;
  int           ndim, width, height, depth, type, sgned, x, y, z,
                sizes[MAX_DIM], y1, frame ;
  char          *dim_names[MAX_DIM] ;
  Real          separations[MAX_DIM], voxel[MAX_DIM], world_vox[MAX_DIM] ;
  
  if (mri->slice_direction != MRI_CORONAL)
    ErrorReturn(NULL, 
                (ERROR_UNSUPPORTED,
                 "MRItoVolume: unsupported slice direction %d", 
                 mri->slice_direction)) ;

  ndim = 4 ;   /* 3 spatial plus one for multiple frames */

  dim_names[0] = MIyspace ;
  dim_names[1] = MIzspace ;
  dim_names[2] = MIxspace ;
  dim_names[3] = MIfspace ;
  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  switch (mri->type)
  {
  case MRI_UCHAR: type = NC_BYTE ;  sgned = 0 ; break ;
  case MRI_FLOAT: type = NC_FLOAT ; sgned = 1 ; break ;
  default:
    ErrorReturn(NULL, 
                (ERROR_UNSUPPORTED, "MRItoVolume: unsupported MRI type %d",
                 mri->type)) ;
    break ;
  }

  volume = create_volume(ndim,dim_names,type, sgned, 0.0, 0.0) ;
  sizes[0] = depth ;
  sizes[1] = height ;
  sizes[2] = width ;
  sizes[3] = mri->nframes ;
  set_volume_sizes(volume, sizes) ;
  alloc_volume_data(volume) ;
#if 0
  if (mri->linear_transform)
    set_voxel_to_world_transform(volume, &mri->transform) ;
#endif
  separations[0] = mri->zsize ;
  separations[1] = (float)mri->yinvert * mri->ysize ;
  separations[2] = mri->xsize ;
  separations[3] = 1 ;

  set_volume_separations(volume, separations) ;

/* 
   set relationship between start of voxel coordinates and
   world space origin.
   */
  voxel[0] = voxel[1] = voxel[2] = voxel[3] = 0.0 ;
  world_vox[0] = mri->xstart ;
  world_vox[1] = mri->zstart ;
  world_vox[2] = (float)mri->yinvert * mri->ystart ;
  world_vox[3] = 0.0 ;
  set_volume_translation(volume, voxel, world_vox) ;
  switch (type)
  {
  case NC_BYTE:
    for (frame = 0 ; frame < mri->nframes ; frame++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        for (y = 0 ; y < height ; y++)
        {
          if (mri->yinvert)
            y1 = height - y - 1 ;
          else
            y1 = y ;
          for (x = 0 ; x < width ; x++)
          {
            GET_MULTIDIM_TYPE_4D(volume->array, BUFTYPE, z, y, x, frame) =
              MRIseq_vox(mri, x, y1, z, frame) ; 
          }
        }
      }
    }
    break ;
  case NC_FLOAT:
    for (frame = 0 ; frame < mri->nframes ; frame++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        for (y = 0 ; y < height ; y++)
        {
          if (mri->yinvert)
            y1 = height - y - 1 ;
          else
            y1 = y ;
          for (x = 0 ; x < width ; x++)
          {
            GET_MULTIDIM_TYPE_4D(volume->array, float, z, y, x, frame) =
              MRIFseq_vox(mri, x, y1, z, frame) ;
          }
        }
      }
    }
    break ;
  default:
    break ;
  }

  return(volume) ;
}
/********************** static utility functions ****************************/


/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Copy the data from an array of rows in the MRI
           data structure to a flat buffer (preparatory to
           writing)
------------------------------------------------------*/
static void
image_to_buffer(BUFTYPE *buf, MRI *mri, int slice)
{
  int           y, width, height ;
  BUFTYPE       *pslice ;
  
  width = mri->width ;
  height = mri->height ;
  for (y=0; y < height ; y++)
  {
    pslice = mri->slices[slice][y] ;
    memcpy(buf, pslice, width*sizeof(BUFTYPE)) ;
    buf += width ;
  }
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Copy the data from a flat buffer to an array
           of rows in the MRI data structure (after reading)
------------------------------------------------------*/
static void
buffer_to_image(BUFTYPE *buf, MRI *mri, int slice, int frame)
{
  int           y, width, height ;
  BUFTYPE       *pslice ;
  
  width = mri->width ;
  height = mri->height ;
  for (y = 0 ; y < height ; y++)
  {
    pslice = &MRIseq_vox(mri, 0, y, slice, frame) ;
    memcpy(pslice, buf, width*sizeof(BUFTYPE)) ;
    buf += width ;
  }
}


static void
int_buffer_to_image(int *buf, MRI *mri, int slice, int frame)
{
  int           y, width, height ;
  int           *pslice ;
  
  width = mri->width ;
  height = mri->height ;
  for (y = 0 ; y < height ; y++)
  {
    pslice = &MRIIseq_vox(mri, 0, y, slice, frame) ;
    memcpy(pslice, buf, width*sizeof(int)) ;
    buf += width ;
  }
}

#if 0
static void
image_to_int_buffer(int *buf, MRI *mri, int slice)
{
  int y, x, width, height, depth ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth;
  for (y=0; y < height ; y++)
  {
    if(mri->type == MRI_UCHAR)
    {
      for(x = 0;x < depth;x++)
        buf[x] = (int)MRIvox(mri, x, y, slice);
    }
    else if(mri->type == MRI_SHORT)
    {
      for(x = 0;x < depth;x++)
        buf[x] = (int)MRISvox(mri, x, y, slice);
    }
    else if(mri->type == MRI_LONG)
    {
      for(x = 0;x < depth;x++)
        buf[x] = (int)MRILvox(mri, x, y, slice);
    }
    else if(mri->type == MRI_FLOAT)
    {
      for(x = 0;x < depth;x++)
        buf[x] = (int)MRIFvox(mri, x, y, slice);
    }
    else
    {
      memcpy(buf, mri->slices[slice][y], width*sizeof(int)) ;
    }

    buf += width ;
  }
}
static void
image_to_long_buffer(long *buf, MRI *mri, int slice)
{
  int y, x, width, height, depth ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth;
  for (y=0; y < height ; y++)
  {
    if(mri->type == MRI_UCHAR)
    {
      for(x = 0;x < depth;x++)
        buf[x] = (long)MRIvox(mri, x, y, slice);
    }
    else if(mri->type == MRI_INT)
    {
      for(x = 0;x < depth;x++)
        buf[x] = (long)MRIIvox(mri, x, y, slice);
    }
    else if(mri->type == MRI_SHORT)
    {
      for(x = 0;x < depth;x++)
        buf[x] = (long)MRISvox(mri, x, y, slice);
    }
    else if(mri->type == MRI_FLOAT)
    {
      for(x = 0;x < depth;x++)
        buf[x] = (long)MRIFvox(mri, x, y, slice);
    }
    else
    {
      memcpy(buf, mri->slices[slice][y], width*sizeof(long)) ;
    }

    buf += width ;
  }
}
static void
image_to_float_buffer(float *buf, MRI *mri, int slice)
{
  int y, x, width, height, depth ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth;
  for (y=0; y < height ; y++)
  {
    if(mri->type == MRI_UCHAR)
    {
      for(x = 0;x < depth;x++)
        buf[x] = (float)MRIvox(mri, x, y, slice);
    }
    else if(mri->type == MRI_INT)
    {
      for(x = 0;x < depth;x++)
        buf[x] = (float)MRIIvox(mri, x, y, slice);
    }
    else if(mri->type == MRI_LONG)
    {
      for(x = 0;x < depth;x++)
        buf[x] = (float)MRILvox(mri, x, y, slice);
    }
    else if(mri->type == MRI_SHORT)
    {
      for(x = 0;x < depth;x++)
        buf[x] = (float)MRISvox(mri, x, y, slice);
    }
    else
    {
      memcpy(buf, mri->slices[slice][y], width*sizeof(float)) ;
    }

    buf += width ;
  }
}
#endif

static void
long_buffer_to_image(long *buf, MRI *mri, int slice, int frame)
{
  int           y, width, height ;
  long          *pslice ;
  
  width = mri->width ;
  height = mri->height ;
  for (y = 0 ; y < height ; y++)
  {
    pslice = &MRILseq_vox(mri, 0, y, slice, frame) ;
    memcpy(pslice, buf, width*sizeof(long)) ;
    buf += width ;
  }
}


static void
float_buffer_to_image(float *buf, MRI *mri, int slice, int frame)
{
  int           y, width, height ;
  float         *pslice ;
  
  width = mri->width ;
  height = mri->height ;
  for (y = 0 ; y < height ; y++)
  {
    pslice = &MRIFseq_vox(mri, 0, y, slice, frame) ;
    memcpy(pslice, buf, width*sizeof(float)) ;
    buf += width ;
  }
}

static void
image_to_short_buffer(short *buf, MRI *mri, int slice)
{
  int y, x, width, height, depth ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth;
  for (y=0; y < height ; y++)
  {
    if(mri->type == MRI_UCHAR)
    {
      for(x = 0;x < depth;x++)
        buf[x] = (short)MRIvox(mri, x, y, slice);
    }
    else if(mri->type == MRI_INT)
    {
      for(x = 0;x < depth;x++)
        buf[x] = (short)MRIIvox(mri, x, y, slice);
    }
    else if(mri->type == MRI_LONG)
    {
      for(x = 0;x < depth;x++)
        buf[x] = (short)MRILvox(mri, x, y, slice);
    }
    else if(mri->type == MRI_FLOAT)
    {
      for(x = 0;x < depth;x++)
        buf[x] = (short)MRIFvox(mri, x, y, slice);
    }
    else
    {
      memcpy(buf, mri->slices[slice][y], width*sizeof(short)) ;
    }

    buf += width ;
  }
}

static void
short_buffer_to_image(short *buf, MRI *mri, int slice, int frame)
{
  int           y, width, height ;
  short         *pslice ;
  
  width = mri->width ;
  height = mri->height ;
  for (y = 0 ; y < height ; y++)
  {
    pslice = &MRISseq_vox(mri, 0, y, slice, frame) ;
    memcpy(pslice, buf, width*sizeof(short)) ;
    buf += width ;
  }
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Read in a volume in MNC format. If read_volume is true,
           then read in the volume data, otherwise only read the
           header.
------------------------------------------------------*/
static MRI *
mncRead(char *fname, int read_volume, int frame)
{
  MRI                 *mri ;
  char                *dim_names[MAX_DIM] ;
  Volume              volume ;
  int                 error, start_frame, end_frame, sizes[MAX_DIM], nframes,
                      ndim ;
  volume_input_struct input_info ;

  dim_names[0] = MIxspace ;
  dim_names[1] = MIyspace ;
  dim_names[2] = MIzspace ;
  dim_names[3] = MIfspace ;

  if (read_volume)
    error = input_volume(fname, 0, dim_names, NC_UNSPECIFIED,
                         TRUE, 0, 0, TRUE, &volume, NULL) ;
  else
    error = start_volume_input(fname, 0, dim_names, NC_UNSPECIFIED,
                         TRUE, 0, 0, TRUE, &volume, NULL, &input_info) ;
  if (error)
    return(NULL) ;

  get_volume_sizes(volume, sizes) ;
  ndim = get_volume_n_dimensions(volume) ;
  if (ndim < 4)
    nframes = 1 ;
  else
    nframes = sizes[3] ;
  if (frame < 0)
  {
    start_frame = 0 ;
    end_frame = nframes-1 ;
  }
  else
  {
    if (frame >= nframes)
    {
      delete_volume(volume) ;
      ErrorReturn(NULL,
                  (ERROR_BADPARM, 
                   "mncRead: specified frame %d out of bounds",frame));
    }
    start_frame = end_frame = frame ;
  }
  mri = MRIfromVolume(volume, start_frame, end_frame) ;
  strncpy(mri->fname, fname, STR_LEN-1) ;
  delete_volume(volume) ;
  return(mri) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
mncWrite(MRI *mri, char *fname, int frame)
{
  Volume        volume ;
  int           type, sgned, error, start                e ;

  if (frame < 0)
  {
    start_frame = 0 ;
    end_frame = mri->nframes-1 ;
  }
  else
    start_frame = end_frame = frame ;

  if (mri->slice_direction != MRI_CORONAL)
    ErrorReturn(ERROR_UNSUPPORTED, 
                (ERROR_UNSUPPORTED,"mncWrite: unsupported slice direction %d", 
                mri->slice_direction)) ;

  switch (mri->type)
  {
  case MRI_UCHAR: type = NC_BYTE ;  sgned = 0 ; break ;
  case MRI_FLOAT: type = NC_FLOAT ; sgned = 1 ; break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED, 
                (ERROR_UNSUPPORTED, "mncWrite: unsupported MRI type %d",
                 mri->type)) ;
    break ;
  }

  volume = MRItoVolume(mri) ;
  if (!volume)
    return(Gerror) ;

  error = output_volume(fname, type, sgned, 0.0, 0.0, volume,"mncWrite", NULL);
  if (error)
  {
    delete_volume(volume) ;
    ErrorReturn(error, (error, "mncWrite: output volume failed")) ;
  }

  delete_volume(volume) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIisValid(MRI *mri)
{
  long   total, bad ;
  float  *fpix ;
  double exponent, val ;
  int    x, y, z, width, height, depth, frame ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  total = bad = 0 ;
  for (frame = 0 ; frame < mri->nframes ; frame++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        switch (mri->type)
        {
        case MRI_FLOAT:
          fpix = &MRIFseq_vox(mri, 0, y, z, frame) ;
          for (x = 0 ; x < width ; x++)
          {
            val = *fpix++ ;
            if (val == 0.0)
              continue ;
            total++ ;
            exponent = log10(fabs(val)) ;
            if (exponent > 10.0)   /* any values this big are indicative */
              return(0) ;
            
            if ((exponent > 6.0) || (exponent < -20))
              bad++ ;
            break ;
          default:
            return(1) ;
            break ;
          }
        }
      }
    }
  }

  if (total)
  {
    float pct ;

    pct = (float)bad / (float)total ;
    return(pct < 0.50f) ;   /* less than 50% of non-zero pixels bad */
  }

  return(1) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIflipByteOrder(MRI *mri_src, MRI *mri_dst)
{
  float  *spix, fval, *dpix ;
  int    x, y, z, width, height, depth, frame ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  for (frame = 0 ; frame < mri_src->nframes ; frame++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        switch (mri_src->type)
        {
        case MRI_FLOAT:
          spix = &MRIFseq_vox(mri_src, 0, y, z, frame) ;
          dpix = &MRIFseq_vox(mri_dst, 0, y, z, frame) ;
          for (x = 0 ; x < width ; x++)
          {
            switch (mri_src->type)
            {
            case PFFLOAT:
              fval = *spix++ ;
              *dpix++ = swapFloat(fval) ;
              break ;
            default:
              ErrorReturn(NULL, 
                          (ERROR_UNSUPPORTED, 
                           "MriFlipBytes: unsupported type %d\n", 
                           mri_src->type)) ;
              break ;
            }
          }
        }
      }
    }
  } return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
analyzeWrite(MRI *mri, char *fname, int frame)
{
  int   type, start_frame, end_frame ;

  if (frame < 0)
  {
    start_frame = 0 ;
    end_frame = mri->nframes-1 ;
  }
  else
    start_frame = end_frame = frame ;

/*
  if (mri->slice_direction != MRI_CORONAL)
    ErrorReturn(ERROR_UNSUPPORTED, 
            (ERROR_UNSUPPORTED,"analyzeWrite: unsupported slice direction %d", 
             mri->slice_direction)) ;
*/

  switch (mri->type)
  {
  case MRI_UCHAR: type = DT_UNSIGNED_CHAR ; break ;
  case MRI_FLOAT: type = DT_FLOAT ;         break ;
  case MRI_INT:   type = DT_SIGNED_INT ;    break ;
  case MRI_SHORT: type = DT_SIGNED_SHORT ;  break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED, 
                (ERROR_UNSUPPORTED, "analyzeWrite: unsupported MRI type %d",
                 mri->type)) ;
    break ;
  }

  write_analyze_header(fname, mri) ;
  write_analyze_image(fname, mri) ;
  write_analyze_mat(fname, mri);

  return(NO_ERROR) ;
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static MRI *
analyzeRead(char *fname, int read_volume, int frame)
{
  MRI    *mri_dst ;
  dsr    hdr;
  int    type, width, height, depth;
  FILE   *fp;
  char   file_name[STRLEN], *c;
  int    swap_flag;
  char   *cbuf;
  short  *sbuf;
  int    *ibuf;
  float  *fbuf;
  double *dbuf;
  int    i, j, k;
  float  fval;
  double dval;
  MATRIX *omat;
  int xmax, ymax, zmax, xsign, ysign, zsign;

  strcpy(file_name, fname);
  c = strrchr(file_name, '.');
  c = (c == NULL ? &file_name[strlen(file_name)] : c);
  sprintf(c, ".hdr");

  if((fp = fopen(file_name, "r")) == NULL)
    ErrorReturn(NULL, (ERROR_BADPARM, "analyzeRead(): can't open file %s", file_name));
  if(!fread(&hdr, sizeof(hdr), 1, fp))
  {
    fclose(fp);
    ErrorReturn(NULL, (ERROR_BADFILE, "analyzeRead(): error reading from file %s", file_name));
  }
  fclose(fp);

  if(hdr.hk.sizeof_hdr == sizeof(hdr))
    swap_flag = 0;
  else
  {
    flipAnalyzeHeader(&hdr);
    if(hdr.hk.sizeof_hdr != sizeof(hdr))
      ErrorReturn(NULL, (ERROR_BADFILE, "analyzeRead(): bad byte ordering in file %s", fname));
    swap_flag = 1;
  }

  if(hdr.dime.datatype == DT_UNSIGNED_CHAR)
    type = MRI_UCHAR;
  else if(hdr.dime.datatype == DT_SIGNED_SHORT)
    type = MRI_SHORT;
  else if(hdr.dime.datatype == DT_SIGNED_INT)
    type = MRI_INT;
  else if(hdr.dime.datatype == DT_FLOAT)
    type = MRI_FLOAT;
  else if(hdr.dime.datatype == DT_DOUBLE)
    type = MRI_FLOAT;
  else
    ErrorReturn(NULL, 
                (ERROR_UNSUPPORTED, "analyzeRead: unsupported data type %d",
                 hdr.dime.datatype)) ;

  width = hdr.dime.dim[1];
  depth = hdr.dime.dim[3];
  height = hdr.dime.dim[2];

  if(!read_volume)
    mri_dst = MRIallocHeader(width, height, depth, type);
  else
  {
    mri_dst = MRIalloc(width, height, depth, type);

    strcpy(file_name, fname);
    c = strrchr(file_name, '.');
    c = (c == NULL ? &file_name[strlen(file_name)] : c);
    sprintf(c, ".img");

    if((fp = fopen(file_name, "r")) == NULL)
      ErrorReturn(NULL, (ERROR_BADPARM, "analyzeRead(): can't open file %s", file_name));

    if(hdr.dime.datatype == DT_UNSIGNED_CHAR)
    {
      cbuf = (char *)malloc(width * height * depth * sizeof(char));

      if(fread(cbuf, sizeof(char), width * height * depth, fp) != width * height * depth)
      {
        fclose(fp);
        ErrorReturn(NULL, (ERROR_BADFILE, "analyzeRead(): error reading from file %s", file_name));
      }

      for(i = 0;i < width;i++)
        for(j = 0;j < height;j++)
          for(k = 0;k < depth;k++)
            MRIvox(mri_dst, i, j, k) = cbuf[k * width * height + j * width + i];

      free(cbuf);
    }

    if(hdr.dime.datatype == DT_SIGNED_SHORT)
    {
      sbuf = (short *)malloc(width * height * depth * sizeof(short));

      if(fread(sbuf, sizeof(short), width * height * depth, fp) != width * height * depth)
      {
        fclose(fp);
        ErrorReturn(NULL, (ERROR_BADFILE, "analyzeRead(): error reading from file %s", file_name));
      }

      if(swap_flag)
      {
        for(i = 0;i < width;i++)
          for(j = 0;j < height;j++)
            for(k = 0;k < depth;k++)
              MRISvox(mri_dst, i, j, k) = swapShort(sbuf[k * width * height + j * width + i]);
      }
      else
      {
        for(i = 0;i < width;i++)
          for(j = 0;j < height;j++)
            for(k = 0;k < depth;k++)
              MRISvox(mri_dst, i, j, k) = sbuf[k * width * height + j * width + i];
      }

      free(sbuf);
    }

    if(hdr.dime.datatype == DT_SIGNED_INT)
    {
      ibuf = (int *)malloc(width * height * depth * sizeof(int));

      if(fread(ibuf, sizeof(int), width * height * depth, fp) != width * height * depth)
      {
        fclose(fp);
        ErrorReturn(NULL, (ERROR_BADFILE, "analyzeRead(): error reading from file %s", file_name));
      }

      if(swap_flag)
      {
        for(i = 0;i < width;i++)
          for(j = 0;j < height;j++)
            for(k = 0;k < depth;k++)
              MRIIvox(mri_dst, i, j, k) = swapInt(ibuf[k * width * height + j * width + i]);
      }
      else
      {
        for(i = 0;i < width;i++)
          for(j = 0;j < height;j++)
            for(k = 0;k < depth;k++)
              MRIIvox(mri_dst, i, j, k) = ibuf[k * width * height + j * width + i];
      }

      free(ibuf);
    }

    if(hdr.dime.datatype == DT_FLOAT)
    {
      fbuf = (float *)malloc(width * height * depth * sizeof(float));

      if(fread(fbuf, sizeof(float), width * height * depth, fp) != width * height * depth)
      {
        fclose(fp);
        ErrorReturn(NULL, (ERROR_BADFILE, "analyzeRead(): error reading from file %s", file_name));
      }

      if(swap_flag)
      {
        for(i = 0;i < width;i++)
          for(j = 0;j < height;j++)
            for(k = 0;k < depth;k++)
            {
              fval = swapFloat(fbuf[k * width * height + j * width + i]);
              MRIFvox(mri_dst, i, j, k) = (fval < 0.0 || fval >= 0.0 ? fval : 0.0);
            }
      }
      else
      {
        for(i = 0;i < width;i++)
          for(j = 0;j < height;j++)
            for(k = 0;k < depth;k++)
            {
              fval = fbuf[k * width * height + j * width + i];
              MRIFvox(mri_dst, i, j, k) = (fval < 0.0 || fval >= 0.0 ? fval : 0.0);
            }
      }

      free(fbuf);
    }

    if(hdr.dime.datatype == DT_DOUBLE)
    {
      dbuf = (double *)malloc(width * height * depth * sizeof(double));

      if(fread(dbuf, sizeof(double), width * height * depth, fp) != width * height * depth)
      {
        fclose(fp);
        ErrorReturn(NULL, (ERROR_BADFILE, "analyzeRead(): error reading from file %s", file_name));
      }

      if(swap_flag)
      {
        for(i = 0;i < width;i++)
          for(j = 0;j < height;j++)
            for(k = 0;k < depth;k++)
            {
              dval = swapDouble(dbuf[k * width * height + j * width + i]);
              MRIFvox(mri_dst, i, j, k) = (float)(dval < 0.0 || dval >= 0.0 ? dval : 0.0);
            }
      }
      else
      {
        for(i = 0;i < width;i++)
          for(j = 0;j < height;j++)
            for(k = 0;k < depth;k++)
            {
              dval = dbuf[k * width * height + j * width + i];
              MRIFvox(mri_dst, i, j, k) = (float)(dval < 0.0 || dval >= 0.0 ? dval : 0.0);
            }
      }

      free(dbuf);
    }

    fclose(fp);

  }

  mri_dst->xsize = hdr.dime.pixdim[1];
  mri_dst->ysize = hdr.dime.pixdim[2];
  mri_dst->zsize = hdr.dime.pixdim[3];

  mri_dst->imnr0 = 0;
  mri_dst->imnr1 = mri_dst->width - 1;

  mri_dst->xstart = -(*(short *)&hdr.hist.originator[0] * mri_dst->xsize);
  mri_dst->xend = mri_dst->xstart + mri_dst->xsize * mri_dst->width;
  mri_dst->ystart = -(*(short *)&hdr.hist.originator[2] * mri_dst->ysize);
  mri_dst->yend = mri_dst->ystart + mri_dst->ysize * mri_dst->height;
  mri_dst->zstart = -(*(short *)&hdr.hist.originator[4] * mri_dst->zsize);
  mri_dst->zend = mri_dst->zstart + mri_dst->zsize * mri_dst->depth;

/* direction */

  mri_dst->xdir = XDIM;
  mri_dst->ydir = ZDIM;
  mri_dst->zdir = -YDIM;
  mri_dst->slice_direction = MRI_UNDEFINED;

  mri_dst->ras_good_flag = 0;

  /* use the mat file for orientation if there is one */
  strcpy(file_name, fname);
  c = strrchr(file_name, '.');
  c = (c == NULL ? &file_name[strlen(file_name)] : c);
  sprintf(c, ".mat");
  if((fp = fopen(file_name, "r")) != NULL)
  {
    fclose(fp);
    omat = MatlabRead(file_name);

    *MATRIX_RELT(omat, 1, 1) = - *MATRIX_RELT(omat, 1, 1);
    *MATRIX_RELT(omat, 1, 2) = - *MATRIX_RELT(omat, 1, 2);
    *MATRIX_RELT(omat, 1, 3) = - *MATRIX_RELT(omat, 1, 3);
    *MATRIX_RELT(omat, 1, 4) = - *MATRIX_RELT(omat, 1, 4);

    xmax = 1 + IMAX3(fabs(*MATRIX_RELT(omat, 1, 1)), fabs(*MATRIX_RELT(omat, 2, 1)), fabs(*MATRIX_RELT(omat, 3, 1)));
    ymax = 1 + IMAX3(fabs(*MATRIX_RELT(omat, 1, 2)), fabs(*MATRIX_RELT(omat, 2, 2)), fabs(*MATRIX_RELT(omat, 3, 2)));
    zmax = 1 + IMAX3(fabs(*MATRIX_RELT(omat, 1, 3)), fabs(*MATRIX_RELT(omat, 2, 3)), fabs(*MATRIX_RELT(omat, 3, 3)));

    xsign = SIGN(*MATRIX_RELT(omat, xmax, 1));
    ysign = SIGN(*MATRIX_RELT(omat, ymax, 2));
    zsign = SIGN(*MATRIX_RELT(omat, zmax, 3));

    if(xmax == 1)
      mri_dst->xdir = XDIM * xsign;
    if(xmax == 2)
      mri_dst->xdir = ZDIM * xsign;
    if(xmax == 3)
      mri_dst->xdir = -YDIM * xsign;

    if(ymax == 1)
      mri_dst->ydir = XDIM * ysign;
    if(ymax == 2)
      mri_dst->ydir = ZDIM * ysign;
    if(ymax == 3)
      mri_dst->ydir = -YDIM * ysign;

    if(zmax == 1)
      mri_dst->zdir = XDIM * zsign;
    if(zmax == 2)
      mri_dst->zdir = ZDIM * zsign;
    if(zmax == 3)
      mri_dst->zdir = -YDIM * zsign;

    /* sanity check -- have all the axes been used? */

    /* must be each of 1, 2, and 3 in some order */
    if(xmax + ymax + zmax != 6 || xmax * ymax * zmax != 6)
      ErrorReturn(NULL, (ERROR_BADPARM, "analyzeRead(%s): error interpreting slice direction", fname));

    mri_dst->x_r = *MATRIX_RELT(omat, 1, 1);
    mri_dst->x_a = *MATRIX_RELT(omat, 1, 2);
    mri_dst->x_s = *MATRIX_RELT(omat, 1, 3);

    mri_dst->y_r = *MATRIX_RELT(omat, 2, 1);
    mri_dst->y_a = *MATRIX_RELT(omat, 2, 2);
    mri_dst->y_s = *MATRIX_RELT(omat, 2, 3);

    mri_dst->z_r = *MATRIX_RELT(omat, 3, 1);
    mri_dst->z_a = *MATRIX_RELT(omat, 3, 2);
    mri_dst->z_s = *MATRIX_RELT(omat, 3, 3);

    mri_dst->c_r = *MATRIX_RELT(omat, 4, 1);
    mri_dst->c_a = *MATRIX_RELT(omat, 4, 2);
    mri_dst->c_s = *MATRIX_RELT(omat, 4, 3);

    mri_dst->ras_good_flag = 1;

  }

  return(mri_dst) ;

} /* end analyzeRead() */

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
write_analyze_header(char *fname, MRI *mri)
{
  FILE  *fp;
  char  hdr_fname[STRLEN], *dot ;
  int   nwritten, i ;
  dsr   hdr ;
  float fmin, fmax ;
#ifdef Linux
  char  c;
#endif

  memset(&hdr, 0, sizeof(hdr)) ;
  hdr.hk.sizeof_hdr = sizeof(hdr) ;
  hdr.hk.extents = 0 ;
  hdr.hk.session_error = 0 ;
  for (i=0;i<8;i++)
    hdr.dime.dim[i] = 0 ;
  hdr.dime.dim[1] = mri->width ; 
  hdr.dime.dim[2] = mri->depth ;
  hdr.dime.dim[3] = mri->height ;
  hdr.dime.dim[4] = mri->nframes ;

  switch (mri->type)
  {
  case MRI_FLOAT:
    hdr.dime.datatype = DT_FLOAT ;
    hdr.dime.bitpix = sizeof(float) * 8 ;
    break ;
  case MRI_UCHAR:
    hdr.dime.datatype = DT_UNSIGNED_CHAR ;
    hdr.dime.bitpix = sizeof(unsigned char) * 8 ;
    break ;
  case MRI_INT:
    hdr.dime.datatype = DT_SIGNED_INT ;
    hdr.dime.bitpix = sizeof(int) * 8 ;
    break;
  case MRI_SHORT:
    hdr.dime.datatype = DT_SIGNED_SHORT ;
    hdr.dime.bitpix = sizeof(short) * 8 ;
    break;
  default:
    ErrorReturn(ERROR_UNSUPPORTED, 
                (ERROR_UNSUPPORTED, 
               "write_analyze_header: unsupported volume type %d", mri->type));
  }
  hdr.dime.dim_un0 = 0 ;
  for (i=0;i<8;i++)
    hdr.dime.pixdim[i] = 0 ;
  hdr.dime.pixdim[4] = 1.0 ;
  hdr.dime.pixdim[3] = mri->xsize ;
  hdr.dime.pixdim[2] = mri->zsize ;
  hdr.dime.pixdim[1] = mri->ysize ;

  hdr.dime.compressed = 0 ;
  hdr.dime.verified = 0 ;
  MRIvalRange(mri, &fmin, &fmax) ;
  hdr.dime.glmax = nint(fmax) ;
  hdr.dime.glmin = nint(fmin) ;
  hdr.hist.views = 0 ;
  hdr.hist.vols_added = 0 ;
  hdr.hist.start_field = 0 ;
  hdr.hist.field_skip = 0 ;
  hdr.hist.omax = 0 ;
  hdr.hist.omin = 0 ;
  hdr.hist.smax = 0 ;
  hdr.hist.smin = 0 ;

  if(mri->ras_good_flag)
  {
    *(short *)&hdr.hist.originator[0] = -(short)floor(mri->xstart / mri->xsize);
    *(short *)&hdr.hist.originator[2] = -(short)floor(mri->zstart / mri->zsize);
    *(short *)&hdr.hist.originator[4] = -(short)floor(mri->ystart / mri->xsize);
  }
  else
  {
    *(short *)&hdr.hist.originator[0] = (short)floor(mri->width / 2);
    *(short *)&hdr.hist.originator[2] = (short)floor(mri->depth / 2);
    *(short *)&hdr.hist.originator[4] = (short)floor(mri->height / 2);
  }

  strcpy(hdr_fname, fname) ;
  dot = strrchr(hdr_fname, '.') ;
  if (dot)
    *dot = 0 ;
  strcat(hdr_fname, ".hdr") ;
  fp = fopen(hdr_fname,"wb");
  if (fp==NULL) 
    ErrorReturn(ERROR_NO_FILE, 
                (ERROR_NO_FILE, "File %s not found\n", hdr_fname)) ;

#ifdef Linux
  hdr.hk.sizeof_hdr = swapInt(hdr.hk.sizeof_hdr) ;
  hdr.hk.extents = swapInt(hdr.hk.extents) ;
  hdr.hk.session_error = swapInt(hdr.hk.session_error) ;
  for (i=0;i<8;i++)
    hdr.dime.dim[i] = swapShort(hdr.dime.dim[i]) ;

  hdr.dime.datatype = swapShort(hdr.dime.datatype) ;
  hdr.dime.bitpix = swapShort(hdr.dime.bitpix) ;
  hdr.dime.dim_un0 = swapShort(hdr.dime.dim_un0) ;
  for (i=0;i<8;i++)
    hdr.dime.pixdim[i] = swapFloat(hdr.dime.pixdim[i]) ;

  hdr.dime.compressed = swapFloat(hdr.dime.compressed) ;
  hdr.dime.verified = swapFloat(hdr.dime.verified) ;
  hdr.dime.glmax = swapInt(hdr.dime.glmax) ;
  hdr.dime.glmin = swapInt(hdr.dime.glmin) ;
  hdr.hist.views = swapInt(hdr.hist.views) ;
  hdr.hist.vols_added = swapInt(hdr.hist.vols_added) ;
  hdr.hist.start_field = swapInt(hdr.hist.start_field) ;
  hdr.hist.field_skip = swapInt(hdr.hist.field_skip) ;
  hdr.hist.omax = swapInt(hdr.hist.omax) ;
  hdr.hist.omin = swapInt(hdr.hist.omin) ;
  hdr.hist.smax = swapInt(hdr.hist.smax) ;
  hdr.hist.smin = swapInt(hdr.hist.smin) ;

  c = hdr.hist.originator[0]; hdr.hist.originator[0] = hdr.hist.originator[1]; hdr.hist.originator[1] = c;
  c = hdr.hist.originator[2]; hdr.hist.originator[2] = hdr.hist.originator[3]; hdr.hist.originator[3] = c;
  c = hdr.hist.originator[4]; hdr.hist.originator[4] = hdr.hist.originator[5]; hdr.hist.originator[5] = c;
  c = hdr.hist.originator[6]; hdr.hist.originator[6] = hdr.hist.originator[7]; hdr.hist.originator[7] = c;
  c = hdr.hist.originator[8]; hdr.hist.originator[8] = hdr.hist.originator[9]; hdr.hist.originator[9] = c;

#endif
  nwritten = fwrite(&hdr, sizeof(char), sizeof(dsr),fp);
  fclose(fp);
  if (nwritten != sizeof(dsr))
    ErrorReturn(ERROR_BADFILE, 
                (ERROR_BADFILE,
                 "write_analyze_header: could write read %d of %d bytes",
                 nwritten, sizeof(dsr))) ;
  return(NO_ERROR) ;
}

static int
write_analyze_mat(char *fname, MRI *mri)
{

  char mat_fname[STRLEN];
  char *dot;
  MATRIX *m, *o, *r;

  strcpy(mat_fname, fname) ;
  dot = strrchr(mat_fname, '.') ;
  if (dot)
    *dot = 0 ;
  strcat(mat_fname, ".mat") ;

  m = MatrixAlloc(4, 4, MATRIX_REAL);

/* fastest: rl (x)
   next:    pa (z)
   slowest: is (y) */

  if(mri->ras_good_flag)
  {

    *MATRIX_RELT(m, 1, 1) = mri->x_r;
    *MATRIX_RELT(m, 2, 1) = mri->x_a;
    *MATRIX_RELT(m, 3, 1) = -mri->x_s;
    *MATRIX_RELT(m, 4, 1) = 0;

    *MATRIX_RELT(m, 1, 2) = mri->z_r;
    *MATRIX_RELT(m, 2, 2) = mri->z_a;
    *MATRIX_RELT(m, 3, 2) = -mri->z_s;
    *MATRIX_RELT(m, 4, 2) = 0;

    *MATRIX_RELT(m, 1, 3) = mri->y_r;
    *MATRIX_RELT(m, 2, 3) = mri->y_a;
    *MATRIX_RELT(m, 3, 3) = -mri->y_s;
    *MATRIX_RELT(m, 4, 3) = 0;

    /* the originator should map to zero */

    *MATRIX_RELT(m, 1, 4) = 0;
    *MATRIX_RELT(m, 2, 4) = 0;
    *MATRIX_RELT(m, 3, 4) = 0;
    *MATRIX_RELT(m, 4, 4) = 1;

    o = MatrixAlloc(4, 1, MATRIX_REAL);
    *MATRIX_RELT(o, 1, 1) = floor(mri->xstart / mri->xsize);
    *MATRIX_RELT(o, 2, 1) = floor(mri->zstart / mri->zsize);
    *MATRIX_RELT(o, 3, 1) = floor(mri->ystart / mri->ysize);
    *MATRIX_RELT(o, 4, 1) = 1;

    r = MatrixMultiply(m, o, NULL);

    *MATRIX_RELT(m, 1, 4) = *MATRIX_RELT(r, 1, 1);
    *MATRIX_RELT(m, 2, 4) = *MATRIX_RELT(r, 2, 1);
    *MATRIX_RELT(m, 3, 4) = *MATRIX_RELT(r, 3, 1);
    *MATRIX_RELT(m, 4, 4) = 1;

    MatrixFree(&r);
    MatrixFree(&o);

  }
  else
  {
    *MATRIX_RELT(m, 1, 1) = -mri->xsize;  *MATRIX_RELT(m, 1, 2) = 0;           *MATRIX_RELT(m, 1, 3) = 0;           *MATRIX_RELT(m, 1, 4) = (float)floor(mri->width / 2);
    *MATRIX_RELT(m, 2, 1) = 0;           *MATRIX_RELT(m, 2, 2) = mri->zsize;  *MATRIX_RELT(m, 2, 3) = 0;           *MATRIX_RELT(m, 2, 4) = -(float)floor(mri->depth / 2);
    *MATRIX_RELT(m, 3, 1) = 0;           *MATRIX_RELT(m, 3, 2) = 0;           *MATRIX_RELT(m, 3, 3) = mri->ysize;  *MATRIX_RELT(m, 3, 4) = -(float)floor(mri->height / 2);
    *MATRIX_RELT(m, 4, 1) = 0;           *MATRIX_RELT(m, 4, 2) = 0;           *MATRIX_RELT(m, 4, 3) = 0;           *MATRIX_RELT(m, 4, 4) = 1;
  }

  MatlabWrite(m, mat_fname, "M");

  MatrixFree(&m);

  return(NO_ERROR);

} /* end write_analyze_mat() */
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
write_analyze_image(char *fname, MRI *mri)
{
  int    x, y, z, bufsize, bytes_per_pix, xd, yd, zd ;
  FILE   *fp;
  int    nwritten ;
  char   *buf ;
  float  f ;

  if(mri->type == MRI_FLOAT)
    bytes_per_pix = sizeof(float);
  else if(mri->type == MRI_UCHAR)
    bytes_per_pix = sizeof(char);
  else if(mri->type == MRI_INT)
    bytes_per_pix = sizeof(int);
  else if(mri->type == MRI_SHORT)
    bytes_per_pix = sizeof(short);
  else
    ErrorReturn(ERROR_UNSUPPORTED, (ERROR_UNSUPPORTED, "write_analyze_image() unsupported data type %d: ", mri->type));

  bufsize =  mri->width * mri->depth;

  fp = fopen(fname,"wb");
  if (fp==NULL) 
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "File %s not found\n",fname)) ;

/* fastest: rl (x)
   next:    pa (z)
   slowest: is (y) */

  buf = (char *)calloc(bufsize, bytes_per_pix) ;
  if (!buf)
  {
    fclose(fp) ;
    ErrorReturn(ERROR_NOMEMORY, 
                (ERROR_NOMEMORY, 
                "write_analyze_image: could not allocate %d x %d x %d buffer",
                 mri->width, mri->height, bytes_per_pix)) ;
  }

  for (y = 0 ; y < mri->height ; y++)
  {
    for (z = 0 ; z < mri->depth ; z++)
    {
      for (x = 0 ; x < mri->width ; x++)
      {
        zd = z  ; yd = (mri->height-y-1) ; xd = x ;
        if (mri->type == MRI_FLOAT)
          f = MRIFvox(mri, xd, yd, zd) ;
        else if(mri->type == MRI_INT)
          f = (float)MRIIvox(mri, xd, yd, zd) ;
        else if(mri->type == MRI_SHORT)
          f = (float)MRISvox(mri, xd, yd, zd) ;
        else
          f = (float)MRIvox(mri, xd, yd, zd) ;
        switch (mri->type)
        {
        default:
          printf("data type %d not supported\n",mri->type);
          exit(1);
          break;
        case MRI_UCHAR:
          *(unsigned char *)(buf+bytes_per_pix*(z*mri->width + x)) = (char)nint(f) ;
          break;
        case MRI_FLOAT:
          *(float *)(buf+bytes_per_pix*(z*mri->width + x)) = orderFloatBytes(f) ;
          break;
        case MRI_INT:
          *(int *)(buf+bytes_per_pix*(z*mri->width + x)) = orderIntBytes((int)f) ;
          break;
        case MRI_SHORT:
          *(short *)(buf+bytes_per_pix*(z*mri->width + x)) = orderShortBytes((short)f) ;
          break;
        }
      }
    }
    nwritten = fwrite(buf, bytes_per_pix, bufsize, fp) ;
    if (nwritten != bufsize)
    {
      free(buf) ;
      fclose(fp) ;
      ErrorReturn(ERROR_BADFILE, 
                  (ERROR_BADFILE, 
                   "write_analyze_image: could not write slice %d "
                   "(%d items read)",
                   z, bufsize)) ;
    }

  }
  fclose(fp);
  free(buf) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIreorder(MRI *mri_src, MRI *mri_dst, int xdim, int ydim, int zdim)
{
  int  width, height, depth, xs, ys, zs, xd, yd, zd, x, y, z ;
  float ras_sign;

  width = mri_src->width ; height = mri_src->height ; depth = mri_src->depth;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  /* check that the source ras coordinates are good and that each direction is used once and only once */
  if(mri_src->ras_good_flag)
    if(abs(xdim) * abs(ydim) * abs(zdim) != 6 || abs(xdim) + abs(ydim) + abs(zdim) != 6)
      mri_dst->ras_good_flag = 0;

  xd = yd = zd = 0 ;

  ras_sign = (xdim < 0 ? -1.0 : 1.0);
  switch (abs(xdim))
  {
  default:
  case XDIM: 
    if (mri_dst->width != width)
      ErrorReturn(NULL,(ERROR_BADPARM, "MRIreorder: incorrect dst width"));
    if(mri_dst->ras_good_flag)
    {
      mri_dst->x_r = mri_src->x_r * ras_sign;
      mri_dst->x_a = mri_src->x_a * ras_sign;
      mri_dst->x_s = mri_src->x_s * ras_sign;
    }
    break ;
  case YDIM: 
    if (mri_dst->height != width)
      ErrorReturn(NULL,(ERROR_BADPARM, "MRIreorder: incorrect dst width"));
    if(mri_dst->ras_good_flag)
    {
      mri_dst->y_r = mri_src->x_r * ras_sign;
      mri_dst->y_a = mri_src->x_a * ras_sign;
      mri_dst->y_s = mri_src->x_s * ras_sign;
    }
    break ;
  case ZDIM: 
    if (mri_dst->depth != width)
      ErrorReturn(NULL,(ERROR_BADPARM, "MRIreorder: incorrect dst width"));
    if(mri_dst->ras_good_flag)
    {
      mri_dst->z_r = mri_src->x_r * ras_sign;
      mri_dst->z_a = mri_src->x_a * ras_sign;
      mri_dst->z_s = mri_src->x_s * ras_sign;
    }
    break ;
  }
  ras_sign = (ydim < 0 ? -1.0 : 1.0);
  switch (abs(ydim))
  {
  default:
  case XDIM: 
    if (mri_dst->width != height)
      ErrorReturn(NULL,(ERROR_BADPARM, "MRIreorder: incorrect dst height"));
    if(mri_dst->ras_good_flag)
    {
      mri_dst->x_r = mri_src->y_r * ras_sign;
      mri_dst->x_a = mri_src->y_a * ras_sign;
      mri_dst->x_s = mri_src->y_s * ras_sign;
    }
    break ;
  case YDIM: 
    if (mri_dst->height != height)
      ErrorReturn(NULL,(ERROR_BADPARM, "MRIreorder: incorrect dst height"));
    if(mri_dst->ras_good_flag)
    {
      mri_dst->y_r = mri_src->y_r * ras_sign;
      mri_dst->y_a = mri_src->y_a * ras_sign;
      mri_dst->y_s = mri_src->y_s * ras_sign;
    }
    break ;
  case ZDIM: 
    if (mri_dst->depth != height)
      ErrorReturn(NULL,(ERROR_BADPARM, "MRIreorder: incorrect dst height"));
    if(mri_dst->ras_good_flag)
    {
      mri_dst->z_r = mri_src->y_r * ras_sign;
      mri_dst->z_a = mri_src->y_a * ras_sign;
      mri_dst->z_s = mri_src->y_s * ras_sign;
    }
    break ;
  }
  ras_sign = (zdim < 0 ? -1.0 : 1.0);
  switch (abs(zdim))
  {
  default:
  case XDIM: 
    if (mri_dst->width != depth)
      ErrorReturn(NULL,(ERROR_BADPARM, "MRIreorder: incorrect dst depth"));
    if(mri_dst->ras_good_flag)
    {
      mri_dst->x_r = mri_src->z_r * ras_sign;
      mri_dst->x_a = mri_src->z_a * ras_sign;
      mri_dst->x_s = mri_src->z_s * ras_sign;
    }
    break ;
  case YDIM: 
    if (mri_dst->height != depth)
      ErrorReturn(NULL,(ERROR_BADPARM, "MRIreorder: incorrect dst depth"));
    if(mri_dst->ras_good_flag)
    {
      mri_dst->y_r = mri_src->z_r * ras_sign;
      mri_dst->y_a = mri_src->z_a * ras_sign;
      mri_dst->y_s = mri_src->z_s * ras_sign;
    }
    break ;
  case ZDIM: 
    if (mri_dst->depth != depth)
      ErrorReturn(NULL,(ERROR_BADPARM, "MRIreorder: incorrect dst depth"));
    if(mri_dst->ras_good_flag)
    {
      mri_dst->z_r = mri_src->z_r * ras_sign;
      mri_dst->z_a = mri_src->z_a * ras_sign;
      mri_dst->z_s = mri_src->z_s * ras_sign;
    }
    break ;
  }

  for (zs = 0 ; zs < depth ; zs++)
  {
    if (zdim < 0)
      z = depth - zs - 1 ;
    else
      z = zs ;
    switch (abs(zdim))
    {
    case XDIM:  xd = z ; break ;
    case YDIM:  yd = z ; break ;
    default:
    case ZDIM:  zd = z ; break ;
    }
    for (ys = 0 ; ys < height ; ys++)
    {
      if (ydim < 0)
        y = height - ys - 1 ;
      else
        y = ys ;
      switch (abs(ydim))
      {
      case XDIM: xd = y ; break ;
      case YDIM: yd = y ; break ;
      default:
      case ZDIM: zd = y ; break ;
      }
      for (xs = 0 ; xs < width ; xs++)
      {
        if (xdim < 0)
          x = width - xs - 1 ;
        else
          x = xs ;
        switch (abs(xdim))
        {
          case XDIM: xd = x ; break ;
          case YDIM: yd = x ; break ;
          default:
          case ZDIM: zd = x ; break ;
        }
        switch (mri_src->type)
        {
        case MRI_SHORT:
          MRISvox(mri_dst, xd, yd, zd) = MRISvox(mri_src, xs, ys, zs) ;
          break ;
        case MRI_FLOAT:
          MRIFvox(mri_dst, xd, yd, zd) = MRIFvox(mri_src, xs, ys, zs) ;
          break ;
        case MRI_INT:
          MRIIvox(mri_dst, xd, yd, zd) = MRIIvox(mri_src, xs, ys, zs) ;
          break ;
        case MRI_UCHAR:
          MRIvox(mri_dst, xd, yd, zd) = MRIvox(mri_src, xs, ys, zs) ;
          break ;
        default:
          ErrorReturn(NULL,
                      (ERROR_UNSUPPORTED, 
                       "MRIreorder: unsupported voxel format %d",
                       mri_src->type)) ;
          break ;
        }
      }
    }
  }
  return(mri_dst) ;
}

static void get_file_limits(char *format_string, int initial, int *min, int *max, int ge_check_flag)
{

  char file_name[STRLEN];
  FILE *fp;
  char volume_string_1[STRLEN], volume_string_2[STRLEN];
  int exists_flag;

  if(ge_check_flag)
  {
    sprintf(file_name, format_string, initial);
    fp = fopen(file_name, "r");
    fseek(fp, 169, 0);
    fread(volume_string_1, 1, 20, fp);
    fclose(fp);
  }

  for(exists_flag = 1, *min = initial - 1;exists_flag;(*min)--)
  {
    sprintf(file_name, format_string, *min);
    fp = fopen(file_name, "r");
    if(fp == NULL)
      exists_flag = 0;
    else
    {
      if(ge_check_flag)
      {
        fseek(fp, 169, 0);
        fread(volume_string_2, 1, 20, fp);
        if(strcmp(volume_string_1, volume_string_2))
          exists_flag = 0;
      }
      fclose(fp);
    }
  }

  (*min) += 2;

  for(exists_flag = 1, *max = initial + 1;exists_flag;(*max)++)
  {
    sprintf(file_name, format_string, *max);
    fp = fopen(file_name, "r");
    if(fp == NULL)
      exists_flag = 0;
    else
    {
      if(ge_check_flag)
      {
        fseek(fp, 169, 0);
        fread(volume_string_2, 1, 20, fp);
        if(strcmp(volume_string_1, volume_string_2))
          exists_flag = 0;
      }
      fclose(fp);
    }
  }

  (*max) -= 2;

} /* end get_file_limits() */


static MRI *
genesisRead(char *fname, int read_volume, int frame)
{
  unsigned char header_header[156], magic_header[4], image_header[116];
  MRI *mri = NULL;
  FILE *fp ;
  long magic;
  int width, height, depth, offset_to_pixel_data, offset_to_image_header;
  int series_header_offset;
  int slice_number, max_slice_number, min_slice_number;
  char fname_format[STRLEN], *i_pos, fname_use[STRLEN], fname_copy[STRLEN];
  float d_fov_x, d_fov_y, d_fov_z;
  short *pixel_data = NULL;
  char fname_mod[STR_LEN];
  float top_left[2][3], top_right[2][3], bottom_right[2][3], normal[2][3];
  int ras_index = 0;
  float coords_in[12];
  float x_vec[3], y_vec[3], z_vec[3];
  int xmax, ymax, zmax;
  int xsign, ysign, zsign;

  frame = 0; /* fix this */

  strcpy(fname_copy, fname);
  if((i_pos = strrchr(fname_copy, 'I')) == NULL)
    ErrorReturn(NULL, (ERROR_BADPARM, "bad file name construction: %s", fname));
  *i_pos = '\0';
  i_pos++;
  if(*i_pos == '.')
    {
    sprintf(fname_format, "%sI.%%.3d", fname_copy);
    i_pos++;
    }
  else
    sprintf(fname_format, "%sI%%d.MR", fname_copy);

  get_file_limits(fname_format, strtol(i_pos, NULL, 10), &min_slice_number, &max_slice_number, 1);

  depth = max_slice_number - min_slice_number + 1;

  if(max_slice_number <= 0)
    ErrorReturn(NULL, (ERROR_BADPARM, "genesisRead(%s): error reading files in sequence", fname));

  if((fp = fopen(fname, "rb")) == NULL)
    ErrorReturn(NULL, (ERROR_BADPARM, "genesisRead(%s, %d): could not open file (1)",
                       fname, frame)) ;

  /* check the magic number of the file */
  fseek(fp, 0L, SEEK_SET) ;
  if(fread(magic_header, 4, 1, fp) < 1)
  {
    fclose(fp) ;
    ErrorReturn(NULL, (ERROR_BADFILE, "genesisRead(%s, %d): error reading file",
                       fname, frame)) ;
  }
  magic = orderLongBytes(*((long *)magic_header)) ; 
  if(magic != GE_MAGIC)
  {
    fclose(fp) ;
    ErrorReturn(NULL, (ERROR_UNSUPPORTED,
                       "genesisRead: corrupt file or unknown file type")) ;
  }

  /* read the control header */
  fseek(fp, 0L, SEEK_SET) ;
  if(fread(header_header, 156, 1, fp) < 1)
  {
    fclose(fp) ;
    ErrorReturn(NULL, (ERROR_BADFILE, "genesisRead(%s, %d): error reading file",
                       fname, frame)) ;
  }

  width = orderIntBytes(*((int *)(&header_header[8])));
  height = orderIntBytes(*((int *)(&header_header[8])));
  series_header_offset = orderIntBytes(*((int *)(&header_header[144])));
  offset_to_image_header = orderIntBytes(*((int *)(&header_header[148]))) ;

  fseek(fp, offset_to_image_header, SEEK_SET);
  if(fread(image_header, 116, 1, fp) < 1)
  {
    fclose(fp) ;
    ErrorReturn(NULL, (ERROR_BADFILE, "genesisRead(%s, %d): error reading file",
                       fname, frame)) ;
  }

  if(read_volume)
    mri = MRIalloc(width, height, depth, MRI_SHORT) ;
  else
    mri = MRIallocHeader(width, height, depth, MRI_SHORT) ;

  mri->zsize = orderFloatBytes(*((float *)(&image_header[26])));
  mri->xsize = orderFloatBytes(*((float *)(&image_header[50])));
  mri->ysize = orderFloatBytes(*((float *)(&image_header[54])));
  mri->thick = mri->xsize;

  d_fov_x = mri->xsize * mri->width;
  d_fov_y = mri->ysize * mri->height;
  d_fov_z = mri->zsize * mri->depth;
  mri->xend = d_fov_x / 2.;
  mri->yend = d_fov_y / 2.;
  mri->zend = d_fov_z / 2.;
  mri->xstart = -mri->xend;
  mri->ystart = -mri->yend;
  mri->zstart = -mri->zend;

  mri->fov = d_fov_x;

  fclose(fp) ;

  mri->imnr0 = 1 ;
  mri->imnr1 = depth;

  mri->ps = 1.0 /*0.001*/;
  mri->tr = 0 ;
  mri->te = 0 ;
  mri->ti = 0 ;
  strcpy(mri->fname, fname) ;

  mri->slice_direction = MRI_UNDEFINED;

  if(read_volume)
  {

    pixel_data = (short *)malloc(2 * mri->width * mri->height);

    for(slice_number = min_slice_number;slice_number <= max_slice_number;slice_number++)
      {
      sprintf(fname_use, fname_format, slice_number);
      if((fp = fopen(fname_use, "rb")) == NULL)
        ErrorReturn(NULL, (ERROR_BADPARM,
                    "genesisRead(%s, %d): could not open file (2)",
                    fname_mod, frame)) ;

      /* check the magic number of the file */
      fseek(fp, 0L, SEEK_SET) ;
      if(fread(header_header, 156, 1, fp) < 1)
      {
        fclose(fp) ;
        ErrorReturn(NULL, (ERROR_BADFILE, "genesisRead(%s, %d): error reading file",
                           fname, frame)) ;
      }
      magic = orderLongBytes(*((long *)header_header)) ; 
      if(magic != GE_MAGIC)
      {
        fclose(fp) ;
        ErrorReturn(NULL, (ERROR_UNSUPPORTED,
                         "genesisRead: corrupt file or unknown file type")) ;
      }
      offset_to_pixel_data = orderIntBytes(*((int *)(&header_header[4]))) ;

      fseek(fp, offset_to_pixel_data, SEEK_SET) ;
      fread(pixel_data, 1, mri->width * mri->height * 2, fp);
#ifdef Linux
      swab(pixel_data, pixel_data, mri->width*mri->height*2);
#endif

      short_buffer_to_image(pixel_data, mri, slice_number - min_slice_number, 0);

      if(ras_index < 2)
      {
        fseek(fp, offset_to_image_header + 142, SEEK_SET);
        fread(coords_in, sizeof(float), 12, fp);
        normal[ras_index][0] = orderFloatBytes(coords_in[0]);
        normal[ras_index][1] = orderFloatBytes(coords_in[1]);
        normal[ras_index][2] = orderFloatBytes(coords_in[2]);
        top_left[ras_index][0] = orderFloatBytes(coords_in[3]);
        top_left[ras_index][1] = orderFloatBytes(coords_in[4]);
        top_left[ras_index][2] = orderFloatBytes(coords_in[5]);
        top_right[ras_index][0] = orderFloatBytes(coords_in[6]);
        top_right[ras_index][1] = orderFloatBytes(coords_in[7]);
        top_right[ras_index][2] = orderFloatBytes(coords_in[8]);
        bottom_right[ras_index][0] = orderFloatBytes(coords_in[9]);
        bottom_right[ras_index][1] = orderFloatBytes(coords_in[10]);
        bottom_right[ras_index][2] = orderFloatBytes(coords_in[11]);
        ras_index++;
      }

      fclose(fp) ;
    }

    free(pixel_data);

    x_vec[0] = (top_right[0][0] - top_left[0][0]);
    x_vec[1] = (top_right[0][1] - top_left[0][1]);
    x_vec[2] = (top_right[0][2] - top_left[0][2]);

    y_vec[0] = (bottom_right[0][0] - top_right[0][0]);
    y_vec[1] = (bottom_right[0][1] - top_right[0][1]);
    y_vec[2] = (bottom_right[0][2] - top_right[0][2]);

    z_vec[0] = top_left[1][0] - top_left[0][0];
    z_vec[1] = top_left[1][1] - top_left[0][1];
    z_vec[2] = top_left[1][2] - top_left[0][2];

    /* determine the axis reordering that's needed */
    /* start with the x axis; is the R, A, or S coordinate the longest? */
/*
printf("%g, %g, %g\n", x_vec[0], x_vec[1], x_vec[2]);
printf("%g, %g, %g\n", y_vec[0], y_vec[1], y_vec[2]);
printf("%g, %g, %g\n", z_vec[0], z_vec[1], z_vec[2]);
*/
    xmax = IMAX3(fabs(x_vec[0]), fabs(x_vec[1]), fabs(x_vec[2]));
    ymax = IMAX3(fabs(y_vec[0]), fabs(y_vec[1]), fabs(y_vec[2]));
    zmax = IMAX3(fabs(z_vec[0]), fabs(z_vec[1]), fabs(z_vec[2]));

    xsign = SIGN(x_vec[xmax]);
    ysign = SIGN(y_vec[ymax]);
    zsign = SIGN(z_vec[zmax]);

    if(xmax == 0)
      mri->xdir = -XDIM * xsign;
    if(xmax == 1)
      mri->xdir = ZDIM * xsign;
    if(xmax == 2)
      mri->xdir = -YDIM * xsign;

    if(ymax == 0)
      mri->ydir = -XDIM * ysign;
    if(ymax == 1)
      mri->ydir = ZDIM * ysign;
    if(ymax == 2)
      mri->ydir = -YDIM * ysign;

    if(zmax == 0)
      mri->zdir = -XDIM * zsign;
    if(zmax == 1)
      mri->zdir = ZDIM * zsign;
    if(zmax == 2)
      mri->zdir = -YDIM * zsign;

    /* sanity check -- have all the axes been used? */
    xmax = abs(xmax + 1);
    ymax = abs(ymax + 1);
    zmax = abs(zmax + 1);

    /* must be each of 1, 2, and 3 in some order */
    if(xmax + ymax + zmax != 6 || xmax * ymax * zmax != 6)
{
printf("eisd: %d %d %d\n", xmax, ymax, zmax);
      ErrorReturn(NULL, (ERROR_BADPARM, "genesisRead(%s): error interpreting slice direction", fname));
}
    /* assign {xyz}_{ras} coordinates */
/*
    mri->x_r = x_vec[0] / mri->width;
    mri->x_a = x_vec[1] / mri->width;
    mri->x_s = x_vec[2] / mri->width;
    mri->y_r = y_vec[0] / mri->depth;
    mri->y_a = y_vec[1] / mri->depth;
    mri->y_s = y_vec[2] / mri->depth;
    mri->z_r = z_vec[0];
    mri->z_a = z_vec[1];
    mri->z_s = z_vec[2];
*/
    mri->x_r = x_vec[0] / mri->width;
    mri->x_a = x_vec[1] / mri->width;
    mri->x_s = x_vec[2] / mri->width;
    mri->y_r = y_vec[0] / mri->height;
    mri->y_a = y_vec[1] / mri->height;
    mri->y_s = y_vec[2] / mri->height;
    mri->z_r = z_vec[0];
    mri->z_a = z_vec[1];
    mri->z_s = z_vec[2];
    mri->c_r = top_left[0][0] + x_vec[0] / 2 + y_vec[0] / 2 + z_vec[0] * mri->depth / 2;
    mri->c_a = top_left[0][1] + x_vec[1] / 2 + y_vec[1] / 2 + z_vec[1] * mri->depth / 2;
    mri->c_s = top_left[0][2] + x_vec[2] / 2 + y_vec[2] / 2 + z_vec[2] * mri->depth / 2;
    mri->ras_good_flag = 1;

  } /* end if(read_volume) */

  return(mri);

} /* end genesisRead() */


static MRI *
gelxRead(char *fname, int read_volume, int frame)
{
  unsigned char header_header[24], magic_header[4], image_header[420];
  MRI *mri = NULL;
  FILE *fp ;
  long magic;
  int width, height, depth, offset_to_pixel_data;
  int slice_number, max_slice_number, min_slice_number;
  char fname_format[STRLEN], *i_pos, fname_use[STRLEN], fname_copy[STRLEN];
  float d_fov_x, d_fov_y, d_fov_z;
  short *pixel_data = NULL;
  char fname_mod[STR_LEN];
  float top_left[2][3], top_right[2][3], bottom_right[2][3], normal[2][3];
  int ras_index = 0;
  float x_vec[3], y_vec[3], z_vec[3];
  int xmax, ymax, zmax;
  int xsign, ysign, zsign;

  frame = 0; /* fix this */

  strcpy(fname_copy, fname);
  if((i_pos = strrchr(fname_copy, 'i')) == NULL)
    ErrorReturn(NULL, (ERROR_BADPARM, "bad file name construction: %s", fname));
  *i_pos = '\0';
  i_pos++;

  sprintf(fname_format, "%si%%d", fname_copy);
  get_file_limits(fname_format, strtol(i_pos, NULL, 10), &min_slice_number, &max_slice_number, 1);

  depth = max_slice_number - min_slice_number + 1;

  if(max_slice_number <= 0)
    ErrorReturn(NULL, (ERROR_BADPARM, "gelxRead(%s): error reading files in sequence", fname));

  if((fp = fopen(fname, "rb")) == NULL)
    ErrorReturn(NULL, (ERROR_BADPARM, "gelxRead(%s, %d): could not open file (1)",
                       fname, frame)) ;

  /* check the magic number of the file */
  fseek(fp, 3228, SEEK_SET) ;
  if(fread(magic_header, 4, 1, fp) < 1)
  {
    fclose(fp) ;
    ErrorReturn(NULL, (ERROR_BADFILE, "gelxRead(%s, %d): error reading file",
                       fname, frame)) ;
  }
  magic = orderLongBytes(*((long *)magic_header)) ; 
  if(magic != GE_MAGIC)
  {
    fclose(fp) ;
    ErrorReturn(NULL, (ERROR_UNSUPPORTED,
                       "gelxRead: corrupt file or unknown file type")) ;
  }

  /* read the control header */
  fseek(fp, 3228, SEEK_SET) ;
  if(fread(header_header, 24, 1, fp) < 1)
  {
    fclose(fp) ;
    ErrorReturn(NULL, (ERROR_BADFILE, "gelxRead(%s, %d): error reading file",
                       fname, frame)) ;
  }

  height = orderIntBytes(*((int *)(&header_header[8])));
  width = orderIntBytes(*((int *)(&header_header[12])));

  fseek(fp, 2184, SEEK_SET);
  if(fread(image_header, 420, 1, fp) < 1)
  {
    fclose(fp) ;
    ErrorReturn(NULL, (ERROR_BADFILE, "gelxRead(%s, %d): error reading file",
                       fname, frame)) ;
  }

  if(read_volume)
    mri = MRIalloc(width, height, depth, MRI_SHORT) ;
  else
    mri = MRIallocHeader(width, height, depth, MRI_SHORT) ;

  mri->zsize = orderFloatBytes(*((float *)(&image_header[28])));
  mri->xsize = orderFloatBytes(*((float *)(&image_header[52])));
  mri->ysize = orderFloatBytes(*((float *)(&image_header[56])));
  mri->thick = mri->xsize;

  d_fov_x = mri->xsize * mri->width;
  d_fov_y = mri->ysize * mri->height;
  d_fov_z = mri->zsize * mri->depth;
  mri->xend = d_fov_x / 2.;
  mri->yend = d_fov_y / 2.;
  mri->zend = d_fov_z / 2.;
  mri->xstart = -mri->xend;
  mri->ystart = -mri->yend;
  mri->zstart = -mri->zend;

  mri->fov = d_fov_x;

  mri->slice_direction = orderShortBytes(*((short *)(&image_header[116])));

  fclose(fp) ;

  mri->imnr0 = 1 ;
  mri->imnr1 = depth;

  mri->ps = 1 /*0.001*/;
  mri->tr = 0 ;
  mri->te = 0 ;
  mri->ti = 0 ;
  strcpy(mri->fname, fname) ;

  pixel_data = NULL;

  if(read_volume)
  {

    pixel_data = (short *)malloc(2 * mri->width * mri->height);

    for(slice_number = min_slice_number;slice_number <= max_slice_number;slice_number++)
      {
      sprintf(fname_use, fname_format, slice_number);
      if((fp = fopen(fname_use, "rb")) == NULL)
        ErrorReturn(NULL, (ERROR_BADPARM,
                    "gelxRead(%s, %d): could not open file (2)",
                    fname_mod, frame)) ;

      /* check the magic number of the file */
      fseek(fp, 3228, SEEK_SET) ;
      if(fread(header_header, 24, 1, fp) < 1)
      {
        fclose(fp) ;
        ErrorReturn(NULL, (ERROR_BADFILE, "gelxRead(%s, %d): error reading file",
                           fname, frame)) ;
      }
      magic = orderLongBytes(*((long *)header_header)) ; 
      if(magic != GE_MAGIC)
      {
        fclose(fp) ;
        ErrorReturn(NULL, (ERROR_UNSUPPORTED,
                         "gelxRead: corrupt file or unknown file type")) ;
      }

      offset_to_pixel_data = 8432;

      fseek(fp, offset_to_pixel_data, SEEK_SET) ;
      fread(pixel_data, 1, mri->width * mri->height * 2, fp);
      fclose(fp) ;

#ifdef Linux
      swab(pixel_data, pixel_data, mri->width*mri->height*2);
#endif

      short_buffer_to_image(pixel_data, mri, slice_number - min_slice_number, 0);

      /* untested... */
      if(ras_index < 2)
      {
        normal[ras_index][0] = orderFloatBytes(image_header[142 + 0*4]);
        normal[ras_index][1] = orderFloatBytes(image_header[142 + 1*4]);
        normal[ras_index][2] = orderFloatBytes(image_header[142 + 2*4]);
        top_left[ras_index][0] = orderFloatBytes(image_header[142 + 3*4]);
        top_left[ras_index][1] = orderFloatBytes(image_header[142 + 4*4]);
        top_left[ras_index][2] = orderFloatBytes(image_header[142 + 5*4]);
        top_right[ras_index][0] = orderFloatBytes(image_header[142 + 6*4]);
        top_right[ras_index][1] = orderFloatBytes(image_header[142 + 7*4]);
        top_right[ras_index][2] = orderFloatBytes(image_header[142 + 8*4]);
        bottom_right[ras_index][0] = orderFloatBytes(image_header[142 + 9*4]);
        bottom_right[ras_index][1] = orderFloatBytes(image_header[142 + 10*4]);
        bottom_right[ras_index][2] = orderFloatBytes(image_header[142 + 11*4]);
        ras_index++;
      }
      /* ...to here, and below */

      fclose(fp) ;
    }

    free(pixel_data);

    /* untested... */
    x_vec[0] = (top_right[0][0] - top_left[0][0]);
    x_vec[1] = (top_right[0][1] - top_left[0][1]);
    x_vec[2] = (top_right[0][2] - top_left[0][2]);

    y_vec[0] = (bottom_right[0][0] - top_right[0][0]);
    y_vec[1] = (bottom_right[0][1] - top_right[0][1]);
    y_vec[2] = (bottom_right[0][2] - top_right[0][2]);

    z_vec[0] = top_left[1][0] - top_left[0][0];
    z_vec[1] = top_left[1][1] - top_left[0][1];
    z_vec[2] = top_left[1][2] - top_left[0][2];

    /* determine the axis reordering that's needed */
    /* start with the x axis; is the R, A, or S coordinate the longest? */

    xmax = IMAX3(fabs(x_vec[0]), fabs(x_vec[1]), fabs(x_vec[2]));
    ymax = IMAX3(fabs(y_vec[0]), fabs(y_vec[1]), fabs(y_vec[2]));
    zmax = IMAX3(fabs(z_vec[0]), fabs(z_vec[1]), fabs(z_vec[2]));

    xsign = SIGN(x_vec[xmax]);
    ysign = SIGN(y_vec[ymax]);
    zsign = SIGN(z_vec[zmax]);

    if(xmax == 0)
      mri->xdir = -XDIM * xsign;
    if(xmax == 1)
      mri->xdir = ZDIM * xsign;
    if(xmax == 2)
      mri->xdir = -YDIM * xsign;

    if(ymax == 0)
      mri->ydir = -XDIM * ysign;
    if(ymax == 1)
      mri->ydir = ZDIM * ysign;
    if(ymax == 2)
      mri->ydir = -YDIM * ysign;

    if(zmax == 0)
      mri->zdir = -XDIM * zsign;
    if(zmax == 1)
      mri->zdir = ZDIM * zsign;
    if(zmax == 2)
      mri->zdir = -YDIM * zsign;

    /* sanity check -- have all the axes been used? */
    xmax = abs(xmax + 1);
    ymax = abs(ymax + 1);
    zmax = abs(zmax + 1);

    if(xmax + ymax + zmax != 6 || xmax * ymax * zmax != 9)
      ErrorReturn(NULL, (ERROR_BADPARM, "gelxRead(%s): error interpreting slice direction", fname));
    /* ...to here */

  } /* end if(read_volume) */

  mri->xdir = XDIM;
  mri->ydir = YDIM;
  mri->zdir = ZDIM;

  return(mri);

} /* end gelxRead() */

static MRI *
siemensRead(char *fname, int read_volume, int frame)
{

  MRI *mri = NULL;
  char buf[4];
  FILE *fp;
  int i;
  int width, height, depth;
  char fname_format[STRLEN], *dot, fname_use[STRLEN];
  int initial_image_number, ino_min, ino_max;
  char slice_direction[7];
  short *pixel_data;
  double center_x, center_y, center_z;
  double normal_x, normal_y, normal_z;
  double row_vec_x, row_vec_y, row_vec_z;
  double col_vec_x, col_vec_y, col_vec_z;
  double x_vec[3], y_vec[3], z_vec[3];
  int xmax, ymax, zmax;
  int xsign, ysign, zsign;
  double d;

  frame = 0;

  if((dot = strrchr(fname, '.')) == NULL)
    ErrorReturn(NULL,
                (ERROR_BADPARM, "siemensRead(%s): bad file name", fname));

  if(strcmp(dot, ".ima") != 0)
    ErrorReturn(NULL,
                (ERROR_BADPARM, "siemensRead(%s): bad file name", fname));

  for(dot--;isdigit((int)*dot) && dot > fname;dot--);

  if(dot != fname)
    {
    dot++;
    strncpy(fname_format, fname, dot - fname);
    fname_format[dot-fname] = '\0';
    strcat(fname_format, "%d.ima");
    }
  else
    {
    if(!isdigit((int)*dot))
      dot++;
    strncpy(fname_format, fname, dot - fname);
    fname_format[dot-fname] = '\0';
    strcat(fname_format, "%d.ima");
    }

  initial_image_number = strtol(dot, NULL, 10);

  get_file_limits(fname_format, initial_image_number, &ino_min, &ino_max, 0);

  if((fp = fopen(fname, "r")) == NULL)
    ErrorReturn(NULL, (ERROR_BADFILE, 
                       "siemensRead(%s): error opening file", fname));

  fseek(fp, 2864, SEEK_SET);
  fread(buf, 4, 1, fp);
  height = width = orderIntBytes(*(int *)buf);
  depth = ino_max - ino_min + 1;

  if(read_volume)
    mri = MRIalloc(width, height, depth, MRI_SHORT) ;
  else
    mri = MRIallocHeader(width, height, depth, MRI_SHORT) ;

  fseek(fp, 1544, SEEK_SET);
  fread(&d, 8, 1, fp);
  d = orderDoubleBytes(d);
  mri->thick = d;
  mri->zsize = mri->thick;
  fseek(fp, 5000, SEEK_SET);
  fread(&d, 8, 1, fp);
  d = orderDoubleBytes(d);
  mri->xsize = d;
  fread(&d, 8, 1, fp);
  d = orderDoubleBytes(d);
  mri->ysize = d;

  mri->xend = mri->width/2;
  mri->yend = mri->height/2;
  mri->zend = mri->depth/2;
  mri->xstart = -mri->xend;
  mri->ystart = -mri->yend;
  mri->zstart = -mri->zend;

  mri->fov = 256;

  mri->imnr0 = 1 ;
  mri->imnr1 = ino_max - ino_min + 1;
  mri->ps = 1 /*0.001*/;
  mri->tr = 0 ;
  mri->te = 0 ;
  mri->ti = 0 ;
  strcpy(mri->fname, fname) ;

  fseek(fp, 5814, SEEK_SET);
  fread(slice_direction, 1, 7, fp);

  fseek(fp, 3768, SEEK_SET);
  fread(&center_x, 8, 1, fp);  center_x = orderDoubleBytes(center_x);
  fread(&center_y, 8, 1, fp);  center_y = orderDoubleBytes(center_y);
  fread(&center_z, 8, 1, fp);  center_z = orderDoubleBytes(center_z);
  fread(&normal_x, 8, 1, fp);  normal_x = orderDoubleBytes(normal_x);
  fread(&normal_y, 8, 1, fp);  normal_y = orderDoubleBytes(normal_y);
  fread(&normal_z, 8, 1, fp);  normal_z = orderDoubleBytes(normal_z);
  fseek(fp, 3832, SEEK_SET);
  fread(&row_vec_x, 8, 1, fp);  row_vec_x = orderDoubleBytes(row_vec_x);
  fread(&row_vec_y, 8, 1, fp);  row_vec_y = orderDoubleBytes(row_vec_y);
  fread(&row_vec_z, 8, 1, fp);  row_vec_z = orderDoubleBytes(row_vec_z);
  fread(&col_vec_x, 8, 1, fp);  col_vec_x = orderDoubleBytes(col_vec_x);
  fread(&col_vec_y, 8, 1, fp);  col_vec_y = orderDoubleBytes(col_vec_y);
  fread(&col_vec_z, 8, 1, fp);  col_vec_z = orderDoubleBytes(col_vec_z);

  x_vec[0] = row_vec_x;
  x_vec[1] = row_vec_y;
  x_vec[2] = row_vec_z;

  y_vec[0] = -col_vec_x;
  y_vec[1] = -col_vec_y;
  y_vec[2] = -col_vec_z;

  z_vec[0] = normal_x;
  z_vec[1] = normal_y;
  z_vec[2] = normal_z;

  /* determine the axis reordering that's needed */
  /* start with the x axis; is the R, A, or S coordinate the longest? */

  xmax = IMAX3(fabs(x_vec[0]), fabs(x_vec[1]), fabs(x_vec[2]));
  ymax = IMAX3(fabs(y_vec[0]), fabs(y_vec[1]), fabs(y_vec[2]));
  zmax = IMAX3(fabs(z_vec[0]), fabs(z_vec[1]), fabs(z_vec[2]));

  xsign = SIGN(x_vec[xmax]);
  ysign = SIGN(y_vec[ymax]);
  zsign = SIGN(z_vec[zmax]);

  if(xmax == 0)
    mri->xdir = -XDIM * xsign;
  if(xmax == 1)
    mri->xdir = ZDIM * xsign;
  if(xmax == 2)
    mri->xdir = -YDIM * xsign;

  if(ymax == 0)
    mri->ydir = -XDIM * ysign;
  if(ymax == 1)
    mri->ydir = ZDIM * ysign;
  if(ymax == 2)
    mri->ydir = -YDIM * ysign;

  if(zmax == 0)
    mri->zdir = -XDIM * zsign;
  if(zmax == 1)
    mri->zdir = ZDIM * zsign;
  if(zmax == 2)
    mri->zdir = -YDIM * zsign;

  /* sanity check -- have all the axes been used? */
  xmax = abs(xmax + 1);
  ymax = abs(ymax + 1);
  zmax = abs(zmax + 1);

  /* must be each of 1, 2, and 3 in some order */

  if(xmax + ymax + zmax != 6 || xmax * ymax * zmax != 6)
    ErrorReturn(NULL, (ERROR_BADPARM, "siemensRead(%s): error interpreting slice direction", fname));

  /* assign {xyz}_{ras} coordinates */
  mri->x_r = x_vec[0] * mri->xsize;
  mri->x_a = x_vec[1] * mri->xsize;
  mri->x_s = x_vec[2] * mri->xsize;
  mri->y_r = y_vec[0] * mri->ysize;
  mri->y_a = y_vec[1] * mri->ysize;
  mri->y_s = y_vec[2] * mri->ysize;
  mri->z_r = z_vec[0] * mri->zsize;
  mri->z_a = z_vec[1] * mri->zsize;
  mri->z_s = z_vec[2] * mri->zsize;
  mri->c_r = center_x + mri->depth / 2 * mri->z_r;
  mri->c_a = center_y + mri->depth / 2 * mri->z_a;
  mri->c_s = center_z + mri->depth / 2 * mri->z_s;

  mri->ras_good_flag = 1;

  mri->slice_direction = MRI_UNDEFINED;

  fclose(fp);

  pixel_data = NULL;

  if(read_volume)
    {
    pixel_data = (short *)malloc(mri->height * mri->width * 2);
    for(i = ino_min;i <= ino_max;i++)
      {
      sprintf(fname_use, fname_format, i);
      if((fp = fopen(fname_use, "r")) == NULL)
        ErrorReturn(NULL, (ERROR_BADPARM, 
                   "siemensRead(%s): bad generated file name", fname_use));
      fseek(fp, 6144, SEEK_SET);

      fread(pixel_data, 2, mri->width*mri->height, fp);
#ifdef Linux
swab(pixel_data, pixel_data, mri->height * mri->width * 2);
#endif
      short_buffer_to_image(pixel_data, mri, i - ino_min, frame) ;
      }

    free(pixel_data);

    }

  return(mri);

} /* end siemensRead() */

#define UNUSED_SPACE_SIZE 256
#define MGH_VERSION       1

static MRI *
mghRead(char *fname, int read_volume, int frame)
{
  MRI  *mri ;
  FILE  *fp ;
  int   start_frame, end_frame, width, height, depth, nframes, type, x, y, z,
        bpv, dof, bytes, version, ival ;
  BUFTYPE *buf ;
  char   unused_buf[UNUSED_SPACE_SIZE+1] ;
  float  fval ;
  short  sval ;

  fp = fopen(fname, "rb") ;
  if (!fp)
    ErrorReturn(NULL, (ERROR_BADPARM,"mghRead(%s, %d): could not open file",
                       fname, frame)) ;
  version = freadInt(fp) ;
  width = freadInt(fp) ;
  height = freadInt(fp) ;
  depth =  freadInt(fp) ;
  nframes = freadInt(fp) ;
  type = freadInt(fp) ;
  dof = freadInt(fp) ;

  /* so stuff can be added to the header in the future */
  fread(unused_buf, sizeof(char), UNUSED_SPACE_SIZE, fp) ;

  switch (type)
  {
  default:
  case MRI_FLOAT:  bpv = sizeof(float) ; break ;
  case MRI_UCHAR:  bpv = sizeof(char)  ; break ;
  case MRI_SHORT:  bpv = sizeof(short) ; break ;
  case MRI_INT:     bpv = sizeof(int) ; break ;
  }
  bytes = width * height * bpv ;  /* bytes per slice */
  if (frame >= 0)
  {
    start_frame = end_frame = frame ;
    fseek(fp, frame*width*height*depth*bpv, SEEK_CUR) ;
    nframes = 1 ;
  }
  else
  {  /* hack - # of frames < -1 means to only read in that
        many frames. Otherwise I would have had to change the whole
        MRIread interface and that was too much of a pain. Sorry.
     */
    if (frame < -1)  
    { nframes = frame*-1 ; } 
    start_frame = 0 ; end_frame = nframes-1 ;
  }
  if (!read_volume)
  {
    mri = MRIallocHeader(width, height, depth, type) ;
    mri->dof = dof ;
  }
  else
  {
    if (type == MRI_UCHAR)
      buf = (BUFTYPE *)calloc(bytes, sizeof(BUFTYPE)) ;
    else
      buf = NULL ;
    mri = MRIallocSequence(width, height, depth, type, nframes) ;
    mri->dof = dof ;
    for (frame = start_frame ; frame <= end_frame ; frame++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        switch (type)
        {
          case MRI_INT:
          for (y = 0 ; y < height ; y++)
          {
            for (x = 0 ; x < width ; x++)
            {
              ival = freadInt(fp) ; 
              MRIIseq_vox(mri,x,y,z,frame-start_frame) = ival ;
            }
          }
          break ;
          case MRI_SHORT:
          for (y = 0 ; y < height ; y++)
          {
            for (x = 0 ; x < width ; x++)
            {
              sval = freadShort(fp) ; 
              MRISseq_vox(mri,x,y,z,frame-start_frame) = sval ;
            }
          }
          break ;
          case MRI_FLOAT:
          for (y = 0 ; y < height ; y++)
          {
            for (x = 0 ; x < width ; x++)
            {
              fval = freadFloat(fp) ; 
              MRIFseq_vox(mri,x,y,z,frame-start_frame) = fval ;
            }
          }
          break ;
        case MRI_UCHAR:
          if (fread(buf, sizeof(BUFTYPE), bytes, fp) != bytes)
            ErrorReturn(NULL,
                        (ERROR_BADFILE, "%s: could not read %dth slice (%d)",
                         Progname, z, bytes)) ;
          buffer_to_image(buf, mri, z, frame-start_frame) ;
          break ;
        default:
          ErrorReturn(NULL, 
                      (ERROR_UNSUPPORTED, "mghRead: unsupported type %d",
                       mri->type)) ;
          break ;
        }
      }
    }
    if (buf)
      free(buf) ;
  }

  fclose(fp) ;
  return(mri) ;
}

static int
mghWrite(MRI *mri, char *fname, int frame)
{
  FILE  *fp ;
  int   ival, start_frame, end_frame, x, y, z, width, height, depth ;
  char  buf[UNUSED_SPACE_SIZE+1] ;
  float fval ;
  short sval ;

  if (frame >= 0)
    start_frame = end_frame = frame ;
  else
  {
    start_frame = 0 ; end_frame = mri->nframes-1 ;
  }
  fp = fopen(fname, "wb") ;
  if (!fp)
    ErrorReturn(ERROR_BADPARM, 
                (ERROR_BADPARM,"mghWrite(%s, %d): could not open file",
                 fname, frame)) ;

  /* WARNING - adding or removing anything before nframes will
     cause mghAppend to fail.
  */
  width = mri->width ; height = mri->height ; depth = mri->depth ; 
  fwriteInt(MGH_VERSION, fp) ;
  fwriteInt(mri->width, fp) ;
  fwriteInt(mri->height, fp) ;
  fwriteInt(mri->depth, fp) ;
  fwriteInt(mri->nframes, fp) ;
  fwriteInt(mri->type, fp) ;
  fwriteInt(mri->dof, fp) ;

  /* so stuff can be added to the header in the future */
  memset(buf, 0, UNUSED_SPACE_SIZE*sizeof(char)) ;
  fwrite(buf, sizeof(char), UNUSED_SPACE_SIZE, fp) ;

  for (frame = start_frame ; frame <= end_frame ; frame++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        switch (mri->type)
        {
        case MRI_SHORT:
          for (x = 0 ; x < width ; x++)
          {
            if (z == 74 && y == 16 && x == 53)
              DiagBreak() ;
            sval = MRISseq_vox(mri,x,y,z,frame) ;
            fwriteShort(sval, fp) ;
          }
          break ;
        case MRI_INT:
          for (x = 0 ; x < width ; x++)
          {
            if (z == 74 && y == 16 && x == 53)
              DiagBreak() ;
            ival = MRIIseq_vox(mri,x,y,z,frame) ;
            fwriteInt(ival, fp) ;
          }
          break ;
        case MRI_FLOAT:
          for (x = 0 ; x < width ; x++)
          {
            if (z == 74 && y == 16 && x == 53)
              DiagBreak() ;
            fval = MRIFseq_vox(mri,x,y,z,frame) ;
            fwriteFloat(fval, fp) ;
          }
          break ;
        case MRI_UCHAR:
          if (fwrite(&MRIseq_vox(mri,0,y,z,frame), sizeof(BUFTYPE), width, fp) 
              != width)
            ErrorReturn(ERROR_BADFILE, 
                        (ERROR_BADFILE, 
                         "mghWrite: could not write %d bytes to %s",
                         width, fname)) ;
          break ;
        default:
          ErrorReturn(ERROR_UNSUPPORTED, 
                      (ERROR_UNSUPPORTED, "mghWrite: unsupported type %d",
                       mri->type)) ;
          break ;
        }
      }
    }
  }

  fclose(fp) ;
  return(NO_ERROR) ;
}
static int
mghAppend(MRI *mri, char *fname, int frame)
{
  FILE  *fp ;
  int   start_frame, end_frame, x, y, z, width, height, depth, nframes ;

  if (frame >= 0)
    start_frame = end_frame = frame ;
  else
  {
    start_frame = 0 ; end_frame = mri->nframes-1 ;
  }
  fp = fopen(fname, "rb") ;
  if (!fp)   /* doesn't exist */
    return(mghWrite(mri, fname, frame)) ;
  fclose(fp) ;
  fp = fopen(fname, "r+b") ;
  if (!fp)
    ErrorReturn(ERROR_BADPARM, 
                (ERROR_BADPARM,"mghAppend(%s, %d): could not open file",
                 fname, frame)) ;

  /* WARNING - this is dependent on the order of writing in mghWrite */
  width = mri->width ; height = mri->height ; depth = mri->depth ;
  fseek(fp, 4*sizeof(int), SEEK_SET) ;
  nframes = freadInt(fp) ;
  fseek(fp, 4*sizeof(int), SEEK_SET) ;
  fwriteInt(nframes+end_frame-start_frame+1, fp) ;
  fseek(fp, 0, SEEK_END) ;

  for (frame = start_frame ; frame <= end_frame ; frame++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        switch (mri->type)
        {
        case MRI_FLOAT:
          for (x = 0 ; x < width ; x++)
          {
            fwriteFloat(MRIFseq_vox(mri,x,y,z,frame), fp) ;
         ght]), 2, mri->wak ;
        case MRI_UCHAR:
          if (fwrite(&MRIseq_vox(mri,0,y,z,frame), sizeof(BUFTYPE), width, fp) 
              != width)
            ErrorReturn(ERROR_BADFILE, 
                        (ERROR_BADFILE, 
                         "mghAppend: could not write %d bytes to %s",
                         width, fname)) ;
          break ;
        default:
          ErrorReturn(ERROR_UNSUPPORTED, 
                      (ERROR_UNSUPPORTED, "mghAppend: unsupported type %d",
                       mri->type)) ;
          break ;
        }
      }
    }
d2440 5
  }

  /* ----- allocate the mri structure ----- */
  if(read_volume)
    mri = MRIalloc(hdr.dime.dim[1], hdr.dime.dim[2], hdr.dime.dim[3], dtype);
  else
    mri = MRIalloc(hdr.dime.dim[1], hdr.dime.dim[2], hdr.dime.dim[3], dtype);
  fclose(fp) ;
  return(NO_ERROR) ;
}

static MRI *
brikRead(char *fname, int read_volume, int frame)
{

  FILE *fin;
  char *dot;
  char header_name[STRLEN];
  char header_type_line[STRLEN], header_name_line[STRLEN], header_count_line[STRLEN];
  char header_line[STRLEN];
  char name[STRLEN], type[STRLEN];
  int count;
  MRI *mri, *mri2;
  char byteorder_string[STRLEN], *tilde;
  int brick_type;
  char *buf;
  int i, j, k;
  char swapbuf[4];
  int orient_vals[3];
  int orient_code[] = { -XDIM, XDIM, ZDIM, -ZDIM, -YDIM, YDIM };
  char *c;

  frame = 0;

  strcpy(header_name, fname);
  if((dot = strrchr(header_name, '.')) == NULL)
    ErrorReturn(NULL, (ERROR_BADPARM, 
                "brikRead(): can't find '.' in filename: %s", fname));

  dot++;
  sprintf(dot, "HEAD");
  if((fin = fopen(header_name, "r")) == NULL)
    ErrorReturn(NULL, (ERROR_BADFILE, "can't open file %s", header_name));

  mri = MRIallocHeader(-1, -1, -1, MRI_SHORT);

  orient_vals[0] = orient_vals[1] = orient_vals[2] = -1;

  mri->xsize = mri->ysize = mri->zsize = 1.0;
  brick_type = 1; /* default to short */

  while(!feof(fin))
  {
    header_type_line[0] = '\0';
    while(!feof(fin) && strncmp(header_type_line, "type", 4))
      fgets(header_type_line, STRLEN, fin);
    if(!feof(fin))
    {
      header_name_line[0] = '\0';
      fgets(header_name_line, STRLEN, fin);
      if(strncmp(header_name_line, "name", 4))
        ErrorReturn(NULL, (ERROR_BADFILE, 
                    "readBrik: error in BRIK header file %s", header_name));
      header_count_line[0] = '\0';
      fgets(header_count_line, STRLEN, fin);
      if(strncmp(header_count_line, "count", 5))
        ErrorReturn(NULL, (ERROR_BADFILE, 
                    "readBrik: error in BRIK header file %s", header_name));
      sscanf(header_type_line, "type = %s", type);
      sscanf(header_name_line, "name = %s", name);
      sscanf(header_count_line, "count = %d", &count);

      if(strcmp(name, "DELTA") == 0)
      {
        fgets(header_line, STRLEN, fin);
        sscanf(header_line, "%f %f %f", &(mri->xsize), &(mri->ysize), &(mri->zsize));
      }
      else if(strcmp(name, "DATASET_DIMENSIONS") == 0)
      {
        fgets(header_line, STRLEN, fin);
        sscanf(header_line, "%d %d %d", &(mri->width), &(mri->height), &(mri->depth));
      }
      else if(strcmp(name, "BRICK_TYPES") == 0)
      {
        fgets(header_line, STRLEN, fin);
        sscanf(header_line, "%d", &brick_type);
      }
      else if(strcmp(name, "BYTEORDER_STRING") == 0)
      {
        fgets(byteorder_string, STRLEN, fin);
        for(tilde = byteorder_string;
            (tilde = strchr(byteorder_string, '~')) != NULL;)
          *tilde = '\0';
      }
      else if(strcmp(name, "ORIENT_SPECIFIC") == 0)
      {
        fgets(header_line, STRLEN, fin);
        sscanf(header_line, "%d %d %d", &orient_vals[0], &orient_vals[1], &orient_vals[2]);
      }

    }
  }

  fclose(fin);

  if(orient_vals[0] == -1 || orient_vals[1] == -1 || orient_vals[2] == -1)
    ErrorReturn(NULL, (ERROR_BADFILE, "brikRead: error getting image orientation"));

  if(mri->width < 0 || mri->height < 0 || mri->depth < 0)
    ErrorReturn(NULL, (ERROR_BADFILE, "brikRead: error settting image dimensions"));

/* from afni source: mrilib.h:
   { "byte" , "short" , "int" , "float" , "double" , "complex" , "rgb" }
*/
  if(brick_type == 0)
    mri->type = MRI_UCHAR;
  else if(brick_type == 1)
    mri->type = MRI_SHORT;
  else if(brick_type == 2)
    mri->type = MRI_INT;
  else if(brick_type == 3)
    mri->type = MRI_FLOAT;
  else
    ErrorReturn(NULL, (ERROR_BADFILE, 
                "bad data type %d in file %s", brick_type, fname));

  mri->imnr0 = 1;
  mri->imnr1 = mri->depth;
  mri->thick = mri->zsize;
  mri->xend = (mri->xsize * mri->width) / 2.;
  mri->yend = (mri->ysize * mri->height) / 2.;
  mri->zend = (mri->zsize * mri->depth) / 2.;
  mri->xstart = -mri->xend;
  mri->ystart = -mri->yend;
  mri->zstart = -mri->zend;
  mri->fov = mri->xend - mri->xstart;
  strcpy(mri->fname, fname);

  buf = NULL;

  if(!read_volume)
  {
    mri2 = MRIallocHeader(mri->width, mri->height, mri->depth, mri->type);
    MRIfree(&mri);
  }
  else
  {
    if((fin = fopen(fname, "r")) == NULL)
      ErrorReturn(NULL, (ERROR_BADPARM, 
                         "can't open file %s", fname));

    buf = (char *)malloc(mri->width * mri->height * mri->depth * data_size[mri->type]);

    if(fread(buf, data_size[mri->type], mri->width * mri->height * mri->depth, fin) < mri->width * mri->height * mri->depth)
      ErrorReturn(NULL, (ERROR_BADFILE, "error reading from file %s", fname));

    fclose(fin);

    for(c = byteorder_string;*c == ' ' || *c == '\'';c++);

#ifdef Linux
    if(strncmp(c, "MSB", 3) == 0)
    {
      if(data_size[mri->type] == 2)
        swab(buf, buf, mri->width * mri->height * mri->depth * data_size[mri->type]);
      else if(data_size[mri->type] == 4)
        for(i = 0;i < mri->width * mri->height * mri->depth * data_size[mri->type];i+=4)
        {
          memcpy(swapbuf, &buf[i], 4);
          buf[i + 0] = swapbuf[3];
          buf[i + 1] = swapbuf[2];
          buf[i + 2] = swapbuf[1];
          buf[i + 3] = swapbuf[0];
        }
    }
#else
    if(strncmp(c, "LSB", 3) == 0)
    {
      if(data_size[mri->type] == 2)
        swab(buf, buf, mri->width * mri->height * mri->depth * data_size[mri->type]);
      else if(data_size[mri->type] == 4)
        for(i = 0;i < mri->width * mri->height * mri->depth * data_size[mri->type];i+=4)
        {
          memcpy(swapbuf, &buf[i], 4);
          buf[i + 0] = swapbuf[3];
          buf[i + 1] = swapbuf[2];
          buf[i + 2] = swapbuf[1];
          buf[i + 3] = swapbuf[0];
        }
    }
#endif

    mri2 = MRIalloc(mri->width, mri->height, mri->depth, mri->type);

    if(mri2->type == MRI_UCHAR)
    {
      for(i = 0;i < mri->depth;i++)
        for(j = 0;j < mri->height;j++)
          for(k = 0;k < mri->width;k++)
            MRIvox(mri2, k, j, i) = buf[k + j*mri->width + i*mri->width*mri->height];
    }
    else if(mri2->type == MRI_SHORT)
    {
      for(i = 0;i < mri->depth;i++)
        for(j = 0;j < mri->height;j++)
          for(k = 0;k < mri->width;k++)
            MRISvox(mri2, k, j, i) = ((short *)buf)[k + j*mri->width + i*mri->width*mri->height];
    }
    else
      ErrorReturn(NULL, (ERROR_UNSUPPORTED, "brikRead(): unsupported type %d", mri->type));

    MRIfree(&mri);

    free(buf);

  }

  mri2->xdir = orient_code[orient_vals[0]];
  mri2->ydir = orient_code[orient_vals[1]];
  mri2->zdir = orient_code[orient_vals[2]];

  strcpy(mri2->fname, fname);

  return(mri2);

} /* end brikRead() */

static int
brikWrite(MRI *mri, char *fname, int frame)
{

  char header_name[STRLEN];
  char *dot;
  FILE *fp;
  int i;
  short *buf;
  int orient_code[3], mri_dir[3];

  strcpy(header_name, fname);
  if((dot = strrchr(header_name, '.')) == NULL)
    ErrorReturn(0, (ERROR_BADPARM, 
                "brikWrite(): can't find '.' in filename: %s", fname));

  dot++;
  sprintf(dot, "HEAD");

  if((fp = fopen(fname, "w")) == NULL)
  {
    ErrorReturn(0, (ERROR_BADFILE, "brikWrite(): can't open file %s", fname));
  }
  buf = (short *)malloc(mri->width * mri->height * sizeof(short));
  for(i = 0;i <= mri->imnr1 - mri->imnr0;i++)
  {
    image_to_short_buffer(buf, mri, i);
#ifdef Linux
    swab(buf, buf, mri->height * mri->width * sizeof(short));
#endif
    fwrite(buf, 2, mri->width * mri->height, fp);
  }
  free(buf);
  fclose(fp);

  if((fp = fopen(header_name, "w")) == NULL)
  {
    ErrorReturn(0, (ERROR_BADFILE, "brikWrite(): can't open file %s", header_name));
  }

  fprintf(fp, "type = float-attribute\n");
  fprintf(fp, "name = DELTA\n");
  fprintf(fp, "count = 3\n");
  fprintf(fp, "  %g      %g      %g\n\n", mri->xsize, mri->ysize, mri->zsize);

  fprintf(fp, "type = integer-attribute\n");
  fprintf(fp, "name = DATASET_DIMENSIONS\n");
  fprintf(fp, "count = 5\n");
  fprintf(fp, "  %d %d %d %d %d\n\n", mri->width, mri->height, mri->depth, 0, 0);

  fprintf(fp, "type = integer-attribute\n");
  fprintf(fp, "name = BRICK_TYPES\n");
  fprintf(fp, "count = 1\n");
  fprintf(fp, "  1\n\n");

  fprintf(fp, "type = string-attribute\n");
  fprintf(fp, "name = BYTEORDER_STRING\n");
  fprintf(fp, "count = 10\n");
  fprintf(fp, "  'MSB_FIRST~\n\n");

  mri_dir[0] = mri->xdir;
  mri_dir[1] = mri->ydir;
  mri_dir[2] = mri->zdir;

  for(i = 0;i < 3;i++)
  {
    if(mri_dir[i] == -XDIM)
      orient_code[i] = 0;
    if(mri_dir[i] == XDIM)
      orient_code[i] = 1;
    if(mri_dir[i] == -ZDIM)
      orient_code[i] = 2;
    if(mri_dir[i] == -ZDIM)
      orient_code[i] = 3;
    if(mri_dir[i] == -YDIM)
      orient_code[i] = 4;
    if(mri_dir[i] == YDIM)
      orient_code[i] = 5;
  }

  fprintf(fp, "type = integer-attribute\n");
  fprintf(fp, "name = ORIENT_SPECIFIC\n");
  fprintf(fp, "count = 3\n");
  fprintf(fp, "  %d %d %d\n\n", orient_code[0], orient_code[1], orient_code[2]);

  fclose(fp);

  return(1);

} /* end brikWrite() */

static MRI *
bshortRead(char *fname, int read_volume, int frame)
{

  FILE *fin;
  char *error, *dot, *num;
  int nstart, nend, norig, i, keep_going;
  char fname2[STRLEN], fname_passed[STRLEN];
  int width, height, frames, swap;
  MRI *mri;
  short *buf;

  frame = 0;

  strcpy(fname_passed, fname);

  if((dot = strrchr(fname, '.')) == NULL)
    ErrorReturn(NULL, (ERROR_BADPARM, "bshortRead: bad bshort file name %s", fname));

  num = dot-3;

  if(num < fname)
    ErrorReturn(NULL, (ERROR_BADPARM, "bshortRead: bad bshort file name %s", fname));

  norig = strtol(num, &error, 10);
  if(error == NULL)
    ErrorReturn(NULL, (ERROR_BADPARM, "bshortRead: bad bshort file name %s", fname));

  *num = '\0';

  /* assume fname is there -- now find the first and last images */

  for(i = norig,keep_going = 1;keep_going;)
  {
    i--;
    sprintf(fname2, "%s%03d.bshort", fname, i);
    if((fin = fopen(fname2, "r")) == NULL)
      keep_going = 0;
    else
      fclose(fin);
  }

  nstart = i + 2;
  
  for(i = norig,keep_going = 1;keep_going;)
  {
    i++;
    sprintf(fname2, "%s%03d.bshort", fname, i);
    if((fin = fopen(fname2, "r")) == NULL)
      keep_going = 0;
    else
      fclose(fin);
  }

  nend = i - 1;

  sprintf(fname2, "%s%03d.hdr", fname, nstart);
  if((fin = fopen(fname2, "r")) == NULL)
    ErrorReturn(NULL, (ERROR_BADPARM, "bshortRead: can't open file %s", fname2));
  fscanf(fin, "%d %d %d %d", &width, &height, &frames, &swap);
  fclose(fin);

  buf = (short *)malloc(height * width * 2);

  /* support multiple images within a bshort only if there's just the one file */
  if(nstart == nend)
  {
    sprintf(fname2, "%s%03d.bshort", fname, nstart);
    if((fin = fopen(fname2, "r")) == NULL)
      ErrorReturn(NULL, (ERROR_BADPARM, "bshortRead: can't open file %s", fname2));

    if(!read_volume)
      mri = MRIalloc(height, width, frames, MRI_SHORT);
    else
    {
      mri = MRIallocHeader(height, width, frames, MRI_SHORT);
      for(i = 0;i < frames;i++)
      {
        fread(buf, height * width, 2, fin);
#ifdef Linux
        if(!swap)
          swab(buf, buf, height * width * 2);
#else
        if(swap)
          swab((const void *)buf, (void *)buf, height * width * 2);
#endif
        short_buffer_to_image((short *)buf, mri, i, frame) ;
      }

    }

    fclose(fin);

  }
  else
  {
    if(!read_volume)
      mri = MRIallocHeader(height, width, nend - nstart + 1, MRI_SHORT);
    else
    {
      mri = MRIalloc(height, width, nend - nstart + 1, MRI_SHORT);

      for(i = nstart;i <= nend;i++)
      {
        sprintf(fname2, "%s%03d.bshort", fname, i);
        if((fin = fopen(fname2, "r")) == NULL)
          ErrorReturn(NULL, (ERROR_BADPARM, "bshortRead: can't open file %s", fname2));

        fread(buf, height * width, 2, fin);
#ifdef Linux
        if(!swap)
          swab(buf, buf, height * width * 2);
#else
        if(swap)
          swab((const void *)buf, (void *)buf, height * width * 2);
#endif
        short_buffer_to_image((short *)buf, mri, i-nstart, frame) ;

        fclose(fin);
      }
    }

  }

  mri->xdir = XDIM;
  mri->ydir = YDIM;
  mri->zdir = ZDIM;

  strcpy(mri->fname, fname_passed);

  return(mri);

} /* end bshortRead() */

static int
bshortWrite(MRI *mri, char *fname, int frame)
{

  FILE *fp;
  int i;
  char *dot, *num, *hnum, *error;
  int norig;
  char header_fname[STRLEN];
  short *buf;

  if((dot = strrchr(fname, '.')) == NULL)
    ErrorReturn(0, (ERROR_BADPARM, "bshortWrite(): bad bshort file name %s", fname));

  if(strcmp(dot+1, "bshort"))
    ErrorReturn(0, (ERROR_BADPARM, "bshortWrite(): bad bshort file name %s", fname));

  strcpy(header_fname, fname);
  sprintf((dot - fname) + header_fname, ".hdr");

  num = dot-3;

  if(num < fname)
    ErrorReturn(0, (ERROR_BADPARM, "bshortWrite(): bad bshort file name %s", fname));

  norig = strtol(num, &error, 10);
  if(error == NULL)
    ErrorReturn(0, (ERROR_BADPARM, "bshortWrite(): bad bshort file name %s", fname));

  hnum = (num - fname) + header_fname;

  buf = (short *)malloc(mri->height * mri->width * sizeof(short));

  for(i = mri->imnr0;i <= mri->imnr1;i++)
  {
    sprintf(num, "%03d", i - mri->imnr0 + norig);
    num[3] = '.';
    if((fp = fopen(fname, "w")) == NULL)
    {
      ErrorReturn(0, (ERROR_BADFILE,
                      "bshortWrite(): can't open file %s for writing", fname));
    }
    image_to_short_buffer(buf, mri, i - mri->imnr0);
#ifdef Linux
    swab(buf, buf, mri->height * mri->width * sizeof(short));
#endif
    fwrite(buf, 2, mri->width * mri->height, fp);
    fclose(fp);

    sprintf(hnum, "%03d", i - mri->imnr0 + norig);
    hnum[3] = '.';
    if((fp = fopen(header_fname, "w")) == NULL)
    {
      ErrorReturn(0, (ERROR_BADFILE,
                      "bshortWrite(): can't open file %s for writing", fname));
    }
    fprintf(fp, "%d %d %d %d\n", mri->height, mri->width, 1, 0);
    fclose(fp);

  }

  free(buf);

  return(1);

} /* end bshortWrite() */

static MRI *sdtRead(char *fname, int read_volume, int frame)
{

  char header_fname[STR_LEN];
  char line[STR_LEN];
  char *colon, *dot;
  FILE *fp;
  MRI *mri;
  int ndim = -1, data_type = -1;
  int dim[4];
  float xsize = 1.0, ysize = 1.0, zsize = 1.0, dummy_size;
  int orientation = MRI_CORONAL;

  dim[0] = -1;
  dim[1] = -1;
  dim[2] = -1;
  dim[3] = -1;

  /* form the header file name */
  strcpy(header_fname, fname);
  if((dot = strrchr(header_fname, '.')))
    sprintf(dot+1, "spr");
  else
    strcat(header_fname, ".spr");

  /* open the header */
  if((fp = fopen(header_fname, "r")) == NULL)
    ErrorReturn(NULL, (ERROR_BADFILE, "sdtRead(%s): could not open header file %s\n", fname, header_fname));

  while(!feof(fp))
  {
    fgets(line, STR_LEN, fp);
    if((colon = strchr(line, ':')))
    {
      *colon = '\0';
      colon++;
      if(strcmp(line, "numDim") == 0)
      {
        sscanf(colon, "%d", &ndim);
        if(ndim < 3 || ndim > 4)
          {
          fclose(fp);
          ErrorReturn(NULL, (ERROR_UNSUPPORTED, "sdtRead(%s): only 3 or 4 dimensions supported (numDim = %d)\n", fname, ndim));
          }
      }
      else if(strcmp(line, "dim") == 0)
      {
        if(ndim == -1)
        {
          fclose(fp);
          ErrorReturn(NULL, (ERROR_BADFILE, "sdtRead(%s): 'dim' before 'numDim' in header file %s\n", fname, header_fname));
        }
        if(ndim == 3)
        {
          sscanf(colon, "%d %d %d", &dim[0], &dim[1], &dim[2]);
          dim[3] = 1;
        }
        else
        {
          sscanf(colon, "%d %d %d %d", &dim[0], &dim[1], &dim[2], &dim[3]);
          if(dim[3] != 1)
          {
            fclose(fp);
            ErrorReturn(NULL, (ERROR_UNSUPPORTED, "sdtRead(%s): nframes != 1 unsupported for sdt (dim(4) = %d)\n", fname, dim[3]));
          }
        }
      }
      else if(strcmp(line, "dataType") == 0)
      {
        while(isspace((int)*colon))
          colon++;
        if(strncmp(colon, "BYTE", 4) == 0)
          data_type = MRI_UCHAR;
        else if(strncmp(colon, "WORD", 4) == 0)
          data_type = MRI_SHORT;
        else if(strncmp(colon, "LWORD", 5) == 0)
          data_type = MRI_INT;
        else if(strncmp(colon, "REAL", 4) == 0)
          data_type = MRI_FLOAT;
        else if(strncmp(colon, "COMPLEX", 7) == 0)
        {
          fclose(fp);
          ErrorReturn(NULL, (ERROR_UNSUPPORTED, "sdtRead(%s): unsupported data type '%s'\n", fname, colon));
        }
        else
        {
          fclose(fp);
          ErrorReturn(NULL, (ERROR_BADFILE, "sdtRead(%s): unknown data type '%s'\n", fname, colon));
        }
      }
      else if(strcmp(line, "interval") == 0)
      {
        if(ndim == 3)
          sscanf(colon, "%f %f %f", &xsize, &ysize, &zsize);
        else
          sscanf(colon, "%f %f %f %f", &xsize, &ysize, &zsize, &dummy_size);
        xsize *= 10.0;
        ysize *= 10.0;
        zsize *= 10.0;
      }
      else if(strcmp(line, "sdtOrient") == 0)
      {
        while(isspace((int)*colon))
          colon++;
        if(strncmp(colon, "sag", 3) == 0)
          orientation = MRI_SAGITTAL;
        else if(strncmp(colon, "ax", 2) == 0)
          orientation = MRI_HORIZONTAL;
        else if(strncmp(colon, "cor", 3) == 0)
          orientation = MRI_CORONAL;
        else
        {
          fclose(fp);
          ErrorReturn(NULL, (ERROR_BADFILE, "sdtRead(%s): unknown orientation %s\n", fname, colon));
        }
      }
      else
      {
      }
    }
  }

  fclose(fp);

  if(data_type == -1)
    ErrorReturn(NULL, (ERROR_BADFILE, "sdtRead(%s): data type undefined\n", fname));
  if(dim[0] == -1 || dim[1] == -1 || dim[2] == -1)
    ErrorReturn(NULL, (ERROR_BADFILE, "sdtRead(%s): one or more dimensions undefined\n", fname));

  if(read_volume)
  {
    if((fp = fopen(fname, "r")) == NULL)
      ErrorReturn(NULL, (ERROR_BADFILE, "sdtRead(%s): error opening data file %s\n", fname, fname));

    mri = MRIreadRaw(fp, dim[0], dim[1], dim[2], data_type);

    if(mri == NULL)
      return(NULL);

    fclose(fp);

  }
  else
  {
    mri = MRIallocHeader(dim[0], dim[1], dim[2], data_type);
    if(mri == NULL)
      return(NULL);
  }

  mri->xsize = xsize;
  mri->ysize = ysize;
  mri->zsize = zsize;

  mri->slice_direction = orientation;

  if(orientation == MRI_CORONAL)
  {
    mri->xdir = XDIM;
    mri->ydir = YDIM;
    mri->zdir = ZDIM;
  }
  if(orientation == MRI_SAGITTAL)
  {
    mri->xdir = ZDIM;
    mri->ydir = YDIM;
    mri->zdir = XDIM;
  }
  if(orientation == MRI_HORIZONTAL)
  {
    mri->xdir = XDIM;
    mri->ydir = -ZDIM;
    mri->zdir = YDIM;

  }

  mri->thick = mri->zsize;
  mri->xend = mri->xsize * mri->width / 2.;
  mri->yend = mri->ysize * mri->height / 2.;
  mri->zend = mri->zsize * mri->depth / 2.;
  mri->xstart = -mri->xend;
  mri->ystart = -mri->yend;
  mri->zstart = -mri->zend;

  mri->imnr0 = 1;
  mri->imnr1 = dim[2];

  mri->ps = 1.0 /*0.001*/;
  mri->tr = 0 ;
  mri->te = 0 ;
  mri->ti = 0 ;

  strcpy(mri->fname, fname) ;

  return(mri);

} /* end sdtRead() */

static void flipAnalyzeHeader(dsr *header)
{

  int i;
  char c;

  header->hk.sizeof_hdr = swapInt(header->hk.sizeof_hdr);
  header->hk.extents = swapInt(header->hk.extents);
  header->hk.session_error = swapShort(header->hk.session_error);

  for(i = 0;i < 8;i++)
  {
    header->dime.dim[i] = swapShort(header->dime.dim[i]);
    header->dime.pixdim[i] = swapFloat(header->dime.pixdim[i]);
  }

  header->dime.unused1 = swapShort(header->dime.unused1);
  header->dime.datatype = swapShort(header->dime.datatype);
  header->dime.bitpix = swapShort(header->dime.bitpix);
  header->dime.dim_un0 = swapShort(header->dime.dim_un0);
  header->dime.vox_offset = swapFloat(header->dime.vox_offset);
  header->dime.roi_scale = swapFloat(header->dime.roi_scale);
  header->dime.funused1 = swapFloat(header->dime.funused1);
  header->dime.funused2 = swapFloat(header->dime.funused2);
  header->dime.cal_max = swapFloat(header->dime.cal_max);
  header->dime.cal_min = swapFloat(header->dime.cal_min);
  header->dime.compressed = swapInt(header->dime.compressed);
  header->dime.verified = swapInt(header->dime.verified);
  header->dime.glmax = swapInt(header->dime.glmax);
  header->dime.glmin = swapInt(header->dime.glmin);

  header->hist.views = swapInt(header->hist.views);
  header->hist.vols_added = swapInt(header->hist.vols_added);
  header->hist.start_field = swapInt(header->hist.start_field);
  header->hist.field_skip = swapInt(header->hist.field_skip);
  header->hist.omax = swapInt(header->hist.omax);
  header->hist.omin = swapInt(header->hist.omin);
  header->hist.smax = swapInt(header->hist.smax);
  header->hist.smin = swapInt(header->hist.smin);

  c = header->hist.originator[0]; header->hist.originator[0] = header->hist.originator[1]; header->hist.originator[1] = c;
  c = header->hist.originator[2]; header->hist.originator[2] = header->hist.originator[3]; header->hist.originator[3] = c;
  c = header->hist.originator[4]; header->hist.originator[4] = header->hist.originator[5]; header->hist.originator[5] = c;
  c = header->hist.originator[6]; header->hist.originator[6] = header->hist.originator[7]; header->hist.originator[7] = c;
  c = header->hist.originator[8]; header->hist.originator[8] = header->hist.originator[9]; header->hist.originator[9] = c;

}  /*  end flipAnalyzeHeader()  */

/* EOF */
