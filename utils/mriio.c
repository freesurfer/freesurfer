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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <memory.h>

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

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/

#define MM_PER_METER  1000.0f
#define INFO_FNAME    "COR-.info"
#define MAX_DIM       4

/*-----------------------------------------------------
                    STATIC DATA
-------------------------------------------------------*/

static char MIfspace[] = "frame space" ;

/*-----------------------------------------------------
                    STATIC PROTOTYPES
-------------------------------------------------------*/

static void buffer_to_image(BUFTYPE *buf, MRI *mri, int slice, int frame) ;
static void image_to_buffer(BUFTYPE *buf, MRI *mri, int slice) ;
static MRI *mncRead(char *fname, int read_volume, int frame) ;
static MRI *mghRead(char *fname, int read_volume, int frame) ;
static int mghWrite(MRI *mri, char *fname, int frame) ;
static int mghAppend(MRI *mri, char *fname, int frame) ;
static int mncWrite(MRI *mri, char *fname, int frame) ;
static MRI *analyzeRead(char *fname, int read_volume, int frame) ;
static int analyzeWrite(MRI *mri, char *fname, int frame) ;
static int read_analyze_header(char *fpref, dsr *bufptr) ;
static int write_analyze_header(char *fpref, MRI *mri) ;
static int write_analyze_image(char *fpref, MRI *mri) ;
static int read_float_analyze_image(char *fname, MRI *mri, dsr *hdr) ;
static int read_byte_analyze_image(char *fname, MRI *mri, dsr *hdr) ;

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
  char *number, *dot, buf[100] ;

  strcpy(outFname, inFname) ;
  number = strrchr(outFname, '#') ;
  dot = strrchr(outFname, '.') ;

  if (number)   /* '#' in filename indicates frame # */
  {
    if (sscanf(number+1, "%d", pframe) < 1)
      *pframe = -1 ;
    *number = 0 ;
  }
  else
    *pframe = -1 ;

  if (dot)
  {
    dot = StrUpper(strcpy(buf, dot+1)) ;
    if (!strcmp(dot, "MNC"))
      *ptype = MRI_MINC_FILE ;
    else if (!strcmp(dot, "MINC"))
      *ptype = MRI_MINC_FILE ;
    else if (!strcmp(dot, "MGH"))
      *ptype = MRI_MGH_FILE ;
    else if (!strcmp(dot, "IMG"))
      *ptype = MRI_ANALYZE_FILE ;
    else
      *ptype = MRI_CORONAL_SLICE_DIRECTORY ;
  }
  else
    *ptype = MRI_CORONAL_SLICE_DIRECTORY ;

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
  int     slice, bytes ;

  mri = MRIalloc(width, height, depth, type) ;
  if (!mri)
    return(NULL) ;

  bytes = width*height ;
  buf = (BUFTYPE *)calloc(bytes, sizeof(BUFTYPE)) ;

  /* every width x height bytes should be another slice */
  for (slice = 0 ; slice < depth ; slice++)
  {
    if (fread(buf, sizeof(BUFTYPE), bytes, fp) != bytes)
      ErrorReturn(NULL,
                  (ERROR_BADFILE, "%s: could not read %dth slice (%d)",
                   Progname, slice, bytes)) ;
    buffer_to_image(buf, mri, slice, 0) ;
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

  MRIunpackFileName(fpref, &frame, &type, fname) ;
  switch (type)
  {
  case MRI_MGH_FILE:
    mri = mghRead(fname, 1, frame) ;
    if (!mri)
      return(NULL) ;
    break ;
  case MRI_MINC_FILE:
    mri = mncRead(fname, 1, frame) ;
    if (!mri)
      return(NULL) ;
    break ;
  case MRI_ANALYZE_FILE:
    mri = analyzeRead(fname, 1, frame) ;
    if (!mri)
      return(NULL) ;
    break ;
  default:   /* coronal slice data */
    mri = MRIreadInfo(fpref) ;
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
  char    cmd[100], fname[200];
  int     imnr0, imnr1, width, height, depth, ptype, type, frame ;


  MRIunpackFileName(fpref, &frame, &type, fname) ;
  switch (type)
  {
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
  default:   /* coronal slice data */
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
    if (fscanf(fp,"%s %s", cmd, mri->transform_fname) == 2)
    {
      if (!stricmp(cmd, "xform") || !stricmp(cmd, "transform"))
      {
        if (*mri->transform_fname != '/') /* relative path, add prefix */
          sprintf(fname, "%s/%s", fpref, mri->transform_fname) ;
        else
          strcpy(fname, mri->transform_fname) ; /* absolute path */
        FileNameAbsolute(fname, mri->transform_fname) ;
        if (!FileExists(mri->transform_fname))  /* try typical location */
          sprintf(mri->transform_fname,"%s/../transforms/talairach.xfm",fpref);

        if (input_transform_file(mri->transform_fname, &mri->transform) != OK)
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
        mri->linear_transform = NULL ;
    }
    else
      mri->linear_transform = NULL ;
    
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
  char     fname[100] ;

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
  char    fname[200];

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
    char fname[100] ;

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
  char *dot, buf[100], *number ;

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
  int           type, sgned, error, start_frame, end_frame ;

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

  if (mri->slice_direction != MRI_CORONAL)
    ErrorReturn(ERROR_UNSUPPORTED, 
            (ERROR_UNSUPPORTED,"analyzeWrite: unsupported slice direction %d", 
             mri->slice_direction)) ;

  switch (mri->type)
  {
  case MRI_UCHAR: type = DT_UNSIGNED_CHAR ; break ;
  case MRI_FLOAT: type = DT_FLOAT ;         break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED, 
                (ERROR_UNSUPPORTED, "analyzeWrite: unsupported MRI type %d",
                 mri->type)) ;
    break ;
  }

  write_analyze_header(fname, mri) ;
  write_analyze_image(fname, mri) ;

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
  MRI    *mri, *mri_dst ;
  dsr    hdr;
  int    type, max_dim, width, height, depth ;

  read_analyze_header(fname, &hdr) ;
  switch (hdr.dime.datatype)
  {
  case DT_SIGNED_SHORT:
  case DT_UNSIGNED_CHAR:
    type = MRI_UCHAR ;
    break ;
  case DT_SIGNED_INT:
  case DT_FLOAT:
    type = MRI_FLOAT ;
    break ;
  default:
    ErrorReturn(NULL, 
                (ERROR_UNSUPPORTED, "analyzeRead: unsupported data type %d",
                 hdr.dime.datatype)) ;
    break ;
  }
  width = hdr.dime.dim[3] ; height = hdr.dime.dim[2] ;depth = hdr.dime.dim[1] ;
  if (width > height)
    max_dim = width > depth ? width : depth ;
  else
    max_dim = height > depth ? height : depth ;

#define SUPPORT_TEXAS 0
#if SUPPORT_TEXAS
  if (hdr.dime.pixdim[0] > 3.5f)/* hack to support FIL and texas */
  {
    width = hdr.dime.dim[1]; depth = hdr.dime.dim[2]; height = hdr.dime.dim[3];
  }
#endif

#if 1
  /* these are the sizes the image 'should' be */
#if SUPPORT_TEXAS
  if (hdr.dime.pixdim[0] > 3.5f)/* hack to support FIL and texas */
  {
    /*    hdr.dime.pixdim[2] *= 1.3 ;*/
#if 1
    width = max_dim / hdr.dime.pixdim[1] ;
    height = max_dim / hdr.dime.pixdim[2] ;
    depth = max_dim / hdr.dime.pixdim[3] ;
#else
    width = max_dim / hdr.dime.pixdim[3] ;
    height = max_dim / hdr.dime.pixdim[2] ;
    depth = max_dim / hdr.dime.pixdim[1] ;
#endif
  }
  else
#endif
  {
    width = max_dim / hdr.dime.pixdim[3] ;
    height = max_dim / hdr.dime.pixdim[2] ;
    depth = max_dim / hdr.dime.pixdim[1] ;
  }
#endif

  mri = 
    MRIallocSequence(width, height, depth, type, hdr.dime.dim[4]);
#if SUPPORT_TEXAS
  if (hdr.dime.pixdim[0] > 3.5f)/* hack to support FIL and texas */
  {
    mri->xsize = hdr.dime.pixdim[1] ;
    mri->ysize = hdr.dime.pixdim[2] ;
    mri->zsize = hdr.dime.pixdim[3] ;
  }
  else
#endif
  {
    mri->xsize = hdr.dime.pixdim[3] ;
    mri->ysize = hdr.dime.pixdim[2] ;
    mri->zsize = hdr.dime.pixdim[1] ;
  }

  if (FEQUAL(mri->xsize, 1))
  {
    mri->zstart = mri->ystart = mri->xstart ;
    mri->zend = mri->yend = mri->xend ;
  }
  else if (FEQUAL(mri->ysize, 1))
  {
    mri->zstart = mri->xstart = mri->ystart ;
    mri->zend = mri->xend = mri->yend ;
  }
  else if (FEQUAL(mri->zsize, 1))
  {
    mri->ystart = mri->xstart = mri->zstart ;
    mri->yend = mri->xend = mri->zend ;
  }
  
  if (type == MRI_FLOAT)
    read_float_analyze_image(fname, mri, &hdr) ;
  else
    read_byte_analyze_image(fname, mri, &hdr) ;

#if 1
  mri_dst = MRIinterpolate(mri, NULL) ;
  MRIfree(&mri) ;
#if SUPPORT_TEXAS
  if (hdr.dime.pixdim[0] > 3.5f)/* hack to support FIL and texas */
  {
    mri = MRIreorder(mri_dst, NULL, ZDIM, -XDIM, -YDIM) ;
    MRIreorder(mri, mri_dst, XDIM, YDIM, -ZDIM) ;
#if 0
    MRIfree(&mri_dst) ;
    mri_dst = mri ;
#endif
  }
#endif
  if (mri_dst->width < 256)
  {
    int  x, y, z ;

    x = (256 - mri_dst->width) / 2 ;
    y = (256 - mri_dst->height) / 2 ;
    z = (256 - mri_dst->depth) / 2 ;
    mri = MRIalloc(256, 256, 256, MRI_UCHAR) ;
    MRIextractInto(mri_dst, mri, 0, 0, 0, width, height, depth, x, y, z) ;
    MRIfree(&mri_dst) ;
    mri_dst = mri ;
  }

  return(mri_dst) ;
#else
  return(mri) ;
#endif
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
read_analyze_header(char *fname, dsr *bufptr)
{
  FILE *fp;
  char hdr_fname[STRLEN], *dot ;
  int  nread, i ;

  strcpy(hdr_fname, fname) ;
  dot = strrchr(hdr_fname, '.') ;
  if (dot)
    *dot = 0 ;
  strcat(hdr_fname, ".hdr") ;
  fp = fopen(hdr_fname,"rb");
  if (fp==NULL) 
    ErrorReturn(ERROR_NO_FILE, 
                (ERROR_NO_FILE, "File %s not found\n", hdr_fname)) ;

  nread = fread(bufptr,sizeof(char), sizeof(dsr),fp);
  fclose(fp);
  if (nread != sizeof(dsr))
    ErrorReturn(ERROR_BADFILE, 
                (ERROR_BADFILE,
                 "read_analyze_header: could only read %d of %d bytes",
                 nread, sizeof(dsr))) ;
#ifdef Linux
  bufptr->hk.sizeof_hdr = swapInt(bufptr->hk.sizeof_hdr) ;
  bufptr->hk.extents = swapInt(bufptr->hk.extents) ;
  bufptr->hk.session_error = swapInt(bufptr->hk.session_error) ;
  for (i=0;i<8;i++)
    bufptr->dime.dim[i] = swapShort(bufptr->dime.dim[i]) ;

  bufptr->dime.datatype = swapShort(bufptr->dime.datatype) ;
  bufptr->dime.bitpix = swapShort(bufptr->dime.bitpix) ;
  bufptr->dime.dim_un0 = swapShort(bufptr->dime.dim_un0) ;
  for (i=0;i<8;i++)
    bufptr->dime.pixdim[i] = swapFloat(bufptr->dime.pixdim[i]) ;

  bufptr->dime.compressed = swapFloat(bufptr->dime.compressed) ;
  bufptr->dime.verified = swapFloat(bufptr->dime.verified) ;
  bufptr->dime.glmax = swapInt(bufptr->dime.glmax) ;
  bufptr->dime.glmin = swapInt(bufptr->dime.glmin) ;
  bufptr->hist.views = swapInt(bufptr->hist.views) ;
  bufptr->hist.vols_added = swapInt(bufptr->hist.vols_added) ;
  bufptr->hist.start_field = swapInt(bufptr->hist.start_field) ;
  bufptr->hist.field_skip = swapInt(bufptr->hist.field_skip) ;
  bufptr->hist.omax = swapInt(bufptr->hist.omax) ;
  bufptr->hist.omin = swapInt(bufptr->hist.omin) ;
  bufptr->hist.smax = swapInt(bufptr->hist.smax) ;
  bufptr->hist.smin = swapInt(bufptr->hist.smin) ;
#endif
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    printf("reading analyze header file %s\n",hdr_fname);
    printf("\nheader_key\n");
    printf("sizeof_hdr    = %d\n",bufptr->hk.sizeof_hdr);
    printf("data_type     = %s\n",bufptr->hk.data_type);
    printf("db_name       = %s\n",bufptr->hk.db_name);
    printf("extents       = %d\n",bufptr->hk.extents);
    printf("session_error = %d\n",bufptr->hk.session_error);
    printf("regular       = %c\n",bufptr->hk.regular);
    printf("hkey_un0      = %c\n",bufptr->hk.hkey_un0);
    printf("\nimage_dimension\n");
    printf("dim           = ");
    for (i=0;i<8;i++)
      printf("%d ",bufptr->dime.dim[i]);
    printf("\n");
    printf("datatype      = %d\n",bufptr->dime.datatype);
    printf("bitpix        = %d\n",bufptr->dime.bitpix);
    printf("dim_un0       = %d\n",bufptr->dime.dim_un0);
    printf("pixdim        = ");
    for (i=0;i<8;i++)
      printf("%f ",bufptr->dime.pixdim[i]);
    printf("\n");
    printf("compressed    = %f\n",bufptr->dime.compressed);
    printf("verified      = %f\n",bufptr->dime.verified);
    printf("glmax         = %d\n",bufptr->dime.glmax);
    printf("glmin         = %d\n",bufptr->dime.glmin);
    printf("\ndata_history\n");
    printf("descrip       = %s\n",bufptr->hist.descrip);
    printf("aux_file      = %s\n",bufptr->hist.aux_file);
    printf("orient        = %c\n",bufptr->hist.orient);
    printf("originator    = %s\n",bufptr->hist.originator);
    printf("generated     = %s\n",bufptr->hist.generated);
    printf("scannum       = %s\n",bufptr->hist.scannum);
    printf("patient_id    = %s\n",bufptr->hist.patient_id);
    printf("exp_date      = %s\n",bufptr->hist.exp_date);
    printf("exp_time      = %s\n",bufptr->hist.exp_time);
    printf("hist_un0      = %s\n",bufptr->hist.hist_un0);
    printf("views         = %d\n",bufptr->hist.views);
    printf("vols_added    = %d\n",bufptr->hist.vols_added);
    printf("start_field   = %d\n",bufptr->hist.start_field);
    printf("field_skip    = %d\n",bufptr->hist.field_skip);
    printf("omax          = %d\n",bufptr->hist.omax);
    printf("omin          = %d\n",bufptr->hist.omin);
    printf("smax          = %d\n",bufptr->hist.smax);
    printf("smin          = %d\n",bufptr->hist.smin);
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
read_float_analyze_image(char *fname, MRI *mri, dsr *hdr)
{
  int    x, y, z, bufsize, bytes_per_pix, width, height, depth ;
  float  f, d;
  FILE   *fp;
  int    nread, datatype, i ;
  char   *buf ;

  bytes_per_pix = hdr->dime.bitpix/8 ;
  width = hdr->dime.dim[1] ;
  height = hdr->dime.dim[2] ;
  depth = hdr->dime.dim[3] ;
  datatype = hdr->dime.datatype ;
  bufsize = width * height ;
  
  fp = fopen(fname,"rb");
  if (fp==NULL) 
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "File %s not found\n",fname)) ;

  buf = (char *)calloc(bufsize, bytes_per_pix) ;
  if (!buf)
  {
    fclose(fp) ;
    ErrorReturn(ERROR_NOMEMORY, 
                (ERROR_NOMEMORY, 
                 "read_analyze_image: could not allocate %d x %d x %d buffer",
                 width, height, bytes_per_pix)) ;
  }

  for (z = 0 ; z < depth ; z++)
  {
    nread = fread(buf, bytes_per_pix, bufsize, fp) ;
    if (nread != bufsize)
    {
      free(buf) ;
      fclose(fp) ;
      ErrorReturn(ERROR_BADFILE, 
                  (ERROR_BADFILE, 
                   "read_analyze_image: could not slice %d (%d items read)",
                   z, bufsize)) ;
    }

    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        switch (datatype)
        {
        default:
          printf("data type %d not supported\n",datatype);
          exit(1);
          break;
        case DT_SIGNED_INT:
          i = *(int *)(buf+bytes_per_pix*(y*width+x)) ;
#ifdef Linux
          i = swapInt(i) ;
#endif
          f = (float)i ;
          break;
        case DT_FLOAT:
          f = *(float *)(buf+bytes_per_pix*(y*width+x));
#ifdef Linux
          f = swapFloat(f) ;
#endif
          break;
        case DT_DOUBLE:
          d = *(double *)(buf+bytes_per_pix*(y*width+x));
#ifdef Linux
          d = swapDouble(d) ;
#endif
          f = (float)d ;
          break;
        }
        MRIFvox(mri, z, height-y-1, width-x-1) = f ;
      }
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
static int
read_byte_analyze_image(char *fname, MRI *mri, dsr *hdr)
{
  int    x, y, z, bufsize, bytes_per_pix, width, height, depth, xd, yd, zd ;
  float  scale ;
  FILE   *fp;
  int    nread, datatype, xoff, yoff, zoff ;
  char   *buf, b ;
  short  s, smin, smax ;

  bytes_per_pix = hdr->dime.bitpix/8 ;
  width = hdr->dime.dim[1] ;
  height = hdr->dime.dim[2] ;
  depth = hdr->dime.dim[3] ;
#if SUPPORT_TEXAS
  if (hdr->dime.pixdim[0] > 3.5f)/* hack to support FIL and texas */
  {
    xoff = (mri->width - width) / 2 ;
    yoff = (mri->height - depth) / 2 ;
    zoff = (mri->depth - height) / 2 ;
  }
  else
#endif
  {
    xoff = (mri->width - depth) / 2 ;
    yoff = (mri->height - height) / 2 ;
    zoff = (mri->depth - width) / 2 ;
  }
  datatype = hdr->dime.datatype ;
  bufsize = width * height ;

  smin = 10000 ; smax = 0 ;
  scale = 255.0f / hdr->dime.glmax ;
  fp = fopen(fname,"rb");
  if (fp==NULL) 
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "File %s not found\n",fname)) ;

  buf = (char *)calloc(bufsize, bytes_per_pix) ;
  if (!buf)
  {
    fclose(fp) ;
    ErrorReturn(ERROR_NOMEMORY, 
                (ERROR_NOMEMORY, 
                "read_analyze_image: could not allocate %d x %d x %d buffer",
                 width, height, bytes_per_pix)) ;
  }

  for (z = 0 ; z < depth ; z++)
  {
    nread = fread(buf, bytes_per_pix, bufsize, fp) ;
    if (nread != bufsize)
    {
      free(buf) ;
      fclose(fp) ;
      ErrorReturn(ERROR_BADFILE, 
                  (ERROR_BADFILE, 
                   "read_analyze_image: could not slice %d (%d items read)",
                   z, bufsize)) ;
    }

    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        switch (datatype)
        {
        default:
          printf("data type %d not supported\n",datatype);
          exit(1);
          break;
        case DT_UNSIGNED_CHAR:
          b = (char)(*(unsigned char *)(buf+bytes_per_pix*(y*width+x)));
          xd = z + xoff ; yd = (height-y-1) + yoff ; zd = (width-x-1) + zoff;
          MRIvox(mri, xd, yd, zd) = b ;
          break;
        case DT_SIGNED_SHORT:
          s = *(short *)(buf+bytes_per_pix*(y*width+x)) ;
#ifdef Linux
          s = swapShort(s) ;
#endif
          if (s > smax)
            smax = s ;
          if (s < smin)
            smin = s ;
          b = (char)(scale * (float)(s-hdr->dime.glmin));
          break;
        }
#if 0
        xd = z + xoff ; yd = (height-y-1) + yoff ; zd = (width-x-1) + zoff ;
        MRIvox(mri, xd, yd, zd) = b ;
#endif
      }
    }
  }
  
  if (datatype == DT_SIGNED_SHORT)
  {
    if (!smax)
      ErrorReturn(ERROR_BADFILE, 
                  (ERROR_BADFILE, "read_analyze_image: 0 image\n")) ;
    scale = 255.0f / smax ;
    if  (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr, "min = %d, max = %d\n", smin, smax) ;
    rewind(fp) ;

    for (z = 0 ; z < depth ; z++)
    {
      nread = fread(buf, bytes_per_pix, bufsize, fp) ;
      if (nread != bufsize)
      {
        free(buf) ;
        fclose(fp) ;
        ErrorReturn(ERROR_BADFILE, 
                    (ERROR_BADFILE, 
                   "read_analyze_image: could not slice %d (%d items read)",
                   z, bufsize)) ;
      }
      
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          s = *(short *)(buf+bytes_per_pix*(y*width+x)) ;
#ifdef Linux
          s = swapShort(s) ;
#endif
          b = (char)(scale * (float)(s-smin));
#if SUPPORT_TEXAS
          if (hdr->dime.pixdim[0] > 3.5f)/* hack to support FIL and texas */
          {
            xd = x + xoff ; yd = z + yoff ; zd = y + zoff;
          }
          else
#endif
          {
            xd = z + xoff ; yd = (height-y-1) + yoff ; zd = (width-x-1) + zoff;
          }
          MRIvox(mri, xd, yd, zd) = b ;
        }
      }
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
static int
write_analyze_header(char *fname, MRI *mri)
{
  FILE  *fp;
  char  hdr_fname[STRLEN], *dot ;
  int   nwritten, i ;
  dsr   hdr ;
  float fmin, fmax ;

  memset(&hdr, 0, sizeof(hdr)) ;
  hdr.hk.sizeof_hdr = sizeof(hdr) ;
  hdr.hk.extents = 0 ;
  hdr.hk.session_error = 0 ;
  for (i=0;i<8;i++)
    hdr.dime.dim[i] = 0 ;
  hdr.dime.dim[3] = mri->width ; 
  hdr.dime.dim[2] = mri->height ;
  hdr.dime.dim[1] = mri->depth ;
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
  hdr.dime.pixdim[2] = mri->ysize ;
  hdr.dime.pixdim[1] = mri->zsize ;

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
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
write_analyze_image(char *fname, MRI *mri)
{
  int    x, y, z, bufsize, bytes_per_pix, width, height, depth, xd, yd, zd ;
  FILE   *fp;
  int    nwritten ;
  char   *buf ;
  float  f ;

  bytes_per_pix = mri->type == MRI_FLOAT ? sizeof(float) : sizeof(char) ;
  width = mri->depth ;
  height = mri->height ;
  depth = mri->width ;
  bufsize = width * height ;

  fp = fopen(fname,"wb");
  if (fp==NULL) 
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "File %s not found\n",fname)) ;

  buf = (char *)calloc(bufsize, bytes_per_pix) ;
  if (!buf)
  {
    fclose(fp) ;
    ErrorReturn(ERROR_NOMEMORY, 
                (ERROR_NOMEMORY, 
                "read_analyze_image: could not allocate %d x %d x %d buffer",
                 width, height, bytes_per_pix)) ;
  }

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        xd = z  ; yd = (height-y-1) ; zd = (width-x-1) ;
        if (mri->type == MRI_FLOAT)
          f = MRIFvox(mri, xd, yd, zd) ;
        else
          f = MRIvox(mri, xd, yd, zd) ;
#ifdef Linux
        f = swapFloat(f) ;
#endif
        switch (mri->type)
        {
        default:
          printf("data type %d not supported\n",mri->type);
          exit(1);
          break;
        case MRI_UCHAR:
          *(unsigned char *)(buf+bytes_per_pix*(y*width+x)) = (char)nint(f) ;
          break;
        case MRI_FLOAT:
          *(float *)(buf+bytes_per_pix*(y*width+x)) = f ;
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
  
  width = mri_src->width ; height = mri_src->height ; depth = mri_src->depth;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  xd = yd = zd = 0 ;
  switch (abs(xdim))
  {
  default:
  case XDIM: 
    if (mri_dst->width != width)
      ErrorReturn(NULL,(ERROR_BADPARM, "MRIreorder: incorrect dst width"));
    break ;
  case YDIM: 
    if (mri_dst->height != width)
      ErrorReturn(NULL,(ERROR_BADPARM, "MRIreorder: incorrect dst width"));
    break ;
  case ZDIM: 
    if (mri_dst->depth != width)
      ErrorReturn(NULL,(ERROR_BADPARM, "MRIreorder: incorrect dst width"));
    break ;
  }
  switch (abs(ydim))
  {
  default:
  case XDIM: 
    if (mri_dst->width != height)
      ErrorReturn(NULL,(ERROR_BADPARM, "MRIreorder: incorrect dst height"));
    break ;
  case YDIM: 
    if (mri_dst->height != height)
      ErrorReturn(NULL,(ERROR_BADPARM, "MRIreorder: incorrect dst height"));
    break ;
  case ZDIM: 
    if (mri_dst->depth != height)
      ErrorReturn(NULL,(ERROR_BADPARM, "MRIreorder: incorrect dst height"));
    break ;
  }
  switch (abs(zdim))
  {
  default:
  case XDIM: 
    if (mri_dst->width != depth)
      ErrorReturn(NULL,(ERROR_BADPARM, "MRIreorder: incorrect dst depth"));
    break ;
  case YDIM: 
    if (mri_dst->height != depth)
      ErrorReturn(NULL,(ERROR_BADPARM, "MRIreorder: incorrect dst depth"));
    break ;
  case ZDIM: 
    if (mri_dst->depth != depth)
      ErrorReturn(NULL,(ERROR_BADPARM, "MRIreorder: incorrect dst depth"));
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
        MRIvox(mri_dst, xd, yd, zd) = MRIvox(mri_src, xs, ys, zs) ;
      }
    }
  }
  return(mri_dst) ;
}

#define UNUSED_SPACE_SIZE 256
#define MGH_VERSION       1

static MRI *
mghRead(char *fname, int read_volume, int frame)
{
  MRI  *mri ;
  FILE  *fp ;
  int   start_frame, end_frame, width, height, depth, nframes, type, x, y, z,
        bpv, dof, bytes, version ;
  BUFTYPE *buf ;
  char   unused_buf[UNUSED_SPACE_SIZE+1] ;

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
    mri = MRIallocHeader(width, height, depth, type) ;
  else
  {
    if (type == MRI_UCHAR)
      buf = (BUFTYPE *)calloc(bytes, sizeof(BUFTYPE)) ;
    else
      buf = NULL ;
    mri = MRIallocSequence(width, height, depth, type, nframes) ;
    for (frame = start_frame ; frame <= end_frame ; frame++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        switch (type)
        {
          case MRI_FLOAT:
          for (y = 0 ; y < height ; y++)
          {
            for (x = 0 ; x < width ; x++)
            {
              MRIFseq_vox(mri,x,y,z,frame-start_frame) = freadFloat(fp) ;
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
  int   start_frame, end_frame, x, y, z, width, height, depth ;
  char        buf[UNUSED_SPACE_SIZE+1] ;

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
        case MRI_FLOAT:
          for (x = 0 ; x < width ; x++)
          {
            fwriteFloat(MRIFseq_vox(mri,x,y,z,frame), fp) ;
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
          }
          break ;
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
  }

  fclose(fp) ;
  return(NO_ERROR) ;
}
  c = header->hist.originator[8]; header->hist.originator[8] = header->hist.originator[9]; header->hist.originator[9] = c;

}  /*  end flipAnalyzeHeader()  */

/* EOF */
