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

static void buffer_to_image(BUFTYPE *buf, MRI *mri, int slice) ;
static void image_to_buffer(BUFTYPE *buf, MRI *mri, int slice) ;
static MRI *mncRead(char *fname, int read_volume, int frame) ;
static MRI *analyzeRead(char *fname, int read_volume, int frame) ;
static int mncWrite(MRI *mri, char *fname, int frame) ;
static int read_analyze_header(char *fpref, dsr *bufptr) ;
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
    buffer_to_image(buf, mri, slice) ;
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
  case MRI_MINC_FILE:
    mri = mncRead(fname, 1, frame) ;
    if (!mri)
      return(NULL) ;
    break ;
  case MRI_ANALYZE_FILE:
    mri = analyzeRead(fname, 1, frame) ;
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
      buffer_to_image(buf, mri, slice) ;
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
buffer_to_image(BUFTYPE *buf, MRI *mri, int slice)
{
  int           y, width, height ;
  BUFTYPE       *pslice ;
  
  width = mri->width ;
  height = mri->height ;
  for (y = 0 ; y < height ; y++)
  {
    pslice = mri->slices[slice][y] ;
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

  /* these are the sizes the image 'should' be */
  width = max_dim / hdr.dime.pixdim[3] ;
  height = max_dim / hdr.dime.pixdim[2] ;
  depth = max_dim / hdr.dime.pixdim[1] ;

  mri = 
    MRIallocSequence(width, height, depth, type, hdr.dime.dim[4]);
  mri->xsize = hdr.dime.pixdim[3] ;
  mri->ysize = hdr.dime.pixdim[2] ;
  mri->zsize = hdr.dime.pixdim[1] ;

  if (type == MRI_FLOAT)
    read_float_analyze_image(fname, mri, &hdr) ;
  else
    read_byte_analyze_image(fname, mri, &hdr) ;

  mri_dst = MRIinterpolate(mri, NULL) ;
  MRIfree(&mri) ;
  return(mri_dst) ;
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

static int
read_byte_analyze_image(char *fname, MRI *mri, dsr *hdr)
{
  int    x, y, z, bufsize, bytes_per_pix, width, height, depth, xd, yd, zd ;
  float  scale ;
  FILE   *fp;
  int    nread, datatype, xoff, yoff, zoff ;
  char   *buf, b ;
  short  s ;

  bytes_per_pix = hdr->dime.bitpix/8 ;
  width = hdr->dime.dim[1] ;
  height = hdr->dime.dim[2] ;
  depth = hdr->dime.dim[3] ;
  xoff = (mri->width - depth) / 2 ;
  yoff = (mri->height - height) / 2 ;
  zoff = (mri->depth - width) / 2 ;
  datatype = hdr->dime.datatype ;
  bufsize = width * height ;

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
          break;
        case DT_SIGNED_SHORT:
          s = *(short *)(buf+bytes_per_pix*(y*width+x)) ;
#ifdef Linux
          s = swapShort(s) ;
#endif
          b = (char)(scale * (float)(s-hdr->dime.glmin));
          break;
        }
        xd = z + xoff ; yd = (height-y-1) + yoff ; zd = (width-x-1) + zoff ;
        MRIvox(mri, xd, yd, zd) = b ;
      }
    }
  }
  fclose(fp);
  free(buf) ;
  return(NO_ERROR) ;
}

  c = header->hist.originator[8]; header->hist.originator[8] = header->hist.originator[9]; header->hist.originator[9] = c;

}  /*  end flipAnalyzeHeader()  */

/* EOF */
