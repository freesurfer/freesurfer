/**
 * @file  gifti_local.c
 * @brief local utilities for GIFTI library
 *
 * This file has some some extra functions for use with the GIFTI
 * utilities. The official utilities reside in gifti_io.c and gifti_xml.c
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2010/03/02 19:38:03 $
 *    $Revision: 1.14 $
 *
 * Copyright (C) 2007-2010,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 *
 */

#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>
#include <pwd.h>
#include <time.h>

#include "error.h" // return codes
#include "utils.h" // strcpyalloc
#include "nifti1.h"
#include "gifti_local.h"

/*
 *
 */
static giiDataArray* gifti_alloc_and_add_darray (gifti_image* image)
{
  if (!image)
  {
    fprintf (stderr,"** gifti_alloc_and_add_darray: NULL image\n");
    return NULL;
  }

  /* Try to add an empty array. */
  if (gifti_add_empty_darray(image,1))
  {
    fprintf (stderr,"** gifti_alloc_and_add_darray: gifti_add_empty_darray "
             "failed\n");
    return NULL;
  }

  /* Return the array we just allocated. */
  return image->darray[image->numDA-1];
}


/*
 *
 */
static double gifti_get_DA_value_2D (giiDataArray* da, int row, int col)
{
  int dim0_index, dim1_index;

  if (!da || !da->data)
  {
    fprintf (stderr,"** gifti_get_DA_value_2D, invalid params: data=%p\n",
             da);
    return 0;
  }

  if (da->num_dim != 2)
  {
    fprintf (stderr,"** gifti_get_DA_value_2D, array dim is %d\n",
             da->num_dim);
    return 0;
  }

  /* Get the dim0 and dims[1] indices based on our order. */
  if (GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord)
  {
    dim0_index = row;
    dim1_index = col;
  }
  else if (GIFTI_IND_ORD_COL_MAJOR == da->ind_ord)
  {
    dim0_index = col;
    dim1_index = row;
  }
  else
  {
    fprintf (stderr,"** gifti_get_DA_value_2D, unknown ind_ord: %d\n",
             da->ind_ord);
    return 0;
  }

  /* Check the indices. */
  if (dim0_index < 0 || dim0_index >= da->dims[0] ||
      dim1_index < 0 || dim1_index >= da->dims[1])
  {
    fprintf(stderr,"** gifti_get_DA_value_2D, invalid params: "
            "dim0_index=%d (max=%d), dim1_index=%d (max=%d)\n",
            dim0_index, da->dims[0], dim1_index, da->dims[1]);
    return 0;
  }

  /* Switch on the data type and return the appropriate
     element. Indexing depends on the data order. */
  switch (da->datatype)
  {
  default :
    fprintf(stderr,"** gifti_get_DA_value_2D, unsupported type %d-"
            "unknown, or can't convert to double\n",da->datatype);
    return 0;
  case NIFTI_TYPE_UINT8:
  {
    if ( GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord )
      return (double)*((unsigned char*)
                       (da->data) + (dim0_index*da->dims[1]) + dim1_index);
    else
      return (double)*((unsigned char*)
                       (da->data) + dim0_index + (dim1_index*da->dims[0]));
    break;
  }
  case NIFTI_TYPE_INT16:
  {
    if ( GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord )
      return (double)*((short*)
                       (da->data) + (dim0_index*da->dims[1]) + dim1_index);
    else
      return (double)*((short*)
                       (da->data) + dim0_index + (dim1_index*da->dims[0]));
    break;
  }
  case NIFTI_TYPE_INT32:
  {
    if ( GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord )
      return (double)*((int*)
                       (da->data) + (dim0_index*da->dims[1]) + dim1_index);
    else
      return (double)*((int*)
                       (da->data) + dim0_index + (dim1_index*da->dims[0]));
    break;
  }
  case NIFTI_TYPE_FLOAT32:
  {
    if ( GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord )
      return (double)*((float*)
                       (da->data) + (dim0_index*da->dims[1]) + dim1_index);
    else
      return (double)*((float*)
                       (da->data) + dim0_index + (dim1_index*da->dims[0]));
    break;
  }
  case NIFTI_TYPE_INT8:
  {
    if ( GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord )
      return (double)*((char*)
                       (da->data) + (dim0_index*da->dims[1]) + dim1_index);
    else
      return (double)*((char*)
                       (da->data) + dim0_index + (dim1_index*da->dims[0]));
    break;
  }
  case NIFTI_TYPE_UINT16:
  {
    if ( GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord )
      return (double)*((unsigned short*)
                       (da->data) + (dim0_index*da->dims[1]) + dim1_index);
    else
      return (double)*((unsigned short*)
                       (da->data) + dim0_index + (dim1_index*da->dims[0]));
    break;
  }
  case NIFTI_TYPE_UINT32:
  {
    if ( GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord )
      return (double)*((unsigned int*)
                       (da->data) + (dim0_index*da->dims[1]) + dim1_index);
    else
      return (double)*((unsigned int*)
                       (da->data) + dim0_index + (dim1_index*da->dims[0]));
    break;
  }
  }

  return 0;
}


/*
 *
 */
static void gifti_set_DA_value_2D (giiDataArray* da,
                                   int row, int col, double value)
{
  int dim0_index, dim1_index;

  if (!da || !da->data)
  {
    fprintf (stderr,"** gifti_set_DA_value_2D, invalid params: data=%p\n",
             da);
    return;
  }

  if (da->num_dim != 2)
  {
    fprintf (stderr,"** gifti_set_DA_value_2D, array dim is %d\n",
             da->num_dim);
    return;
  }

  /* Get the dim0 and dims[1] indices based on our order. */
  if (GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord)
  {
    dim0_index = row;
    dim1_index = col;
  }
  else
  {
    dim0_index = col;
    dim1_index = row;
  }

  /* Check the indices. */
  if (dim0_index < 0 || dim0_index >= da->dims[0] ||
      dim1_index < 0 || dim1_index >= da->dims[1])
  {
    fprintf(stderr,"** gifti_set_DA_value_2D, invalid params: "
            "dim0_index=%d (max=%d), dim1_index=%d (max=%d)\n",
            dim0_index, da->dims[0], dim1_index, da->dims[1]);
    return;
  }

  /* Switch on the data type and write the appropriate
     element. Indexing depends on the data order. */
  switch (da->datatype)
  {
  default :
    fprintf(stderr,"** gifti_set_DA_value_2D, unsupported type %d-"
            "unknown, or can't convert to double\n",da->datatype);
    return;
  case NIFTI_TYPE_UINT8:
  {
    if ( GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord )
      *((unsigned char*)(da->data) + (dim0_index*da->dims[1]) + dim1_index) =
        (unsigned char)value;
    else
      *((unsigned char*)(da->data) + dim0_index + (dim1_index*da->dims[0])) =
        (unsigned char)value;
    break;
  }
  case NIFTI_TYPE_INT16:
  {
    if ( GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord )
      *((short*)(da->data) + (dim0_index*da->dims[1]) + dim1_index) =
        (short)value;
    else
      *((short*)(da->data) + dim0_index + (dim1_index*da->dims[0])) =
        (short)value;
    break;
  }
  case NIFTI_TYPE_INT32:
  {
    if ( GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord )
      *((int*)(da->data) + (dim0_index*da->dims[1]) + dim1_index) =
        (int)value;
    else
      *((int*)(da->data) + dim0_index + (dim1_index*da->dims[0])) =
        (int)value;
    break;
  }
  case NIFTI_TYPE_FLOAT32:
  {
    if ( GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord )
      *((float*)(da->data) + (dim0_index*da->dims[1]) + dim1_index) =
        (float)value;
    else
      *((float*)(da->data) + dim0_index + (dim1_index*da->dims[0])) =
        (float)value;
    break;
  }
  case NIFTI_TYPE_INT8:
  {
    if ( GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord )
      *((char*)(da->data) + (dim0_index*da->dims[1]) + dim1_index) =
        (char)value;
    else
      *((char*)(da->data) + dim0_index + (dim1_index*da->dims[0])) =
        (char)value;
    break;
  }
  case NIFTI_TYPE_UINT16:
  {
    if ( GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord )
      *((unsigned short*)(da->data) + (dim0_index*da->dims[1]) + dim1_index) =
        (unsigned short)value;
    else
      *((unsigned short*)(da->data) + dim0_index + (dim1_index*da->dims[0])) =
        (unsigned short)value;
    break;
  }
  case NIFTI_TYPE_UINT32:
  {
    if ( GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord )
      *((unsigned int*)(da->data) + (dim0_index*da->dims[1]) + dim1_index) =
        (unsigned int)value;
    else
      *((unsigned int*)(da->data) + dim0_index + (dim1_index*da->dims[0])) =
        (unsigned int)value;
    break;
  }
  }

  return;
}



/*-----------------------------------------------------
  Parameters:    input file name of GIFTI file

  Returns value: freesurfer surface structure

  Description:   reads a GIFTI file, putting vertices, 
                 and faces into an MRIS_SURFACE structure.
  ------------------------------------------------------*/
MRI_SURFACE * mrisReadGIFTIfile(const char *fname)
{
  /* Attempt to read the file. */
  gifti_image* image = gifti_read_image (fname, 1);
  if (NULL == image)
  {
    fprintf (stderr,"mrisReadGIFTIfile: gifti_read_image() returned NULL\n");
    return NULL;
  }

  /* check for compliance */
  int valid = gifti_valid_gifti_image (image, 1);
  if (valid == 0)
  {
    fprintf (stderr,"mrisReadGIFTIfile: GIFTI file %s is invalid!\n", fname);
    gifti_free_image (image);
    return NULL;
  }

  /* Make sure we got at least coords and faces. */
  giiDataArray* coords = gifti_find_DA (image, NIFTI_INTENT_POINTSET, 0);
  if (NULL == coords)
  {
    fprintf (stderr,"mrisReadGIFTIfile: no coordinates in file %s\n", fname);
    gifti_free_image (image);
    return NULL;
  }
  giiDataArray* faces = gifti_find_DA (image,  NIFTI_INTENT_TRIANGLE, 0);
  if (NULL == faces)
  {
    fprintf (stderr,"mrisReadGIFTIfile: no topology in file %s\n", fname);
    gifti_free_image (image);
    return NULL;
  }

  /* Check the number of vertices and faces. */
  long long num_vertices = 0;
  long long num_cols = 0;
  gifti_DA_rows_cols (coords, &num_vertices, &num_cols);
  if (num_vertices <= 0 || num_cols != 3)
  {
    fprintf (stderr,"mrisReadGIFTIfile: malformed coords data array in file "
             "%s: num_vertices=%d num_cols=%d\n",
             fname, (int)num_vertices, (int)num_cols);
    gifti_free_image (image);
    return NULL;
  }
  long long num_faces = 0;
  num_cols = 0;
  gifti_DA_rows_cols (faces, &num_faces, &num_cols);
  if (num_faces <= 0 || num_cols != 3)
  {
    fprintf (stderr,"mrisReadGIFTIfile: malformed faces data array in file "
             "%s: num_faces=%d num_cols=%d\n",
             fname, (int)num_faces, (int)num_cols);
    gifti_free_image (image);
    return NULL;
  }

  /* Try to allocate a surface. */
  MRIS* mris = MRISalloc (num_vertices, num_faces);
  if (NULL == mris)
  {
    fprintf (stderr,"mrisReadGIFTIfile: failed to allocate an MRIS with "
             "%d vertices and %d faces\n",
             (int)num_vertices,(int) num_faces);
    gifti_free_image (image);
    return NULL;
  }

  /* Set some meta data in the mris. */
  strcpy(mris->fname,fname);
  mris->type = MRIS_TRIANGULAR_SURFACE;
  char* hemi = gifti_get_meta_value(&coords->meta,
                                    "AnatomicalStructurePrimary'");
  if (hemi && (strcmp(hemi,"CortexRight")==0))
    mris->hemisphere = RIGHT_HEMISPHERE ;
  else if (hemi && (strcmp(hemi,"CortexLeft")==0))
    mris->hemisphere = LEFT_HEMISPHERE ;

  /* Copy in the vertices. */
  float x, y, z, xhi, xlo, yhi, ylo, zhi, zlo ;
  xhi=yhi=zhi= -10000;
  xlo=ylo=zlo= 10000;
  int vertex_index;
  for (vertex_index = 0; vertex_index < num_vertices; vertex_index++)
  {
    mris->vertices[vertex_index].x =
      (float) gifti_get_DA_value_2D (coords, vertex_index, 0);
    mris->vertices[vertex_index].y =
      (float) gifti_get_DA_value_2D (coords, vertex_index, 1);
    mris->vertices[vertex_index].z =
      (float) gifti_get_DA_value_2D (coords, vertex_index, 2);
    mris->vertices[vertex_index].num = 0;
    mris->vertices[vertex_index].origarea = -1;
    x = mris->vertices[vertex_index].x;
    y = mris->vertices[vertex_index].y;
    z = mris->vertices[vertex_index].z;
    if (x>xhi) xhi=x;
    if (x<xlo) xlo=x;
    if (y>yhi) yhi=y;
    if (y<ylo) ylo=y;
    if (z>zhi) zhi=z;
    if (z<zlo) zlo=z;
  }
  mris->xlo = xlo ;
  mris->ylo = ylo ;
  mris->zlo = zlo ;
  mris->xhi = xhi ;
  mris->yhi = yhi ;
  mris->zhi = zhi ;
  mris->xctr = (xhi+xlo)/2;
  mris->yctr = (yhi+ylo)/2;
  mris->zctr = (zhi+zlo)/2;

  /* Copy in the faces. */
  int face_index;
  for (face_index = 0; face_index < num_faces; face_index++)
  {
    int face_vertex_index;
    for (face_vertex_index = 0;
         face_vertex_index < VERTICES_PER_FACE;
         face_vertex_index++)
    {
      vertex_index =
        (int) gifti_get_DA_value_2D (faces, face_index, face_vertex_index);
      mris->faces[face_index].v[face_vertex_index] = vertex_index;
      mris->vertices[vertex_index].num++;
    }
  }
  // each vertex has a face list (faster that face list in some operations)
  for (vertex_index = 0; vertex_index < num_vertices; vertex_index++)
  {
    mris->vertices[vertex_index].f =
      (int *)calloc(mris->vertices[vertex_index].num,sizeof(int));
    mris->vertices[vertex_index].n =
      (uchar *)calloc(mris->vertices[vertex_index].num,sizeof(uchar));
    mris->vertices[vertex_index].num = 0; // this gets re-calc'd next...
  }
  for (face_index = 0 ; face_index < mris->nfaces ; face_index++)
  {
    FACE *face = &mris->faces[face_index] ;
    int n;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
      mris->vertices[face->v[n]].f[mris->vertices[face->v[n]].num++]
        = face_index; // note that .num is auto-incremented!
  }
  for (vertex_index = 0; vertex_index < num_vertices; vertex_index++)
  {
    int n,m;
    for (n=0;n<mris->vertices[vertex_index].num;n++)
    {
      for (m=0;m<VERTICES_PER_FACE;m++)
      {
        if (mris->faces[mris->vertices[vertex_index].f[n]].v[m] == 
            vertex_index)
          mris->vertices[vertex_index].n[n] = m;
      }
    }
  }

  /* other data structure essentials, namely:
   *  mrisFindNeighbors(mris);
   *  mrisComputeVertexDistances(mris);
   *  mrisReadTransform(mris, fname) ;
   *  mris->radius = MRISaverageRadius(mris) ;
   *  MRIScomputeMetricProperties(mris) ;
   *  MRISstoreCurrentPositions(mris) ;
   */
  MRIScomputeNormals(mris);
  UpdateMRIS(mris,fname);

  /* And we're done. */
  gifti_free_image (image);

  return mris;
}


/*
 *
 */
int mrisReadScalarGIFTIfile(MRI_SURFACE *mris,const  char *fname)
{
  /* Attempt to read the file. */
  gifti_image* image = gifti_read_image (fname, 1);
  if (NULL == image)
  {
    fprintf 
      (stderr,
       "mrisReadScalarGIFTIfile: gifti_read_image() returned NULL\n");
    return ERROR_BADFILE;
  }

  /* check for compliance */
  int valid = gifti_valid_gifti_image (image, 1);
  if (valid == 0)
  {
    fprintf 
      (stderr,
       "mrisReadScalarGIFTIfile: GIFTI file %s is invalid!\n", fname);
    gifti_free_image (image);
    return ERROR_BADFILE;
  }

  /* check for 'shape' data */
  giiDataArray* scalars = gifti_find_DA (image, NIFTI_INTENT_SHAPE, 0);
  if (NULL == scalars)
  {
    fprintf 
      (stderr,
       "mrisReadScalarGIFTIfile: no shape data found in file %s\n", 
       fname);
    gifti_free_image (image);
    return ERROR_BADFILE;
  }

  /* Check the number of vertices */
  long long num_vertices = 0;
  long long num_cols = 0;
  gifti_DA_rows_cols (scalars, &num_vertices, &num_cols);
  if (num_vertices <= 0 ||
      num_vertices != mris->nvertices ||
      num_cols != 1)
  {
    fprintf 
      (stderr,
       "mrisReadScalarGIFTIfile: malformed coords data array in file "
       "%s: num_vertices=%d num_cols=%d expected nvertices=%d\n",
       fname, (int)num_vertices, (int)num_cols, mris->nvertices);
    gifti_free_image (image);
    return ERROR_BADFILE;
  }

  /* Copy in all our data. */
  int scalar_index;
  for (scalar_index = 0; scalar_index < mris->nvertices; scalar_index++)
  {
    if (mris->vertices[scalar_index].ripflag) continue;
    mris->vertices[scalar_index].curv =
      (float) gifti_get_DA_value_2D (scalars, scalar_index, 0);
  }

  /* And we're done. */
  gifti_free_image (image);
  return(NO_ERROR) ;
}



/*-----------------------------------------------------------
  MRISreadGiftiAsMRI() - reads GIFTI functional frames into
  an MRI volume struct, which is a retro-fit usage to store
  multiple frames of data (where in this case, a frame is one
  complete vector of vertices).
  -----------------------------------------------------------*/
MRI *MRISreadGiftiAsMRI(const char *fname, int read_volume)
{
  /* Attempt to read the file. */
  gifti_image* image = gifti_read_image (fname, 1);
  if (NULL == image)
  {
    fprintf 
      (stderr,
       "MRISreadGiftiAsMRI: gifti_read_image() returned NULL\n");
    return NULL;
  }

  /* check for compliance */
  int valid = gifti_valid_gifti_image (image, 1);
  if (valid == 0)
  {
    fprintf 
      (stderr,
       "MRISreadGiftiAsMRI: GIFTI file %s is invalid!\n", fname);
    gifti_free_image (image);
    return NULL;
  }

  /* check for overlay data */
  giiDataArray* scalars = NULL;
  int frame_count = 0;
  long long num_vertices = 0;
  long long num_cols = 0;
  do
  {
    scalars = gifti_find_DA (image, NIFTI_INTENT_NONE, frame_count);
    if (NULL == scalars)
    {
      scalars = gifti_find_DA (image, NIFTI_INTENT_TIME_SERIES, frame_count);
      if (NULL == scalars)
      {
        scalars = gifti_find_DA (image, NIFTI_INTENT_SHAPE, frame_count);
        if (NULL == scalars) break;
        else
        {
          //printf("Found NIFTI_INTENT_SHAPE data array #%d\n", frame_count);
        }
      }
      else
      {
        //printf("Found NIFTI_INTENT_TIME_SERIES data array #%d\n", 
        //     frame_count);
      }
    }
    else
    {
      //printf("Found NIFTI_INTENT_NONE data array #%d\n", frame_count);
    }
    long long nvertices, ncols;
    gifti_DA_rows_cols (scalars, &nvertices, &ncols);
    if (frame_count == 0)
    {
      num_vertices=nvertices;
      num_cols=ncols;
    }
    else
    {
      if (num_vertices <= 0 ||
          num_vertices != nvertices ||
          ncols != 1)
      {
        fprintf 
          (stderr,
           "MRISreadGiftiAsMRI: malformed data array in file "
           "%s: nvertices=%d ncols=%d expected num_vertices=%d\n",
           fname, (int)nvertices, (int)num_cols, (int)num_vertices);
        gifti_free_image (image);
        return NULL;
      }
    }
    frame_count++;
  } while ( scalars );

  if (frame_count == 0)
  {
    fprintf 
      (stderr,
       "MRISreadGiftiAsMRI: no overlay data found in file %s\n", 
       fname);
    gifti_free_image (image);
    return NULL;
  }

  /* if we don't need to read the volume, just return a header */
  MRI *mri;
  if (!read_volume)
  {
    mri = MRIallocHeader(num_vertices,1,1,MRI_FLOAT);
    mri->nframes = frame_count;
    return(mri);
  }

  /* Copy in each scalar frame to 'volume' frame. */
  mri =  MRIallocSequence(num_vertices,1,1,MRI_FLOAT,frame_count);
  int frame_counter;
  for (frame_counter = 0; frame_counter < frame_count; frame_counter++)
  {
    scalars = gifti_find_DA (image, NIFTI_INTENT_NONE, frame_counter);
    if (NULL == scalars)
    {
      scalars = gifti_find_DA (image, NIFTI_INTENT_TIME_SERIES, frame_counter);
      if (NULL == scalars)
      {
        scalars = gifti_find_DA (image, NIFTI_INTENT_SHAPE, frame_counter);
      }
      if (NULL == scalars) break;
    }
    int scalar_index;
    for (scalar_index = 0; scalar_index < num_vertices; scalar_index++)
    {
      float val = (float) gifti_get_DA_value_2D (scalars, scalar_index, 0);
      MRIsetVoxVal(mri,scalar_index,0,0,frame_counter,val);
    }
  }

  /* And we're done. */
  gifti_free_image (image);
  return(mri) ;
}


/*
 * Extracts label table data from a gifti file:
 * puts the LabelTable in the freesurfer colortable struct
 * and puts the DataArray containing an 'index' for each vertex into
 * the vertex's 'annotation' element (which should be contained in
 * the colortable)
 */
int mrisReadLabelTableGIFTIfile(MRI_SURFACE *mris, const char *fname)
{
  MRISclearAnnotations(mris);

  /* Attempt to read the file. */
  gifti_image* image = gifti_read_image (fname, 1);
  if (NULL == image)
  {
    fprintf 
      (stderr,
       "mrisReadLabelTableGIFTIfile: gifti_read_image() returned NULL\n");
    return ERROR_BADFILE;
  }

  /* check for compliance */
  int valid = gifti_valid_gifti_image (image, 1);
  if (valid == 0)
  {
    fprintf 
      (stderr,
       "mrisReadLabelTableGIFTIfile: GIFTI file %s is invalid!\n", fname);
    gifti_free_image (image);
    return ERROR_BADFILE;
  }

  /* check for 'label table' data */
  if (image->labeltable.length == 0)
  {
    fprintf 
      (stderr,
       "mrisReadLabelTableGIFTIfile: no labeltable data found in file %s\n", 
       fname);
    gifti_free_image (image);
    return ERROR_BADFILE;
  }

  /* check validity of labeltable data */
  if (!gifti_valid_LabelTable(&image->labeltable,1))
  {
    fprintf 
      (stderr,
       "mrisReadLabelTableGIFTIfile: invalid labeltable found in file %s\n", 
       fname);
    gifti_free_image (image);
    return ERROR_BADFILE;
  }

  /* copy label table contents to our color_table struct */
  COLOR_TABLE* ct = (COLOR_TABLE *)calloc(1, sizeof(COLOR_TABLE));
  if (ct == NULL)
  {
    fprintf 
      (stderr,
       "mrisReadLabelTableGIFTIfile: could not alloc colortable memory\n");
    gifti_free_image (image);
    return ERROR_BADFILE;
  }
  memset(ct,0,sizeof(COLOR_TABLE));
  ct->nentries = image->labeltable.length + 1;
  ct->version = 2;
  ct->entries = (COLOR_TABLE_ENTRY**)
                calloc(ct->nentries, sizeof(COLOR_TABLE_ENTRY*));
  if (ct->entries == NULL)
  {
    fprintf 
      (stderr,
       "mrisReadLabelTableGIFTIfile: could not alloc colortable entries\n");
    gifti_free_image (image);
    return ERROR_BADFILE;
  }
  memset(ct->entries,0,sizeof(ct->entries));
  strncpy(ct->fname, fname, sizeof(ct->fname));
  int label_index;
  float* rgba = image->labeltable.rgba;
  for (label_index = 0; label_index < image->labeltable.length; label_index++)
  {
    ct->entries[label_index] = (CTE*) malloc(sizeof(CTE));
    if (ct->entries[label_index] == NULL)
    {
      fprintf 
        (stderr,
         "mrisReadLabelTableGIFTIfile: could not alloc colortable entry\n");
      gifti_free_image (image);
      return ERROR_BADFILE;
    }
    strncpy(ct->entries[label_index]->name,
            image->labeltable.label[label_index],
            sizeof(ct->entries[label_index]->name));
    if( rgba ) 
    {
      ct->entries[label_index]->rf = rgba[0];
      ct->entries[label_index]->ri = floor((rgba[0])*256);
      if (ct->entries[label_index]->ri > 255) ct->entries[label_index]->ri=255;
      ct->entries[label_index]->gf = rgba[1];
      ct->entries[label_index]->gi = floor((rgba[1])*256);
      if (ct->entries[label_index]->gi > 255) ct->entries[label_index]->gi=255;
      ct->entries[label_index]->bf = rgba[2];
      ct->entries[label_index]->bi = floor((rgba[2])*256);
      if (ct->entries[label_index]->bi > 255) ct->entries[label_index]->bi=255;
      ct->entries[label_index]->af = rgba[3];
      ct->entries[label_index]->ai = floor((rgba[3])*256);
      if (ct->entries[label_index]->ai > 255) ct->entries[label_index]->ai=255;
      rgba += 4;
      /*
      printf("RGBA: %d %d %d %d %f %f %f %f\n",
             ct->entries[label_index]->ri,
             ct->entries[label_index]->gi,
             ct->entries[label_index]->bi,
             ct->entries[label_index]->ai,
             ct->entries[label_index]->rf,
             ct->entries[label_index]->gf,
             ct->entries[label_index]->bf,
             ct->entries[label_index]->af);
      */
    }
    else
    {
      fprintf 
        (stderr,
         "mrisReadLabelTableGIFTIfile: missing rgba floats!\n");
      gifti_free_image (image);
      return ERROR_BADFILE;
    }

    // sanity-cross-check:
    int annotation = CTABrgb2Annotation(ct->entries[label_index]->ri,
                                        ct->entries[label_index]->gi,
                                        ct->entries[label_index]->bi);
    if (annotation != image->labeltable.index[label_index])
    {
      printf("WARNING: rgb2annot:%8.8X != labeltable.index:%8.8X\n",
             annotation,
             image->labeltable.index[label_index]);
    }
  }
  ct->entries[label_index] = NULL;
  CTABfindDuplicateNames(ct);
  // and of course store this colortable
  mris->ct = ct;

  // sanity-check
  int numEntries=0;
  CTABgetNumberOfValidEntries(mris->ct,&numEntries);
  if (numEntries != image->labeltable.length)
  {
    fprintf 
      (stderr,
       "mrisReadLabelTableGIFTIfile: "
       "ct_entries:%d != labeltable_entries:%d\n", 
       numEntries, image->labeltable.length);
    gifti_free_image (image);
    return ERROR_BADFILE;
  }

  /* check for 'label' data */
  giiDataArray* labels = gifti_find_DA (image, NIFTI_INTENT_LABEL, 0);
  if (NULL == labels)
  {
    labels = gifti_find_DA (image, NIFTI_INTENT_NONE, 0);
    if (NULL == labels)
    {
      fprintf 
        (stderr,
         "mrisReadLabelTableGIFTIfile: no label data found in file %s\n", 
         fname);
      gifti_free_image (image);
      return ERROR_BADFILE;
    }
  }

  /* Check the number of vertices */
  long long num_vertices = 0;
  long long num_cols = 0;
  gifti_DA_rows_cols (labels, &num_vertices, &num_cols);
  if (num_vertices <= 0 ||
      num_vertices != mris->nvertices ||
      num_cols != 1)
  {
    fprintf 
      (stderr,
       "mrisReadLabelTableGIFTIfile: malformed labels data array in file "
       "%s: num_vertices=%d num_cols=%d expected nvertices=%d\n",
       fname, (int)num_vertices, (int)num_cols, mris->nvertices);
    gifti_free_image (image);
    return ERROR_BADFILE;
  }

  /* Copy label node data (annotation for each vertex) */
  unsigned int* label_data = labels->data;
  for (label_index = 0; label_index < mris->nvertices; label_index++)
  {
    if (mris->vertices[label_index].ripflag) continue;
    mris->vertices[label_index].annotation = *(label_data + label_index);
    // cross-check:
    int index = -1;
    int result = CTABfindAnnotation(mris->ct,
                                    mris->vertices[label_index].annotation,
                                    &index);
    if ((result != NO_ERROR) || 
        (index < 0) || 
        (index > image->labeltable.length))
    {
      fprintf 
        (stderr,
         "mrisReadLabelTableGIFTIfile: "
         "label node data not found in colortable! "
         "vno: %d, annot: %8.8X\n",
         label_index, mris->vertices[label_index].annotation);
      gifti_free_image (image);
      return ERROR_BADFILE;
    }
  }

  /* And we're done. */
  gifti_free_image (image);
  return(NO_ERROR) ;
}


/*
 *
 */
static void insertMetaData(MRIS* mris, giiDataArray* dataArray)
{

  if (mris && mris->fname)
  {
    const char *primary=NULL, *secondary=NULL, *geotype=NULL;
    char *name = mris->fname;
    char *topotype="Closed";
    if (strstr(name, "lh.")) primary = "CortexLeft";
    if (strstr(name, "rh.")) primary = "CortexRight";
    if (strstr(name, ".orig"))     secondary = "GrayWhite";
    if (strstr(name, ".smoothwm")) secondary = "GrayWhite";
    if (strstr(name, ".white"))    secondary = "GrayWhite";
    if (strstr(name, ".pial"))     secondary = "Pial";
    if (strstr(name, ".orig"))     geotype = "Reconstruction";
    if (strstr(name, ".smoothwm")) geotype = "Anatomical";
    if (strstr(name, ".white"))    geotype = "Anatomical";
    if (strstr(name, ".pial"))     geotype = "Anatomical";
    if (strstr(name, ".inflated")) geotype = "Inflated";
    if (strstr(name, ".sphere"))   geotype = "Sphere";
    if (strstr(name, ".qsphere"))  geotype = "Sphere";
    if (strstr(name,"pial-outer")) geotype = "Hull";
    if (mris->patch)
    {
      geotype = "Flat";
      topotype = "Cut";
    }

    if (primary) gifti_add_to_meta( &dataArray->meta,
                                    "AnatomicalStructurePrimary",
                                    primary,
                                    1 );
    if (secondary) gifti_add_to_meta( &dataArray->meta,
                                      "AnatomicalStructureSecondary",
                                      secondary,
                                      1 );
    if (geotype) gifti_add_to_meta( &dataArray->meta,
                                    "GeometricType",
                                    geotype,
                                    1 );
    gifti_add_to_meta( &dataArray->meta, "TopologicalType", topotype, 1 );
    gifti_add_to_meta( &dataArray->meta, "Name", name, 1 );
  }

  if (mris && strlen(mris->subject_name))
  {
    gifti_add_to_meta( &dataArray->meta, "SubjectID", mris->subject_name, 1 );
  }

#if 0
#include <uuid/uuid.h>
  uuid_t uuid;
  char uuidstr[2048];
  uuid_generate(uuid);
  uuid_unparse(uuid, uuidstr);
  gifti_add_to_meta( &dataArray->meta, "UniqueID", uuidstr, 1 );
#endif

  struct passwd *pw = getpwuid(geteuid());
  if ((pw != NULL) && (pw->pw_name != NULL))
  {
    gifti_add_to_meta( &dataArray->meta, "UserName", pw->pw_name, 1 );
  }

  time_t tyme = time(NULL);
  struct tm *lt = localtime(&tyme);
  char *date = asctime(lt);
  char *chr = strchr(date,'\r');
  if (chr) *chr = 0; // remove carriage return
  chr = strchr(date,'\n');
  if (chr) *chr = 0; // remove linefeed
  gifti_add_to_meta( &dataArray->meta, "Date", date, 1 );
}


/*-----------------------------------------------------
  Parameters:    MRIS_SURFACE structure,
                 output file name of GIFTI file

  Returns value: 0 if passed, else error code

  Description:   writes a GIFTI file, putting vertices, 
                 and faces from input MRIS_SURFACE structure.
  ------------------------------------------------------*/
int MRISwriteGIFTI(MRIS* mris,const  char *fname)
{
  if (NULL == mris || NULL == fname)
  {
    fprintf (stderr,"MRISwriteGIFTI: invalid parameter\n");
    return ERROR_BADPARM;
  }

  gifti_image* image = (gifti_image *)calloc(1,sizeof(gifti_image));
  if (NULL == image)
  {
    fprintf (stderr,"MRISwriteGIFTI: couldn't allocate image\n");
    return ERROR_NOMEMORY;
  }
  image->version = strcpyalloc(GIFTI_XML_VERSION);

  /* -------------------------------------------------------
   * Coordinates array.
   */
  giiDataArray* coords = gifti_alloc_and_add_darray (image);
  if (NULL == coords)
  {
    fprintf (stderr,"MRISwriteGIFTI: couldn't allocate giiDataArray\n");
    gifti_free_image (image);
    return ERROR_NOMEMORY;
  }

  /* Set its attributes. */
  coords->intent = NIFTI_INTENT_POINTSET;
  coords->datatype = NIFTI_TYPE_FLOAT32;
  coords->ind_ord = GIFTI_IND_ORD_ROW_MAJOR;
  coords->num_dim = 2;
  coords->dims[0] = mris->nvertices; /* In highest first, dim0 = rows */
  coords->dims[1] = 3;               /* In highest first, dim1 = cols */
  coords->encoding = GIFTI_ENCODING_B64GZ; // data stored in gzip'd base64
#if (BYTE_ORDER == LITTLE_ENDIAN)
  coords->endian = GIFTI_ENDIAN_LITTLE;
#else
  coords->endian = GIFTI_ENDIAN_BIG;
#endif

  coords->coordsys = NULL; // empty, unless we find something here...
  MRISreadTransform(mris, fname) ; // tries to get xform based on fname
  if (mris->SRASToTalSRAS_ &&
      mris->SRASToTalSRAS_->rows==4 &&
      mris->SRASToTalSRAS_->cols==4)
  { // found a valid xform, so use it...
    gifti_add_empty_CS( coords );
    coords->coordsys[0]->dataspace = strcpyalloc("NIFTI_XFORM_UNKNOWN");
    coords->coordsys[0]->xformspace = strcpyalloc("NIFTI_XFORM_TALAIRACH");
    MATRIX *xform = mris->SRASToTalSRAS_;
    int r,c;
    for (r=1; r <= 4; r++)
      for (c=1; c <= 4; c++)
        coords->coordsys[0]->xform[r-1][c-1] = xform->rptr[r][c];
  }

  insertMetaData (mris, coords); /* standard meta data */

  coords->nvals = gifti_darray_nvals (coords);
  gifti_datatype_sizes (coords->datatype, &coords->nbyper, NULL);

  /* Allocate the data array. */
  coords->data = NULL;
  coords->data = (void*) calloc (coords->nvals, coords->nbyper);
  if (NULL == coords->data)
  {
    fprintf (stderr,"MRISwriteGIFTI: couldn't allocate coords data of "
             "length %d, element size %d\n",
             (int)coords->nvals, coords->nbyper);
    gifti_free_image (image);
    return ERROR_NOMEMORY;
  }

  /* Copy in all our data. */
  int vertex_index;
  for (vertex_index = 0; vertex_index < mris->nvertices; vertex_index++)
  {
    if (mris->vertices[vertex_index].ripflag) continue;
    gifti_set_DA_value_2D (coords, vertex_index, 0,
                           mris->vertices[vertex_index].x);
    gifti_set_DA_value_2D (coords, vertex_index, 1,
                           mris->vertices[vertex_index].y);
    gifti_set_DA_value_2D (coords, vertex_index, 2,
                           mris->vertices[vertex_index].z);
  }

  /* -------------------------------------------------------
   * Faces array. 
   */
  giiDataArray* faces = gifti_alloc_and_add_darray (image);
  if (NULL == faces)
  {
    fprintf (stderr,"MRISwriteGIFTI: couldn't allocate giiDataArray\n");
    gifti_free_image (image);
    return ERROR_NOMEMORY;
  }

  /* count the real number of faces (the ones that dont have a vertex
     with a ripflag set) */
  int numFaces = 0;
  int face_index;
  for (face_index = 0; face_index < mris->nfaces; face_index++)
  {
    if (mris->vertices[mris->faces[face_index].v[0]].ripflag) continue;
    if (mris->vertices[mris->faces[face_index].v[1]].ripflag) continue;
    if (mris->vertices[mris->faces[face_index].v[2]].ripflag) continue;
    numFaces++;
  }

  /* Set its attributes. */
  faces->intent = NIFTI_INTENT_TRIANGLE;
  faces->datatype = NIFTI_TYPE_INT32;
  faces->ind_ord = GIFTI_IND_ORD_ROW_MAJOR;
  faces->num_dim = 2;
  faces->dims[0] = numFaces;    /* In highest first, dim0 = rows */
  faces->dims[1] = 3;               /* In highest first, dim1 = cols */
  faces->encoding = GIFTI_ENCODING_B64GZ; // data stored in gzip'd base64
#if (BYTE_ORDER == LITTLE_ENDIAN)
  faces->endian = GIFTI_ENDIAN_LITTLE;
#else
  faces->endian = GIFTI_ENDIAN_BIG;
#endif
  faces->coordsys = NULL;
  faces->nvals = gifti_darray_nvals (faces);
  gifti_datatype_sizes (faces->datatype, &faces->nbyper, NULL);

  /* Allocate the data array. */
  faces->data = NULL;
  faces->data = (void*) calloc (faces->nvals, faces->nbyper);
  if (NULL == faces->data)
  {
    fprintf (stderr,"MRISwriteGIFTI: couldn't allocate faces data of "
             "length %d, element size %d\n", (int)faces->nvals, faces->nbyper);
    gifti_free_image (image);
    return ERROR_NOMEMORY;
  }

  /* Copy in all our face data (remembering to ignore faces which
     have a vertex with the ripflag set). */
  int faceNum = 0;
  for (face_index = 0; face_index < mris->nfaces; face_index++)
  {
    if (mris->vertices[mris->faces[face_index].v[0]].ripflag) continue;
    if (mris->vertices[mris->faces[face_index].v[1]].ripflag) continue;
    if (mris->vertices[mris->faces[face_index].v[2]].ripflag) continue;

    gifti_set_DA_value_2D (faces, faceNum, 0,
                           mris->faces[face_index].v[0]);
    gifti_set_DA_value_2D (faces, faceNum, 1,
                           mris->faces[face_index].v[1]);
    gifti_set_DA_value_2D (faces, faceNum, 2,
                           mris->faces[face_index].v[2]);
    faceNum++;
  }

  /* check for compliance */
  int valid = gifti_valid_gifti_image (image, 1);
  if (valid == 0)
  {
    fprintf (stderr,"MRISwriteGIFTI: GIFTI file %s is invalid!\n", fname);
    gifti_free_image (image);
    return ERROR_BADFILE;
  }

  /* Write the file. */
  if (gifti_write_image (image, fname, 1))
  {
    fprintf (stderr,"MRISwriteGIFTI: couldn't write image\n");
    gifti_free_image (image);
    return ERROR_BADFILE;
  }

  gifti_free_image (image);

  return ERROR_NONE;
}


/*-----------------------------------------------------
  Parameters:    MRIS_SURFACE structure,
                 output file name of GIFTI file,
                 input scalar data file

  Returns value: 0 if passed, else error code

  Description:   writes a GIFTI file containing 'shape'
                 data, ie. thickness, curv, sulc...
  ------------------------------------------------------*/
int MRISwriteScalarGIFTI(MRIS* mris,
                         const char *fname,
                         const char *scalar_fname)
{
  if (NULL == mris || NULL == fname)
  {
    fprintf (stderr,"MRISwriteScalarGIFTI: invalid input parameters\n");
    return ERROR_BADPARM;
  }

  gifti_image* image = (gifti_image *)calloc(1,sizeof(gifti_image));
  if (NULL == image)
  {
    fprintf (stderr,"MRISwriteScalarGIFTI: couldn't allocate gifti_image\n");
    return ERROR_NOMEMORY;
  }
  image->version = strcpyalloc(GIFTI_XML_VERSION);

  /* -------------------------------------------------------
   * Scalars array
   */
  if (MRISreadCurvatureFile(mris, scalar_fname))
  {
    fprintf (stderr,"MRISwriteScalarGIFTI: couldn't read %s\n",scalar_fname);
    gifti_free_image (image);
    return ERROR_BADFILE;
  }

  giiDataArray* scalars = gifti_alloc_and_add_darray (image);
  if (NULL == scalars)
  {
    fprintf (stderr,"MRISwriteScalarGIFTI: couldn't allocate giiDataArray\n");
    gifti_free_image (image);
    return ERROR_NOMEMORY;
  }

  /* Set its attributes. */
  scalars->intent = NIFTI_INTENT_SHAPE;
  scalars->datatype = NIFTI_TYPE_FLOAT32;
  scalars->ind_ord = GIFTI_IND_ORD_ROW_MAJOR;
  scalars->num_dim = 2;
  scalars->dims[0] = mris->nvertices;
  scalars->dims[1] = 1;
  scalars->encoding = GIFTI_ENCODING_B64GZ; // data stored in gzip'd base64
#if (BYTE_ORDER == LITTLE_ENDIAN)
  scalars->endian = GIFTI_ENDIAN_LITTLE;
#else
  scalars->endian = GIFTI_ENDIAN_BIG;
#endif
  scalars->coordsys = NULL;
  scalars->nvals = gifti_darray_nvals (scalars);
  gifti_datatype_sizes (scalars->datatype, &scalars->nbyper, NULL);

  /* include some metadata describing this thing */
  insertMetaData (mris, scalars); /* standard meta data */
  gifti_add_to_meta( &scalars->meta, "Name", scalar_fname, 1 );
  char *meta=NULL;
  if (strstr(scalar_fname, ".thickness")) meta = "Thickness";
  if (strstr(scalar_fname, ".curv"))      meta = "CurvatureRadial";
  if (strstr(scalar_fname, ".sulc"))      meta = "SulcalDepth";
  if (strstr(scalar_fname, ".area"))      meta = "Area";
  if (strstr(scalar_fname, ".volume"))    meta = "Volume";
  if (strstr(scalar_fname, ".jacobian"))  meta = "Jacobian";
  if (meta) gifti_add_to_meta( &scalars->meta, "ScalarDataType", meta, 1 );

  /* Allocate the data array. */
  scalars->data = NULL;
  scalars->data = (void*) calloc (scalars->nvals, scalars->nbyper);
  if (NULL == scalars->data)
  {
    fprintf (stderr,"MRISwriteScalarGIFTI: couldn't allocate scalars data of "
             "length %d, element size %d\n",
             (int)scalars->nvals,scalars->nbyper);
    gifti_free_image (image);
    return ERROR_NOMEMORY;
  }

  /* Copy in all our data. */
  int scalar_index;
  for (scalar_index = 0; scalar_index < mris->nvertices; scalar_index++)
  {
    if (mris->vertices[scalar_index].ripflag) continue;
    gifti_set_DA_value_2D (scalars, scalar_index, 0,
                           mris->vertices[scalar_index].curv);
  }

  /* check for compliance */
  int valid = gifti_valid_gifti_image (image, 1);
  if (valid == 0)
  {
    fprintf (stderr,"MRISwriteScalarGIFTI: GIFTI file %s is invalid!\n", 
             fname);
    gifti_free_image (image);
    return ERROR_BADFILE;
  }

  /* Write the file. */
  if (gifti_write_image (image, fname, 1))
  {
    fprintf (stderr,"MRISwriteScalarGIFTI: couldn't write image\n");
    gifti_free_image (image);
    return ERROR_BADFILE;
  }

  gifti_free_image (image);

  return ERROR_NONE;
}


/*-----------------------------------------------------------
  Parameters:    MRI structure (surface-encoded volume),
                 output file name of GIFTI file,

  Returns value: 0 if passed, else error code

  Description:   writes a GIFTI file containing functional or
                 timeseries data
  -----------------------------------------------------------*/
int mriWriteGifti(MRI* mri, const char *fname)
{
  if (NULL == mri || NULL == fname)
  {
    fprintf (stderr,"mriWriteGifti: invalid input parameters\n");
    return ERROR_BADPARM;
  }

  gifti_image* image = (gifti_image *)calloc(1,sizeof(gifti_image));
  if (NULL == image)
  {
    fprintf (stderr,"mriWriteGifti: couldn't allocate gifti_image\n");
    return ERROR_NOMEMORY;
  }
  image->version = strcpyalloc(GIFTI_XML_VERSION);

  /* -------------------------------------------------------
   * One DataArray for each 'frame' in the 'volume' data
   */
  int frame;
  for (frame=0; frame < mri->nframes; frame++)
  {
    giiDataArray* scalars = gifti_alloc_and_add_darray (image);
    if (NULL == scalars)
    {
      fprintf (stderr,"mriWriteGifti: couldn't allocate giiDataArray\n");
      gifti_free_image (image);
      return ERROR_NOMEMORY;
    }

    /* Set its attributes. */
    scalars->intent = NIFTI_INTENT_NONE;
    scalars->datatype = NIFTI_TYPE_FLOAT32;
    scalars->ind_ord = GIFTI_IND_ORD_ROW_MAJOR;
    scalars->num_dim = 2;
    scalars->dims[0] = mri->width;
    scalars->dims[1] = 1;
    scalars->encoding = GIFTI_ENCODING_B64GZ; // data stored in gzip'd base64
#if (BYTE_ORDER == LITTLE_ENDIAN)
    scalars->endian = GIFTI_ENDIAN_LITTLE;
#else
    scalars->endian = GIFTI_ENDIAN_BIG;
#endif
    scalars->coordsys = NULL;
    scalars->nvals = gifti_darray_nvals (scalars);
    gifti_datatype_sizes (scalars->datatype, &scalars->nbyper, NULL);

    /* include some metadata describing this thing */
    insertMetaData (NULL, scalars); /* standard meta data */

    /* Allocate the data array. */
    scalars->data = NULL;
    scalars->data = (void*) calloc (scalars->nvals, scalars->nbyper);
    if (NULL == scalars->data)
    {
      fprintf (stderr,"mriWriteGifti: couldn't allocate scalars data of "
               "length %d, element size %d\n",
               (int)scalars->nvals,scalars->nbyper);
      gifti_free_image (image);
      return ERROR_NOMEMORY;
    }

    /* Copy in all our data. */
    int scalar_index;
    for (scalar_index = 0; scalar_index < mri->width; scalar_index++)
    {
      float val = MRIgetVoxVal(mri,scalar_index,0,0,frame);
      gifti_set_DA_value_2D (scalars, scalar_index, 0, val);
    }

    // next frame
  }

  /* check for compliance */
  int valid = gifti_valid_gifti_image (image, 1);
  if (valid == 0)
  {
    fprintf (stderr,"mriWriteGifti: GIFTI file %s is invalid!\n", 
             fname);
    gifti_free_image (image);
    return ERROR_BADFILE;
  }

  /* Write the file. */
  if (gifti_write_image (image, fname, 1))
  {
    fprintf (stderr,"mriWriteGifti: couldn't write image\n");
    gifti_free_image (image);
    return ERROR_BADFILE;
  }

  gifti_free_image (image);

  return ERROR_NONE;
}


/*
 * Writes .annot data to a label table data gifti file:
 * puts the freesurfer colortable struct into a LabelTable,
 * and puts the .annotation field from each vertex into a 
 * DataArray
 */
int MRISwriteLabelTableGIFTI(MRI_SURFACE *mris, const char *fname)
{
  if (NULL == mris || NULL == fname)
  {
    fprintf (stderr,"MRISwriteLabelTableGIFTI: invalid input parameters\n");
    return ERROR_BADPARM;
  }

  gifti_image* image = (gifti_image *)calloc(1,sizeof(gifti_image));
  if (NULL == image)
  {
    fprintf (stderr,"MRISwriteLabelTableGIFTI: "
             "couldn't allocate gifti_image\n");
    return ERROR_NOMEMORY;
  }
  image->version = strcpyalloc(GIFTI_XML_VERSION);

  /* -------------------------------------------------------
   * LabelTable struct, fill it in with our colortable stuff
   */
  giiLabelTable labeltable;
  labeltable.length = mris->ct->nentries;
  if (labeltable.length == 0)
  {
    fprintf(stderr, "MRISwriteLabelTableGIFTI: "
             "colortable is empty!\n");
    return ERROR_BADFILE;
  }
  labeltable.index = (int *)calloc(labeltable.length, sizeof(int *));
  labeltable.label = (char **)calloc(labeltable.length, sizeof(char *));
  labeltable.rgba = (float *)calloc(labeltable.length, 4*sizeof(float *));
  if ((NULL == labeltable.index) ||
      (NULL == labeltable.label) ||
      (NULL == labeltable.rgba))
  {
    fprintf (stderr,"MRISwriteLabelTableGIFTI: "
             "couldn't allocate giftiLabelTable\n");
    return ERROR_NOMEMORY;
  }
  float* rgba = labeltable.rgba;
  int idx;
  for (idx=0; idx < labeltable.length; idx++)
  {
    labeltable.index[idx] = CTABrgb2Annotation(mris->ct->entries[idx]->ri,
                                               mris->ct->entries[idx]->gi,
                                               mris->ct->entries[idx]->bi);
    //printf("%8.8X\n",labeltable.index[idx]);

    labeltable.label[idx] = strcpyalloc(mris->ct->entries[idx]->name);

    rgba[0] = mris->ct->entries[idx]->rf;
    rgba[1] = mris->ct->entries[idx]->gf;
    rgba[2] = mris->ct->entries[idx]->bf;
    rgba[3] = 1.0f;
    rgba += 4;
    /*
    printf("RGBA: %d %d %d %d %f %f %f %f\n",
           mris->ct->entries[idx]->ri,
           mris->ct->entries[idx]->gi,
           mris->ct->entries[idx]->bi,
           mris->ct->entries[idx]->ai,
           mris->ct->entries[idx]->rf,
           mris->ct->entries[idx]->gf,
           mris->ct->entries[idx]->bf,
           mris->ct->entries[idx]->af);
    */
  }
  // dont forget to stick us in the image
  image->labeltable = labeltable;
  //gifti_disp_LabelTable(NULL,&image->labeltable);

  /* -------------------------------------------------------
   * LabelTables array
   */
  giiDataArray* labels = gifti_alloc_and_add_darray (image);
  if (NULL == labels)
  {
    fprintf (stderr,"MRISwriteLabelTableGIFTI: "
             "couldn't allocate giiDataArray\n");
    gifti_free_image (image);
    return ERROR_NOMEMORY;
  }

  /* Set its attributes. */
  labels->intent = NIFTI_INTENT_LABEL;
  labels->datatype = NIFTI_TYPE_UINT32;
  labels->ind_ord = GIFTI_IND_ORD_ROW_MAJOR;
  labels->num_dim = 1;
  labels->dims[0] = mris->nvertices;
  labels->dims[1] = 1;
  labels->encoding = GIFTI_ENCODING_B64GZ; // data stored in gzip'd base64
#if (BYTE_ORDER == LITTLE_ENDIAN)
  labels->endian = GIFTI_ENDIAN_LITTLE;
#else
  labels->endian = GIFTI_ENDIAN_BIG;
#endif
  labels->coordsys = NULL;
  labels->nvals = gifti_darray_nvals (labels);
  gifti_datatype_sizes (labels->datatype, &labels->nbyper, NULL);

  /* include some metadata describing this thing */
  insertMetaData (mris, labels); /* standard meta data */
  gifti_add_to_meta( &labels->meta, "Name", "node label", 1 );

  /* Allocate the data array. */
  labels->data = NULL;
  labels->data = (void*) calloc (labels->nvals, labels->nbyper);
  if (NULL == labels->data)
  {
    fprintf (stderr,"MRISwriteLabelTableGIFTI: "
             "couldn't allocate labels data of "
             "length %d, element size %d\n",
             (int)labels->nvals,labels->nbyper);
    gifti_free_image (image);
    return ERROR_NOMEMORY;
  }

  /* Copy our 'annotation' data for each vertex */
  unsigned int* label_data = labels->data;
  int label_index;
  for (label_index = 0; label_index < mris->nvertices; label_index++)
  {
    if (mris->vertices[label_index].ripflag) continue;
    *(label_data + label_index) = mris->vertices[label_index].annotation;
    //printf("%8.8X ", *(label_data + label_index));
  }

  /* check for compliance */
  int valid = gifti_valid_gifti_image (image, 1);
  if (valid == 0)
  {
    fprintf (stderr,"MRISwriteLabelTableGIFTI: GIFTI file %s is invalid!\n", 
             fname);
    gifti_free_image (image);
    return ERROR_BADFILE;
  }

  /* Write the file. */
  if (gifti_write_image (image, fname, 1))
  {
    fprintf (stderr,"MRISwriteLabelTableGIFTI: couldn't write image\n");
    gifti_free_image (image);
    return ERROR_BADFILE;
  }

  gifti_free_image (image);

  return ERROR_NONE;
}


