/**
 * @file  gifti_local.c
 * @brief local utilities for GIFTI library
 *
 * This file has some some extra functions for use with the GIFTI
 * utilities. The official utilities reside in gifti.c and gifti_xml.c
 * 
 */
/*
 * Original Author: Kevin Teich 
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/08/02 21:07:27 $
 *    $Revision: 1.1 $
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
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>

#include "gifti_local.h"

DataArray* gifti_alloc_and_add_darray (gifti_image* image)
{
  if (!image)
    {
      fprintf (stderr,"** gifti_alloc_and_add_darray: NULL image\n");
      return NULL;
    }

  /* Try to add an empty array. */
  if (gifti_add_empty_darray(image)) 
    {
      fprintf (stderr,"** gifti_alloc_and_add_darray: gifti_add_empty_darray "
	       "failed\n");
      return NULL;
    }

  /* Return the array we just allocated. */
  return image->darray[image->numDA-1]; 
}


double gifti_get_DA_value_2D (DataArray* da, int row, int col)
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
  if (GIFTI_IND_ORD_HIGH2LOW == da->ind_ord)
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
    case 2: {       /* NIFTI_TYPE_UINT8 */
      if( GIFTI_IND_ORD_HIGH2LOW == da->ind_ord )
	return (double)*((unsigned char*)
			 (da->data) + (dim0_index*da->dims[1]) + dim1_index);
      else
	return (double)*((unsigned char*)
			 (da->data) + dim0_index + (dim1_index*da->dims[0]));
      break;
    }
    case 4: {       /* NIFTI_TYPE_INT16 */
      if( GIFTI_IND_ORD_HIGH2LOW == da->ind_ord )
	return (double)*((short*)
			 (da->data) + (dim0_index*da->dims[1]) + dim1_index);
      else
	return (double)*((short*)
			 (da->data) + dim0_index + (dim1_index*da->dims[0]));
      break;
    }
    case 8: {       /* NIFTI_TYPE_INT32 */
      if( GIFTI_IND_ORD_HIGH2LOW == da->ind_ord )
	return (double)*((int*)
			 (da->data) + (dim0_index*da->dims[1]) + dim1_index);
      else
	return (double)*((int*)
			 (da->data) + dim0_index + (dim1_index*da->dims[0]));
      break;
    }
    case 16: {      /* NIFTI_TYPE_FLOAT32 */
      if( GIFTI_IND_ORD_HIGH2LOW == da->ind_ord )
	return (double)*((float*)
			 (da->data) + (dim0_index*da->dims[1]) + dim1_index);
      else
	return (double)*((float*)
			 (da->data) + dim0_index + (dim1_index*da->dims[0]));
      break;
    }
    case 256: {     /* NIFTI_TYPE_INT8 */
      if( GIFTI_IND_ORD_HIGH2LOW == da->ind_ord )
	return (double)*((char*)
			 (da->data) + (dim0_index*da->dims[1]) + dim1_index);
      else
	return (double)*((char*)
			 (da->data) + dim0_index + (dim1_index*da->dims[0]));
      break;
    }
    case 512: {     /* NIFTI_TYPE_UINT16 */
      if( GIFTI_IND_ORD_HIGH2LOW == da->ind_ord )
	return (double)*((unsigned short*)
			 (da->data) + (dim0_index*da->dims[1]) + dim1_index);
      else
	return (double)*((unsigned short*)
			 (da->data) + dim0_index + (dim1_index*da->dims[0]));
      break;
    }
    case 768: {     /* NIFTI_TYPE_UINT32 */
      if( GIFTI_IND_ORD_HIGH2LOW == da->ind_ord )
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

void gifti_set_DA_value_2D (DataArray* da, int row, int col, double value)
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
  if (GIFTI_IND_ORD_HIGH2LOW == da->ind_ord)
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
    case 2: {       /* NIFTI_TYPE_UINT8 */
      if( GIFTI_IND_ORD_HIGH2LOW == da->ind_ord )
	*((unsigned char*)(da->data) + (dim0_index*da->dims[1]) + dim1_index) =
	  (unsigned char)value;
      else
	*((unsigned char*)(da->data) + dim0_index + (dim1_index*da->dims[0])) =
	  (unsigned char)value;
      break;
    }
    case 4: {       /* NIFTI_TYPE_INT16 */
      if( GIFTI_IND_ORD_HIGH2LOW == da->ind_ord )
	*((short*)(da->data) + (dim0_index*da->dims[1]) + dim1_index) =
	  (short)value;
      else
	*((short*)(da->data) + dim0_index + (dim1_index*da->dims[0])) =
	  (short)value;
      break;
    }
    case 8: {       /* NIFTI_TYPE_INT32 */
      if( GIFTI_IND_ORD_HIGH2LOW == da->ind_ord )
	*((int*)(da->data) + (dim0_index*da->dims[1]) + dim1_index) =
	  (int)value;
      else
	*((int*)(da->data) + dim0_index + (dim1_index*da->dims[0])) =
	  (int)value;
      break;
    }
    case 16: {      /* NIFTI_TYPE_FLOAT32 */
      if( GIFTI_IND_ORD_HIGH2LOW == da->ind_ord )
	*((float*)(da->data) + (dim0_index*da->dims[1]) + dim1_index) =
	  (float)value;
      else
	*((float*)(da->data) + dim0_index + (dim1_index*da->dims[0])) =
	  (float)value;
      break;
    }
    case 256: {     /* NIFTI_TYPE_INT8 */
      if( GIFTI_IND_ORD_HIGH2LOW == da->ind_ord )
	*((char*)(da->data) + (dim0_index*da->dims[1]) + dim1_index) =
	  (char)value;
      else
	*((char*)(da->data) + dim0_index + (dim1_index*da->dims[0])) =
	  (char)value;
      break;
    }
    case 512: {     /* NIFTI_TYPE_UINT16 */
      if( GIFTI_IND_ORD_HIGH2LOW == da->ind_ord )
	*((unsigned short*)(da->data) + (dim0_index*da->dims[1]) + dim1_index) =
	  (unsigned short)value;
      else
	*((unsigned short*)(da->data) + dim0_index + (dim1_index*da->dims[0])) =
	  (unsigned short)value;
      break;
    }
    case 768: {     /* NIFTI_TYPE_UINT32 */
      if( GIFTI_IND_ORD_HIGH2LOW == da->ind_ord )
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
