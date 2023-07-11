 #define COMPILING_MRISURF_TOPOLOGY_FRIEND_CHECKED
/**
 * @brief local utilities for GIFTI library
 *
 * This file has some some extra functions for use with the GIFTI
 * utilities. The official utilities reside in gifti_io.c and gifti_xml.c
 *
 */
/*
 * Original Authors: Kevin Teich and Nick Schmansky
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

#include <pwd.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#include "error.h"  // return codes
#include "gifti.h"
#include "nifti1.h"
#include "timer.h"
#include "utils.h"  // strcpyalloc
#include "libgen.h"

#define TAG_CMDLINE_LEN   1024
#define NUM_VOLGEOM_META  19

/*
 *
 */
static giiDataArray *gifti_alloc_and_add_darray(gifti_image *image)
{
  if (!image) {
    fprintf(stderr, "** gifti_alloc_and_add_darray: NULL image\n");
    return NULL;
  }

  /* Try to add an empty array. */
  if (gifti_add_empty_darray(image, 1)) {
    fprintf(stderr,
            "** gifti_alloc_and_add_darray: gifti_add_empty_darray "
            "failed\n");
    return NULL;
  }

  /* Return the array we just allocated. */
  return image->darray[image->numDA - 1];
}

/*
 *
 */
static double gifti_get_DA_value_2D(giiDataArray *da, int row, int col)
{
  int dim0_index, dim1_index;
  int dims_0 = 0, dims_1 = 0;

  if (!da || !da->data) {
    fprintf(stderr, "** gifti_get_DA_value_2D, invalid params: data=%p\n", da);
    exit(1);
  }

  if (da->num_dim == 1) {
    // support for using this routine to read 1D data, under one condition...
    if (col != 0) {
      fprintf(stderr,
              "** gifti_get_DA_value_2D, array dim is 1 "
              "but trying to access 2D data element (col=%d)\n",
              col);
      exit(1);
    }
    dims_0 = da->dims[0];
    dims_1 = 1;  // 1D data
  }
  else if (da->num_dim != 2) {
    fprintf(stderr, "** gifti_get_DA_value_2D, array dim is %d\n", da->num_dim);
    exit(1);
  }
  else {
    dims_0 = da->dims[0];
    dims_1 = da->dims[1];
  }

  /* Get the dim0 and dims[1] indices based on our order. */
  if (GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord) {
    dim0_index = row;
    dim1_index = col;
  }
  else if (GIFTI_IND_ORD_COL_MAJOR == da->ind_ord) {
    // NJS NOTE: notice that order is treated as row/col, so that the
    // calling sequence can just assume row major
    dim0_index = row;  // col;
    dim1_index = col;  // row;
  }
  else {
    fprintf(stderr, "** gifti_get_DA_value_2D, unknown ind_ord: %d\n", da->ind_ord);
    exit(1);
  }
  if (da->num_dim == 1)  // support for using this routine to read 1D data
  {
    dim0_index = row;
    dim1_index = col;
  }

  /* Check the indices. */
  if (dim0_index < 0 || dim0_index >= dims_0 || dim1_index < 0 || dim1_index >= dims_1) {
    fprintf(stderr,
            "** gifti_get_DA_value_2D, invalid params: "
            "dim0_index=%d (max=%d), dim1_index=%d (max=%d)\n",
            dim0_index,
            dims_0,
            dim1_index,
            dims_1);
    exit(1);
  }

  /* Switch on the data type and return the appropriate
     element. Indexing depends on the data order. */
  switch (da->datatype) {
    default:
      fprintf(stderr,
              "** gifti_get_DA_value_2D, unsupported type %d-"
              "unknown, or can't convert to double\n",
              da->datatype);
      exit(1);
    case NIFTI_TYPE_UINT8: {
      if (GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord)
        return (double)*((unsigned char *)(da->data) + (dim0_index * dims_1) + dim1_index);
      else
        return (double)*((unsigned char *)(da->data) + dim0_index + (dim1_index * dims_0));
      break;
    }
    case NIFTI_TYPE_INT16: {
      if (GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord)
        return (double)*((short *)(da->data) + (dim0_index * dims_1) + dim1_index);
      else
        return (double)*((short *)(da->data) + dim0_index + (dim1_index * dims_0));
      break;
    }
    case NIFTI_TYPE_INT32: {
      if (GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord)
        return (double)*((int *)(da->data) + (dim0_index * dims_1) + dim1_index);
      else
        return (double)*((int *)(da->data) + dim0_index + (dim1_index * dims_0));
      break;
    }
    case NIFTI_TYPE_FLOAT32: {
      if (GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord)
        return (double)*((float *)(da->data) + (dim0_index * dims_1) + dim1_index);
      else
        return (double)*((float *)(da->data) + dim0_index + (dim1_index * dims_0));
      break;
    }
    case NIFTI_TYPE_FLOAT64: {
      if (GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord)
        return (double)*((double *)(da->data) + (dim0_index * dims_1) + dim1_index);
      else
        return (double)*((double *)(da->data) + dim0_index + (dim1_index * dims_0));
      break;
    }
    case NIFTI_TYPE_COMPLEX64: {
      if (GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord)
        return (double)*((double *)(da->data) + (dim0_index * dims_1) + dim1_index);
      else
        return (double)*((double *)(da->data) + dim0_index + (dim1_index * dims_0));
      break;
    }
    case NIFTI_TYPE_INT8: {
      if (GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord)
        return (double)*((char *)(da->data) + (dim0_index * dims_1) + dim1_index);
      else
        return (double)*((char *)(da->data) + dim0_index + (dim1_index * dims_0));
      break;
    }
    case NIFTI_TYPE_UINT16: {
      if (GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord)
        return (double)*((unsigned short *)(da->data) + (dim0_index * dims_1) + dim1_index);
      else
        return (double)*((unsigned short *)(da->data) + dim0_index + (dim1_index * dims_0));
      break;
    }
    case NIFTI_TYPE_UINT32: {
      if (GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord)
        return (double)*((unsigned int *)(da->data) + (dim0_index * dims_1) + dim1_index);
      else
        return (double)*((unsigned int *)(da->data) + dim0_index + (dim1_index * dims_0));
      break;
    }
  }

  exit(1);
}

/*
 *
 */
static void gifti_set_DA_value_2D(giiDataArray *da, int row, int col, double value)
{
  int dim0_index, dim1_index;
  int dims_0 = 0, dims_1 = 0;

  if (!da || !da->data) {
    fprintf(stderr, "** gifti_set_DA_value_2D, invalid params: data=%p\n", da);
    exit(1);
  }

  if (da->num_dim == 1) {
    // support for using this routine to write 1D data, under one condition...
    if (col != 0) {
      fprintf(stderr,
              "** gifti_set_DA_value_2D, array dim is 1 "
              "but trying to access 2D data element (col=%d)\n",
              col);
      exit(1);
    }
    dims_0 = da->dims[0];
    dims_1 = 1;  // 1D data
  }
  else if (da->num_dim != 2) {
    fprintf(stderr, "** gifti_set_DA_value_2D, array dim is %d\n", da->num_dim);
    exit(1);
  }
  else {
    dims_0 = da->dims[0];
    dims_1 = da->dims[1];
  }

  /* Get the dim0 and dims[1] indices based on our order. */
  if (GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord) {
    dim0_index = row;
    dim1_index = col;
  }
  else {
    dim0_index = col;
    dim1_index = row;
  }
  if (da->num_dim == 1)  // support for using this routine to read 1D data
  {
    dim0_index = row;
    dim1_index = col;
  }

  /* Check the indices. */
  if (dim0_index < 0 || dim0_index >= dims_0 || dim1_index < 0 || dim1_index >= dims_1) {
    fprintf(stderr,
            "** gifti_set_DA_value_2D, invalid params: "
            "dim0_index=%d (max=%d), dim1_index=%d (max=%d)\n",
            dim0_index,
            dims_0,
            dim1_index,
            dims_1);
    return;
  }

  /* Switch on the data type and write the appropriate
     element. Indexing depends on the data order. */
  switch (da->datatype) {
    default:
      fprintf(stderr,
              "** gifti_set_DA_value_2D, unsupported type %d-"
              "unknown, or can't convert to double\n",
              da->datatype);
      return;
    case NIFTI_TYPE_UINT8: {
      if (GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord)
        *((unsigned char *)(da->data) + (dim0_index * dims_1) + dim1_index) = (unsigned char)value;
      else
        *((unsigned char *)(da->data) + dim0_index + (dim1_index * dims_0)) = (unsigned char)value;
      break;
    }
    case NIFTI_TYPE_INT16: {
      if (GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord)
        *((short *)(da->data) + (dim0_index * dims_1) + dim1_index) = (short)value;
      else
        *((short *)(da->data) + dim0_index + (dim1_index * dims_0)) = (short)value;
      break;
    }
    case NIFTI_TYPE_INT32: {
      if (GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord)
        *((int *)(da->data) + (dim0_index * dims_1) + dim1_index) = (int)value;
      else
        *((int *)(da->data) + dim0_index + (dim1_index * dims_0)) = (int)value;
      break;
    }
    case NIFTI_TYPE_FLOAT32: {
      if (GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord)
        *((float *)(da->data) + (dim0_index * dims_1) + dim1_index) = (float)value;
      else
        *((float *)(da->data) + dim0_index + (dim1_index * dims_0)) = (float)value;
      break;
    }
    case NIFTI_TYPE_INT8: {
      if (GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord)
        *((char *)(da->data) + (dim0_index * dims_1) + dim1_index) = (char)value;
      else
        *((char *)(da->data) + dim0_index + (dim1_index * dims_0)) = (char)value;
      break;
    }
    case NIFTI_TYPE_UINT16: {
      if (GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord)
        *((unsigned short *)(da->data) + (dim0_index * dims_1) + dim1_index) = (unsigned short)value;
      else
        *((unsigned short *)(da->data) + dim0_index + (dim1_index * dims_0)) = (unsigned short)value;
      break;
    }
    case NIFTI_TYPE_UINT32: {
      if (GIFTI_IND_ORD_ROW_MAJOR == da->ind_ord)
        *((unsigned int *)(da->data) + (dim0_index * dims_1) + dim1_index) = (unsigned int)value;
      else
        *((unsigned int *)(da->data) + dim0_index + (dim1_index * dims_0)) = (unsigned int)value;
      break;
    }
  }

  return;
}


/*
 *
 */
static COLOR_TABLE *makeColorTable(std::map<int, float*> &unique_annot_map, const char *fname)
{
    /* copy label table contents to our color_table struct */
    COLOR_TABLE *ct = (COLOR_TABLE *)calloc(1, sizeof(COLOR_TABLE));
    if (ct == NULL) {
      fprintf(stderr, "makeColorTable(): could not alloc colortable memory\n");
      return NULL;
    }
    memset(ct, 0, sizeof(COLOR_TABLE));
    ct->nentries = unique_annot_map.size();
    ct->version = CTAB_VERSION_TO_WRITE;  //2;
    ct->entries = (COLOR_TABLE_ENTRY **)calloc(ct->nentries + 1, sizeof(COLOR_TABLE_ENTRY *));
    if (ct->entries == NULL) {
      fprintf(stderr, "makeColorTable(): could not alloc colortable entries\n");
      return NULL;
    }
    // memset(ct->entries,0,sizeof(ct->entries)); // original
    memset(ct->entries, 0, sizeof(*ct->entries));  // changed by dng
    strncpy(ct->fname, fname, sizeof(ct->fname)-1);

    std::map<int, float*>::iterator it = unique_annot_map.begin();
    int label_index = 0;
    while (it != unique_annot_map.end())
    {
      ct->entries[label_index] = (CTE *)malloc(sizeof(CTE));
      if (ct->entries[label_index] == NULL) {
        fprintf(stderr, "makeColorTable(): could not alloc colortable entry\n");
        return NULL;
      }
      char unknown_label[256] = {'\0'};
      sprintf(unknown_label, "Unknown_Label_%d", label_index);
      strncpy(
          ct->entries[label_index]->name,
	  unknown_label,
	  sizeof(ct->entries[label_index]->name)-1);

      float *rgba = it->second;
      ct->entries[label_index]->rf = rgba[0];
      ct->entries[label_index]->ri = floor((rgba[0]) * 256);
      if (ct->entries[label_index]->ri > 255) {
        ct->entries[label_index]->ri = 255;
      }
      ct->entries[label_index]->gf = rgba[1];
      ct->entries[label_index]->gi = floor((rgba[1]) * 256);
      if (ct->entries[label_index]->gi > 255) {
        ct->entries[label_index]->gi = 255;
      }
      ct->entries[label_index]->bf = rgba[2];
      ct->entries[label_index]->bi = floor((rgba[2]) * 256);
      if (ct->entries[label_index]->bi > 255) {
        ct->entries[label_index]->bi = 255;
      }
      ct->entries[label_index]->af = rgba[3];
      ct->entries[label_index]->ai = floor((rgba[3]) * 256);
      if (ct->entries[label_index]->ai > 255) {
        ct->entries[label_index]->ai = 255;
      }

      printf("%d %s %d (%d %d %d %d)\n",
             label_index, ct->entries[label_index]->name, it->first,
             ct->entries[label_index]->ri, ct->entries[label_index]->gi, ct->entries[label_index]->bi, ct->entries[label_index]->ai);

      it++;
      label_index++;
    }
    ct->entries[label_index] = NULL;
    CTABfindDuplicateNames(ct);

    return ct;
}


/*-------------------------------------------------------------------
  Parameters:    input file name of GIFTI file
                 optional mris structure to store data found
                 optional data array number to read

  Returns value: freesurfer surface structure

  Description:   reads a GIFTI file, putting vertices,
                 and faces into an MRIS_SURFACE structure,
                 along with any other data, like labels,
                 colors, curv data, stats or values.
                 if daNum is not -1, then read only the
                 data in data array number daNum
  -------------------------------------------------------------------*/
MRIS *mrisReadGIFTIdanum(const char *fname, MRIS *mris, int daNum, std::vector<OverlayInfoStruct> *poverlayinfo)
{
  /*
   * attempt to read the file
   */
  gifti_image *image = gifti_read_image(fname, 1);
  if (NULL == image) {
    fprintf(stderr, "mrisReadGIFTIdanum: gifti_read_image() returned NULL\n");
    return NULL;
  }

  // make sure version is recoded before validation
  if (!strcmp(image->version, "1")) {
    free(image->version);
    image->version = strcpyalloc(GIFTI_XML_VERSION);
  }

  /*
   * check for compliance
   */
  int valid = gifti_valid_gifti_image(image, 1);
  if (valid == 0) {
    fprintf(stderr, "mrisReadGIFTIdanum: GIFTI file %s is invalid!\n", fname);
    gifti_free_image(image);
    return NULL;
  }

  /*
   * check for 'LabelTable' data and read into our colortable if exists
   */
  COLOR_TABLE *ct = NULL;
  int maxkey = -1;
  if (image->labeltable.length > 0) {
    /* check validity of labeltable data */
    if (!gifti_valid_LabelTable(&image->labeltable, 1)) {
      fprintf(stderr, "mrisReadGIFTIdanum: invalid labeltable found in file %s\n", fname);
      gifti_free_image(image);
      return NULL;
    }

    
    /* 1. GIFTI LabelTable is a LUT. 
     *    a. It is used by DataArrays whose values are an key into the LabelTable’s labels. 
     *       A file should contain at most one LabelTable and it must be located in the file prior to any DataArray elements.
     *    b. The label keys are non-negative integers. They are not necessarily ordered or sequential.
     *    c. The label keys serve as index to COLOR_TABLE.  Number of COLOR_TABLE entries to create is max(label key) + 1.
     *       Missing numbers in label keys will leave empty COLOR_TABLE entries.
     *
     * 2. LUT is read into COLOR_TABLE in CTABreadASCII2():
     *    a. first colum in LUT becomes COLOR_TABLE index 
     *    b. number of COLOR_TABLE entries to create is max(1st LUT column) + 1
     *    c. COLOR_TABLE index is 0 .. n
     *    d. w/ or w/o '0  Unknown 0    0      0      0' in the first line, index 0 will be there
     *       the difference: w/  the line, CTABfindAnnotation() returns 0  for annotation=0;
     *                       w/o the line, CTABfindAnnotation() returns -1 for annotation=0
     *    e. skipped numbers in LUT will leave holes (empty entries) in COLOR_TABLE,
     *       their corresponding index will have unknown annotations
     */

    /* Scan through the file and see what our max label key is. 
     * Create COLOR_TABLE with (maxkey + 1) entries.
     */
    for (int nlabel = 0; nlabel < image->labeltable.length; nlabel++)
      maxkey = (image->labeltable.key[nlabel] > maxkey) ? image->labeltable.key[nlabel] : maxkey;

    /* copy label table contents to our color_table struct */
    ct = (COLOR_TABLE *)calloc(1, sizeof(COLOR_TABLE));
    if (ct == NULL) {
      fprintf(stderr, "mrisReadGIFTIdanum: could not alloc colortable memory\n");
      gifti_free_image(image);
      return NULL;
    }
    memset(ct, 0, sizeof(COLOR_TABLE));
    ct->nentries = maxkey + 1;  //image->labeltable.length;
    ct->version = CTAB_VERSION_TO_WRITE;  //2;
    ct->entries = (COLOR_TABLE_ENTRY **)calloc(ct->nentries, sizeof(COLOR_TABLE_ENTRY *));
    if (ct->entries == NULL) {
      fprintf(stderr, "mrisReadGIFTIdanum: could not alloc colortable entries\n");
      gifti_free_image(image);
      return NULL;
    }
    // memset(ct->entries,0,sizeof(ct->entries)); // original
    memset(ct->entries, 0, sizeof(*ct->entries));  // changed by dng
    strncpy(ct->fname, fname, sizeof(ct->fname)-1);

    float *rgba = image->labeltable.rgba;
    if (NULL == rgba) {
      // optional rgba values are missing, so we must create colors for the labels
      image->labeltable.rgba = (float *)calloc(image->labeltable.length, 4 * sizeof(float *));
      if (NULL == image->labeltable.rgba) {
        fprintf(stderr,
                "mrisReadGIFTIdanum: "
                "couldn't allocate memory for labeltable.rgba\n");
        return NULL;
      }
      rgba = image->labeltable.rgba;
      setRandomSeed(12);  // so that color generation is consistent
      int label_index;
      for (label_index = 0; label_index < image->labeltable.length; label_index++) {
        rgba[0] = (float)randomNumber(0.0f, 1.0f);
        rgba[1] = (float)randomNumber(0.0f, 1.0f);
        rgba[2] = (float)randomNumber(0.0f, 1.0f);
        rgba[3] = 1.0f;
        rgba += 4;
      }
    }

    if (Gdiag & DIAG_SHOW)
      printf("[DEBUG] mrisReadGIFTIdanum(): ct->nentries=%d, num_entries_to_read=%d\n", ct->nentries, image->labeltable.length);

    rgba = image->labeltable.rgba;
    int label_index;
    // loop through LabelTable
    for (label_index = 0; label_index < image->labeltable.length; label_index++) {
      // retrieve the label key, use it as COLOR_TABLE index
      int labelkey = image->labeltable.key[label_index];
      if (ct->entries[labelkey] != NULL)
      {
        printf("mrisReadGIFTIdanum(%s): Duplicate labelkey %d:%s, was %s\n",
               fname, labelkey, image->labeltable.label[label_index], ct->entries[labelkey]->name);
      }
      else
      {
        ct->entries[labelkey] = (CTE *)malloc(sizeof(CTE));
        if (ct->entries[labelkey] == NULL) {
          fprintf(stderr, "mrisReadGIFTIdanum: could not alloc colortable entry\n");
          gifti_free_image(image);
          return NULL;
        }

        if (Gdiag & DIAG_SHOW)
          printf("[DEBUG] mrisReadGIFTIdanum(): created ct->entries[%d]\n", labelkey);

        strncpy(
            ct->entries[labelkey]->name,
	    image->labeltable.label[label_index],
	    sizeof(ct->entries[labelkey]->name)-1);

        ct->entries[labelkey]->rf = rgba[0];
        ct->entries[labelkey]->ri = floor((rgba[0]) * 256);
        if (ct->entries[labelkey]->ri > 255) {
          ct->entries[labelkey]->ri = 255;
        }

        ct->entries[labelkey]->gf = rgba[1];
        ct->entries[labelkey]->gi = floor((rgba[1]) * 256);
        if (ct->entries[labelkey]->gi > 255) {
          ct->entries[labelkey]->gi = 255;
        }

        ct->entries[labelkey]->bf = rgba[2];
        ct->entries[labelkey]->bi = floor((rgba[2]) * 256);
        if (ct->entries[labelkey]->bi > 255) {
          ct->entries[labelkey]->bi = 255;
        }

        ct->entries[labelkey]->af = rgba[3];
        ct->entries[labelkey]->ai = floor((rgba[3]) * 256);
        if (ct->entries[labelkey]->ai > 255) {
          ct->entries[labelkey]->ai = 255;
        }

        rgba += 4;
      }
      /*
        printf("RGBA: %d %d %d %d %f %f %f %f\n",
           ct->entries[labelkey]->ri,
           ct->entries[labelkey]->gi,
           ct->entries[labelkey]->bi,
           ct->entries[labelkey]->ai,
           ct->entries[labelkey]->rf,
           ct->entries[labelkey]->gf,
           ct->entries[labelkey]->bf,
           ct->entries[labelkey]->af);
      */
    }
    //ct->entries[label_index] = NULL;
    CTABfindDuplicateNames(ct);

    // we're done, except that the colortable struct 'ct' will get stored
    // in the mris structure element 'ct' at the end of this routine, after
    // mris is known to exist
  }
  // end of LabelTable parsing (into colortable)

  /*
   * Now parse the DataArrays looking for coordinate and face data arrays,
   * so that we can create our mris structure
   */
  giiDataArray *coords = NULL;
  giiDataArray *faces = NULL;
  int numDA;
  for (numDA = 0; numDA < image->numDA; numDA++) {
    if (image->darray[numDA]->intent == NIFTI_INTENT_POINTSET) {
      coords = image->darray[numDA];
    }
    else if (image->darray[numDA]->intent == NIFTI_INTENT_TRIANGLE) {
      faces = image->darray[numDA];
    }
  }

  /*
   * if we found coordinate and face data....create mris struct and fill it
   */
  if (coords && faces) {
    /* Check the number of vertices and faces. */
    long long num_vertices = 0;
    long long num_cols = 0;
    if (coords->ind_ord == GIFTI_IND_ORD_ROW_MAJOR) {
      // RowMajorOrder
      gifti_DA_rows_cols(coords, &num_vertices, &num_cols);
    }
    else {
      // ColumnMajorOrder
      gifti_DA_rows_cols(coords, &num_cols, &num_vertices);
    }
    if (num_vertices <= 0 || num_cols != 3) {
      fprintf(stderr,
              "mrisReadGIFTIdanum: malformed coords data array in file "
              "%s: num_vertices=%d num_cols=%d\n",
              fname,
              (int)num_vertices,
              (int)num_cols);
      gifti_free_image(image);
      return NULL;
    }
    long long num_faces = 0;
    num_cols = 0;
    if (faces->ind_ord == GIFTI_IND_ORD_ROW_MAJOR) {
      // RowMajorOrder
      gifti_DA_rows_cols(faces, &num_faces, &num_cols);
    }
    else {
      // ColumnMajorOrder
      gifti_DA_rows_cols(faces, &num_cols, &num_faces);
    }
    if (num_faces <= 0 || num_cols != 3) {
      fprintf(stderr,
              "mrisReadGIFTIdanum: malformed faces data array in file "
              "%s: num_faces=%d num_cols=%d\n",
              fname,
              (int)num_faces,
              (int)num_cols);
      gifti_free_image(image);
      return NULL;
    }

    /* Try to allocate a surface. */
    mris = MRISalloc(num_vertices, num_faces);
    if (NULL == mris) {
      fprintf(stderr,
              "mrisReadGIFTIdanum: failed to allocate an MRIS with "
              "%d vertices and %d faces\n",
              (int)num_vertices,
              (int)num_faces);
      gifti_free_image(image);
      return NULL;
    }

    /* Set some meta data in the mris. */
    strcpy(mris->fname, fname);
    mris->type = MRIS_TRIANGULAR_SURFACE;  // ??? is this correct ???
    char *hemi = gifti_get_meta_value(&coords->meta, "AnatomicalStructurePrimary");
    if (hemi && (strcmp(hemi, "CortexRight") == 0)) {
      mris->hemisphere = RIGHT_HEMISPHERE;
    }
    else if (hemi && (strcmp(hemi, "CortexLeft") == 0)) {
      mris->hemisphere = LEFT_HEMISPHERE;
    }
    else {
      mris->hemisphere = NO_HEMISPHERE;
    }

    /* gifti uses real RAS by default */
    // This is not correct. As of 12/19/2022, the default surface XYZ is in tkregister space.
    // mris->useRealRAS will be set accordingly based on the value in <dataspace> tag later.
    mris->useRealRAS = 0;

    /* retrieve volume geometry info */
    {
      int vgvalid = 0;  // there are a total of 19 values (NUM_VOLGEOM_META)
      char *stmp = gifti_get_meta_value(&coords->meta, "VolGeomWidth");
      if (stmp && (1 == sscanf(stmp, "%d", &mris->vg.width))) {
        vgvalid++;  // track valid volgeom values found
      }
      stmp = gifti_get_meta_value(&coords->meta, "VolGeomHeight");
      if (stmp && (1 == sscanf(stmp, "%d", &mris->vg.height))) {
        vgvalid++;
      }
      stmp = gifti_get_meta_value(&coords->meta, "VolGeomDepth");
      if (stmp && (1 == sscanf(stmp, "%d", &mris->vg.depth))) {
        vgvalid++;
      }

      stmp = gifti_get_meta_value(&coords->meta, "VolGeomXsize");
      if (stmp && (1 == sscanf(stmp, "%f", &mris->vg.xsize))) {
        vgvalid++;
      }
      stmp = gifti_get_meta_value(&coords->meta, "VolGeomYsize");
      if (stmp && (1 == sscanf(stmp, "%f", &mris->vg.ysize))) {
        vgvalid++;
      }
      stmp = gifti_get_meta_value(&coords->meta, "VolGeomZsize");
      if (stmp && (1 == sscanf(stmp, "%f", &mris->vg.zsize))) {
        vgvalid++;
      }

      stmp = gifti_get_meta_value(&coords->meta, "VolGeomX_R");
      if (stmp && (1 == sscanf(stmp, "%f", &mris->vg.x_r))) {
        vgvalid++;
      }
      stmp = gifti_get_meta_value(&coords->meta, "VolGeomX_A");
      if (stmp && (1 == sscanf(stmp, "%f", &mris->vg.x_a))) {
        vgvalid++;
      }
      stmp = gifti_get_meta_value(&coords->meta, "VolGeomX_S");
      if (stmp && (1 == sscanf(stmp, "%f", &mris->vg.x_s))) {
        vgvalid++;
      }

      stmp = gifti_get_meta_value(&coords->meta, "VolGeomY_R");
      if (stmp && (1 == sscanf(stmp, "%f", &mris->vg.y_r))) {
        vgvalid++;
      }
      stmp = gifti_get_meta_value(&coords->meta, "VolGeomY_A");
      if (stmp && (1 == sscanf(stmp, "%f", &mris->vg.y_a))) {
        vgvalid++;
      }
      stmp = gifti_get_meta_value(&coords->meta, "VolGeomY_S");
      if (stmp && (1 == sscanf(stmp, "%f", &mris->vg.y_s))) {
        vgvalid++;
      }

      stmp = gifti_get_meta_value(&coords->meta, "VolGeomZ_R");
      if (stmp && (1 == sscanf(stmp, "%f", &mris->vg.z_r))) {
        vgvalid++;
      }
      stmp = gifti_get_meta_value(&coords->meta, "VolGeomZ_A");
      if (stmp && (1 == sscanf(stmp, "%f", &mris->vg.z_a))) {
        vgvalid++;
      }
      stmp = gifti_get_meta_value(&coords->meta, "VolGeomZ_S");
      if (stmp && (1 == sscanf(stmp, "%f", &mris->vg.z_s))) {
        vgvalid++;
      }

      stmp = gifti_get_meta_value(&coords->meta, "VolGeomC_R");
      if (stmp && (1 == sscanf(stmp, "%f", &mris->vg.c_r))) {
        vgvalid++;
      }
      stmp = gifti_get_meta_value(&coords->meta, "VolGeomC_A");
      if (stmp && (1 == sscanf(stmp, "%f", &mris->vg.c_a))) {
        vgvalid++;
      }
      stmp = gifti_get_meta_value(&coords->meta, "VolGeomC_S");
      if (stmp && (1 == sscanf(stmp, "%f", &mris->vg.c_s))) {
        vgvalid++;
      }
      stmp = gifti_get_meta_value(&coords->meta, "VolGeomFname");
      if (stmp) {
        vgvalid++;
        memcpy(mris->vg.fname, stmp, sizeof(mris->vg.fname));
      }


      // we got all the values
      if (vgvalid == NUM_VOLGEOM_META) {
        mris->vg.valid = 1;  // finally we can say its valid data
      }

      stmp = gifti_get_meta_value(&coords->meta, "SurfaceCenterX");
      if (stmp && (1 == sscanf(stmp, "%f", &mris->xctr))) {
        vgvalid++;
      }
      stmp = gifti_get_meta_value(&coords->meta, "SurfaceCenterY");
      if (stmp && (1 == sscanf(stmp, "%f", &mris->yctr))) {
        vgvalid++;
      }
      stmp = gifti_get_meta_value(&coords->meta, "SurfaceCenterZ");
      if (stmp && (1 == sscanf(stmp, "%f", &mris->zctr))) {
        vgvalid++;
      }
    }

    /* retrieve TAG_GROUP_AVG_SURFACE_AREA info */
    {
      char *group_avg_surface_area = gifti_get_meta_value(&coords->meta, "TAG_GROUP_AVG_SURFACE_AREA");
      if (group_avg_surface_area)
        sscanf(group_avg_surface_area, "%f", &mris->group_avg_surface_area);
    }

    /* retrieve TAG_CMDLINE info */
    {
      char *ncmds = gifti_get_meta_value(&coords->meta, "NUM_TAG_CMDLINE");
      int numcmds = 0;
      if (ncmds)
        sscanf(ncmds, "%d", &numcmds);

      if (numcmds > MAX_CMDS)
      {
        printf("[WARN] mrisReadGIFTIdanum():  too many commands (%d) in file. Only last %d will be saved!\n", numcmds, MAX_CMDS);
      }

      int toskip = (numcmds > MAX_CMDS) ? (numcmds - MAX_CMDS) : 0;
      if (toskip)
        numcmds = MAX_CMDS;

      while (toskip)
        gifti_get_meta_value(&coords->meta, "TAG_CMDLINE");

      for (int ncmd = 0; ncmd < numcmds; ncmd++)
      {
        char tag[24] = {'\0'};
        sprintf(tag, "TAG_CMDLINE#%d", ncmd);

        char *cmdline = gifti_get_meta_value(&coords->meta, tag);
        if (cmdline == NULL)
	{
          printf("[ERROR] TAG_CMDLINE out of sync! No value found for %s\n", tag);
          break;
	}

        // mris->ncmds will be increased
        MRISaddCommandLine(mris, cmdline);
      }
    }

    /* Copy in the vertices. */
    int vertex_index;
    for (vertex_index = 0; vertex_index < num_vertices; vertex_index++) {
      mris->vertices_topology[vertex_index].num = 0;
      float x = (float)gifti_get_DA_value_2D(coords, vertex_index, 0);
      float y = (float)gifti_get_DA_value_2D(coords, vertex_index, 1);
      float z = (float)gifti_get_DA_value_2D(coords, vertex_index, 2);
      MRISsetXYZ(mris,vertex_index,x,y,z);
      mris->vertices[vertex_index].origarea = -1;
    }
    mrisComputeSurfaceDimensions(mris);
    
    /* Copy in the faces. */
    int face_index;
    for (face_index = 0; face_index < num_faces; face_index++) {
      int face_vertex_index;
      for (face_vertex_index = 0; face_vertex_index < VERTICES_PER_FACE; face_vertex_index++) {
        vertex_index = (int)gifti_get_DA_value_2D(faces, face_index, face_vertex_index);
        mris->faces[face_index].v[face_vertex_index] = vertex_index;
        mris->vertices_topology[vertex_index].num++;
      }
    }
    // each vertex has a face list (faster than face list in some operations)
    for (vertex_index = 0; vertex_index < num_vertices; vertex_index++) {
      mris->vertices_topology[vertex_index].f = (int   *)calloc(mris->vertices_topology[vertex_index].num, sizeof(int));
      mris->vertices_topology[vertex_index].n = (uchar *)calloc(mris->vertices_topology[vertex_index].num, sizeof(uchar));
      mris->vertices_topology[vertex_index].num = 0;  // this gets re-calc'd next...
    }
    for (face_index = 0; face_index < mris->nfaces; face_index++) {
      FACE *face = &mris->faces[face_index];
      int n;
      for (n = 0; n < VERTICES_PER_FACE; n++)
        mris->vertices_topology[face->v[n]].f[mris->vertices_topology[face->v[n]].num++] =
            face_index;  // note that .num is auto-incremented!
    }
    for (vertex_index = 0; vertex_index < num_vertices; vertex_index++) {
      int n, m;
      for (n = 0; n < mris->vertices_topology[vertex_index].num; n++) {
        for (m = 0; m < VERTICES_PER_FACE; m++) {
          if (mris->faces[mris->vertices_topology[vertex_index].f[n]].v[m] == vertex_index) {
            mris->vertices_topology[vertex_index].n[n] = m;
          }
        }
      }
    }

    mrisCompleteTopology(mris);
    
    // check-for and read coordsys struct for talairach xform
    if (coords->coordsys && (coords->numCS > 0)) {
      int idx;
      for (idx = 0; idx < coords->numCS; idx++) {
        if (0 == strcmp(coords->coordsys[idx]->dataspace, "NIFTI_XFORM_SCANNER_ANAT")) {
          mris->useRealRAS = 1;
	}
        else if (0 == strcmp(coords->coordsys[idx]->dataspace, "NIFTI_XFORM_UNKNOWN")) {
          if (0 == strcmp(coords->coordsys[idx]->xformspace, "NIFTI_XFORM_TALAIRACH")) {
            int r, c;
            mris->SRASToTalSRAS_ = MatrixAlloc(4, 4, MATRIX_REAL);
            for (r = 1; r <= 4; r++) {
              for (c = 1; c <= 4; c++) {
                mris->SRASToTalSRAS_->rptr[r][c] = coords->coordsys[idx]->xform[r - 1][c - 1];
              }
            }
          }
        }
      }
    }

    /* other data structure essentials, namely:
     *  mrisComputeVertexDistances(mris);
     *  mrisReadTransform(mris, fname) ;
     *  mris->radius = MRISaverageRadius(mris) ;
     *  MRIScomputeMetricProperties(mris) ;
     *  MRISstoreCurrentPositions(mris) ;
     */
    MRIScomputeNormals(mris);
    UpdateMRIS(mris, fname);
  }
  // completed parsing of coordinate and face data

  // sanity-check, we ought to have an mris struct (either passed-in as a
  // parameter, or created when we found coord and face data)
  if (NULL == mris) {
    fprintf(stderr, "mriseadGIFTIfile: mris is NULL! found when parsing file %s\n", fname);
    gifti_free_image(image);
    return NULL;
  }

  /*
   * and dont forget to store the colortable (if one was found)
   */
  if (ct) {
    mris->ct = ct;
    // sanity-check
    int numEntries = 0;
    CTABgetNumberOfValidEntries(mris->ct, &numEntries);
    if (numEntries != image->labeltable.length) {
      fprintf(
          stderr, "mrisReadGIFTIdanum: ct_entries:%d != labeltable_entries:%d\n", numEntries, image->labeltable.length);
      gifti_free_image(image);
      return NULL;
    }
  }

  /*
   * Now re-parse the DataArrays looking for all the other data type (except
   * coordinate and face data arrays) and fill-in mris structure as needed.
   */
  int nOverlay = 0;                 // track if multiple shape, stat/val data arrays exits
  MRI *overlayMRI = NULL;           // allocated for each overlay
  OverlayInfoStruct statInfo;
  int nStatIntentFrame = 0;         // track number of MRI frames belonging to same stat

  int found_curv_data = 0;          // track if multiple shape data arrays exist
  int found_statval_data = 0;       // track if multiple stat/val data arrays exist
  int prev_stat_intent = -1;
  giiDataArray *node_index = NULL;  // support for sparse data storage
  long long num_index_nodes = 0;    // support for sparse data storage
  int startDAnum = 0;
  int endDAnum = image->numDA;
  if (daNum != -1)  // support for extracting one particular data array
  {
    startDAnum = daNum;
    endDAnum = daNum + 1;
  }
  for (numDA = startDAnum; numDA < endDAnum; numDA++) {
    giiDataArray *darray = image->darray[numDA];

    // did these already
    if ((darray->intent == NIFTI_INTENT_POINTSET) || (darray->intent == NIFTI_INTENT_TRIANGLE)) {
      continue;
    }

    /* support for sparse data storage.  this array contains a list of node
       numbers and it should be the first data array in the file. The remaining
       data arrays in the file that contain data assigned to nodes must contain
       the same number of elements as the NIFTI_INTENT_NODE_INDEX array. */
    if (darray->intent == NIFTI_INTENT_NODE_INDEX) {
      if (numDA != 0) {
        fprintf(stderr,
                "mrisReadGIFTIdanum: NODE_INDEX data array found but its not the "
                "first data array in file %s\n",
                fname);
        gifti_free_image(image);
        return NULL;
      }
      long long num_cols = 0;

      if (darray->ind_ord == GIFTI_IND_ORD_ROW_MAJOR) {
        gifti_DA_rows_cols(darray, &num_index_nodes, &num_cols);
      }
      else {
        gifti_DA_rows_cols(darray, &num_cols, &num_index_nodes);
      }

      if (num_index_nodes <= 0 || num_index_nodes > mris->nvertices || num_cols > 1) {
        fprintf(stderr,
                "mrisReadGIFTIdanum: malformed NODE_INDEX data array in file %s: "
                "num_index_nodes=%d num_cols=%d max nvertices=%d, num_cols>1\n",
                fname,
                (int)num_index_nodes,
                (int)num_cols,
                mris->nvertices);
        gifti_free_image(image);
        return NULL;
      }
      // else good to do, so store this node index info
      node_index = darray;
      continue;
    }
    else {
      /* Check the number of vertices, so we dont trounce the mris struct */
      long long num_vertices = 0;
      long long num_cols = 0;
      long long expected_num_cols = 1;
      if (darray->intent == NIFTI_INTENT_VECTOR) {
        expected_num_cols = 3;
      }
      else if (darray->intent == NIFTI_INTENT_RGB_VECTOR) {
        expected_num_cols = 3;
      }
      else if (darray->intent == NIFTI_INTENT_RGBA_VECTOR) {
        expected_num_cols = 4;
      }
      else if (darray->intent == NIFTI_INTENT_GENMATRIX) {
        expected_num_cols = 9;
      }

      if (darray->ind_ord == GIFTI_IND_ORD_ROW_MAJOR) {
        gifti_DA_rows_cols(darray, &num_vertices, &num_cols);
      }
      else {
        gifti_DA_rows_cols(darray, &num_cols, &num_vertices);
      }

      if (num_vertices <= 0 || num_vertices != mris->nvertices || num_cols > expected_num_cols) {
        fprintf(stderr,
                "mrisReadGIFTIdanum: malformed data array [%d] in file %s: "
                "num_vertices=%d num_cols=%d expected nvertices=%d, num_cols=%d\n",
                numDA,
                fname,
                (int)num_vertices,
                (int)num_cols,
                mris->nvertices,
                (int)expected_num_cols);
        gifti_free_image(image);
        return NULL;
      }
    }

    /*
     * parse each intent type
     */
    if (darray->intent == NIFTI_INTENT_SHAPE) {
      // 'shape' data goes in our 'curv' data element of mris
      if (found_curv_data && poverlayinfo == NULL) {
        fprintf(stderr,
                "WARNING: a prior data array of shape data has already "
                "been read!  Skipping data in array #%d in file %s\n",
                numDA,
                fname);
      }
      else {
        found_curv_data++;

        if (poverlayinfo != NULL)
	{
          fprintf(stderr, "INFO: %s data in array #%d in file %s saved as MRI\n",
                          gifti_intent_to_string(darray->intent), numDA, fname);

          std::vector<int> shape{mris->nvertices, 1, 1, 1};
          overlayMRI = new MRI(shape, MRI_FLOAT);
	}

        if (node_index)  // sparse data storage
        {
          int nindex;
          for (nindex = 0; nindex < num_index_nodes; nindex++) {
            int vno = gifti_get_DA_value_2D(node_index, nindex, 0);
            if (mris->vertices[vno].ripflag) {
              continue;
            }
            mris->vertices[vno].curv = (float)gifti_get_DA_value_2D(darray, nindex, 0);
            if (overlayMRI != NULL)
              MRIsetVoxVal(overlayMRI, vno, 0, 0, 0, mris->vertices[vno].curv);
          }
        }
        else  // regular indexing
        {
          int vno;
          for (vno = 0; vno < mris->nvertices; vno++) {
            if (mris->vertices[vno].ripflag) {
              continue;
            }
            mris->vertices[vno].curv = (float)gifti_get_DA_value_2D(darray, vno, 0);
            if (overlayMRI != NULL)
              MRIsetVoxVal(overlayMRI, vno, 0, 0, 0, mris->vertices[vno].curv);
          }
        }

        if (poverlayinfo == NULL)
          continue;

        
        // save SHAPE information
        char *shapename_base = NULL;
        char *metashapename = gifti_get_meta_value(&darray->meta, "Name");
        if (metashapename != NULL)
	{
          char shapename[1024] = {'\0'};
          // make a copy because basename() may modify the contents of path passed into it
          memcpy(shapename, metashapename, sizeof(shapename));
          shapename_base = basename(shapename);

          // remove the .gii extension
          char *tmpptr = strrchr(shapename_base, '.');
          if (tmpptr != NULL && strcmp(tmpptr+1, "gii") == 0)
            *tmpptr = '\0';
	}

        OverlayInfoStruct overlayInfo;
        overlayInfo.__foverlay = NULL;
        if (shapename_base != NULL)
	{
          overlayInfo.__foverlay = new char[strlen(shapename_base)+1];
          strcpy(overlayInfo.__foverlay, shapename_base);
	}
        overlayInfo.__type = FS_MRISURFOVERLAY_SHAPE;
        overlayInfo.__giftiIntent = NIFTI_INTENT_SHAPE;
        char *metashapedatatype = gifti_get_meta_value(&darray->meta, "ShapeDataType");
        if (metashapedatatype != NULL)
          strcpy(overlayInfo.__shapedatatype, metashapedatatype);
        overlayInfo.__format = GIFTI_FILE;
        overlayInfo.__overlaymri = overlayMRI;

        (*poverlayinfo).push_back(overlayInfo);

        nOverlay++;
      }
    } // NIFTI_INTENT_SHAPE
    else if (darray->intent == NIFTI_INTENT_LABEL) {
      // 'label' data goes into the 'annotation' data element of mris
      if ((NULL == mris->ct) || (NULL == ct))  // sanity-check
      {
        fprintf(stderr, "mrisReadGIFTIdanum: NULL colortable\n");
        gifti_free_image(image);
        return NULL;
      }
      unsigned int *label_data = (unsigned int *)darray->data;
      int nindex = 0;    // index into node_index (if sparse data storage is used)
      int da_index = 0;  // index into the data array at hand
      int vno = 0;       // index into the mris struct (vertex number)
      while (vno < mris->nvertices) {
        if (node_index)  // sparse data storage support
        {
          vno = gifti_get_DA_value_2D(node_index, nindex, 0);
          da_index = nindex;
        }
        else  // regular indexing
        {
          da_index = vno;
        }

        if (mris->vertices[vno].ripflag) {
          continue;
        }

	/*
         * The value at each node is a key into the LabelTable. 
         * Nodes that are “unassigned” should be assigned a label whose color 
         * components’ alpha value is zero. An RGBA color with an alpha value of 
         * zero is fully transparent. Also note that the keys of the Labels are not 
         * necessarily sequential.
	 */
        int table_key = *(label_data + da_index);
        int table_index = 0;
        for (table_index = 0; table_index < image->labeltable.length; table_index++) {
          if (table_key == image->labeltable.key[table_index]) {
            // found the label key for this node
            break;
          }
        }

        int annotation = 0;  // default to no label found
        if ((table_index < image->labeltable.length) && (table_index >= 0)) {
          // table_index pass the tests, table_key (label key) is valid.
          // invalid table_key is getting the default annotation = 0
          if (Gdiag & DIAG_SHOW)
            printf("mrisReadGIFTIdanum(): vno: %d, tkey: %d, tidx: %d, name: %s\n",
		   vno, table_key, table_index, ct->entries[table_key]->name);
          annotation = CTABrgb2Annotation(
              ct->entries[table_key]->ri, ct->entries[table_key]->gi, ct->entries[table_key]->bi);
        }
        mris->vertices[vno].annotation = annotation;

#if 0   // the check below will fail because not every node is assigned 
        // cross-check:
        int index = -1;
        int result = CTABfindAnnotation(mris->ct, mris->vertices[vno].annotation, &index);
        if ((result != NO_ERROR) || (index < 0) || (index >= maxkey)) {
          fprintf(stderr,
                  "mrisReadGIFTIdanum: label node data not found in colortable! "
                  "vno: %d, annot: %8.8X, index: %d\n",
                  vno,
                  mris->vertices[vno].annotation,
                  index);
          gifti_free_image(image);
          return NULL;
        }
#endif

        if (node_index)  // sparse data storage support
        {
          if (++nindex >= num_index_nodes) {
            break;
          }
        }
        else  // regular indexing
        {
          vno++;
        }
      }
    } // NIFTI_INTENT_LABEL
    else if (darray->intent == NIFTI_INTENT_VECTOR) {
      // 'vector' data goes in our 'dx,dy,dz' data element of mris
      int vno;
      for (vno = 0; vno < mris->nvertices; vno++) {
        if (mris->vertices[vno].ripflag) {
          continue;
        }
        mris->vertices[vno].dx = (float)gifti_get_DA_value_2D(darray, vno, 0);
        mris->vertices[vno].dy = (float)gifti_get_DA_value_2D(darray, vno, 1);
        mris->vertices[vno].dz = (float)gifti_get_DA_value_2D(darray, vno, 2);
      }
    } // NIFTI_INTENT_VECTOR
    else if ((darray->intent == NIFTI_INTENT_RGB_VECTOR) || (darray->intent == NIFTI_INTENT_RGBA_VECTOR)) {
      std::map<int, float*> unique_annot_map;

      // 'rgba' data goes in our 'annotation' data element of mris
      for (int vno = 0; vno < mris->nvertices; vno++) {
        if (mris->vertices[vno].ripflag) {
          continue;
        }
        int r, g, b;

        float red = (float)gifti_get_DA_value_2D(darray, vno, 0);
        float green = (float)gifti_get_DA_value_2D(darray, vno, 1);
        float blue = (float)gifti_get_DA_value_2D(darray, vno, 2);

        if (red > 1) {
          r = (int)red;
        }
        else {
          r = (int)floor(red * 256);
        }
        if (r > 255) {
          r = 255;
        }
        if (green > 1) {
          g = (int)green;
        }
        else {
          g = (int)floor(green * 256);
        }
        if (g > 255) {
          g = 255;
        }
        if (blue > 1) {
          b = (int)blue;
        }
        else {
          b = (int)floor(blue * 256);
        }
        if (b > 255) {
          b = 255;
        }

        MRISRGBToAnnot(r, g, b, mris->vertices[vno].annotation);
        std::map<int, float*>::iterator it = unique_annot_map.find(mris->vertices[vno].annotation);
        if (it == unique_annot_map.end())
	{
          // map key/value pair: annotation/(r, g, b)
          float *rgb = new float[3];
          rgb[0] = red; rgb[1] = green; rgb[2] = blue;
          unique_annot_map[mris->vertices[vno].annotation] = rgb;
	}
      }

      // make a color lookup table with unique annotation from dataarray
      mris->ct = makeColorTable(unique_annot_map, fname);
    } // NIFTI_INTENT_RGB_VECTOR || NIFTI_INTENT_RGBA_VECTOR
    else if (darray->intent == NIFTI_INTENT_GENMATRIX) {
      fprintf(stderr,
              "WARNING: ignoring unsupported data array NIFTI_INTENT_GENMATRIX"
              " in file %s\n",
              fname);
    } // NIFTI_INTENT_GENMATRIX
    else {
      // 'statistics' and all other kinds of data we'll put in both our
      // 'stat' and 'val' data elements of the mris structure
      if (found_statval_data && poverlayinfo == NULL) {
        fprintf(stderr,
                "WARNING: a prior data array of stat/val data has already "
                "been read!  Skipping data in array #%d in file %s\n",
                numDA,
                fname);
      }
      else {
        found_statval_data++;

        if (poverlayinfo != NULL)
	{
          fprintf(stderr, "INFO: %s data in array #%d in file %s saved as MRI\n",
                          gifti_intent_to_string(darray->intent), numDA, fname);

          std::vector<int> shape{mris->nvertices, 1, 1, 1};
          overlayMRI = new MRI(shape, MRI_FLOAT);
	}

        if (node_index)  // sparse data storage
        {
          int nindex;
          for (nindex = 0; nindex < num_index_nodes; nindex++) {
            int vno = gifti_get_DA_value_2D(node_index, nindex, 0);
            if (mris->vertices[vno].ripflag) {
              continue;
            }
            mris->vertices[vno].val = (float)gifti_get_DA_value_2D(darray, nindex, 0);
            mris->vertices[vno].stat = (float)gifti_get_DA_value_2D(darray, nindex, 0);
            if (overlayMRI != NULL)
              MRIsetVoxVal(overlayMRI, vno, 0, 0, nStatIntentFrame, mris->vertices[vno].stat);
          }
        }
        else  // regular indexing
        {
          int vno;
          for (vno = 0; vno < mris->nvertices; vno++) {
            if (mris->vertices[vno].ripflag) {
              continue;
            }
            mris->vertices[vno].val = (float)gifti_get_DA_value_2D(darray, vno, 0);
            mris->vertices[vno].stat = (float)gifti_get_DA_value_2D(darray, vno, 0);
            if (overlayMRI != NULL)
              MRIsetVoxVal(overlayMRI, vno, 0, nStatIntentFrame, 0, mris->vertices[vno].stat);
          }
        }

        if (poverlayinfo == NULL)
          continue;

        if (prev_stat_intent != darray->intent)
	{
          nStatIntentFrame = 0;   // reset

          // save STAT information
          if (prev_stat_intent != -1)
            (*poverlayinfo).push_back(statInfo);

          statInfo.__foverlay = new char[strlen(fname)+1];
          strcpy(statInfo.__foverlay, fname);
          statInfo.__type = FS_MRISURFOVERLAY_STATS;
          statInfo.__giftiIntent = darray->intent;
          memset(statInfo.__shapedatatype, 0, sizeof(statInfo.__shapedatatype));
          statInfo.__format = GIFTI_FILE;
          statInfo.__overlaymri = overlayMRI;

          prev_stat_intent = darray->intent;
	}
        else
          nStatIntentFrame++;

        nOverlay++;
      }
    } // 'statistics' and all other kinds of data
  } // for each DAnum

  /*
   * And we're done.
   */
  gifti_free_image(image);

  return mris;
} // end of mrisReadGIFTIdanum()


 /*-----------------------------------------------------------
   mrisReadGIFTIfile()
                 reads a GIFTI file, putting vertices,
                 and faces into an MRIS_SURFACE structure,
                 along with any other data, like labels,
                 colors, curv data, stats or values.
                
   after read, 
     first SHAPE is saved in mris->curv; 
     first <STATS> is saved in mris->val and mris->stat;
     all SHAPE and <STATS> data arrays are saved as multi-frame MRI
  -----------------------------------------------------------*/
MRI_SURFACE *mrisReadGIFTIfile(const char *fname, MRI_SURFACE *mris, std::vector<OverlayInfoStruct> *poverlayinfo)
{
  // default read routine (read all data arrays)
  return mrisReadGIFTIdanum(fname, mris, -1, poverlayinfo);
}

/*-----------------------------------------------------------
  MRISreadGiftiAsMRI() - reads GIFTI functional frames into
  an MRI volume struct, which is a retro-fit usage to store
  multiple frames of data (where in this case, a frame is one
  complete vector of vertices).
  This routine will only read NIFTI_INTENT_TIME_SERIES data
  arrays.
  -----------------------------------------------------------*/
MRI *MRISreadGiftiAsMRI(const char *fname, int read_volume)
{
  /* Attempt to read the file. */
  gifti_image *image = gifti_read_image(fname, 1);
  if (NULL == image) {
    fprintf(stderr, "MRISreadGiftiAsMRI: gifti_read_image() returned NULL\n");
    return NULL;
  }

  // make sure version is recoded before validation
  if (!strcmp(image->version, "1")) {
    free(image->version);
    image->version = strcpyalloc(GIFTI_XML_VERSION);
  }

  /* check for compliance */
  int valid = gifti_valid_gifti_image(image, 1);
  if (valid == 0) {
    fprintf(stderr, "MRISreadGiftiAsMRI: GIFTI file %s is invalid!\n", fname);
    gifti_free_image(image);
    return NULL;
  }

  /* check for overlay data */
  giiDataArray *scalars = NULL;
  int frame_count = 0;
  long long num_vertices = -1;
  long long num_cols = 0;
#define INTENT_CODE_MAX_IDX 4
  int intent_code[INTENT_CODE_MAX_IDX] = {
      NIFTI_INTENT_TIME_SERIES, NIFTI_INTENT_SHAPE, NIFTI_INTENT_NONE, NIFTI_INTENT_NORMAL};
  int intent_code_idx = 0;
  // search all DAs for time series, then shape, then none, then normal.
  // if time series found, check all DAs to make sure all the same size.
  for (intent_code_idx = 0; intent_code_idx < INTENT_CODE_MAX_IDX; intent_code_idx++) {
    int da_num = 0;
    do {
      scalars = gifti_find_DA(image, intent_code[intent_code_idx], da_num);
      if (NULL == scalars) {
        if (++da_num >= image->numDA) {
          break;
        }
        else {
          continue;
        }
      }
      frame_count++;
      long long nvertices = 0, ncols = 0;

      if (scalars->ind_ord == GIFTI_IND_ORD_ROW_MAJOR) {
        gifti_DA_rows_cols(scalars, &nvertices, &ncols);
      }
      else {
        gifti_DA_rows_cols(scalars, &ncols, &nvertices);
      }

      if (num_vertices == -1) {
        num_vertices = nvertices;
        num_cols = ncols;
      }
      else {
        if (num_vertices <= 0 || num_vertices != nvertices || ncols != 1) {
          fprintf(stderr,
                  "MRISreadGiftiAsMRI: malformed time-series data array in file "
                  "%s: nvertices=%d ncols=%d expected num_vertices=%d\n",
                  fname,
                  (int)nvertices,
                  (int)num_cols,
                  (int)num_vertices);
          gifti_free_image(image);
          return NULL;
        }
      }
      if (++da_num >= image->numDA) {
        break;
      }
      // ??? why the check ???
      if ((intent_code[intent_code_idx] != NIFTI_INTENT_TIME_SERIES) &&
          (intent_code[intent_code_idx] != NIFTI_INTENT_SHAPE) && (intent_code[intent_code_idx] != NIFTI_INTENT_NONE) &&
          (intent_code[intent_code_idx] != NIFTI_INTENT_NORMAL)) {
        break;
      }
    } while (scalars);

    if (scalars) {
      break;  // found some data, no need to check other intents
    }
  }

  if (frame_count == 0) {
    fprintf(stderr, "MRISreadGiftiAsMRI: no overlay data found in file %s\n", fname);
    gifti_free_image(image);
    return NULL;
  }

  /* if we don't need to read the volume, just return a header */
  MRI *mri;
  if (!read_volume) {
    mri = MRIallocHeader(num_vertices, 1, 1, MRI_FLOAT, frame_count);
    mri->nframes = frame_count;
    // not sure this is the best way to do this (dng, 4/4/17)
    if (image->numDA > 0) {
      char *stmp = gifti_get_meta_value(&image->darray[0]->meta, "TimeStep");
      if (stmp) sscanf(stmp, "%f", &mri->tr);
    }
    return (mri);
  }

  /* Copy in each scalar frame to 'volume' frame. */
  mri = MRIallocSequence(num_vertices, 1, 1, MRI_FLOAT, frame_count);
  frame_count = 0;
  int da_num;
  // ??? intent_code_idx = 4 here ???
  // ??? need outer loop for through INTENT_CODE_MAX_IDX ???
  for (da_num = 0; da_num < image->numDA; da_num++) {
    scalars = gifti_find_DA(image, intent_code[intent_code_idx], da_num);
    if (NULL == scalars) {
      continue;
    }
    int vno;
    for (vno = 0; vno < num_vertices; vno++) {
      float val = (float)gifti_get_DA_value_2D(scalars, vno, 0);
      MRIsetVoxVal(mri, vno, 0, 0, frame_count, val);
    }
    // printf("frame #%d\n",frame_count);
    frame_count++;
  }

  // not sure this is the best way to do this (dng, 4/4/17)
  if (image->numDA > 0) {
    char *stmp = gifti_get_meta_value(&image->darray[0]->meta, "TimeStep");
    if (stmp) sscanf(stmp, "%f", &mri->tr);
  }

  /* And we're done. */
  gifti_free_image(image);
  return (mri);
} // end of MRISreadGiftiAsMRI()

/*
 * insert username and current date into meta data
 */
static void insertCommonMetaData(giiMetaData *md)
{
  struct passwd *pw = getpwuid(geteuid());
  if ((pw != NULL) && (pw->pw_name != NULL)) {
    gifti_add_to_meta(md, "UserName", pw->pw_name, 1);
  }

  gifti_add_to_meta(md, "Date", currentDateTime().c_str(), 1);
}

/*-------------------------------------------------------------------------
  Parameters:    mris - MRIS_SURFACE structure
                 intent_code - NIFTI_INTENT_* indicating
                   type of data in mris to write
                 out_fname - output file name of GIFTI file
                 curv_fname - if intent code is NIFTI_INTENT_SHAPE,
                   then this file contains the data to load into .curv field
                   prior to writing shape data to gifti output

  Returns value: 0 if passed, else error code

  Description:   writes a GIFTI file, if intent code is...
                 NIFTI_INTENT_POINTSET or _TRIANGLE, then write vertices and
                 face data
                 if NIFTI_INTENT_LABEL, then write LabelTable and annotations
                 if NIFTI_INTENT_<statistic>, then write .stats data
  ------------------------------------------------------------------------*/
int MRISwriteGIFTI(MRIS *mris, int intent_code, const char *out_fname, const char *curv_fname)
{
  if (NULL == mris || NULL == out_fname) {
    fprintf(stderr, "MRISwriteGIFTI: invalid parameter, surf or fname is NULL\n");
    return ERROR_BADPARM;
  }

  if (intent_code == NIFTI_INTENT_SHAPE && NULL == curv_fname) {
    printf("MRISwriteGIFTI: invalid parameter: curv_fname is NULL %s:%d\n", __FILE__, __LINE__);
    return ERROR_BADPARM;
  }

  gifti_image *image = (gifti_image *)calloc(1, sizeof(gifti_image));
  if (NULL == image) {
    fprintf(stderr, "MRISwriteGIFTI: couldn't allocate image\n");
    return ERROR_NOMEMORY;
  }
  image->version = strcpyalloc(GIFTI_XML_VERSION);

  insertCommonMetaData(&image->meta);
  if (strlen(mris->subject_name)) {
    gifti_add_to_meta(&image->meta, "SubjectID", mris->subject_name, 1);
  }

  int error = MRISwriteGIFTIIntent(mris, intent_code, image, out_fname, curv_fname);
  if (error != NO_ERROR)
      return error;

  // make sure version is recoded before validation
  if (!strcmp(image->version, "1")) {
    free(image->version);
    image->version = strcpyalloc(GIFTI_XML_VERSION);
  }

  /* check for compliance */
  int valid = gifti_valid_gifti_image(image, 1);
  if (valid == 0) {
    fprintf(stderr, "MRISwriteGIFTI: GIFTI file %s is invalid!\n", out_fname);
    gifti_free_image(image);
    return ERROR_BADFILE;
  }

  /* Write the file. */
  if (gifti_write_image(image, out_fname, 1)) {
    fprintf(stderr, "MRISwriteGIFTI: couldn't write image\n");
    gifti_free_image(image);
    return ERROR_BADFILE;
  }

  gifti_free_image(image);

  return ERROR_NONE;
} // end of MRISwriteGIFTI()

/*-----------------------------------------------------------
  Parameters:    MRI structure (surface-encoded volume),
                 output file name of GIFTI file,

  Returns value: 0 if passed, else error code

  Description:   writes a GIFTI file containing functional or
                 timeseries data
  -----------------------------------------------------------*/
int mriWriteGifti(MRI *mri, const char *out_fname)
{
  if (NULL == mri || NULL == out_fname) {
    fprintf(stderr, "mriWriteGifti: invalid input parameters\n");
    return ERROR_BADPARM;
  }

  gifti_image *image = (gifti_image *)calloc(1, sizeof(gifti_image));
  if (NULL == image) {
    fprintf(stderr, "mriWriteGifti: couldn't allocate gifti_image\n");
    return ERROR_NOMEMORY;
  }
  image->version = strcpyalloc(GIFTI_XML_VERSION);

  /* include some metadata describing this thing */
  insertCommonMetaData(&image->meta);
  /* -------------------------------------------------------
   * One DataArray for each 'frame' in the 'volume' data
   */
  int frame;
  for (frame = 0; frame < mri->nframes; frame++) {
    giiDataArray *scalars = gifti_alloc_and_add_darray(image);
    if (NULL == scalars) {
      fprintf(stderr, "mriWriteGifti: couldn't allocate giiDataArray\n");
      gifti_free_image(image);
      return ERROR_NOMEMORY;
    }

    /* Set its attributes. */
    scalars->intent = NIFTI_INTENT_NONE;
    if (mri->nframes > 1) {
      scalars->intent = NIFTI_INTENT_TIME_SERIES;
    }
    if (scalars->intent == NIFTI_INTENT_TIME_SERIES) {
      // add TR (repetition time) to metadata:
      char buf[STRLEN];
      sprintf(buf, "%f", mri->tr);
      gifti_add_to_meta(&scalars->meta, "TimeStep", buf, 1);
    }
    scalars->datatype = NIFTI_TYPE_FLOAT32;
    scalars->ind_ord = GIFTI_IND_ORD_ROW_MAJOR;
    scalars->num_dim = 1;
    scalars->dims[0] = mri->width;
    scalars->dims[1] = 0;
    scalars->encoding = GIFTI_ENCODING_B64GZ;  // data stored in gzip'd base64
#if (BYTE_ORDER == LITTLE_ENDIAN)
    scalars->endian = GIFTI_ENDIAN_LITTLE;
#else
    scalars->endian = GIFTI_ENDIAN_BIG;
#endif
    scalars->coordsys = NULL;
    scalars->nvals = gifti_darray_nvals(scalars);
    gifti_datatype_sizes(scalars->datatype, &scalars->nbyper, NULL);

    /* Allocate the data array. */
    scalars->data = NULL;
    scalars->data = (void *)calloc(scalars->nvals, scalars->nbyper);
    if (NULL == scalars->data) {
      fprintf(stderr,
              "mriWriteGifti: couldn't allocate scalars data of "
              "length %d, element size %d\n",
              (int)scalars->nvals,
              scalars->nbyper);
      gifti_free_image(image);
      return ERROR_NOMEMORY;
    }

    /* Copy in all our data. */
    int scalar_index;
    for (scalar_index = 0; scalar_index < mri->width; scalar_index++) {
      float val = MRIgetVoxVal(mri, scalar_index, 0, 0, frame);
      gifti_set_DA_value_2D(scalars, scalar_index, 0, val);
    }

    // next frame
  }

  // make sure version is recoded before validation
  if (!strcmp(image->version, "1")) {
    free(image->version);
    image->version = strcpyalloc(GIFTI_XML_VERSION);
  }

  /* check for compliance */
  int valid = gifti_valid_gifti_image(image, 1);
  if (valid == 0) {
    fprintf(stderr, "mriWriteGifti: GIFTI file %s is invalid!\n", out_fname);
    gifti_free_image(image);
    return ERROR_BADFILE;
  }

  /* Write the file. */
  if (gifti_write_image(image, out_fname, 1)) {
    fprintf(stderr, "mriWriteGifti: couldn't write image\n");
    gifti_free_image(image);
    return ERROR_BADFILE;
  }

  gifti_free_image(image);

  return ERROR_NONE;
} // end of mriWriteGifti()



// write single intent
int MRISwriteGIFTIIntent(MRIS *mris, int intent_code, gifti_image *image, const char *out_fname, const char *curv_fname)
{
  /* -------------------------------------------------------
   * Surface file
   */
  if (intent_code == NIFTI_INTENT_POINTSET || intent_code == NIFTI_INTENT_TRIANGLE) {
    int error = MRISwriteGIFTISurface(mris, image, out_fname); 
    if (error != NO_ERROR)
      return error;  
  }  // end of if NIFTI_INTENT_POINTSET or NIFTI_INTENT_TRIANGLE

  /* -------------------------------------------------------
   * Shape file
   */
  if (intent_code == NIFTI_INTENT_SHAPE) {
    int error = MRISwriteGIFTIShape(mris, image, intent_code, curv_fname);
    if (error != NO_ERROR)
      return error;
  }  // end of if NIFTI_INTENT_SHAPE

  /* -------------------------------------------------------
   * Label file
   */
  if (intent_code == NIFTI_INTENT_LABEL) {
    int error = MRISwriteGIFTILabel(mris, image, intent_code); 
    if (error != NO_ERROR)
      return error;    
  }  // end of if NIFTI_INTENT_LABEL

  /* -------------------------------------------------------
   * RGBA file
   */
  if (intent_code == NIFTI_INTENT_RGBA_VECTOR) {
    int error = MRISwriteGIFTIRGBAVector(mris, image, intent_code); 
    if (error != NO_ERROR)
      return error;    
  }  // end of if NIFTI_INTENT_RGBA_VECTOR

  /* -------------------------------------------------------
   * Statistics file
   */
  if (intent_code == NIFTI_INTENT_CORREL || intent_code == NIFTI_INTENT_TTEST || intent_code == NIFTI_INTENT_FTEST ||
      intent_code == NIFTI_INTENT_ZSCORE || intent_code == NIFTI_INTENT_CHISQ || intent_code == NIFTI_INTENT_BETA ||
      intent_code == NIFTI_INTENT_BINOM || intent_code == NIFTI_INTENT_GAMMA || intent_code == NIFTI_INTENT_POISSON ||
      intent_code == NIFTI_INTENT_NORMAL || intent_code == NIFTI_INTENT_FTEST_NONC ||
      intent_code == NIFTI_INTENT_CHISQ_NONC || intent_code == NIFTI_INTENT_LOGISTIC ||
      intent_code == NIFTI_INTENT_LAPLACE || intent_code == NIFTI_INTENT_UNIFORM ||
      intent_code == NIFTI_INTENT_TTEST_NONC || intent_code == NIFTI_INTENT_WEIBULL ||
      intent_code == NIFTI_INTENT_CHI || intent_code == NIFTI_INTENT_INVGAUSS || intent_code == NIFTI_INTENT_EXTVAL ||
      intent_code == NIFTI_INTENT_PVAL || intent_code == NIFTI_INTENT_LOGPVAL ||
      intent_code == NIFTI_INTENT_LOG10PVAL || intent_code == NIFTI_INTENT_ESTIMATE) {
    int error = MRISwriteGIFTIStats(mris, image, intent_code); 
    if (error != NO_ERROR)
      return error;
  }  // end of if NIFTI_INTENT_<stats>

  return NO_ERROR;
} // end of MRISwriteGIFTIIntent(MRIS *mris, ...) 


/*
 * Shape file
 *       intent_code = NIFTI_INTENT_SHAPE
 */
int MRISwriteGIFTIShape(MRIS *mris, gifti_image *image, int intent_code, const char *curv_fname)
{
    // ??? input mris should have curvature information ???
    /// ??? why read again ??? 
    //if (MRISreadCurvatureFile(mris, curv_fname)) {
    //  fprintf(stderr, "MRISwriteGIFTIShape: couldn't read %s\n", curv_fname);
    //  gifti_free_image(image);
    //  return ERROR_BADFILE;
    //}

    giiDataArray *shape = gifti_alloc_and_add_darray(image);
    if (NULL == shape) {
      fprintf(stderr, "MRISwriteGIFTIShape: couldn't allocate giiDataArray\n");
      gifti_free_image(image);
      return ERROR_NOMEMORY;
    }

    /* Set its attributes. */
    shape->intent = NIFTI_INTENT_SHAPE;
    shape->datatype = NIFTI_TYPE_FLOAT32;
    shape->ind_ord = GIFTI_IND_ORD_ROW_MAJOR;
    shape->num_dim = 1;
    shape->dims[0] = mris->nvertices;
    shape->dims[1] = 0;
    shape->encoding = GIFTI_ENCODING_B64GZ;  // data stored in gzip'd base64
#if (BYTE_ORDER == LITTLE_ENDIAN)
    shape->endian = GIFTI_ENDIAN_LITTLE;
#else
    shape->endian = GIFTI_ENDIAN_BIG;
#endif
    shape->coordsys = NULL;
    shape->nvals = gifti_darray_nvals(shape);
    gifti_datatype_sizes(shape->datatype, &shape->nbyper, NULL);

    /* include some metadata describing this shape */
    gifti_add_to_meta(&shape->meta, "Name", curv_fname, 1);
    const char *meta = NULL;
    if (strstr(curv_fname, ".thickness")) {
      meta = "Thickness";
    }
    if (strstr(curv_fname, ".curv")) {
      meta = "CurvatureRadial";
    }
    if (strstr(curv_fname, ".sulc")) {
      meta = "SulcalDepth";
    }
    if (strstr(curv_fname, ".area")) {
      meta = "Area";
    }
    if (strstr(curv_fname, ".volume")) {
      meta = "Volume";
    }
    if (strstr(curv_fname, ".jacobian")) {
      meta = "Jacobian";
    }
    if (meta) {
      gifti_add_to_meta(&shape->meta, "ShapeDataType", meta, 1);
    }

    /* Allocate the data array. */
    shape->data = NULL;
    shape->data = (void *)calloc(shape->nvals, shape->nbyper);
    if (NULL == shape->data) {
      fprintf(stderr,
              "MRISwriteGIFTIShape: couldn't allocate shape data of "
              "length %d, element size %d\n",
              (int)shape->nvals,
              shape->nbyper);
      gifti_free_image(image);
      return ERROR_NOMEMORY;
    }

    /* Copy in all our data. */
    int vno;
    for (vno = 0; vno < mris->nvertices; vno++) {
      if (mris->vertices[vno].ripflag) {
        continue;
      }
      gifti_set_DA_value_2D(shape, vno, 0, mris->vertices[vno].curv);
    }

    return NO_ERROR;
} // end of MRISwriteGIFTIShape(MRIS *mris, ...)


/*
 * Statistics file
 *       intent_code = NIFTI_INTENT_<stats>
 */
int MRISwriteGIFTIStats(MRIS *mris, gifti_image *image, int intent_code)
{
    giiDataArray *stats = gifti_alloc_and_add_darray(image);
    if (NULL == stats) {
      fprintf(stderr, "MRISwriteGIFTIStats: couldn't allocate giiDataArray\n");
      gifti_free_image(image);
      return ERROR_NOMEMORY;
    }

    /* Set its attributes. */
    stats->intent = intent_code;
    stats->datatype = NIFTI_TYPE_FLOAT32;
    stats->ind_ord = GIFTI_IND_ORD_ROW_MAJOR;
    stats->num_dim = 1;
    stats->dims[0] = mris->nvertices;
    stats->dims[1] = 0;
    stats->encoding = GIFTI_ENCODING_B64GZ;  // data stored in gzip'd base64
#if (BYTE_ORDER == LITTLE_ENDIAN)
    stats->endian = GIFTI_ENDIAN_LITTLE;
#else
    stats->endian = GIFTI_ENDIAN_BIG;
#endif
    stats->coordsys = NULL;
    stats->nvals = gifti_darray_nvals(stats);
    gifti_datatype_sizes(stats->datatype, &stats->nbyper, NULL);

    /* include some metadata describing this thing */
    gifti_add_to_meta(&stats->meta, "Intent_code", gifti_intent_to_string(intent_code), 1);
    if (intent_code == NIFTI_INTENT_UNIFORM) {
      gifti_add_to_meta(&stats->meta, "Intent_p1", "0", 1);  // lower end
      gifti_add_to_meta(&stats->meta, "Intent_p2", "1", 1);  // upper end
    }

    /* Allocate the data array. */
    stats->data = NULL;
    stats->data = (void *)calloc(stats->nvals, stats->nbyper);
    if (NULL == stats->data) {
      fprintf(stderr,
              "MRISwriteGIFTIStats: couldn't allocate stats data of "
              "length %d, element size %d\n",
              (int)stats->nvals,
              stats->nbyper);
      gifti_free_image(image);
      return ERROR_NOMEMORY;
    }

    /* Copy in all our data. */
    int vno;
    for (vno = 0; vno < mris->nvertices; vno++) {
      if (mris->vertices[vno].ripflag) {
        continue;
      }
      gifti_set_DA_value_2D(stats, vno, 0, mris->vertices[vno].stat);
    }

    return NO_ERROR;
} // end of MRISwriteGIFTIStats(MRIS *mris, ...)


/*
 * Label file
 *       intent_code = NIFTI_INTENT_LABEL
 */
int MRISwriteGIFTILabel(MRIS *mris, gifti_image *image, int intent_code)
{
    /*
     * Writes .annot data to a label table data gifti file:
     * puts the freesurfer colortable struct into a LabelTable,
     * and puts the .annotation field from each vertex into a
     * DataArray
     */
    /*
     * LabelTable struct, fill it in with our colortable stuff
     */

    /* We have to run through our table and count our non-null entries.
     * The count will be LabelTable length.
     */
    int num_entries_to_write = 0;
    for (int i = 0; i < mris->ct->nentries; i++)
      if (NULL != mris->ct->entries[i]) num_entries_to_write++;

    if (Gdiag & DIAG_SHOW)
      printf("[DEBUG] MRISwriteGIFTILabel(): mris->ct->entries=%d, to_write=%d\n", mris->ct->nentries, num_entries_to_write);

    giiLabelTable labeltable;
    labeltable.length = num_entries_to_write;  //mris->ct->nentries;
    if (labeltable.length == 0) {
      fprintf(stderr, "MRISwriteGIFTILabel: colortable is empty!\n");
      return ERROR_BADFILE;
    }
    labeltable.key = (int *)calloc(labeltable.length, sizeof(int *));
    labeltable.label = (char **)calloc(labeltable.length, sizeof(char *));
    labeltable.rgba = (float *)calloc(labeltable.length, 4 * sizeof(float *));
    if ((NULL == labeltable.key) || (NULL == labeltable.label) || (NULL == labeltable.rgba)) {
      fprintf(stderr, "MRISwriteGIFTILabel: couldn't allocate giftiLabelTable\n");
      return ERROR_NOMEMORY;
    }
    float *rgba = labeltable.rgba;
    int idx = 0;
    for (int n = 0; n < mris->ct->nentries; n++) {
      // the key could be the freesurfer 'annotation' value, which is
      // supposed to be unique to the FreeSurferColorLUT, but for gifti
      // purposes, it is more intutive and obvious to use the index.
      // also, a display application might choose to interpret the
      // label data at each vertex as indicies rather than keys (which
      // i think ignores the gifti spec, but is reasonable to do so).

      // output only the non-null entries
      if (mris->ct->entries[n] != NULL && strlen(mris->ct->entries[n]->name) != 0) {
        labeltable.key[idx] = n;
        // labeltable.key[idx] = CTABrgb2Annotation(mris->ct->entries[idx]->ri,
        //                                       mris->ct->entries[idx]->gi,
        //                                       mris->ct->entries[idx]->bi);
        // printf("%8.8X\n",labeltable.key[idx]);

        printf("idx=%d, name=%s\n",idx,mris->ct->entries[n]->name);
        labeltable.label[idx] = strcpyalloc(mris->ct->entries[n]->name);

        rgba[0] = mris->ct->entries[n]->rf;
        rgba[1] = mris->ct->entries[n]->gf;
        rgba[2] = mris->ct->entries[n]->bf;
        rgba[3] = 1.0f;

        idx++;      // next label
        rgba += 4;  // next color
      }
#if 0
      else {
        char tmpname[30];
        sprintf(tmpname, "unknown_%d", idx);
        printf("idx=%d, name=NULL, assigned as %s (is the colortable correct?)\n", idx, tmpname);
        labeltable.label[idx] = strcpyalloc(tmpname);
      }

      if (mris->ct->entries[idx] == NULL ||
          (strlen(mris->ct->entries[idx]->name) == 0) ||
          (strcmp(labeltable.label[idx], "unknown") == 0) ||
          (strcmp(labeltable.label[idx], "Unknown") == 0)) {
        // make certain unknown region is completely empty, invisible
        printf("idx=%d, make certain unknown region is completely empty, invisible\n", idx);
        rgba[0] = rgba[1] = rgba[2] = rgba[3] = 0.0f;
      }
      else {
        rgba[0] = mris->ct->entries[idx]->rf;
        rgba[1] = mris->ct->entries[idx]->gf;
        rgba[2] = mris->ct->entries[idx]->bf;
        rgba[3] = 1.0f;
      }
      rgba += 4;  // next color
#endif
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
    // gifti_disp_LabelTable(NULL,&image->labeltable);

    /*
     * Labels array
     */
    giiDataArray *labels = gifti_alloc_and_add_darray(image);
    if (NULL == labels) {
      fprintf(stderr, "MRISwriteGIFTILabel: couldn't allocate giiDataArray\n");
      gifti_free_image(image);
      return ERROR_NOMEMORY;
    }

    /* Set its attributes. */
    labels->intent = NIFTI_INTENT_LABEL;
    labels->datatype = NIFTI_TYPE_INT32;
    labels->ind_ord = GIFTI_IND_ORD_ROW_MAJOR;
    labels->num_dim = 1;
    labels->dims[0] = mris->nvertices;
    labels->dims[1] = 0;
    labels->encoding = GIFTI_ENCODING_B64GZ;  // data stored in gzip'd base64
#if (BYTE_ORDER == LITTLE_ENDIAN)
    labels->endian = GIFTI_ENDIAN_LITTLE;
#else
    labels->endian = GIFTI_ENDIAN_BIG;
#endif
    labels->coordsys = NULL;
    labels->nvals = gifti_darray_nvals(labels);
    gifti_datatype_sizes(labels->datatype, &labels->nbyper, NULL);

    /* include some metadata describing this as a label */
    gifti_add_to_meta(&labels->meta, "Name", "node label", 1);
    if (labeltable.length == 2) {
      // in the special case of a label table consisting of one label
      // (assuming the first label is 'unknown') use this one label as name,
      // for instance in the case of the V1 label
      gifti_add_to_meta(&labels->meta, "Name", labeltable.label[1], 1);
    }

    /* Allocate the data array. */
    labels->data = NULL;
    labels->data = (void *)calloc(labels->nvals, labels->nbyper);
    if (NULL == labels->data) {
      fprintf(stderr,
              "MRISwriteGIFTILabel: couldn't allocate labels data of "
              "length %d, element size %d\n",
              (int)labels->nvals,
              labels->nbyper);
      gifti_free_image(image);
      return ERROR_NOMEMORY;
    }

    /* Copy our 'annotation' data for each vertex (actually an index) */
    unsigned int *label_data = (unsigned int *)labels->data;
    int label_index, theIdx, result;
    for (label_index = 0; label_index < mris->nvertices; label_index++) {
      if (mris->vertices[label_index].ripflag) {
        continue;
      }
      result = CTABfindAnnotation(mris->ct, mris->vertices[label_index].annotation, &theIdx);
      if (result) {
        return ERROR_BADFILE;
      }

      *(label_data + label_index) = theIdx;
      // printf("%8.8X ", *(label_data + label_index));
    }

    return NO_ERROR;
} // end of MRISwriteGIFTILabel()


/*
 * RGBA file
 *       output NIFTI_INTENT_RGBA_VECTOR
 */
int MRISwriteGIFTIRGBAVector(MRIS *mris, gifti_image *image, int intent_code)
{
    /*
     * RGBA
     */
    giiDataArray *rgba = gifti_alloc_and_add_darray(image);
    if (NULL == rgba) {
      fprintf(stderr, "MRISwriteGIFTIRGBAVector(): couldn't allocate giiDataArray\n");
      gifti_free_image(image);
      return ERROR_NOMEMORY;
    }

    /* Set its attributes. */
    rgba->intent = NIFTI_INTENT_RGBA_VECTOR;
    rgba->datatype = NIFTI_TYPE_FLOAT32;
    rgba->ind_ord = GIFTI_IND_ORD_ROW_MAJOR;
    rgba->num_dim = 2;
    rgba->dims[0] = mris->nvertices;        /* In highest first, dim0 = rows */
    rgba->dims[1] = 4;                      /* In highest first, dim1 = cols */
    rgba->encoding = GIFTI_ENCODING_B64GZ;  // data stored in gzip'd base64
#if (BYTE_ORDER == LITTLE_ENDIAN)
    rgba->endian = GIFTI_ENDIAN_LITTLE;
#else
    rgba->endian = GIFTI_ENDIAN_BIG;
#endif

    rgba->coordsys = NULL;             // empty, unless we find something here...
    rgba->nvals = gifti_darray_nvals(rgba);
    gifti_datatype_sizes(rgba->datatype, &rgba->nbyper, NULL);

    /* Allocate the data array. */
    rgba->data = NULL;
    rgba->data = (void *)calloc(rgba->nvals, rgba->nbyper);
    if (NULL == rgba->data) {
      fprintf(stderr,
              "MRISwriteGIFTIRGBAVector(): couldn't allocate rgba data of "
              "length %d, element size %d\n",
              (int)rgba->nvals,
              rgba->nbyper);
      gifti_free_image(image);
      return ERROR_NOMEMORY;
    }

    /* Copy in all our data. */
    int vertex_index;
    for (vertex_index = 0; vertex_index < mris->nvertices; vertex_index++) {
      if (mris->vertices[vertex_index].ripflag)
        continue;

      int annot = mris->vertices[vertex_index].annotation;

      int r, g, b;
      MRISAnnotToRGB(annot, r, g, b);

      gifti_set_DA_value_2D(rgba, vertex_index, 0, (float)r/255.0);
      gifti_set_DA_value_2D(rgba, vertex_index, 1, (float)g/255.0);
      gifti_set_DA_value_2D(rgba, vertex_index, 2, (float)b/255.0);
      gifti_set_DA_value_2D(rgba, vertex_index, 3, 1.0f);
    }

    return NO_ERROR;
} // end of MRISwriteGIFTIRGBAVector()



/*
 * Surface file
 *       output NIFTI_INTENT_POINTSET or NIFTI_INTENT_TRIANGLE
 */
int MRISwriteGIFTISurface(MRIS *mris, gifti_image *image, const char *out_fname)
{
    /*
     * Coordinates
     */
    giiDataArray *coords = gifti_alloc_and_add_darray(image);
    if (NULL == coords) {
      fprintf(stderr, "MRISwriteGIFTISurface: couldn't allocate giiDataArray\n");
      gifti_free_image(image);
      return ERROR_NOMEMORY;
    }

    /* Set its attributes. */
    coords->intent = NIFTI_INTENT_POINTSET;
    coords->datatype = NIFTI_TYPE_FLOAT32;
    coords->ind_ord = GIFTI_IND_ORD_ROW_MAJOR;
    coords->num_dim = 2;
    coords->dims[0] = mris->nvertices;        /* In highest first, dim0 = rows */
    coords->dims[1] = 3;                      /* In highest first, dim1 = cols */
    coords->encoding = GIFTI_ENCODING_B64GZ;  // data stored in gzip'd base64
#if (BYTE_ORDER == LITTLE_ENDIAN)
    coords->endian = GIFTI_ENDIAN_LITTLE;
#else
    coords->endian = GIFTI_ENDIAN_BIG;
#endif

    /* The implementation as of 12/19/2022 set
     * <DataSpace> = NIFTI_XFORM_UNKNOWN
     * <MatrixData> = mris->SRASToTalSRAS_
     * <TransformedSpace> = NIFTI_XFORM_TALAIRACH
     */
    coords->coordsys = NULL;             // empty, unless we find something here...

    if (mris->useRealRAS)
    {
      // surface XYZ coordinates are in scanner space
      if (mris->vg.valid)
      {
        MATRIX *S = vg_i_to_r(&mris->vg);
        MATRIX *T = TkrVox2RASfromVolGeom(&mris->vg);
        MATRIX *Sinv = MatrixInverse(S, NULL);
        MATRIX *xform = MatrixMultiply(T, Sinv, NULL);

        gifti_add_empty_CS(coords);
        int idx = coords->numCS - 1;

        //  <DataSpace> = NIFTI_XFORM_SCANNER_ANAT
        //  <MatrixData> = transform matrix go from scanner space to Freesurfer tkregister space
        //  <TransformedSpace> = NIFTI_XFORM_UNKNOWN (Freesurfer tkregister space)
        coords->coordsys[idx]->dataspace = strcpyalloc("NIFTI_XFORM_SCANNER_ANAT");
        coords->coordsys[idx]->xformspace = strcpyalloc("NIFTI_XFORM_UNKNOWN");

        for (int r = 1; r <= 4; r++)
          for (int c = 1; c <= 4; c++)
            coords->coordsys[idx]->xform[r - 1][c - 1] = xform->rptr[r][c];

        MatrixFree(&S);
        MatrixFree(&T);
        MatrixFree(&Sinv);
        MatrixFree(&xform);
      }
      else
      {
        gifti_add_empty_CS(coords);
        int idx = coords->numCS - 1;

        coords->coordsys[idx]->dataspace = strcpyalloc("NIFTI_XFORM_SCANNER_ANAT");
        coords->coordsys[idx]->xformspace = strcpyalloc("NIFTI_XFORM_SCANNER_ANAT");

        MATRIX *xform = MatrixIdentity(4, NULL);
        for (int r = 1; r <= 4; r++)
          for (int c = 1; c <= 4; c++)
            coords->coordsys[idx]->xform[r - 1][c - 1] = xform->rptr[r][c];

        MatrixFree(&xform);
      }
    }
    else
    {
      // surface XYZ coordinates are in tkregister space
      if (mris->vg.valid)
      {
        MATRIX *S = vg_i_to_r(&mris->vg);
        MATRIX *T = TkrVox2RASfromVolGeom(&mris->vg);
        MATRIX *Tinv = MatrixInverse(T, NULL);
        MATRIX *xform = MatrixMultiply(S, Tinv, NULL);

        gifti_add_empty_CS(coords);
        int idx = coords->numCS - 1;

        //  <DataSpace> = NIFTI_XFORM_UNKNOWN (Freesurfer tkregister space)
        //  <MatrixData> = transform matrix go from Freesurfer tkregister space to scanner space
        //  <TransformedSpace> = NIFTI_XFORM_SCANNER_ANAT
        coords->coordsys[idx]->dataspace = strcpyalloc("NIFTI_XFORM_UNKNOWN");
        coords->coordsys[idx]->xformspace = strcpyalloc("NIFTI_XFORM_SCANNER_ANAT");

        for (int r = 1; r <= 4; r++)
          for (int c = 1; c <= 4; c++)
            coords->coordsys[idx]->xform[r - 1][c - 1] = xform->rptr[r][c];

        MatrixFree(&S);
        MatrixFree(&T);
        MatrixFree(&Tinv);
        MatrixFree(&xform);
      }
      else
      {
        // ??? read into a local MRIS ???
        MRISreadTransform(mris, out_fname);  // tries to get xform from out_fname
        if (mris->SRASToTalSRAS_ && mris->SRASToTalSRAS_->rows == 4 && mris->SRASToTalSRAS_->cols == 4) {
          gifti_add_empty_CS(coords);
          int idx = coords->numCS - 1;
          // found a valid xform, so use it...
          coords->coordsys[idx]->dataspace = strcpyalloc("NIFTI_XFORM_UNKNOWN");
          coords->coordsys[idx]->xformspace = strcpyalloc("NIFTI_XFORM_TALAIRACH");
          MATRIX *xform = mris->SRASToTalSRAS_;
          int r, c;
          for (r = 1; r <= 4; r++)
            for (c = 1; c <= 4; c++) {
              coords->coordsys[idx]->xform[r - 1][c - 1] = xform->rptr[r][c];
            }
          }
      }
    }

    coords->nvals = gifti_darray_nvals(coords);
    gifti_datatype_sizes(coords->datatype, &coords->nbyper, NULL);

    /* Allocate the data array. */
    coords->data = NULL;
    coords->data = (void *)calloc(coords->nvals, coords->nbyper);
    if (NULL == coords->data) {
      fprintf(stderr,
              "MRISwriteGIFTISurface: couldn't allocate coords data of "
              "length %d, element size %d\n",
              (int)coords->nvals,
              coords->nbyper);
      gifti_free_image(image);
      return ERROR_NOMEMORY;
    }

    /* Copy in all our data. */
    int vertex_index;
    for (vertex_index = 0; vertex_index < mris->nvertices; vertex_index++) {
      if (mris->vertices[vertex_index].ripflag) {
        continue;
      }
      gifti_set_DA_value_2D(coords, vertex_index, 0, mris->vertices[vertex_index].x);
      gifti_set_DA_value_2D(coords, vertex_index, 1, mris->vertices[vertex_index].y);
      gifti_set_DA_value_2D(coords, vertex_index, 2, mris->vertices[vertex_index].z);
    }

    /*
     * Faces
     */
    giiDataArray *faces = gifti_alloc_and_add_darray(image);
    if (NULL == faces) {
      fprintf(stderr, "MRISwriteGIFTISurface: couldn't allocate giiDataArray\n");
      gifti_free_image(image);
      return ERROR_NOMEMORY;
    }

    /* count the real number of faces (the ones that dont have a vertex
       with a ripflag set) */
    int numFaces = 0;
    int face_index;
    for (face_index = 0; face_index < mris->nfaces; face_index++) {
      if (mris->vertices[mris->faces[face_index].v[0]].ripflag) {
        continue;
      }
      if (mris->vertices[mris->faces[face_index].v[1]].ripflag) {
        continue;
      }
      if (mris->vertices[mris->faces[face_index].v[2]].ripflag) {
        continue;
      }
      numFaces++;
    }

    /* Set its attributes. */
    faces->intent = NIFTI_INTENT_TRIANGLE;
    faces->datatype = NIFTI_TYPE_INT32;
    faces->ind_ord = GIFTI_IND_ORD_ROW_MAJOR;
    faces->num_dim = 2;
    faces->dims[0] = numFaces;               /* In highest first, dim0 = rows */
    faces->dims[1] = 3;                      /* In highest first, dim1 = cols */
    faces->encoding = GIFTI_ENCODING_B64GZ;  // data stored in gzip'd base64
#if (BYTE_ORDER == LITTLE_ENDIAN)
    faces->endian = GIFTI_ENDIAN_LITTLE;
#else
    faces->endian = GIFTI_ENDIAN_BIG;
#endif
    faces->coordsys = NULL;
    faces->nvals = gifti_darray_nvals(faces);
    gifti_datatype_sizes(faces->datatype, &faces->nbyper, NULL);

    /* Allocate the data array. */
    faces->data = NULL;
    faces->data = (void *)calloc(faces->nvals, faces->nbyper);
    if (NULL == faces->data) {
      fprintf(stderr,
              "MRISwriteGIFTISurface: couldn't allocate faces data of "
              "length %d, element size %d\n",
              (int)faces->nvals,
              faces->nbyper);
      gifti_free_image(image);
      return ERROR_NOMEMORY;
    }

    /* Copy in all our face data (remembering to ignore faces which
       have a vertex with the ripflag set). */
    int faceNum = 0;
    for (face_index = 0; face_index < mris->nfaces; face_index++) {
      if (mris->vertices[mris->faces[face_index].v[0]].ripflag) {
        continue;
      }
      if (mris->vertices[mris->faces[face_index].v[1]].ripflag) {
        continue;
      }
      if (mris->vertices[mris->faces[face_index].v[2]].ripflag) {
        continue;
      }

      gifti_set_DA_value_2D(faces, faceNum, 0, mris->faces[face_index].v[0]);
      gifti_set_DA_value_2D(faces, faceNum, 1, mris->faces[face_index].v[1]);
      gifti_set_DA_value_2D(faces, faceNum, 2, mris->faces[face_index].v[2]);
      faceNum++;
    }

    /* standard meta data for surfaces */
    //printf("[DEBUG] MRISwriteGIFTISurface(): mris->fname=%s\n", (strlen(mris->fname) != 0) ? mris->fname : "null");
    if (strlen(mris->fname) != 0) {
      //printf("[DEBUG] MRISwriteGIFTISurface(): writing standard metadata ...\n");
      const char *primary = NULL, *secondary = NULL, *geotype = NULL;
      char *name = mris->fname;
      if (strstr(name, "lh.")) {
        primary = "CortexLeft";
      }
      if (strstr(name, "rh.")) {
        primary = "CortexRight";
      }
      if (strstr(name, ".orig")) {
        secondary = "GrayWhite";
      }
      if (strstr(name, ".smoothwm")) {
        secondary = "GrayWhite";
      }
      if (strstr(name, ".white")) {
        secondary = "GrayWhite";
      }
      if (strstr(name, ".graymid")) {
        secondary = "MidThickness";
      }
      if (strstr(name, ".gray")) {
        secondary = "Pial";
      }
      if (strstr(name, ".pial")) {
        secondary = "Pial";
      }
      if (strstr(name, ".orig")) {
        geotype = "Reconstruction";
      }
      if (strstr(name, ".smoothwm")) {
        geotype = "Reconstruction";  //"Anatomical";
      }
      if (strstr(name, ".white")) {
        geotype = "Anatomical";
      }
      if (strstr(name, ".gray")) {
        geotype = "Anatomical";
      }
      if (strstr(name, ".graymid")) {
        geotype = "Anatomical";
      }
      if (strstr(name, ".pial")) {
        geotype = "Anatomical";
      }
      if (strstr(name, ".inflated")) {
        geotype = "Inflated";
      }
      if (strstr(name, ".sphere")) {
        geotype = "Spherical";
      }
      if (strstr(name, ".qsphere")) {
        geotype = "Spherical";
      }
      if (strstr(name, "pial-outer")) {
        geotype = "Hull";
      }
      const char *topotype;
      if (mris->patch) {
        //geotype = "Flat";
        topotype = "Cut";
      } else {
        topotype = "Closed";
      }

      if (primary) gifti_add_to_meta(&coords->meta, "AnatomicalStructurePrimary", primary, 1);
      if (secondary) gifti_add_to_meta(&coords->meta, "AnatomicalStructureSecondary", secondary, 1);
      if (geotype) gifti_add_to_meta(&coords->meta, "GeometricType", geotype, 1);
      gifti_add_to_meta(&faces->meta, "TopologicalType", topotype, 1);
      gifti_add_to_meta(&coords->meta, "Name", name, 1);
      gifti_add_to_meta(&faces->meta, "Name", name, 1);
    } // end writing standard gifti metadata

    // add volume geometry info if valid, and surface center-coords
    if (mris->vg.valid) {
      char stmp[100];

      gifti_add_to_meta(&coords->meta, "VolGeomFname", mris->vg.fname, 1);
      sprintf(stmp, "%d", mris->vg.width);
      gifti_add_to_meta(&coords->meta, "VolGeomWidth", stmp, 1);
      sprintf(stmp, "%d", mris->vg.height);
      gifti_add_to_meta(&coords->meta, "VolGeomHeight", stmp, 1);
      sprintf(stmp, "%d", mris->vg.depth);
      gifti_add_to_meta(&coords->meta, "VolGeomDepth", stmp, 1);

      sprintf(stmp, "%f", mris->vg.xsize);
      gifti_add_to_meta(&coords->meta, "VolGeomXsize", stmp, 1);
      sprintf(stmp, "%f", mris->vg.ysize);
      gifti_add_to_meta(&coords->meta, "VolGeomYsize", stmp, 1);
      sprintf(stmp, "%f", mris->vg.zsize);
      gifti_add_to_meta(&coords->meta, "VolGeomZsize", stmp, 1);

      sprintf(stmp, "%f", mris->vg.x_r);
      gifti_add_to_meta(&coords->meta, "VolGeomX_R", stmp, 1);
      sprintf(stmp, "%f", mris->vg.x_a);
      gifti_add_to_meta(&coords->meta, "VolGeomX_A", stmp, 1);
      sprintf(stmp, "%f", mris->vg.x_s);
      gifti_add_to_meta(&coords->meta, "VolGeomX_S", stmp, 1);

      sprintf(stmp, "%f", mris->vg.y_r);
      gifti_add_to_meta(&coords->meta, "VolGeomY_R", stmp, 1);
      sprintf(stmp, "%f", mris->vg.y_a);
      gifti_add_to_meta(&coords->meta, "VolGeomY_A", stmp, 1);
      sprintf(stmp, "%f", mris->vg.y_s);
      gifti_add_to_meta(&coords->meta, "VolGeomY_S", stmp, 1);

      sprintf(stmp, "%f", mris->vg.z_r);
      gifti_add_to_meta(&coords->meta, "VolGeomZ_R", stmp, 1);
      sprintf(stmp, "%f", mris->vg.z_a);
      gifti_add_to_meta(&coords->meta, "VolGeomZ_A", stmp, 1);
      sprintf(stmp, "%f", mris->vg.z_s);
      gifti_add_to_meta(&coords->meta, "VolGeomZ_S", stmp, 1);

      sprintf(stmp, "%f", mris->vg.c_r);
      gifti_add_to_meta(&coords->meta, "VolGeomC_R", stmp, 1);
      sprintf(stmp, "%f", mris->vg.c_a);
      gifti_add_to_meta(&coords->meta, "VolGeomC_A", stmp, 1);
      sprintf(stmp, "%f", mris->vg.c_s);
      gifti_add_to_meta(&coords->meta, "VolGeomC_S", stmp, 1);

      sprintf(stmp, "%f", mris->xctr);
      gifti_add_to_meta(&coords->meta, "SurfaceCenterX", stmp, 1);
      sprintf(stmp, "%f", mris->yctr);
      gifti_add_to_meta(&coords->meta, "SurfaceCenterY", stmp, 1);
      sprintf(stmp, "%f", mris->zctr);
      gifti_add_to_meta(&coords->meta, "SurfaceCenterZ", stmp, 1);
    }

    // group avg surface area, TAG_GROUP_AVG_SURFACE_AREA
    if (!FZERO(mris->group_avg_surface_area))
    {
      char group_avg_surface_area[100] = {'\0'};

      sprintf(group_avg_surface_area, "%.20f", mris->group_avg_surface_area);
      gifti_add_to_meta(&coords->meta, "TAG_GROUP_AVG_SURFACE_AREA", group_avg_surface_area, 1);
    }

    // TAG_CMDLINE
    if (mris->ncmds > 0)
    {
      char ncmds[20] = {'\0'};
      sprintf(ncmds, "%d", mris->ncmds);
      gifti_add_to_meta(&coords->meta, "NUM_TAG_CMDLINE", ncmds, 1);

      for (int ncmd = 0; ncmd < mris->ncmds; ncmd++)
      {
        char cmdline[TAG_CMDLINE_LEN] = {'\0'};
        snprintf(cmdline, sizeof(cmdline), "%s", mris->cmdlines[ncmd]);

        char tag[24] = {'\0'};
        sprintf(tag, "TAG_CMDLINE#%d", ncmd);
        gifti_add_to_meta(&coords->meta, tag, cmdline, 1);
      }
    }

    return NO_ERROR;
} // end of MRISwriteGIFTISurface()



int MRISwriteGIFTI(MRIS* mris, const MRI *mri, int intent_code, const char *out_fname, const char *curv_fname, const char *datatype)
{
  if (NULL == mris || NULL == out_fname) {
    fprintf(stderr, "MRISwriteGIFTI: invalid parameter\n");
    return ERROR_BADPARM;
  }

  gifti_image *image = (gifti_image *)calloc(1, sizeof(gifti_image));
  if (NULL == image) {
    fprintf(stderr, "MRISwriteGIFTI: couldn't allocate image\n");
    return ERROR_NOMEMORY;
  }
  image->version = strcpyalloc(GIFTI_XML_VERSION);

  insertCommonMetaData(&image->meta);
  if (strlen(mris->subject_name)) {
    gifti_add_to_meta(&image->meta, "SubjectID", mris->subject_name, 1);
  }

  // 
  int error = MRISwriteGIFTIIntent(mris, mri, image, intent_code, out_fname, curv_fname, datatype);
  if (error != NO_ERROR)
    return error;

  // make sure version is recoded before validation
  if (!strcmp(image->version, "1")) {
    free(image->version);
    image->version = strcpyalloc(GIFTI_XML_VERSION);
  }

  /* check for compliance */
  int valid = gifti_valid_gifti_image(image, 1);
  if (valid == 0) {
    fprintf(stderr, "MRISwriteGIFTI: GIFTI file %s is invalid!\n", out_fname);
    gifti_free_image(image);
    return ERROR_BADFILE;
  }

  /* Write the file. */
  if (gifti_write_image(image, out_fname, 1)) {
    fprintf(stderr, "MRISwriteGIFTI: couldn't write image\n");
    gifti_free_image(image);
    return ERROR_BADFILE;
  }

  gifti_free_image(image);

  return ERROR_NONE;
} // end of MRISwriteGIFTI(MRIS *mris, const MRI *mri, ...)


int MRISwriteGIFTICombined(MRIS *mris, std::vector<OverlayInfoStruct> *poverlays, const char *out_fname)
{
  if (NULL == mris || NULL == out_fname) {
    fprintf(stderr, "MRISwriteGIFTICombined: invalid parameter\n");
    return ERROR_BADPARM;
  }

  gifti_image *image = (gifti_image *)calloc(1, sizeof(gifti_image));
  if (NULL == image) {
    fprintf(stderr, "MRISwriteGIFTICombined: couldn't allocate image\n");
    return ERROR_NOMEMORY;
  }
  image->version = strcpyalloc(GIFTI_XML_VERSION);

  insertCommonMetaData(&image->meta);
  if (strlen(mris->subject_name)) {
    gifti_add_to_meta(&image->meta, "SubjectID", mris->subject_name, 1);
  }

  // write surface
  int error = MRISwriteGIFTISurface(mris, image, out_fname); 
  if (error != NO_ERROR)
    return error;  

  // write overlays
  int noverlay = (*poverlays).size();
  for (int n = 0; n < noverlay; n++)
  {
    MRI *overlaymri = (*poverlays)[n].__overlaymri;

    int overlaytype = (*poverlays)[n].__type;
    int giftiintent = MRISurfOverlay::getGIFTIIntent(overlaytype);

    const char *shapedatatype = (*poverlays)[n].__shapedatatype;
    int error = MRISwriteGIFTIIntent(mris, overlaymri, image, giftiintent, out_fname, (*poverlays)[n].__foverlay, shapedatatype);
    if (error != NO_ERROR)
      return error;
  }

  // make sure version is recoded before validation
  if (!strcmp(image->version, "1")) {
    free(image->version);
    image->version = strcpyalloc(GIFTI_XML_VERSION);
  }

  /* check for compliance */
  int valid = gifti_valid_gifti_image(image, 1);
  if (valid == 0) {
    fprintf(stderr, "MRISwriteGIFTICombined: GIFTI file %s is invalid!\n", out_fname);
    gifti_free_image(image);
    return ERROR_BADFILE;
  }

  /* Write the file. */
  if (gifti_write_image(image, out_fname, 1)) {
    fprintf(stderr, "MRISwriteGIFTICombined: couldn't write image\n");
    gifti_free_image(image);
    return ERROR_BADFILE;
  }

  gifti_free_image(image);

  return ERROR_NONE;
} // end of MRISwriteGIFTICombined()


int MRISwriteGIFTIIntent(MRIS *mris, const MRI *mri, gifti_image *image, int intent_code, const char *out_fname, const char *curv_fname, const char *datatype)
{
  /* -------------------------------------------------------
   * Surface file
   */
  if (intent_code == NIFTI_INTENT_POINTSET || intent_code == NIFTI_INTENT_TRIANGLE) {
    return MRISwriteGIFTISurface(mris, image, out_fname); 
  }  // end of if NIFTI_INTENT_POINTSET or NIFTI_INTENT_TRIANGLE

  /* -------------------------------------------------------
   * Shape file
   */
  if (intent_code == NIFTI_INTENT_SHAPE) {
    return MRISwriteGIFTIShape(mris, mri, image, intent_code, curv_fname, datatype);
  }  // end of if NIFTI_INTENT_SHAPE

  /* -------------------------------------------------------
   * Label file
   */
  if (intent_code == NIFTI_INTENT_LABEL) {
    return MRISwriteGIFTILabel(mris, image, intent_code); 
  }  // end of if NIFTI_INTENT_LABEL

  /* -------------------------------------------------------
   * Statistics file
   */
  if (intent_code == NIFTI_INTENT_CORREL || intent_code == NIFTI_INTENT_TTEST || intent_code == NIFTI_INTENT_FTEST ||
      intent_code == NIFTI_INTENT_ZSCORE || intent_code == NIFTI_INTENT_CHISQ || intent_code == NIFTI_INTENT_BETA ||
      intent_code == NIFTI_INTENT_BINOM || intent_code == NIFTI_INTENT_GAMMA || intent_code == NIFTI_INTENT_POISSON ||
      intent_code == NIFTI_INTENT_NORMAL || intent_code == NIFTI_INTENT_FTEST_NONC ||
      intent_code == NIFTI_INTENT_CHISQ_NONC || intent_code == NIFTI_INTENT_LOGISTIC ||
      intent_code == NIFTI_INTENT_LAPLACE || intent_code == NIFTI_INTENT_UNIFORM ||
      intent_code == NIFTI_INTENT_TTEST_NONC || intent_code == NIFTI_INTENT_WEIBULL ||
      intent_code == NIFTI_INTENT_CHI || intent_code == NIFTI_INTENT_INVGAUSS || intent_code == NIFTI_INTENT_EXTVAL ||
      intent_code == NIFTI_INTENT_PVAL || intent_code == NIFTI_INTENT_LOGPVAL ||
      intent_code == NIFTI_INTENT_LOG10PVAL || intent_code == NIFTI_INTENT_ESTIMATE) {
    return MRISwriteGIFTIStats(mris, mri, image, intent_code, curv_fname, datatype); 
  }  // end of if NIFTI_INTENT_<stats>

  return NO_ERROR;
} // end of MRISwriteGIFTIIntent(MRIS *mris, const MRI *mri, ...)


/*
 * Shape file
 *       intent_code = NIFTI_INTENT_SHAPE
 */
int MRISwriteGIFTIShape(MRIS *mris, const MRI *mri, gifti_image *image, int intent_code, const char *curv_fname, const char *shapedatatype)
{
#if 0
    // data is in mri
    // ??? input mris should have curvature information ???
    /// ??? why read again ??? 
    if (MRISreadCurvatureFile(mris, curv_fname)) {
      fprintf(stderr, "MRISwriteGIFTIShape: couldn't read %s\n", curv_fname);
      gifti_free_image(image);
      return ERROR_BADFILE;
    }
#endif

    giiDataArray *shape = gifti_alloc_and_add_darray(image);
    if (NULL == shape) {
      fprintf(stderr, "MRISwriteGIFTIShape: couldn't allocate giiDataArray\n");
      gifti_free_image(image);
      return ERROR_NOMEMORY;
    }

    /* Set its attributes. */
    shape->intent = NIFTI_INTENT_SHAPE;
    shape->datatype = NIFTI_TYPE_FLOAT32;
    shape->ind_ord = GIFTI_IND_ORD_ROW_MAJOR;
    shape->num_dim = 1;
    shape->dims[0] = mris->nvertices;
    shape->dims[1] = 0;
    shape->encoding = GIFTI_ENCODING_B64GZ;  // data stored in gzip'd base64
#if (BYTE_ORDER == LITTLE_ENDIAN)
    shape->endian = GIFTI_ENDIAN_LITTLE;
#else
    shape->endian = GIFTI_ENDIAN_BIG;
#endif
    shape->coordsys = NULL;
    shape->nvals = gifti_darray_nvals(shape);
    gifti_datatype_sizes(shape->datatype, &shape->nbyper, NULL);

    /* include some metadata describing this shape */
    gifti_add_to_meta(&shape->meta, "Name", curv_fname, 1);
    if (shapedatatype) {
      gifti_add_to_meta(&shape->meta, "ShapeDataType", shapedatatype, 1);
    }

    /* Allocate the data array. */
    shape->data = NULL;
    shape->data = (void *)calloc(shape->nvals, shape->nbyper);
    if (NULL == shape->data) {
      fprintf(stderr,
              "MRISwriteGIFTIShape: couldn't allocate shape data of "
              "length %d, element size %d\n",
              (int)shape->nvals,
              shape->nbyper);
      gifti_free_image(image);
      return ERROR_NOMEMORY;
    }

    /* Copy in all our data. */
    // loop through MRI crs
    for (int s = 0; s < mri->depth; s++)
    {
      for (int r = 0; r < mri->height; r++)
      {
        for (int c = 0; c < mri->width; c++)
        {
          if (mris->vertices[c].ripflag)
            continue;

          float curv = MRIgetVoxVal(mri, c, r, s, 0);
          gifti_set_DA_value_2D(shape, c, 0, curv);
        }
      }
    }
#if 0
    int vno;
    for (vno = 0; vno < mris->nvertices; vno++) {
      if (mris->vertices[vno].ripflag) {
        continue;
      }
      gifti_set_DA_value_2D(shape, vno, 0, mris->vertices[vno].curv);
    }
#endif

    return NO_ERROR;
} // end of MRISwriteGIFTIShape(MRIS *mris, const MRI *mri, ...)


/*
 * Statistics file
 *       intent_code = NIFTI_INTENT_<stats>
 */
int MRISwriteGIFTIStats(MRIS *mris, const MRI *mri, gifti_image *image, int intent_code, const char *curv_fname, const char *statsdatatype)
{
    giiDataArray *stats = gifti_alloc_and_add_darray(image);
    if (NULL == stats) {
      fprintf(stderr, "MRISwriteGIFTIStats: couldn't allocate giiDataArray\n");
      gifti_free_image(image);
      return ERROR_NOMEMORY;
    }

    for (int f = 0; f < mri->nframes; f++)
    {
      /* Set its attributes. */
      stats->intent = intent_code;
      stats->datatype = NIFTI_TYPE_FLOAT32;
      stats->ind_ord = GIFTI_IND_ORD_ROW_MAJOR;
      stats->num_dim = 1;
      stats->dims[0] = mris->nvertices;
      stats->dims[1] = 0;
      stats->encoding = GIFTI_ENCODING_B64GZ;  // data stored in gzip'd base64
#if (BYTE_ORDER == LITTLE_ENDIAN)
      stats->endian = GIFTI_ENDIAN_LITTLE;
#else
      stats->endian = GIFTI_ENDIAN_BIG;
#endif
      stats->coordsys = NULL;
      stats->nvals = gifti_darray_nvals(stats);
      gifti_datatype_sizes(stats->datatype, &stats->nbyper, NULL);

      /* include some metadata describing this thing */
      gifti_add_to_meta(&stats->meta, "Intent_code", gifti_intent_to_string(intent_code), 1);
      if (intent_code == NIFTI_INTENT_UNIFORM) {
        gifti_add_to_meta(&stats->meta, "Intent_p1", "0", 1);  // lower end
        gifti_add_to_meta(&stats->meta, "Intent_p2", "1", 1);  // upper end
      }

      /* Allocate the data array. */
      stats->data = NULL;
      stats->data = (void *)calloc(stats->nvals, stats->nbyper);
      if (NULL == stats->data) {
        fprintf(stderr,
                "MRISwriteGIFTIStats: couldn't allocate stats data of "
                "length %d, element size %d\n",
                (int)stats->nvals,
                stats->nbyper);
        gifti_free_image(image);
        return ERROR_NOMEMORY;
     }

     /* Copy in all our data. */
     // loop through MRI crs
      for (int s = 0; s < mri->depth; s++)
      {
        for (int r = 0; r < mri->height; r++)
        {
          for (int c = 0; c < mri->width; c++)
	  {
            if (mris->vertices[c].ripflag)
              continue;

            float stat = MRIgetVoxVal(mri, c, r, s, f);
            gifti_set_DA_value_2D(stats, c, 0, stat);
          }
        }
      }
#if 0
      int vno;
      for (vno = 0; vno < mris->nvertices; vno++) {
        if (mris->vertices[vno].ripflag) {
          continue;
        }
        gifti_set_DA_value_2D(stats, vno, 0, mris->vertices[vno].stat);
      }
#endif
    }

    return NO_ERROR;
} // end of MRISwriteGIFTIStats(MRIS *mris, const MRI *mri, ...)


int getShapeStatIntentCount(const char *fgifti, int *nVertices, int *nFaces)
{
  /*
   * attempt to read the file
   */
  gifti_image *image = gifti_read_image(fgifti, 1);
  if (NULL == image) {
    fprintf(stderr, "getShapeStatIntentCount(): gifti_read_image() returned NULL\n");
    return 0;
  }

  // make sure version is recoded before validation
  if (!strcmp(image->version, "1")) {
    free(image->version);
    image->version = strcpyalloc(GIFTI_XML_VERSION);
  }

  /*
   * check for compliance
   */
  int valid = gifti_valid_gifti_image(image, 1);
  if (valid == 0) {
    fprintf(stderr, "getShapeStatIntentCount(): GIFTI file %s is invalid!\n", fgifti);
    gifti_free_image(image);
    return 0;
  }

  /*
   * Now parse the DataArrays, count NIFTI_INTENT_SHAPE and NIFTI_INTENT_<stat>
   */
  int count = 0;
  int endDAnum = image->numDA;
  for (int numDA = 0; numDA < endDAnum; numDA++) {
    giiDataArray *darray = image->darray[numDA];

    if (darray->intent == NIFTI_INTENT_POINTSET)
    {
      // get number of vertices
      long long num_vertices = 0;
      long long  num_cols = 0;
      if (darray->ind_ord == GIFTI_IND_ORD_ROW_MAJOR) {
        // RowMajorOrder
        gifti_DA_rows_cols(darray, &num_vertices, &num_cols);
      }
      else {
        // ColumnMajorOrder
        gifti_DA_rows_cols(darray, &num_cols, &num_vertices);
      }
      if (num_vertices <= 0 || num_cols != 3) {
        fprintf(stderr,
                "getShapeStatIntentCount(): malformed coords data array in file "
                "%s: num_vertices=%d num_cols=%d\n",
                fgifti,
                (int)num_vertices,
                (int)num_cols);
        gifti_free_image(image);
        return 0;
      }

      if (nVertices != NULL)
        *nVertices = num_vertices;
      continue;
    }
    else if (darray->intent == NIFTI_INTENT_TRIANGLE)
    {
      // get number of triangle faces
      long long num_faces = 0;
      long long num_cols = 0;
      if (darray->ind_ord == GIFTI_IND_ORD_ROW_MAJOR) {
        // RowMajorOrder
        gifti_DA_rows_cols(darray, &num_faces, &num_cols);
      }
      else {
        // ColumnMajorOrder
        gifti_DA_rows_cols(darray, &num_cols, &num_faces);
      }
      if (num_faces <= 0 || num_cols != 3) {
        fprintf(stderr,
                "getShapeStatIntentCount(): malformed faces data array in file "
                "%s: num_faces=%d num_cols=%d\n",
                fgifti,
                (int)num_faces,
                (int)num_cols);
        gifti_free_image(image);
        return 0;
      }

      if (nFaces != NULL)
        *nFaces = num_faces;
      continue;
    }
    else if ((darray->intent == NIFTI_INTENT_LABEL)      || 
             (darray->intent == NIFTI_INTENT_GENMATRIX)  ||
             (darray->intent == NIFTI_INTENT_VECTOR)     || 
             (darray->intent == NIFTI_INTENT_RGB_VECTOR) ||
             (darray->intent == NIFTI_INTENT_RGBA_VECTOR))
      // skip these intents
      continue;

    count++;
  } 

  /*
   * And we're done.
   */
  gifti_free_image(image);

  return count;
} // end of getShapeStatIntentCount()
