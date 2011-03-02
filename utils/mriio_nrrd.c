/**
 * @file  mriio_nrrd.c
 * @brief Provides Nrrd IO to Freesurfer
 *
 * Implements mriNrrdRead/Write using NrrdIO lib.
 */
/*
 * Original Author: Nick Schmansky
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:46 $
 *    $Revision: 1.3 $
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

#include <ctype.h>
#include <unistd.h>
#include <memory.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <ctype.h>
#include <dirent.h>
#include <time.h>
#include <errno.h>
#include <fcntl.h>

#include "utils.h"
#include "error.h"
#include "proto.h"
#include "mri.h"
#include "NrrdIO/NrrdIO.h"
extern MRI *mriNrrdReadDiffusion(char *fname, int read_volume);

MRI *mriNrrdRead(char *fname, int read_volume)
{
  if (!nrrdSanity())
  {
    fprintf(stderr, "\n");
    fprintf(stderr, "!!! nrrd sanity check FAILED: fix and re-compile\n");
    char *err = biffGet(NRRD);
    fprintf(stderr, "%s\n", err);
    free(err);
    return NULL;
  }

  /* create a nrrd; at this point this is just an empty container */
  Nrrd *nin = nrrdNew();
  NrrdIoState *nio = nrrdIoStateNew();

  /* read in the nrrd from file */
  if (nrrdLoad(nin, fname, nio))
  {
    char *err = biffGetDone(NRRD);
    fprintf(stderr, "mriNrrdRead: trouble reading \"%s\":\n%s", fname, err);
    free(err);
    return NULL;
  }

  /* if it has more than 3 dimensions, then maybe its diffusion data that
     the ITK reader might understand (assuming we're building with ITK)*/
  if ( ((nin->dim != 3) && (nin->dim != 4)) || (nin->spaceDim != 3) )
  {
    /* say something about the array */
    printf("mriNrrdRead: \"%s\" is a %d-dimensional nrrd of type %d (%s)\n",
           fname, nin->dim, nin->type,
           airEnumStr(nrrdType, nin->type));
    printf("mriNrrdRead: the array contains %d elements, %d bytes in size\n",
           (int)nrrdElementNumber(nin), (int)nrrdElementSize(nin));
    if (nin->content) printf("mriNrrdRead: content: %s\n", nin->content);

    return mriNrrdReadDiffusion(fname, read_volume);
  }

  /* print out the key/value pairs present */
  int kvn = nrrdKeyValueSize(nin);
  if (kvn)
  {
    int kvi;
    for (kvi=0; kvi<kvn; kvi++)
    {
      char *val, *key;
      nrrdKeyValueIndex(nin, &key, &val, kvi);
      printf("mriNrrdRead: key:value %d = %s:%s\n", kvi, key, val);
      free(key);
      free(val);
      key = val = NULL;
    }
  }

  // Get the component type.
  int type = MRI_UCHAR;
  switch (nin->type)
  {
  case nrrdTypeChar:
  case nrrdTypeUChar: type = MRI_UCHAR; break;
  case nrrdTypeShort:
  case nrrdTypeUShort: type = MRI_SHORT; break;
  case nrrdTypeInt:
  case nrrdTypeUInt: type = MRI_INT; break;
  case nrrdTypeLLong:
  case nrrdTypeULLong: type = MRI_LONG; break;
  case nrrdTypeFloat: type = MRI_FLOAT; break;
  default:
    printf("mriNrrdRead: Unsupported type: %d (%s)\n",
           nin->type, airEnumStr(nrrdType, nin->type));
    return NULL;
  }

  // alloc mri struct with the correct dimensions.
  int width  = nin->axis[0].size;
  int height = nin->axis[1].size;
  int depth  = nin->axis[2].size;
  int nframes = 1; // default nin->dim = 3
  if (nin->dim == 4) nframes = nin->axis[3].size; // multiple frames found
  MRI *mri = MRIallocSequence( width, height, depth, type, nframes );
  if( NULL == mri )
  {
    printf("mriNrrdRead: Couldn't allocate MRI of size %d %d %d %d\n",
           width, height, depth, nframes);
    return NULL;
  }

  // Copy all the pixel data.
  int x,y,z,f;
  if (type == MRI_UCHAR)
  {
    for (f = 0; f < nframes; f++)
      for (z = 0; z < depth; z++)
        for (y = 0; y < height; y++)
          for (x = 0; x < width; x++)
          {
            int index = x + (y * width) + 
              (z * width * height) + (f * width * height * depth);
            unsigned char *_uc = (unsigned char*)nin->data;
            MRIseq_vox(mri, x, y, z, f) = (BUFTYPE)_uc[index];
          }
  }
  else if (type == MRI_SHORT)
  {
    for (f = 0; f < nframes; f++)
      for (z = 0; z < depth; z++)
        for (y = 0; y < height; y++)
          for (x = 0; x < width; x++)
          {
            int index = x + (y * width) + 
              (z * width * height) + (f * width * height * depth);
            short *_s = (short*)nin->data;
            MRISseq_vox( mri, x, y, z, f ) = (short)_s[index];
          }
  }
  else if (type == MRI_INT)
  {
    for (f = 0; f < nframes; f++)
      for (z = 0; z < depth; z++)
        for (y = 0; y < height; y++)
          for (x = 0; x < width; x++)
          {
            int index = x + (y * width) + 
              (z * width * height) + (f * width * height * depth);
            int *_i = (int*)nin->data;
            MRIIseq_vox( mri, x, y, z, f ) = (int)_i[index];
          }
  }
  else if (type == MRI_LONG)
  {
    for (f = 0; f < nframes; f++)
      for (z = 0; z < depth; z++)
        for (y = 0; y < height; y++)
          for (x = 0; x < width; x++)
          {
            int index = x + (y * width) + 
              (z * width * height) + (f * width * height * depth);
            long *_l = (long*)nin->data;
            MRILseq_vox( mri, x, y, z, f ) = (long)_l[index];
          }
  }
  else if (type == MRI_FLOAT)
  {
    for (f = 0; f < nframes; f++)
      for (z = 0; z < depth; z++)
        for (y = 0; y < height; y++)
          for (x = 0; x < width; x++)
          {
            int index = x + (y * width) + 
              (z * width * height) + (f * width * height * depth);
            float *_f = (float*)nin->data;
            MRIFseq_vox( mri, x, y, z, f ) = (float)_f[index];
          }
  }
  else
  {
    printf("mriNrrdRead: Unsupported type=%d\n", type);
    return NULL;
  }

  // get and set the origin
  mri->c_r = (float)nin->spaceOrigin[0];
  mri->c_a = (float)nin->spaceOrigin[1];
  mri->c_s = (float)nin->spaceOrigin[2];
  mri->ras_good_flag = 1;

  // get and set the spacing
  mri->xsize = (float)nin->axis[0].spaceDirection[0];
  mri->ysize = (float)nin->axis[1].spaceDirection[1];
  mri->zsize = (float)nin->axis[2].spaceDirection[2];

  // set orientation string
  switch (nin->space)
  {
  case nrrdSpaceRightAnteriorSuperior:
  {
    MRIorientationStringToDircos(mri, "RAS");
    break;
  }
  case nrrdSpaceLeftAnteriorSuperior:
  {
    MRIorientationStringToDircos(mri, "LAS");
    break;
  }
  case nrrdSpaceLeftPosteriorSuperior:
  {
    MRIorientationStringToDircos(mri, "LPS");
    break;
  }
  }

  return mri;
}



int mriNrrdWrite(MRI *mri, char *fname)
{
  ErrorReturn
    (ERROR_UNSUPPORTED,
     (ERROR_UNSUPPORTED,
      "mriNrrdWrite(): Nrrd output not supported!"));
}

