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
 *    $Date: 2008/04/17 03:33:17 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2008,
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

  /* if it has more than 3 dimensions, then maybe its diffusion data */
  if (nin->dim != 3) 
  {
    /* say something about the array */
    printf("mriNrrdRead: \"%s\" is a %d-dimensional nrrd of type %d (%s)\n",
           fname, nin->dim, nin->type,
           airEnumStr(nrrdType, nin->type));
    printf("mriNrrRead: the array contains %d elements, %d bytes in size\n",
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
  // recall that we know that nin->dim == 3
  int width  = nin->axis[0].size;
  int height = nin->axis[1].size;
  int depth  = nin->axis[2].size;
  MRI *mri = MRIalloc( width, height, depth, type );
  if( NULL == mri )
  {
    printf("mriNrrdRead: Couldn't allocate MRI of size %d %d %d\n",
           width, height, depth);
    return NULL;
  }

  // Copy all the pixel data.
  int x,y,z;
  for( z = 0; z < depth; z++ )
  {
    for( y = 0; y < height; y++ )
    {
      for( x = 0; x < width; x++ )
      {
        int index = (x) + (y * width) + (z * width * height);
        switch( type )
        {
        case MRI_UCHAR:
        {
          char *_uc = (char*)nin->data;
          MRIseq_vox( mri, x, y, z, 0 ) = (BUFTYPE)_uc[index];
          break;
        }
        case MRI_SHORT:
        {
          short *_s = (short*)nin->data;
          MRISseq_vox( mri, x, y, z, 0 ) = (short)_s[index];
          break;
        }
        case MRI_INT:
        {
          int *_i = (int*)nin->data;
          MRIIseq_vox( mri, x, y, z, 0 ) = (int)_i[index];
          break;
        }
        case MRI_LONG:
        {
          long *_l = (long*)nin->data;
          MRILseq_vox( mri, x, y, z, 0 ) = (long)_l[index];
          break;
        }
        case MRI_FLOAT:
        {
          float *_f = (float*)nin->data;
          MRIFseq_vox( mri, x, y, z, 0 ) = (float)_f[index];
          break;
        }
        default:
          return NULL;
        }
      }
    }
  }

  // Get and set the spacing.
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

