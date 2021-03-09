/**
 * @brief Tests a call to MRItoImageView, ImageWrite, and TiffWriteImage
 *
 * This is meant to generate a TIFF from a volume. It can be used to
 * test if the output of TiffWriteImage in imageio.c produces valid
 * TIFFs. It calls MRItoImageView to make an Image, then ImageWrite
 * with a .tif filename to call TiffWriteImage.
 */
/*
 * Original Author: Kevin Teich
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
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>

#include "mri.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "version.h"
#include "mrimorph.h"
#include "mri_circulars.h"

const char *Progname;

int main ( int argc, char** argv ) 
{

  int nargs;
  char mri_fname[STRLEN];
  char tiff_fname[STRLEN];
  MRI* mri;
  int slice;
  IMAGE* image;

  Progname = argv[0];

  nargs = handleVersionOption(argc, argv, "tiff_write_image");
  argc -= nargs ;
  if (1 == argc)
    exit (0);

  if (3 != argc) {
    printf ("usage: %s <in vol> <out tiff>\n", Progname);
    printf ("  <in vol> is a MRIread-able volume file, and <out tiff> is the\n"
	    "  name of a TIFF file to write.");
    printf ("\n");
    exit (0);
  }

  strncpy (mri_fname, argv[1], sizeof(mri_fname));
  strncpy (tiff_fname, argv[2], sizeof(tiff_fname));

  mri = MRIread (mri_fname);
  if (!mri)
    ErrorExit(ERROR_BADPARM, "%s: could not read source volume %s",
              Progname, mri_fname) ;

  slice = mri->width / 2;

  image = MRItoImageView (mri, NULL, slice, MRI_CORONAL, 0);
  if (!image)
    ErrorExit(Gerror, "MRItoImageView failed");

  ImageWrite (image, tiff_fname);
  ImageFree (&image);

  exit (0);
}
