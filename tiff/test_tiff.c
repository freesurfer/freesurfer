/**
 * @file  test_tiff.c
 * @brief simple test of libtiff
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/02/26 00:49:04 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2006-2008,
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
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <unistd.h>
#include "tiffio.h"

const char *Progname = NULL;

int main ( int argc, char *argv[] ) {

  TIFF* tiff;
  uint32 zImageWidth, zImageHeight;
  size_t cPixels;
  uint32* data;

  Progname = argv[0] ;

  if ( argc != 2 ) {
    printf( "\n" );
    printf( "Usage: %s FILENAME\n", Progname );
    printf( "Try opening a TIFF file.\n" );
    printf( "\n" );
    printf( "FILENAME is a TIFF image file to attempt to open.\n" );
    printf( "\n" );
    printf( "\n" );
  }

  tiff = NULL;
  tiff = TIFFOpen( argv[1], "r" );
  if ( !tiff ) {
    printf( "%s: ERROR: File %s could not be opened.\n", Progname, argv[1] );
    exit( 1 );
  }

  zImageWidth = zImageHeight = -1;
  cPixels = -1;
  data = NULL;

  TIFFGetField( tiff, TIFFTAG_IMAGEWIDTH, &zImageWidth );
  TIFFGetField( tiff, TIFFTAG_IMAGELENGTH, &zImageHeight );

  cPixels = zImageWidth * zImageHeight;

  data = (uint32*)_TIFFmalloc( cPixels * sizeof(uint32) );
  if ( data ) {
    if ( TIFFReadRGBAImage( tiff, zImageWidth, zImageHeight, data, 1 ) ) {

      if ( !data ) {
        printf( "%s: TIFFReadRGBAImageOriented returned 1, but data "
                "was NULL.\n", Progname );
        exit( 1 );
      }

      if ( zImageWidth == -1 || zImageHeight == -1 ) {
        printf( "%s: TIFFReadRGBAImageOriented returned 1, but "
                "zImageWidth was %d and zImageHeight was %d.\n",
                Progname, (int)zImageWidth, (int)zImageHeight );
        exit( 1 );
      }

    } else {
      printf( "%s: TIFFReadRGBAImageOriented failed.\n", Progname );
      exit( 1 );
    }
  } else {
    printf( "%s: _TIFFmalloc returned NULL.\n", Progname );
    exit( 1 );
  }

  TIFFClose( tiff );

  printf( "%s: %s successfully opened, decoded, and closed.\n",
          Progname, argv[1] );

  return 0;
}
