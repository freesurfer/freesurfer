/*
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


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mri.h"
#include "error.h"

#define SET_METHOD_XYZ 0
#define SET_METHOD_RANDOM 1
#define SET_METHOD_CONSTANT 2

char sSetMethods[][256] = {
                            "XYZ", "Random", "Constant"
                          };

const char* Progname="makevol";

void PrintUsage ( const char* isError ) {

  if ( isError ) {
    printf( "\nERROR: %s\n\n", isError );
  }

  printf( "Usage: makevol [OPTION]...\n" );
  printf( "Create a volume with given parameters.\n" );
  printf( "\n" );
  printf( "Options:\n" );
  printf( "-f, --filename FILENAME   Write volume to the given file name,\n" );
  printf( "                          implying type. Default=new_volume.mgz\n");
  printf( "\n" );
  printf( "-x, --width WIDTH         Use integer WIDTH as the x dimension. Default=256\n" );
  printf( "-y, --height HEIGHT       Use integer HEIGHT as the y dimension. Default=256\n" );
  printf( "-z, --depth DEPTH         Use integer DEPTH as the z dimension. Default=256\n" );
  printf( "\n" );
  printf( "--sizex SIZEX             Use float SIZEX as the x resolution. Default=1.0\n");
  printf( "--sizey SIZEY             Use float SIZEY as the y resolution. Default=1.0\n");
  printf( "--sizez SIZEZ             Use float SIZEZ as the z resolution. Default=1.0\n");
  printf( "\n" );
  printf( "--set-method METHOD [VALUE]\n" );
  printf( "                          Use METHOD to fill the values. Default=xyz. METHOD\n" );
  printf( "                          can be:\n");
  printf( "                            xyz: Value is set to its x,yz, coords\n" );
  printf( "                            random: Random values from 0-255\n" );
  printf( "                            constant: Set all values to VALUE\n" );
  printf( "\n" );
}

int main ( int argc, char** argv ) {

  MRI* mri = NULL;
  int zX = 256;
  int zY = 256;
  int zZ = 256;
  int err = NO_ERROR;
  int nX = 0;
  int nY = 0;
  int nZ = 0;
  float sizeX = 1.0;
  float sizeY = 1.0;
  float sizeZ = 1.0;
  int setMethod = SET_METHOD_XYZ;
  int setValue = 0;
  char fnVol[256] = "new_volume.mgz";
  int i;
  char* arg = NULL;

  for ( i = 1; i < argc; i++ ) {

    arg = argv[i];

    if ( argv[i] && *argv[i] != '-' ) {
      printf( "ERROR: Unrecognized argument %s\n", argv[i] );
      PrintUsage( NULL );
      exit( 1 );
    }

    while ( arg[0] == '-' )
      arg = arg+1;

    if ( strlen(arg) <= 0 )
      continue;

    if ( strcmp(arg,"h") == 0 ||
         strcmp(arg,"help") == 0 ) {
      PrintUsage( NULL );
      exit( 0 );
    }

    if ( strcmp(arg,"f") == 0 ||
         strcmp(arg,"filename") == 0 ) {
      if ( i+1 >= argc ) {
        PrintUsage( "No argument to filename option." );
        exit( 1 );
      }
      strcpy( fnVol, argv[i+1] );
      i++;
    }

    if ( strcmp(arg,"x") == 0 ||
         strcmp(arg,"width") == 0 ) {
      if ( i+1 >= argc ) {
        PrintUsage( "No argument to width option." );
        exit( 1 );
      }
      zX = atoi(argv[i+1]);
      i++;
    }
    if ( strcmp(arg,"y") == 0 ||
         strcmp(arg,"height") == 0 ) {
      if ( i+1 >= argc ) {
        PrintUsage( "No argument to height option." );
        exit( 1 );
      }
      zY = atoi(argv[i+1]);
      i++;
    }
    if ( strcmp(arg,"z") == 0 ||
         strcmp(arg,"depth") == 0 ) {
      if ( i+1 >= argc ) {
        PrintUsage( "No argument to depth option." );
        exit( 1 );
      }
      zZ = atoi(argv[i+1]);
      i ++;
    }

    if ( strcmp(arg,"sizex") == 0 ) {
      if ( i+1 >= argc ) {
        PrintUsage( "No argument to sizex option." );
        exit( 1 );
      }
      sizeX = atof(argv[i+1]);
      i++;
    }
    if ( strcmp(arg,"sizey") == 0 ) {
      if ( i+1 >= argc ) {
        PrintUsage( "No argument to sizey option." );
        exit( 1 );
      }
      sizeY = atof(argv[i+1]);
      i++;
    }
    if ( strcmp(arg,"sizez") == 0 ) {
      if ( i+1 >= argc ) {
        PrintUsage( "No argument to sizez optoin." );
        exit( 1 );
      }
      sizeZ = atof(argv[i+1]);
      i++;
    }

    if ( strcmp(arg,"set-method") == 0 ) {
      if ( i+1 >= argc ) {
        PrintUsage( "No argument to set-method option." );
        exit( 1 );
      }
      if ( strcmp( argv[i+1], "xyz" ) == 0 ) {
        setMethod = SET_METHOD_XYZ;
        i++;
        printf( "set_method is xyz\n" );
      } else if ( strcmp( argv[i+1], "random" ) == 0 ) {
        setMethod = SET_METHOD_RANDOM;
        i++;
        printf( "set_method is random\n" );
      } else if ( strncmp( argv[i+1], "constant", 9 ) == 0 ) {
        if ( i+2 >= argc ) {
          PrintUsage( "No value argument to constant method option." );
          exit( 1 );
        }
        setMethod = SET_METHOD_CONSTANT;
        setValue = atoi( argv[i+2] );
        i+=2;
        printf( "set_method is constant, %d\n", setValue );
      } else {
        PrintUsage( "Unrecognized argument to set-method option" );
        exit( 1 );
      }
    }
  }

  printf( "Creating volume %s\n"
          "  width = %d height = %d depth = %d\n"
          "  xsize = %f ysize = %f zsize = %f\n"
          "  set method = %s, constant = %d\n",
          fnVol, zX, zY, zZ,
          sizeX, sizeY, sizeZ,
          sSetMethods[setMethod], setValue );

  mri = MRIalloc( zX, zY, zZ, MRI_UCHAR );
  if ( NULL == mri ) {
    fprintf( stderr, "Couldn't create volume.\n" );
    return 1;
  }
  MRIsetResolution( mri, sizeX, sizeY, sizeZ );

  switch ( setMethod ) {
  case SET_METHOD_CONSTANT:
    MRIvalueFill( mri, 0 );
    break;
  case SET_METHOD_RANDOM:
    for ( nZ = 0; nZ < zZ; nZ++ ) {
      for ( nY = 0; nY < zY; nY++ ) {
        for ( nX = 0; nX < zX; nX++ ) {
          MRIvox( mri, nX, nY, nZ ) = (int)((random()/(float)RAND_MAX)*255.0);
        }
      }
    }
    break;
  case SET_METHOD_XYZ:
    for ( nZ = 0; nZ < zZ; nZ++ ) {
      for ( nY = 0; nY < zY; nY++ ) {
        for ( nX = 0; nX < zX; nX++ ) {
          MRIvox( mri, nX, nY, nZ ) =
            (((float)nZ/(float)zZ)*255.0/3.0) +
            (((float)nY/(float)zY)*255.0/3.0) +
            (((float)nX/(float)zX)*255.0/3.0) ;
        }
      }
    }
    break;
  default:
    for ( nZ = (zZ/2) - (zZ/4); nZ < (zZ/2) + (zZ/4); nZ++ ) {
      for ( nY = (zY/2) - (zY/4); nY < (zY/2) + (zY/4); nY++ ) {
        for ( nX = (zX/2) - (zX/4); nX < (zX/2) + (zX/4); nX++ ) {
          MRIvox( mri, nX, nY, nZ ) = 255;
        }
      }
    }
  }

  err = MRIwrite( mri, fnVol );
  if ( NO_ERROR != err ) {
    fprintf( stderr, "Couldn't write volume.\n" );
    return 1;
  }

  MRIfree( &mri );

  return 0;
}
