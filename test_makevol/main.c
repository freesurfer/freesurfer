#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mri.h"
#include "error.h"

#define SET_METHOD_XYZ 0
#define SET_METHOD_RANDOM 1
#define SET_METHOD_CONSTANT 2 

char* Progname="makevol";

int main ( int argc, char** argv ) {

  MRI* mri = NULL;
  int zX = 256;
  int zY = 256;
  int zZ = 256;
  int err = NO_ERROR;
  int nX = 0;
  int nY = 0;
  int nZ = 0;
  int setMethod = SET_METHOD_XYZ;
  int setValue = 0;
  char fnVol[256] = "new_volume.mgz";
  int i;
  char* arg = NULL;

  for( i = 1; i < argc; i++ ) {
    
    if( argv[i] && *argv[i] != '-' ) {
      printf( "Unrecognized argument %s\n", argv[i] );
      continue;
    }

    arg = argv[i];
    while( arg[0] == '-' )
      arg = arg+1;

    if( strlen(arg) <= 0 ) 
      continue;
    
    if( strcmp(arg,"f") == 0 ||
	strcmp(arg,"filename") == 0 ) {
      strcpy( fnVol, argv[i+1] );
      i+=2;
    }

    if( strcmp(arg,"x") == 0 ||
	strcmp(arg,"width") == 0 ) {
      zX = atoi(argv[i+1]);
      i+=2;
    }
    if( strcmp(arg,"y") == 0 ||
	strcmp(arg,"height") == 0 ) {
      zY = atoi(argv[i+1]);
      i+=2;
    }
    if( strcmp(arg,"z") == 0 ||
	strcmp(arg,"depth") == 0 ) {
      zZ = atoi(argv[i+1]);
      i += 2;
    }

    if( strcmp(arg,"set-method") == 0 ) {
      if( strcmp( argv[i+1], "xyz" ) == 0 ) {
	setMethod = SET_METHOD_XYZ;
	i += 2;
	printf( "set_method is xyz\n" );
      } else if( strcmp( argv[i+1], "random" ) == 0 ) {
	setMethod = SET_METHOD_RANDOM;
	i += 2;
	printf( "set_method is random\n" );
      } else if( strncmp( argv[i+1], "constant", 9 ) == 0 ) {
	setMethod = SET_METHOD_CONSTANT;
	setValue = atoi( argv[i+2] );
	i += 3;
	printf( "set_method is constant, %d\n", setValue );
      } else {
	printf( "set_method not good, should be xyz, random, or constant\n" );
	exit( 1 );
      }
    }
  }

  printf( "Creating volume %s\n"
	  "  width = %d height = %d depth = %d\n"
	  "  set method = %d, constant = %d\n",
	  fnVol, zX, zY, zZ, setMethod, setValue );
  
  mri = MRIalloc( zX, zY, zZ, MRI_UCHAR );
  if( NULL == mri ) {
    fprintf( stderr, "Couldn't create volume.\n" );
    return 1;
  }

  switch( setMethod ) {
  case SET_METHOD_CONSTANT:
    MRIvalueFill( mri, 0 );
    break;
  case SET_METHOD_RANDOM:
    for( nZ = 0; nZ < zZ; nZ++ ) {
      for( nY = 0; nY < zY; nY++ ) {
	for( nX = 0; nX < zX; nX++ ) {
	  MRIvox( mri, nX, nY, nZ ) = (int)((random()/(float)RAND_MAX)*255.0);
	}
      }
    }
    break;
  case SET_METHOD_XYZ:
    for( nZ = 0; nZ < zZ; nZ++ ) {
      for( nY = 0; nY < zY; nY++ ) {
	for( nX = 0; nX < zX; nX++ ) {
	  MRIvox( mri, nX, nY, nZ ) = 
	    (((float)nZ/(float)zZ)*255.0/3.0) + 
	    (((float)nY/(float)zY)*255.0/3.0) + 
	    (((float)nX/(float)zX)*255.0/3.0) ;
	}
      }
    }
    break;
  default:
    for( nZ = (zZ/2) - (zZ/4); nZ < (zZ/2) + (zZ/4); nZ++ ) {
      for( nY = (zY/2) - (zY/4); nY < (zY/2) + (zY/4); nY++ ) {
	for( nX = (zX/2) - (zX/4); nX < (zX/2) + (zX/4); nX++ ) {
	  MRIvox( mri, nX, nY, nZ ) = 255;
	}
      }
    }
  }

  err = MRIwrite( mri, fnVol );
  if( NO_ERROR != err ) {
    fprintf( stderr, "Couldn't write volume.\n" );
    return 1;
  }


  return 0;
}
