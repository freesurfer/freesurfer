#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>
#include "xUtilities.h"

char xUtil_ksaErrorString [xUtil_tErr_knNumErrorCodes][256] = {
  "No error.",
  "Invalid error code."
};

xUtil_tErr xUtil_BreakStringIntoPathAndStem ( char* isPathAndStem,
                char* osPath,
                char* osStem ) {

  xUtil_tErr eResult              = xUtil_tErr_NoError;
  char*      sSection             = "";
  char       sPathSection[20][60];
  int        nSection             = 0;
  int        nNumSections         = 0;

  sSection = strtok( isPathAndStem, "/" );
  if( NULL == sSection ) {

    /* no path, for some reason */
    strcpy( osPath, "" );
    strcpy( osStem, isPathAndStem );
    goto cleanup;

  } else {

    /* got the first section of the path, now get the rest. */
    strcpy( sPathSection[0], sSection );
    nSection = 1;
    while( NULL != (sSection = strtok( NULL, "/" ))) {
      strcpy( sPathSection[nSection++], sSection );
    }
    
    /* the last part was the stem */
    strcpy( osStem, sPathSection[--nSection] );
    
    /* now build the path into a single string */
    nNumSections = nSection;
    strcpy( osPath, "" );
    for( nSection = 0; nSection < nNumSections; nSection++ ) {
      sprintf( osPath, "%s/%s", osPath, sPathSection[nSection] );
    }
  }

  goto cleanup;

  goto error;
 error:

  if( xUtil_tErr_NoError != eResult ) {
    DebugPrint "Error %d in xUtil_BreakStringIntoPathAndStem: %s\n",
      eResult, xUtil_GetErrorString(eResult) EndDebugPrint;
  }

 cleanup:

  return eResult;
}

char* xUtil_GetErrorString ( xUtil_tErr ieCode ) {

  xUtil_tErr eCode = ieCode;

  if ( ieCode    < 0
       || ieCode >= xUtil_tErr_knNumErrorCodes ) {
    eCode = xUtil_tErr_InvalidErrorCode;
  }

  return xUtil_ksaErrorString [eCode];
}

struct timeval sStartTime;
struct timeval sEndTime;

void xUtil_StartTimer () {

  gettimeofday( &sStartTime, NULL );
}

void xUtil_StopTimer ( char* isMessage ) {

  gettimeofday( &sEndTime, NULL );

  DebugPrint "%s: %lu usec\n", isMessage,
    (sEndTime.tv_sec*1000000 + sEndTime.tv_usec) -
    (sStartTime.tv_sec*1000000 + sStartTime.tv_usec)
    EndDebugPrint;
}




