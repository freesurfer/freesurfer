#include "xTypes.h"
#include "xDebug.h"

typedef enum {
  xUtil_tErr_NoError = 0,
  xUtil_tErr_InvalidErrorCode,
  xUtil_tErr_knNumErrorCodes
} xUtil_tErr;

xUtil_tErr xUtil_BreakStringIntoPathAndStem ( char* isPathAndStem,
                char* osPath,
                char* osStem );

void xUtil_StartTimer ();
void xUtil_StopTimer ( char* isMessage );

char* xUtil_GetErrorString( xUtil_tErr iCode );
