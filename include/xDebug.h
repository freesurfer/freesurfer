#ifndef xDebug_H
#define xDebug_H

#include <stdio.h>
#include "xTypes.h"

#ifndef kDebugging
#define kDebugging          1 
#endif

#define xDebug_Nothing 0
#define xDebug_Print   1
#define xDebug_File    2

void xDbg_Init ();
void xDbg_ShutDown ();
void xDbg_PrintStatus ();

#if kDebugging

extern tBoolean xDbg_gbOutput;
extern FILE*    xDbg_gStream;
extern int      xDbg_gType;
extern char*    xDbg_gsRequest;

#define InitDebugging     xDbg_Init();
#define DeleteDebugging   xDbg_ShutDown();
   
#define DisableDebuggingOutput  xDbg_gbOutput = FALSE;
#define EnableDebuggingOutput   if ( xDbg_gStream != xDebug_Nothing)         \
                                  xDbg_gbOutput = TRUE;
#define DebugCode               
#define EndDebugCode            
#define DebugPrint              if( xDbg_gbOutput ) { fprintf( xDbg_gStream,
#define EndDebugPrint           ); }
#define Here(n)                 if( xDbg_gbOutput )                          \
                                   fprintf( xDbg_gStream, "--> here %d\n", n );
#define IsDebugging        (xDebug_Print == xDbg_gType)

#else

/* to do: find a way to make these not print any debugging output. */
#define InitDebugging       
#define DisableDebuggingOutput  
#define EnableDebuggingOutput
#define DeleteDebugging
#define DebugCode               
#define EndDebugCode            
#define DebugPrint              
#define EndDebugPrint           
#define Here(n)
#endif

void xDbg_PrintStatus ();

#endif
