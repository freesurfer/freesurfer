#ifndef xDebug_H
#define xDebug_H

#include <stdio.h>

// define true and false
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

#ifndef kDebugging
#define kDebugging          1 
#endif

#if kDebugging

extern char gDebuggingOn;

/* print debugging stuff to stderr */
#define InitDebugging           if (getenv("XDEBUG")) gDebuggingOn = TRUE;
#define DeleteDebugging
#define DisableDebuggingOutput  gDebuggingOn = FALSE;
#define EnableDebuggingOutput   if (getenv("XDEBUG")) gDebuggingOn = TRUE;
#define DebugCode               
#define EndDebugCode            
#define DebugPrint              if(gDebuggingOn) { fprintf ( stderr,
#define EndDebugPrint           ); }
#define Here(n)                 if(gDebuggingOn){fprintf(stderr,"--> here %d\n",\
                                n);}

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

#endif
