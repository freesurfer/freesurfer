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

#define kDebugging          1 

/* these can be redefined to do anything, or nothing at all. they can be used
   to print to an error log file, to print with a specfic prefix or suffix,
   print line numbers, etc. debugging output can be enabled or disabled in
   the middle of the code. use InitDebugging to do any one-time set up work
   and DeleteDebugging to take it down. */
#ifdef kDebugging

extern char gDebuggingOn;

// regular cerr debugging output.
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

// no debugging output.
#define InitDebugging       
#define DisableDebuggingOutput  
#define EnableDebuggingOutput
#define DeleteDebugging
#define DebugCode               /*
#define EndDebugCode            */
#define DebugPrint              /*
#define EndDebugPrint           */
#define Here(n)
#endif

#endif
