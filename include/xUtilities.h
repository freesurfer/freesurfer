/**
 * @file  xUtilities.h
 * @brief general purpose utils
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/02/02 19:25:19 $
 *    $Revision: 1.6 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 *
 */


#include "xTypes.h"
#include "xDebug.h"

typedef enum
{
  xUtil_tErr_Invalid = -1,
  xUtil_tErr_NoError = 0,
  xUtil_tErr_InvalidErrorCode,
  xUtil_tErr_knNumErrorCodes
} xUtil_tErr;

xUtil_tErr xUtil_BreakStringIntoPathAndStem ( char* isPathAndStem,
    char* osPath,
    char* osStem );

/*
 * Poor man's profiler. Use these to time code execution length. Pass a
 * message top be printed to the debug stream when the time is stopped.
 */

void xUtil_StartTimer ();
void xUtil_StopTimer ( char* isMessage );

/*
 * These are 'safe' replacements for the normal string routines. At least,
 * I think they are safe. The functions without n use sizeof on the buffer
 * and then call the n sized functions. After each n function copies in
 * normally, it sets the last character in the buffer to null. If it wasn't
 * null, it was print out a debugging note for an overflow alert.
 */

void xUtil_strcpy   ( char* ipDest, char* ipSrc );
void xUtil_strncpy  ( char* ipDest, char* ipSrc, int inSize );
void xUtil_snprintf ( char* ipDest, int inSize, char* isFormat, ... );
void xUtil_sprintf  ( char* ipDest, char* isFormat, ... );



char* xUtil_GetErrorString( xUtil_tErr iCode );


/*
 * For doing ctrl-c to cancel long operations but not quit the program.
 */
void xUtil_InitializeUserCancel ();
void xUtil_StartListeningForUserCancel ();
void xUtil_StopListeningForUserCancel ();
int xUtil_DidUserCancel ();
void xUtil_HandleUserCancelCallback (int signal);
