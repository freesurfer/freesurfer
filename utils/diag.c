/*
 *       FILE NAME:   diag.c
 *
 *       DESCRIPTION: diagnostic routines
 *
 *       AUTHOR:      Bruce Fischl
 *       DATE:        2/5/96
 *
*/

/*-----------------------------------------------------
                    INCLUDE FILES
-------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>

#include "image.h"
#include "diag.h"
#include "window.h"
#include "proto.h"

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/

/*-----------------------------------------------------
                    GLOBAL VARIABLES
-------------------------------------------------------*/

unsigned long  Gdiag = 0 ;

/*-----------------------------------------------------
                    GLOBAL FUNCTIONS
-------------------------------------------------------*/

/*------------------------------------------------------------------------
       Parameters:

      Description:
  
    Return Values:
        diag bits.

------------------------------------------------------------------------*/
unsigned long
DiagInit(void)
{
  char *cp ;

  cp = getenv("diag") ;

  if (!cp) cp = getenv("DIAG") ;
  if (cp) 
  {
    sscanf(cp, "0x%lx", &Gdiag) ;
    fprintf(stderr, "diagnostics set to 0x%lx\n", Gdiag) ;
  }

#if 0
  if (getenv("logging") != NULL)
    DebugOpenLogFile("trace.log") ;
  
  if (getenv("flushing") != NULL)
    flushing = 1 ;
#endif
  return(Gdiag) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
DiagShowImage(unsigned long diag_bits, int win, int which, IMAGE *I, 
              char *fmt, ...)
{
  char    name[100] ;
  va_list args ;
  
  if (!(diag_bits & Gdiag))
    return(-1) ;

  va_start(args, fmt) ;
/*  fmt = va_arg(args, char *) ;*/
  vsprintf(name, fmt, args) ;

  if (win == DIAG_NEW_WINDOW)    /* create a new window for this image */
  {
  }

  WinShowImage(win, I, which) ;
  WinSetName(win, which, name) ;
  if (diag_bits & DIAG_WAIT)     /* wait for a keystroke before continuing */
    fgetc(stdin) ;    

  return(win) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
DiagDrawBox(unsigned long diag_bits, int win, int row0, int col, 
                 int rows, int cols, int color)
{
   int i;
   i=win; i=row0; i=col; i=rows; i=cols; i=color; /* prevent warning (dng) */

  if (!(diag_bits & Gdiag))
    return(-1) ;

  if (diag_bits & DIAG_WAIT)     /* wait for a keystroke before continuing */
    fgetc(stdin) ;    

  return(0) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
DiagCloseWindow(unsigned long diag_bits, int win)
{
   int i;
   i=(int)diag_bits; i=win;

  return(0) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
DiagCreateWindow(unsigned long diag_bits, int wrows, int wcols,
                 int rows,int cols)
{
  int win ;
   int i;
   i=wcols; i=(int)diag_bits;

/*
  create a set of rows x cols windows, each one of which is wrows x wcols
  in size.
*/
#if 0
  win = WinAlloc("window", 50,50, wcols,wrows) ;
  WinShow(win) ;
#else
  win = WinCreate("window", 1, wrows, rows, cols) ;
#endif

  return(win) ;
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
DiagPrintf(unsigned long diag_bits, char *fmt, ...)
{
  va_list args ;
  
  if (!(diag_bits & Gdiag))
    return(-1) ;

  va_start(args, fmt) ;
  return(vprintf(fmt, args)) ;
}
