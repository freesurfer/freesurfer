/**
 * @brief diagnostic routines
 *
 */
/*
 * Original Author: Bruce Fischl
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

/*-----------------------------------------------------
                    INCLUDE FILES
-------------------------------------------------------*/
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "const.h"
#include "diag.h"
#include "error.h"
#include "fsinit.h"
#include "image.h"
#include "proto.h"
#include "windiag.h"

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/

/*-----------------------------------------------------
                    GLOBAL VARIABLES
-------------------------------------------------------*/

FILE *Gstdout;
FILE *Gstderr;
FILE *Gstdin;
FILE *Gdiag_fp = NULL;

int Gprofile = 0;
int Gvx = -1;
int Gvy = -1;
int Gvz = -1;

int Gsx = -1;
int Gsy = -1;
int Gsz = -1;

unsigned long Gdiag = 0;
int Gdiag_no = -1;
int Gdiag_no2 = -1;
int Gx = -1;
int Gy = -1;
int Gz = -1;
int Gx2 = -1;
int Gy2 = -1;
int Gz2 = -1;
#define DEFAULT_IMAGE_SIZE 512
int IMAGE_SIZE = DEFAULT_IMAGE_SIZE;

/*-----------------------------------------------------
                     STATIC DATA
-------------------------------------------------------*/

static char diag_fname[STRLEN] = "diag.log";
static int (*diag_vprintf)(const char *fmt, va_list args) = vprintf;
static int (*diag_vfprintf)(FILE *fp, const char *fmt, va_list args) = vfprintf;

/*-----------------------------------------------------
                    GLOBAL FUNCTIONS
-------------------------------------------------------*/

void assertFailed(const char* file, int line, const char* tst) {
    fprintf(stdout, "ASSERTION FAILED: %s:%d %s\n", file, line, tst);
    fprintf(stderr, "ASSERTION FAILED: %s:%d %s\n", file, line, tst);
    fflush(stdout);
    fflush(stderr);
    
    while (!!getenv("FREESURFER_HANG_ON_ASSERT")) {
      fprintf(stdout, "FREESURFER_HANG_ON_ASSERT\n");
      fprintf(stderr, "FREESURFER_HANG_ON_ASSERT\n");
      fflush(stdout);
      fflush(stderr);
      static volatile long count;
      count = 0;
      while (count < 4000000000L) count++;
    }
    
    *(int*)(-1) = 0;
}


/*------------------------------------------------------------------------
       Parameters:

      Description:

    Return Values:
        diag bits.

------------------------------------------------------------------------*/
unsigned long DiagInit(char *fname,
                       int (*vfprint)(FILE *fp, const char *fmt, va_list args),
                       int (*vprint)(const char *fmt, va_list args))
{
  char *cp = 0;
  unsigned long diag = 0L;

  FSinit();
  Gstdout = stdout;
  Gstdin = stdin;
  Gstderr = stderr;  // for use in gdb
  if (fname) strcpy(diag_fname, fname);
  if (vfprint) diag_vfprintf = vfprint;
  if (vprint) diag_vprintf = vprint;

  cp = getenv("DIAG_NO");
  if (cp) Gdiag_no = atoi(cp);

  cp = getenv("IMAGE_SIZE");
  if (cp) IMAGE_SIZE = atoi(cp);

  cp = getenv("DIAGX");
  if (cp) Gx = atoi(cp);
  cp = getenv("DIAGY");
  if (cp) Gy = atoi(cp);
  cp = getenv("DIAGZ");
  if (cp) Gz = atoi(cp);

  cp = getenv("DIAGVX");
  if (cp) Gvx = atoi(cp);
  cp = getenv("DIAGVY");
  if (cp) Gvy = atoi(cp);
  cp = getenv("DIAGVZ");
  if (cp) Gvz = atoi(cp);

  cp = getenv("PROFILE");
  if (cp) {
    Gprofile = atof(cp);
    printf("turning profiling diagnostics on (%d)...\n", Gprofile);
  }
  cp = getenv("diag");
  if (!cp) cp = getenv("DIAG");
  if (cp) {
    sscanf(cp, "0x%lx", &diag);
    Gdiag |= diag;
  }

  if (getenv("DIAG_VERBOSE")) Gdiag |= DIAG_VERBOSE;
  if (Gdiag & DIAG_VERBOSE) printf("verbose diagnostics enabled...\n");

#if 0
  if (Gdiag)
    fprintf(stderr, "diagnostics set to 0x%lx\n", Gdiag) ;
#endif

#if 0
  if (getenv("logging") != NULL)
    DebugOpenLogFile("trace.log") ;

  if (getenv("flushing") != NULL)
    flushing = 1 ;
#endif
  return (Gdiag);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int DiagShowImage(unsigned long diag_bits, int win, int which, IMAGE *I, char *fmt, ...)
{
#if 0
  char    name[100] ;
  va_list args ;

  if (diag_bits && !(diag_bits & Gdiag))
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
#endif
  return (win);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int DiagDrawBox(unsigned long diag_bits, int win, int row0, int col, int rows, int cols, int color)
{
  // int i;
  // i = win;
  // i = row0;
  // i = col;
  // i = rows;
  // i = cols;
  // i = color; /* prevent warning (dng) */



  if (diag_bits && !(diag_bits & Gdiag)) return (-1);

  if (diag_bits & DIAG_WAIT) /* wait for a keystroke before continuing */
    fgetc(stdin);

  return (0);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int DiagCloseWindow(unsigned long diag_bits, int win)
{
  if (diag_bits && !(diag_bits & Gdiag)) return (-1);

  win = win + (int)diag_bits; /* to remove warning */
  return (0);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int DiagCreateWindow(unsigned long diag_bits, int wrows, int wcols, int rows, int cols)
{
  int win = -1;
#if 0
  int i;
  i=wcols;

  if (diag_bits && !(diag_bits & Gdiag))
    return(-1) ;

  /*
    create a set of rows x cols windows, each one of which is wrows x wcols
    in size.
  */
#if 0
  win = WinAlloc("window", 50,50, wcols,wrows) ;
  WinShow(win) ;
#else
  win = WinCreate("window", 1, wrows, wcols, rows, cols) ;
#endif
#endif
  return (win);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int DiagFprintf(unsigned long diag_bits, char *fmt, ...)
{
  static int first = 1;
  va_list args;
  FILE *fp;
  int len;

  if (diag_bits && !(diag_bits & Gdiag)) return (-1);

  if (first)
    fp = fopen(diag_fname, "w");
  else
    fp = fopen(diag_fname, "a+");
  first = 0;
  va_start(args, fmt);
  len = vfprintf(fp, fmt, args);
  fclose(fp);
  return (len);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int DiagPrintf(unsigned long diag_bits, const char *fmt, ...)
{
  va_list args;

  if (diag_bits && !(diag_bits & Gdiag)) return (-1);

  va_start(args, fmt);
  return ((*diag_vprintf)(fmt, args));
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          dummy for break points in debugger
------------------------------------------------------*/
void DiagBreak(void) {}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          dummy for break points in debugger
------------------------------------------------------*/
void DiagHeartbeat(float pct_done)
{
  static float old_pct = -10.0f;

#if 0
  if ((old_pct < 0.0f) || (old_pct > pct_done) || (pct_done < 0.00001f))
    fprintf(stderr, "\n") ;
  old_pct = pct_done ;
  if (Gdiag & DIAG_HEARTBEAT)
    fprintf(stderr, "%2.1f%% finished     \r",100.0f*pct_done);
  if (pct_done >= 0.999f)
    fprintf(stderr, "\n") ;
#else
  if (pct_done < old_pct) old_pct = -20;
  if (pct_done - old_pct > .10) {
    old_pct = pct_done;
    if (Gdiag & DIAG_HEARTBEAT) fprintf(stderr, "%2.1f%% finished\n", 100.0f * pct_done);
  }
#endif
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          dummy for break points in debugger
------------------------------------------------------*/
void DiagShowPctDone(float pct_done, int nprints)
{
  static float old_pct = -10.0f;

  if (pct_done < old_pct) old_pct = -20;
  if (pct_done - old_pct > 1 / (float)nprints) {
    old_pct = pct_done;
    if (Gdiag & DIAG_HEARTBEAT) fprintf(stderr, "%2.1f%% finished\n", 100.0f * pct_done);
  }
}

int check_finite(const char *where, double what)
{
  if (!std::isfinite(what)) {
    ErrorPrintf(ERROR_BADPARM, "%s not finite!\n", where);
    DiagBreak();
    return (0);
  }
  return (1);
}
