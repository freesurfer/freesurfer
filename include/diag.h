/*
 *       FILE NAME:   diag.h
 *
 *       DESCRIPTION: diagnostic routine prototypes
 *
 *       AUTHOR:      Bruce Fischl
 *       DATE:        2/5/96
 *
*/

#ifndef DIAG_H
#define DIAG_H

#include <stdarg.h>
#include "image.h"

int  DiagShowImage(unsigned long diag_bits, int win, int which,
                   IMAGE *I, char *fmt,...) ;
int  DiagDrawBox(unsigned long diag_bits, int win, int row0, int col, 
                 int rows, int cols, int color) ;
int  DiagCloseWindow(unsigned long diag_bits, int win) ;
int  DiagCreateWindow(unsigned long diag_bits, 
                      int wrows, int wcols, int rows, int cols) ;
unsigned long  DiagInit(char *fname, 
                  int (*vfprint)(FILE *fp, const char *fmt, va_list args),
                  int (*vprint)(const char *fmt, va_list args)) ;

int  DiagPrintf(unsigned long diag_bits, char *fmt, ...) ;
int  DiagFprintf(unsigned long diag_bits, char *fmt, ...) ;
void DiagBreak(void) ;  /* dummy for break points in debugger */
void DiagHeartbeat(float pct_done) ;

/* diagnostic codes */
#define DIAG_LP         0x00000001L   /* license plate detection diagnostics */
#define DIAG_EXM        0x00000002L
#define DIAG_TIMER      0x00000004L
#define DIAG_WRITE      0x00000008L
#define DIAG_MATRIX     0x00000010L
#define DIAG_KERNELS    0x00000020L
#define DIAG_SHOW       0x00000040L
#define DIAG_INPUT      0x00000080L
#define DIAG_MOVIE      0x00000100L
#define DIAG_LOGMAP     0x00000200L
#define DIAG_LOG_QUAD   0x00000400L
#define DIAG_CONGRAPH   0x00000800L
#define DIAG_BACKPROP   0x00001000L
#define DIAG_LOGDIFF    0x00002000L
#define DIAG_HEARTBEAT  0x00004000L

#define DIAG_WAIT    0x08000000L
#define DIAG_VERBOSE 0x10000000L   /* allows 2 levels for each type */

/* supported colors */
#define DIAG_BLACK   -1
#define DIAG_RED     0
#define DIAG_WHITE   1
#define DIAG_GREEN   2
#define DIAG_BLUE    3
#define DIAG_CYAN    4
#define DIAG_YELLOW  5
#define DIAG_PURPLE  6

/* misc. stuff */
#define DIAG_NEW_WINDOW   -1


extern unsigned long Gdiag ;    /* global diagnostic flag */

#define PT_NONE           0
#define PT_CIRCLE         1

#define PLOT_LINES_SOLID          0
#define PLOT_LINES_DASHED         1
#define PLOT_LINES_DOUBLE_DASHED  2
#define PLOT_IMPULSES             3
#define PLOT_CIRCLES              4

#define PLOT_INT        0
#define PLOT_DOUBLE     1
#define PLOT_CHAR       2
#define PLOT_SEQUENTIAL 3
#define NEW_WINDOW     -1

#define PLOT_HORIZONTAL 0
#define PLOT_VERTICAL   1

#define BLACK  0
#define WHITE  1
#define FLIP   2
#define RED    3
#define GREEN  4



#endif





