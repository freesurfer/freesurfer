/**
 * @file  utils.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:00 $
 *    $Revision: 1.26 $
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
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


/*
   @(#)utils.h  1.9
   10/16/95
*/
/*------------------------------------------------------------------------
      File Name:  utils.h

         Author:  Bruce Fischl

        Created:  Jan. 1994

    Description:

------------------------------------------------------------------------*/
#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>

#define MATLAB_FILE   0
#define HIPS_FILE     1
#define LIST_FILE     2
#define UNKNOWN_FILE  3

double randomNumber(double low, double hi) ;
int    setRandomSeed(long seed) ;
double normAngle(double angle) ;
float deltaAngle(float angle1, float angle2) ;
double calcDeltaPhi(double phi1, double phi2) ;
#if 1
double latan2(double y, double x) ;
#else
#define latan2(y,x)  atan2(y,x)
#endif
float  angleDistance(float theta1, float theta2) ;
int    QuadEqual(double a1, double a2) ;
void   fComplementCode(double *pdIn, double *pdOut, int iLen) ;
#ifndef _HOME_
char *fgetl(char *s, int n, FILE *fp) ;
#endif

int IntSqrt(int n) ;

char *StrRemoveSpaces(char *str) ;
char *StrUpper(char *str) ;
char *StrLower(char *str) ;
char *StrSkipNumber(char *str) ;
char *StrReplace(char *src, char *dst, char csrc, int cdst) ;

char *FileNameOnly(char *str, char *fname) ;
char *FileNameFromWildcard(char *inStr, char *outStr) ;
int  FileExists(char *fname) ;
int  FileType(char *fname) ;
int  FileNumber(char *fname) ;
int  FileNumberOfEntries(char *fname) ;
char *FileName(char *full_name) ;
char *FileFullName(char *full_name) ;
char *FileTmpName(char *base) ;
char *FileTmpName(char *basename) ;
void FileRename(char *inName, char *outName) ;
char *FileNameAbsolute(char *fname, char *absFname) ;
char *FileNamePath(char *fname, char *pathName) ;
char *FileNameRemoveExtension(char *in_fname, char *out_fname) ;
char *FileNameExtension(char *fname, char *ext) ;
char *AppendString(char *src, char *app);

int devIsinf(float value);
int devIsnan(float value);
int devFinite(float value);

int getMemoryUsed(); // return total virtual memory used by Progname in Kbytes. works only under Linux /proc system
void printMemoryUsed(); // print function of the above.
char *strcpyalloc(char *str);
int ItemsInString(char *str);

#endif
