/**
 * @file  matfile.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:08:59 $
 *    $Revision: 1.9 $
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


#ifndef MATFILE_H
#define MATFILE_H

#ifdef Darwin
typedef unsigned char  Byte;
#endif

#include "znzlib.h"
#include "matrix.h"
#include "machine.h"


/*----------------------------------------------------------
  MATHD - Header for version 4 matrix data element (not for the file).
  ----------------------------------------------------------*/
typedef struct
{
  long32  type ;
  long32  mrows ;
  long32  ncols ;
  long32  imagf ;
  long32  namlen ;
}
MATHD ;

/*----------------------------------------------------------
  MATFILE - structure to store information about a matrix
  in a matfile (not the matfile itself).
  ----------------------------------------------------------*/
typedef struct
{
  /*common -- for version 4 and 5*/
  long32  type;   // precision (int, float)
  long32  mrows ; // number of rows of matrix
  long32  ncols ; // number of cols of matrix
  long32  imagf;  // 1 if imaginary
  long32  namlen ; // length of the matlab variable name
  int     version ; // 4 or 5
  char  *data ;     // buffer for real data
  char  *idata ;    // buffer for imaginary data
  char endian; // only for 5 (endian for 4 is in type)
}
MATFILE ;


/*----------------------------------------------------------
  MLFC -
  ----------------------------------------------------------*/
typedef struct
{
  char *mfile;    // Name of matfile
  int nvars;      // Number of variables in matfile
  char *varname[1000]; // Variable name (all variables)
  MATRIX *varmtx[1000];  // Matrices ???
}
MATFILECONTENTS, MLFC;

char *MatReadHeader0(FILE *fp, MATFILE *mf);
char    *MatReadHeader(FILE *fp, MATFILE *mf, long32 *compressed) ;
char    *znzMatReadHeader(FILE *fp, MATFILE *mf, char **data) ;
MATFILE *MatFileRead(const char *fname, int type) ;
MATRIX  *MatlabRead(const char *fname) ;
MATRIX  *MatlabRead2(const char *fname) ;
int     MatlabWrite(MATRIX *mat, const char *fname, char *name) ;
int     MatFileWrite(const char *fname,
                     float *data, int rows, int cols, char *name) ;
int Matlab_Install_printf( int (*new_printf)(const char *szFormat, ...) );
MLFC *ReadMatlabFileContents(const char *fname);
int   MLFCprint(FILE *fp, MLFC *mlfc);
int MLFCfree(MLFC **ppmlfc);
MATRIX *ReadMatlabFileVariable(char *fname, char *varname);


#define MAT_BYTE     0
#define MAT_DOUBLE   1
#define MAT_INT      2
#define MAT_SHORT    3
#define MAT_FLOAT    4

// Endianness for version 4
#define MATFILE_PC      0000
#define MATFILE_SPARC   1000

// Endianness for version 5
#define MATFILE_PC5     'I'
#define MATFILE_SPARC5  'M'

#define MATFILE_DOUBLE  00
#define MATFILE_FLOAT   10
#define MATFILE_LONG    20
#define MATFILE_SHORT   30
#define MATFILE_USHORT  40
#define MATFILE_BYTE    50

// ?????
#define MATFILE_FULL    0
#define MATFILE_TEXT    1
#define MATFILE_SPARSE  2


#endif

