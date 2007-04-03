/**
 * @file  colortab.h
 * @brief color table utilities
 *
 * An entry in a color table has:
 *   1. string name
 *   2. rgb (in both int and float)
 * An annotation is is an int packed with the int values of the
 * rgb in the first 3 bytes of the annotation int.
 */
/*
 * Original Authors: Kevin Teich, Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2007/04/03 21:17:52 $
 *    $Revision: 1.16 $
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


#ifndef COLORTAB_H
#define COLORTAB_H

#include "const.h"
#include "stdio.h"

/* A color table entry. The index of the entry in the table itself is
   the structure index. */
typedef struct
{
  char  name[STRLEN];  /* Structure name. */
  int   ri, gi, bi, ai;  /* 0-255 range. */
  float rf, gf, bf, af;   /* 0-1 range.  */
}
COLOR_TABLE_ENTRY, CTE;

/* The color table. Has an arry of entries. A structure index maps to
   a color table index. The table can be sparse, there need not be an
   entry for every structure. An entry may be NULL if the entry
   doesn't exist in the table. */
typedef struct
{
  CTE   **entries;  /* Array of CTE ptrs; some may be NULL */
  int   nentries;  /* Size of entries array. */
  char  fname[STRLEN];  /* Original file name. */
  int   version;  /* Version number, if read from binary */
}
COLOR_TABLE, CT ;


/* Reads a color table from an ASCII file. In this format, each line
   contains a structure index, an rgb triple in 0-255 format, and then
   a transparency integer, which is 255-alpha. Allocates the entire
   table. */
COLOR_TABLE *CTABreadASCII(char *fname);

/* Reads and writes a table to and from a binary stream. In this
   format, a negative version number is written (see special case in
   code) and then binary information for each entry. */
COLOR_TABLE *CTABreadFromBinary(FILE *fp);
int         CTABwriteIntoBinary(COLOR_TABLE *ct, FILE *fp);

/* Allocates an empty table. Here, all entries are allocated and
   filled with random colors. */
COLOR_TABLE *CTABalloc(int nentries) ;

/* Delete a CTAB. */
int         CTABfree(COLOR_TABLE **pct) ;

/* Returns a deep copy of the table. */
COLOR_TABLE *CTABdeepCopy(COLOR_TABLE *ct);

/* Copy the file name. */
int CTABcopyFileName(COLOR_TABLE *ct, char *name, size_t name_len);

/* Returns the number of non-NULL, valid entries. */
int CTABgetNumberOfValidEntries(COLOR_TABLE *ct, int *num);

/* Returns the total number of entries, including NULL ones. */
int CTABgetNumberOfTotalEntries(COLOR_TABLE *ct, int *num);

/* Returns whether or not an entry is valid. */
int CTABisEntryValid(COLOR_TABLE *ct, int index, int *valid);

/* Returns the integer or floating point color values or color + alpha
   at a given index. Returns an error code if the index is out of
   bounds. */
int CTABrgbAtIndexi(COLOR_TABLE *ct, int index,
                    int*r, int*g, int*b);
int CTABrgbAtIndexf(COLOR_TABLE *ct, int index,
                    float*r, float*g, float*b);
int CTABrgbaAtIndexi(COLOR_TABLE *ct, int index,
                     int*r, int*g, int*b, int*a);
int CTABrgbaAtIndexf(COLOR_TABLE *ct, int index,
                     float*r, float*g, float*b, float*a);

/* Gets the name of a structure. name should be an allocated string
   and name_len should be its size. */
int CTABcopyName(COLOR_TABLE *ct, int index, char *name, size_t name_len);

/* Returns the special annotation format version of an index. This is
   used in surface fields as a compressed color value. The r, g, and b
   values are shifted into a single integer. */
int CTABannotationAtIndex(COLOR_TABLE *ct, int index, int *annot);
int CTABfindAnnotation(COLOR_TABLE *ctab, int annotation, int *index);

/* Searches the table for duplicate annotation values, returns the count. */
int CTABfindDuplicateAnnotations(COLOR_TABLE *ct);

/* Searches through the table and finds the first entry with the given
   information. Sets index to the found index. Sets it to -1 if the
   color wasn't found. */
int CTABfindRGBi(COLOR_TABLE *ct, int r, int g, int b, int*index);
int CTABfindName(COLOR_TABLE *ct, char *name, int*index);

/* Print the color table to a file stream. This is a valid format for
   the ASCII file, but is not guaranteed to be the same as the
   original file. */
int CTABprintASCII(COLOR_TABLE *ct, FILE *fp);

/* Print the color table to a file. */
int CTABwriteFileASCII(COLOR_TABLE *ct, char *fname);

#endif
