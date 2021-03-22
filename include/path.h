/*
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


#ifndef PATHS_INCLUDED
#define PATHS_INCLUDED

#include <stdio.h>
#ifdef Darwin
// Mac OS X doesnt have gnulib, so it doesn't have getline,
// so the source for getline is found here in utils and built on the Mac
#include "getline.h"
#endif

#include "label.h"

/* One of these for every point in the path. */
typedef struct
{
  float x;
  float y;
  float z;
  int vno;
}
PATH_POINT;

/* Main path struct. */
typedef struct
{
  int n_points;         /* number of points in this path. */
  char name[100];       /* original file name */
  PATH_POINT *points;   /* array of size n_points */
}
PATH;


/* Read in multiple paths from a path file. On return, num_read will
   point to the number of paths in the paths variable. Pass in the
   address to a PATH** variable. On return it will point to a list of
   allocated PATH objects. Returns an error code. */
int    PathReadMany (char *fname, int *num_read, PATH ***paths);

/* Write multiple paths to a file. paths should be an array of
   pointers to PATH structs, each one with a valid PATH
   object. num_paths should be the number of paths in that
   array. Returns an error code. */
int    PathWriteMany (char *fname, int num_paths, PATH **paths);

/* Allocate and free paths. PathAlloc returns the new path or NULL if
   there was an error. PathFree returns an error code.*/
PATH*  PathAlloc (int n_points, const char* name);
int    PathFree (PATH** path);

/* Returns whether or not a file is a path file. If an error occurs,
   it will just return 0. The file name version will open and close a
   file, and stream version will act on an open file and leave it
   open, and will change the filepos. */
int    PathIsPathFile (char* fname);     /* Will open and close file. */
int    PathIsPathFileStream (FILE* fp);  /* Will not close. */

/* Converts a single path into a single label. Unlike the label files
   created by mri_path2label that encode multiple paths into a single
   label file, this function just copies the coordinate positions into
   a new label. This will create a new label and return it; the caller
   is responsible for freeing it. */
int    PathConvertToLabel (PATH* path, LABEL** label);

/* Convert a label into a single path. Unlike the label files
   converted by mri_path2label that can encode multiple paths in a
   single label, this function just copies the coordinates in a label
   object into a new path. This will create a new path and return it;
   the caller is responsible for freeing it. */
int    PathCreateFromLabel (LABEL* label, PATH** path);

#endif
