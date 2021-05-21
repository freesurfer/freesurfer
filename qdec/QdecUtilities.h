/**
 * @brief Misc utilities
 *
 * A bunch of misc utilities for Qdec.
 */
/*
 * Original Author: Kevin Teich
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

#ifndef QdecUtilities_h
#define QdecUtilities_h

#include <stdio.h>
#include <unistd.h>
#ifndef MAXPATHLEN
#include <sys/param.h>
#endif

#include <string>

using namespace std;

class QdecUtilities
{

public:

  // Calls IsFileReadable and throws an error if it fails.
  static void AssertFileIsReadable ( string const& ifn );

  // Returns true if a file exists and is openable with read
  // permissions.
  static bool IsFileReadable ( string const& ifn );

  // extract the path name from a file name and return a pointer to it
  static const char *FileNamePath(const char *fname, const char *pathName);

  // read the value of specified key from resource file .Qdecrc
  static const char *GetQdecrcResourceString(const char *key);

 private:

  static const char *GetResourceString(const char *key, 
                                       const char *filename);
};

#endif
