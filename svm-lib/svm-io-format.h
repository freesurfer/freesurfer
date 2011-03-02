/**
 * @file  svm-io-format.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.4 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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


////SVM-LIB////////////////////////////////////////////////////////////////
//
// Name: Vector
//
// This file defines a simple vector class: a 1D array of double elements
// in contiguous memory with binary and text i/o.
//
//  Polina Golland polina@ai.mit.edu
//
///////////////////////////////////////////////////////////////////////////



#ifndef __SVM_IO_FORMAT_H__
#define __SVM_IO_FORMAT_H__

#include <stdio.h>
#include <iostream>
#include <string.h>

inline const char* readFormat(int a) {
  return "%d ";
}

inline const char* writeFormat(int a) {
  return "%d ";
}


inline const char* readFormat(float a) {
  return "%f ";
}

inline const char* writeFormat(float a) {
  return "%12.8e ";
}


inline const char* readFormat(double a) {
  return "%le ";
}

inline const char* writeFormat(double a) {
  return "%20.16e ";
}


template <class P>
bool read(P& p, const char *fileName, bool binary = false) {
  FILE *f = ((strcmp(fileName,"-"))?fopen(fileName,"r"):stdin);
  if ( f == NULL ) {
    std::cerr << "I/O Error: could not open file " << fileName << "\n";
    return false;
  }

  bool retCode = p.read(f,binary);
  fclose(f);
  return retCode;
}

template <class P>
bool write(const P& p, const char *fileName,
           bool binary = false, bool appendFlag = false) {
  FILE *f = ((strcmp(fileName,"-"))?
             fopen(fileName,((appendFlag)?"a":"w")):stdout);
  if ( f == NULL ) {
    std::cerr << "I/O Error: could not open file " << fileName << "\n";
    return false;
  }

  bool retCode = p.write(f,binary);
  fclose(f);
  return retCode;
}

#endif //  __SVM_IO_FORMAT_H__



