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
#include <iostream.h>
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
bool read(P& p, const char *fileName, bool binary = false) 
{
  FILE *f = ((strcmp(fileName,"-"))?fopen(fileName,"r"):stdin);
  if ( f == NULL ) {
    cerr << "I/O Error: could not open file " << fileName << "\n";
    return false;
  }

  bool retCode = p.read(f,binary);
  fclose(f);
  return retCode;
}

template <class P>  
bool write(const P& p, const char *fileName, 
	   bool binary = false, bool appendFlag = false) 
{
  FILE *f = ((strcmp(fileName,"-"))?
	     fopen(fileName,((appendFlag)?"a":"w")):stdout);
  if ( f == NULL ) {
    cerr << "I/O Error: could not open file " << fileName << "\n";
    return false;
  }

  bool retCode = p.write(f,binary);
  fclose(f);
  return retCode; 
}

#endif //  __SVM_IO_FORMAT_H__



