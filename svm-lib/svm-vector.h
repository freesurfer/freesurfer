/**
 * @file  svm-vector.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:17 $
 *    $Revision: 1.3 $
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



#ifndef __SVM_VECTOR_H__
#define __SVM_VECTOR_H__

#include <stdio.h>
#include <iostream>

template <class E>
class Vector {
  int _size;
  E *_ptr;

  // Memory manipulation for vector re-sizing. These two functions
  // are private, as we don't want the caller to access them directry,
  // but rather use the public init function.

  void deallocate() {
    if ( _size > 0 ) {
      _size = 0;
      delete [] _ptr;
    }
  }

  void allocate(int size) {
    _size = size;
    _ptr = new E[_size];
  }


public:

  // Constructors and destructors
  Vector() :
      _size(0), _ptr(NULL) {}

  Vector ( int size ) {
    allocate(size);
  }

  Vector(int size, const E *ptr) {
    init(size,ptr);
  }

  ~Vector() {
    deallocate();
  }



  // Public memory allocation function
  void init(int size) {
    if ( _size != size ) {
      deallocate();
      allocate(size);
    }
  }
  void init(int size, const E* ptr) {
    init(size);
    for ( int i = 0; i < _size; i++ )
      _ptr[i] = ptr[i];
  }


  // Assignment
  Vector& operator = ( const Vector& v) {
    if ( _ptr == v._ptr )
      return *this;

    init(v.size());
    for ( int i = 0; i < _size; i++ )
      _ptr[i] = v._ptr[i];
    return *this;
  }


  // Accessors
  int size() const {
    return _size;
  };
  E* const data() const {
    return _ptr;
  };

  // Element accessors
  E& operator [] (int i) {
    return _ptr[i];
  };
  const E& operator [] (int i) const {
    return _ptr[i];
  };


  // I/O

  bool read(FILE *f, bool binary = false) {
    if ( binary )
      return ( fread(_ptr,sizeof(E),_size,f) - _size == 0 );

    for ( int i = 0; i < _size; i++ )
      if ( fscanf(f, readFormat(E(0)), &(_ptr[i])) != 1 )
        return false;
    return true;
  }

  bool write (FILE *f, bool binary = false) const {
    if ( binary )
      return ( fwrite(_ptr,sizeof(E),_size,f) - _size == 0 );

    for ( int i = 0; i < _size; i++ )
      fprintf(f, writeFormat(E(0)), _ptr[i]);
    fprintf(f,"\n");
    return true;
  }

};


#endif //  __SVM_VECTOR_H__


