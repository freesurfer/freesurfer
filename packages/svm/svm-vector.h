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
  E* data() const {
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


