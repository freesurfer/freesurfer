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
// Name: Matrix
//
// This file defines a simple matrix class: a 2D array of double elements
// ech row stored in a contiguous memory segment.
//
//
//  Polina Golland polina@ai.mit.edu
//
///////////////////////////////////////////////////////////////////////////


#ifndef __SVM_MATRIX_H__
#define __SVM_MATRIX_H__

#include <stdio.h>

#include "svm-vector.h"


template <class E>
class Matrix {
  int _rows, _cols;
  E **_ptr;

  // Private memory manipulation functions. The caller should use
  // init to resize the matrix.
  void allocate (int rowCount, int colCount) {
    _rows = rowCount;
    _cols = colCount;
    _ptr = new E*[_rows];
    for ( int i = 0; i < _rows; i++ )
      _ptr[i] = new E[_cols];
  }

  void deallocate () {
    if ( _rows > 0 ) {
      for ( int i = 0; i < _rows; i++ )
        delete [] _ptr[i];
      delete [] _ptr;
      _rows = _cols = 0;
    }
  }



public:

  // Constructors and desctructors
Matrix() : _rows(0), _cols(0), _ptr(NULL) {}
  Matrix(int rowCount, int colCount) {
    allocate(rowCount,colCount);
  }
  Matrix(int rowCount, int colCount, const E* const* ptr) {
    init(rowCount,colCount,ptr);
  }

  ~Matrix() {
    deallocate();
  }

  // Public memory allocation function
  void init(int rowCount, int colCount) {
    if ( _rows != rowCount || _cols != colCount ) {
      deallocate();
      allocate(rowCount,colCount);
    }
  }

  void init(int rowCount, int colCount, const E* const* ptr) {
    init(rowCount,colCount);
    for ( int i = 0; i < _rows; i++ )
      for ( int j = 0; j < _cols; j++ )
        _ptr[i][j] = ptr[i][j];
  }


  // Assignment
  Matrix& operator = (const Matrix &m) {
    if ( _ptr == m._ptr )
      return *this;

    init(m.rows(),m.cols());
    for ( int i = 0; i < _rows; i++ )
      for ( int j = 0; j < _cols; j++ )
        _ptr[i][j] = m._ptr[i][j];

    return *this;
  }


  // Row access
  E* operator [] (int i ) {
    return _ptr[i];
  }

  E* operator [] (int i) const {
    return _ptr[i];
  }

  E** data() const {
    return _ptr;
  };


  void setRow(int i, const E* v) {
    for ( int j = 0; j < _cols; j++ )
      _ptr[i][j] = v[j];
  }

  void setRow(int i, const Vector<E>& v) {
    for ( int j = 0; j < _cols; j++ )
      _ptr[i][j] = v[j];
  }

  void setCol(int i, const E* v) {
    for ( int j = 0; j < _rows; j++ )
      _ptr[j][i] = v[j];
  }

  void setCol(int i, const Vector<E>& v) {
    for ( int j = 0; j < _rows; j++ )
      _ptr[j][i] = v[j];
  }


  // Size accessors
  int rows() const {
    return _rows;
  }
  int cols() const {
    return _cols;
  }


  // I/O

  bool read (FILE *f, bool binary = false) {
    if ( binary ) {
      for ( int i = 0; i < _rows; i++ )
        if ( fread(_ptr[i],sizeof(E),_cols,f) - _cols != 0 )
          return false;
    } else {
      for ( int i = 0; i < _rows; i++ )
        for ( int j = 0; j < _cols; j++ )
          if ( fscanf(f, readFormat(E(0)), &(_ptr[i][j])) != 1 )
            return false;
    }
    return true;
  }

  bool write (FILE *f, bool binary = false) const {
    if ( binary ) {
      for ( int i = 0; i < _rows; i++ )
        if ( fwrite(_ptr[i],sizeof(E),_cols,f) - _cols != 0 )
          return false;
    } else {
      for ( int i = 0; i < _rows; i++ ) {
        for ( int j = 0; j < _cols; j++ )
          fprintf(f,writeFormat(E(0)), _ptr[i][j]);
        fprintf(f,"\n");
      }
    }
    return true;
  }
};




#endif //  __SVM_MATRIX_H__




