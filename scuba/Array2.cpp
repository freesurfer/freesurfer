/**
 * @file  Array2.cpp
 * @brief Simple class for a 2D array
 *
 * This is a simple class representing a 2D array, accessible by x,y
 * pairs or by a Point2<int>.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:36 $
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


#include "Array2.h"

template <class T>
Array2<T>::Array2 ( int izX, int izY ) {

  this->mzX = izX;
  this->mzY = izY;
  this->mArray = new T[this->mzX * this->mzY];
}

template <class T>
Array2<T>::Array2 ( int izX, int izY, T iInitValue ) {

  this->mzX = izX;
  this->mzY = izY;
  this->mArray = new T[this->mzX * this->mzY];

  // Set all values to the values we got.
  for ( int n = 0; n < this->mzX * this->mzY; n++)
    this->mArray[n]= iInitValue;
}

template <class T>
Array2<T>::~Array2 () {

  delete[] this->mArray;
};

template <class T>
void
Array2<T>::SetAll ( T const iValue ) {

  // Set all values to the value we got.
  for ( int n = 0; n < this->mzX * this->mzY; n++)
    this->mArray[n]= iValue;
}

