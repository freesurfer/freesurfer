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
 *    $Author: kteich $
 *    $Date: 2007/10/11 21:45:43 $
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

