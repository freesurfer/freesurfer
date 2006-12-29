/**
 * @file  Array2.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:13 $
 *    $Revision: 1.2 $
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

  for ( int n = 0; n < this->mzX * this->mzY; n++) {
    this->mArray[n]= iInitValue;
  }
}

template <class T>
void
Array2<T>::SetAll ( T const iValue ) {
  for ( int n = 0; n < this->mzX * this->mzY; n++) {
    this->mArray[n]= iValue;
  }
}
