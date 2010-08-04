/**
 * @file  minmax.hpp
 * @brief Min/Max monitoring class
 *
 */
/*
 * Original Authors: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/08/04 19:01:36 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2002-2010,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 *
 */

#ifndef MINMAX_HPP
#define MINMAX_HPP

#include <limits>
#include <iostream>

namespace Freesurfer {

  template<typename T>
  class MinMax {
  public:
    //! Default constructor
    MinMax( void ) : minVal(std::numeric_limits<T>::max()),
		     maxVal(std::numeric_limits<T>::min()) {};

    //! Accumulator
    void Accumulate( const T val ) {
      if( val > this->maxVal ) {
	this->maxVal = val;
      }
      if( val < this->minVal ) {
	this->minVal = val;
      }
    }

    //! Accessor for maximum
    T GetMax( void ) const {
      return( this->maxVal );
    }

    //! Accessor for minimum
    T GetMin( void ) const {
      return( this->minVal );
    }
			

  private:
    T minVal;
    T maxVal;

  };

}

template<typename T>
std::ostream& operator<<( std::ostream& os, const Freesurfer::MinMax<T>& me ) {

  os << "[ "
     << me.GetMin()
     << ", "
     << me.GetMax()
     << " ]";

  return( os );
}


#endif
