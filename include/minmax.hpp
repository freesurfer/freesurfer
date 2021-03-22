/**
 * @brief Min/Max monitoring class
 *
 */
/*
 * Original Authors: Richard Edgar
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

#ifndef MINMAX_HPP
#define MINMAX_HPP

#include <limits>
#include <iostream>

namespace Freesurfer
{

template<typename T>
class MinMax
{
public:
  //! Default constructor
  MinMax( void ) : minVal(std::numeric_limits<T>::max()),
    maxVal(std::numeric_limits<T>::min()),
    nAcc(0), intTotal(0), doubleTotal(0) {};

  //! Accumulator
  void Accumulate( const T val )
  {
    if( val > this->maxVal )
    {
      this->maxVal = val;
    }
    if( val < this->minVal )
    {
      this->minVal = val;
    }

    nAcc++;
    intTotal += val;
    doubleTotal += val;
  }

  //! Accessor for maximum
  T GetMax( void ) const
  {
    return( this->maxVal );
  }

  //! Accessor for minimum
  T GetMin( void ) const
  {
    return( this->minVal );
  }

  //! Accessor for number accumulated
  size_t GetN( void ) const
  {
    return( this->nAcc );
  }

  //! Accessor for integer running total
  long long GetIntTotal( void ) const
  {
    return( this->intTotal );
  }

  //! Accessor for double running total
  double GetDoubleTotal( void ) const
  {
    return( this->doubleTotal );
  }


private:
  T minVal;
  T maxVal;

  size_t nAcc;
  long long intTotal;
  double doubleTotal;
};

}

template<typename T>
std::ostream& operator<<( std::ostream& os, const Freesurfer::MinMax<T>& me )
{

  os << "[ "
     << me.GetMin()
     << ", "
     << me.GetMax()
     << " ]";

  return( os );
}


#endif
