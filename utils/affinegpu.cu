/**
 * @file  affinegpu.cu
 * @brief Holds affine transformation class for the GPU
 *
 * Holds an affine transformation type for the GPU
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:24 $
 *    $Revision: 1.6 $
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

#include <cstdlib>


#include <iostream>



#include "affinegpu.hpp"

// ====================================================

namespace GPU
{
namespace Classes
{

AffineTransformation::AffineTransformation( void )
{
  /*!
  Default constructor zeros the matrix
  */

  // Have a little sanity check
  if( kVectorSize != 4 )
  {
    std::cerr << __FUNCTION__
              << ": Incompatible universe detected"
              << std::endl
              << "Please adjust number of spatial dimensions "
              << "to 3 and try again"
              << std::endl;
    exit( EXIT_FAILURE );
  }

  for( unsigned int i=0; i<kMatrixSize; i++ )
  {
    this->matrix[i] = 0;
  }
}

AffineTransformation::AffineTransformation( const MATRIX* src )
{
  /*!
  Constructor to take values from a real 4x4 matrix.
  Double checks everything before performing assignment
  */

  if( src->type != MATRIX_REAL )
  {
    std::cerr << __FUNCTION__ << ": Invalid matrix type " <<
              src->type << std::endl;
    exit( EXIT_FAILURE );
  }

  if( static_cast<unsigned int>(src->rows) != kVectorSize )
  {
    std::cerr << __FUNCTION__ << ": Invalid number of rows " <<
              src->rows << std:: endl;
    exit( EXIT_FAILURE );
  }

  if( static_cast<unsigned int>(src->cols) != kVectorSize )
  {
    std::cerr << __FUNCTION__ << ": Invalid number of cols " <<
              src->cols << std:: endl;
    exit( EXIT_FAILURE );
  }

  // Do the copy
  for( unsigned int i=0; i<kVectorSize; i++ )
  {
    for( unsigned int j=0; j<kVectorSize; j++ )
    {
      this->operator()( i, j ) = *MATRIX_RELT( src, i+1, j+1 );
    }
  }
}



void AffineTransformation::SetTransform( const MATRIX* src )
{
  /*!
  Duplicates the functionality of the constructor from
  a real 4x4 matrix
  */
  if( src->type != MATRIX_REAL )
  {
    std::cerr << __FUNCTION__ << ": Invalid matrix type " <<
              src->type << std::endl;
    exit( EXIT_FAILURE );
  }

  if( static_cast<unsigned int>(src->rows) != kVectorSize )
  {
    std::cerr << __FUNCTION__ << ": Invalid number of rows " <<
              src->rows << std:: endl;
    exit( EXIT_FAILURE );
  }

  if( static_cast<unsigned int>(src->cols) != kVectorSize )
  {
    std::cerr << __FUNCTION__ << ": Invalid number of cols " <<
              src->cols << std:: endl;
    exit( EXIT_FAILURE );
  }

  // Do the copy
  for( unsigned int i=0; i<kVectorSize; i++ )
  {
    for( unsigned int j=0; j<kVectorSize; j++ )
    {
      this->operator()( i, j ) = *MATRIX_RELT( src, i+1, j+1 );
    }
  }
}

void AffineTransformation::SetTransform( const AffineMatrix* src )
{
  memcpy( this->matrix, src->mat, kMatrixSize*sizeof(float) );
}


}
}
