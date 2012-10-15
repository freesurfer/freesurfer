/**
 * @file  kvlImageSmoother.cxx
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:39 $
 *    $Revision: 1.3 $
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
#include "kvlImageSmoother.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "vnl/algo/vnl_svd.h"
#if 0
#include "vnl/algo/vnl_qr.h"
#else
#include "vnl_kvl_qr.h"
#endif


namespace kvl
{

//
//
//
ImageSmoother
::ImageSmoother()
{
  m_MaskImage = 0;
  m_Image = 0;
  m_WeightImage = 0;
  m_PolynomialOrder = 2;
  m_PolynomialOrderUsed = m_PolynomialOrder;

}




//
//
//
ImageSmoother
::~ImageSmoother()
{

}



//
//
//
void
ImageSmoother
::Initialize() const
{

  // Check if we need to initialize
  if ( !m_BasisFunctions.empty() )
  {
    return;
  }


  // Sanity check on input
  if ( !m_MaskImage )
  {
    itkExceptionMacro( << "No mask image set; can't initialize" );
  }

  std::cout << "Initalizing" << std::endl;

  // Get standard polynomial basis functions
  vnl_matrix< double >  nonOrthogonalizedBasisFunctions = ImageSmoother::GetNonOrthogonalizedBasisFunctions( m_MaskImage, m_PolynomialOrder );


  // Orthogonalize basis functions
  if ( m_PolynomialOrder > 0 )
  {
    std::cout << "Orthogonalizing basis functions..." << std::flush;
    vnl_qr< double >  qr( nonOrthogonalizedBasisFunctions );
    //std::cout << nonOrthogonalizedBasisFunctions.rows() << std::endl;
    //std::cout << nonOrthogonalizedBasisFunctions.cols() << std::endl;
    //std::cout << nonOrthogonalizedBasisFunctions << std::endl;
    m_BasisFunctions = qr.Q( true );
    m_Orthogonalizer = qr.R( true );
  }
  else
  {
    // Dunno why, but VXL's qr decomposition segfaults for nx1 matrix
    m_BasisFunctions = nonOrthogonalizedBasisFunctions;
    m_Orthogonalizer = vnl_matrix< double >( 1, 1, 1.0f );
  }

  std::cout << "done!"<< std::endl;

}




//
//
//
vnl_matrix< double >
ImageSmoother
::GetNonOrthogonalizedBasisFunctions( const MaskImageType*  maskImage, int polynomialOrder,
                                      const int* indicesWhereUnityIsReached )
{

  const int  numberOfBasisFunctions = ( polynomialOrder + 1 ) *
                                      ( polynomialOrder + 2 ) / 2 *
                                      ( polynomialOrder + 3 ) / 3;
  std::cout << "numberOfBasisFunctions: " << numberOfBasisFunctions << std::endl;
  int  numberOfVoxelsInMask = 0;
  for ( itk::ImageRegionConstIterator< MaskImageType >  it( maskImage, maskImage->GetBufferedRegion() );
        !it.IsAtEnd(); ++it )
  {
    if ( it.Value() )
    {
      numberOfVoxelsInMask++;
    }
  }
  std::cout << "numberOfVoxelsInMask: " << numberOfVoxelsInMask << std::endl;


  // If no indicesWhereUnityIsReached is given, (x,y,z) will run between -1 (first index)
  // and 1 (last index) of the image grid. However, we can modify the latter behavior by
  // specifying explicitly providing indicesWhereUnityIsReached, specifying at what index
  // 1 is reached
  itk::Index< 3 >  indexWhereUnityIsReached;
  for ( int i = 0; i < 3; i++ )
  {
    indexWhereUnityIsReached[ i ] = maskImage->GetBufferedRegion().GetIndex()[ i ] +
                                    ( maskImage->GetBufferedRegion().GetSize()[ i ] - 1 );
  }
  if ( indicesWhereUnityIsReached )
  {
    for ( int i = 0; i < 3; i++ )
    {
      indexWhereUnityIsReached[ i ] = indicesWhereUnityIsReached[ i ];
    }
  }
  std::cout << "indexWhereUnityIsReached " << indexWhereUnityIsReached << std::endl;


  //
  std::vector< double >  x( numberOfVoxelsInMask );
  std::vector< double >  y( numberOfVoxelsInMask );
  std::vector< double >  z( numberOfVoxelsInMask );
  std::vector< double >::iterator  xIt = x.begin();
  std::vector< double >::iterator  yIt = y.begin();
  std::vector< double >::iterator  zIt = z.begin();
  const double  xOffset = maskImage->GetBufferedRegion().GetIndex()[ 0 ] +
                          ( indexWhereUnityIsReached[ 0 ] - maskImage->GetBufferedRegion().GetIndex()[ 0 ] ) / 2.0f;
  const double  xSlope = 2.0f / ( indexWhereUnityIsReached[ 0 ] - maskImage->GetBufferedRegion().GetIndex()[ 0 ] );
  const double  yOffset = maskImage->GetBufferedRegion().GetIndex()[ 1 ] +
                          ( indexWhereUnityIsReached[ 1 ] - maskImage->GetBufferedRegion().GetIndex()[ 1 ] ) / 2.0f;
  const double  ySlope = 2.0f / ( indexWhereUnityIsReached[ 1 ] - maskImage->GetBufferedRegion().GetIndex()[ 1 ] );
  const double  zOffset = maskImage->GetBufferedRegion().GetIndex()[ 2 ] +
                          ( indexWhereUnityIsReached[ 2 ] - maskImage->GetBufferedRegion().GetIndex()[ 2 ] ) / 2.0f;
  const double  zSlope = 2.0f / ( indexWhereUnityIsReached[ 2 ] - maskImage->GetBufferedRegion().GetIndex()[ 2 ] );
  for ( itk::ImageRegionConstIteratorWithIndex< MaskImageType >  it( maskImage, maskImage->GetBufferedRegion() );
        !it.IsAtEnd(); ++it )
  {
    if ( !it.Value() )
    {
      continue;
    }

    *xIt = xSlope * ( it.GetIndex()[ 0 ] - xOffset );
    *yIt = ySlope * ( it.GetIndex()[ 1 ] - yOffset );
    *zIt = zSlope * ( it.GetIndex()[ 2 ] - zOffset );

    ++xIt;
    ++yIt;
    ++zIt;
  }

  std::cout << "Successfully put up x, y, and z" << std::endl;


  std::cout << "Constructing basis functions" << std::endl;
  vnl_matrix< double >  nonOrthogonalizedBasisFunctions( numberOfVoxelsInMask, numberOfBasisFunctions );
  int  basisFunctionNumber = 0;
  for ( int order = 0; order <= polynomialOrder; order++ )
  {
    for ( int xOrder = 0; xOrder <= order; xOrder++ )
    {
      for ( int yOrder = 0; yOrder <= order-xOrder; yOrder++ )
      {
        const int zOrder = order - yOrder - xOrder;

        std::cout << "      (xOrder, yOrder, zOrder): (" << xOrder << ", "
                  << yOrder << ", "
                  << zOrder << ")" << std::endl;

        // Loop over all voxels
        std::vector< double >::const_iterator  xIt = x.begin();
        std::vector< double >::const_iterator  yIt = y.begin();
        std::vector< double >::const_iterator  zIt = z.begin();
        for ( int i = 0; i < numberOfVoxelsInMask; ++i, ++xIt, ++yIt, ++zIt )
        {
          nonOrthogonalizedBasisFunctions( i, basisFunctionNumber ) = pow( *xIt, xOrder ) *
              pow( *yIt, yOrder ) *
              pow( *zIt, zOrder );
        }

        basisFunctionNumber++;
      }
    }

  } // End loop over all orders
  //std::cout << "nonOrthogonalizedBasisFunctions:\n" << nonOrthogonalizedBasisFunctions << std::endl;

  return nonOrthogonalizedBasisFunctions;
}




//
//
//
ImageSmoother::ImageType::Pointer
ImageSmoother::GetSmoothedImage( const MaskImageType* maskImage, int voxelSizeFactor ) const
{

  // Sanity check on input
  if ( !m_Image || !m_MaskImage )
  {
    itkExceptionMacro( << "No image or mask image set" );
  }


  // Estimate the parameters
  this->EstimateParameters();

  // Expand into image
  ImageType::Pointer  smoothedImage = 0;
  if ( !maskImage )
  {
    smoothedImage = this->ExpandPolynomialToImage( m_BasisFunctions, m_MaskImage );
  }
  else
  {
#if 0
    // Create a uniform mask image from scratch, making sure that the outer voxels coincide
    // exactly with the outer voxels in the original mask image. This is important so that
    // the polynomials later on are exactly in the same reference frame
    MaskImageType::SizeType  size;
    for ( int i = 0; i < 3; i++ )
    {
      size[ i ] = upsamplingFactor * ( m_MaskImage->GetBufferedRegion().GetSize()[ i ] - 1 ) + 1;
    }
    MaskImageType::Pointer  maskImage = MaskImageType::New();
    maskImage->SetRegions( size );
    maskImage->Allocate();
    maskImage->FillBuffer( true );
#endif

    // Tricky stuff here: we need to make sure (x,y,z) reach the value 1 at exactly the
    // correct spot...
    int indicesWhereUnityIsReached[ 3 ];
    for ( int i = 0; i < 3; i++ )
    {
      indicesWhereUnityIsReached[ i ] = voxelSizeFactor * ( m_MaskImage->GetBufferedRegion().GetIndex()[ i ] +
                                        ( m_MaskImage->GetBufferedRegion().GetSize()[ i ] - 1 ) );
    }

    // Create non-orthogonized basis functions
    vnl_matrix< double >  nonOrthogonalizedBasisFunctions = ImageSmoother::GetNonOrthogonalizedBasisFunctions( maskImage, m_PolynomialOrder,
        indicesWhereUnityIsReached );

    // Orthogonalize using the same m_Orthogonalizer as used for m_BasisFunctions
    vnl_matrix< double >  basisFunctions = nonOrthogonalizedBasisFunctions * vnl_matrix_inverse< double >( m_Orthogonalizer );

    // Now expand into image
    smoothedImage = this->ExpandPolynomialToImage( basisFunctions, maskImage );
  }

  return smoothedImage;
}




//
//
//
ImageSmoother::ImageType::Pointer
ImageSmoother
::ExpandPolynomialToImage( const vnl_matrix< double >& basisFunctions, const MaskImageType* maskImage ) const
{

  // Now that we have the parameters, fill in the smoothed image
  // by combining the basis functions accordinlgy
  vnl_vector< double >  smoothedR = basisFunctions.extract( basisFunctions.rows(), m_Parameters.size() ) * m_Parameters;

  // Create an empty image and fill in the appropriate voxels
  ImageType::Pointer  smoothedImage = ImageType::New();
  smoothedImage->SetRegions( maskImage->GetBufferedRegion() );
  smoothedImage->Allocate();
  smoothedImage->FillBuffer( 0 );

  itk::ImageRegionConstIterator< MaskImageType >  maskIt( maskImage, maskImage->GetBufferedRegion() );
  itk::ImageRegionIterator< ImageType >  smoothedIt( smoothedImage, smoothedImage->GetBufferedRegion() );
  for ( int rowNumber = 0; !maskIt.IsAtEnd(); ++maskIt, ++smoothedIt )
  {
    if ( !maskIt.Value() )
    {
      continue;
    }

    smoothedIt.Value() = smoothedR( rowNumber );
    rowNumber++;
  }

  return smoothedImage;
}


//
//
//
void
ImageSmoother
::EstimateParameters() const
{

  // Make sure we have basis functions constructed
  this->Initialize();


  // Make sure we always have a weight image. If none is provided
  // by the user, create a default one
  ImageType::ConstPointer  weightImage = m_WeightImage;
  if ( !weightImage )
  {
    ImageType::Pointer  defaultWeightImage = ImageType::New();
    defaultWeightImage->SetRegions( m_MaskImage->GetBufferedRegion() );
    defaultWeightImage->Allocate();
    defaultWeightImage->FillBuffer( 1 );

    weightImage = defaultWeightImage;
  }


  // Determine how many of the basis functions we're actually going to use
  const int  numberOfBasisFunctionsUsed = ( m_PolynomialOrderUsed + 1 ) *
                                          ( m_PolynomialOrderUsed + 2 ) / 2 *
                                          ( m_PolynomialOrderUsed + 3 ) / 3;
  std::cout << "estimating polynomial parameters up to order " << m_PolynomialOrderUsed << std::endl;
  std::cout << "       numberOfBasisFunctionsUsed: " << numberOfBasisFunctionsUsed << std::endl;


  // Set up the system of equations LHS * parameters = RHS
  // where LHS = ( A' * W * A ) and RHS = A' * W * r
  // with A matrix with basis functions
  //      W diagonal matrix with weights
  //      r vector with intensities to be smoothed
  //std::cout << "Setting up linear system" << std::endl;
  vnl_matrix< double >  A = m_BasisFunctions.extract( m_BasisFunctions.rows(), numberOfBasisFunctionsUsed );
  vnl_matrix< double >  WtimesA = A;
  vnl_vector< double >  r( m_BasisFunctions.rows() );

  itk::ImageRegionConstIterator< MaskImageType >  maskIt( m_MaskImage, m_MaskImage->GetBufferedRegion() );
  itk::ImageRegionConstIterator< ImageType >  weightIt( weightImage, weightImage->GetBufferedRegion() );
  itk::ImageRegionConstIterator< ImageType >  it( m_Image, m_Image->GetBufferedRegion() );
  for ( int rowNumber = 0; !maskIt.IsAtEnd(); ++maskIt, ++weightIt, ++it )
  {
    if ( !maskIt.Value() )
    {
      continue;
    }

    WtimesA.scale_row( rowNumber, weightIt.Value() );
    r( rowNumber ) = it.Value();
    rowNumber++;
  }

  vnl_matrix< double >  LHS = A.transpose() * WtimesA;
  vnl_vector< double >  RHS = WtimesA.transpose() * r;
  //std::cout << "LHS:\n" << LHS << std::endl;
  //std::cout << "RHS:\n" << RHS << std::endl;


  // To solve M * x = y for x, do
  //
  // vnl_matrix< double >  M;
  // vnl_svd< double >  svd( M )
  // vnl_vector<T>  y;
  // vnl_vector< double > x = svd.solve( y );
  //
  // A short cut to this is
  //
  // x = vnl_matrix_inverse< double >( M ) * y;
  //
  // which under the hood does exactly the same
  vnl_svd< double >  svd( LHS );
  m_Parameters = svd.solve( RHS );
  //std::cout << "m_Parameters: " << m_Parameters << std::endl;

}


} // end namespace kvl
