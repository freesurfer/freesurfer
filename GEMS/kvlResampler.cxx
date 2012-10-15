/**
 * @file  kvlResampler.cxx
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:40 $
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
#include "kvlRegisterer.h"

#include "itkAffineTransform.h"


namespace kvl
{

//
//
//
Registerer
::Registerer()
{
  m_FixedImage = 0;
  m_MovingImage = 0;

  m_Transform = TransformType::New();

}




//
//
//
Registerer
::~Registerer()
{

}



//
//
//
void
Registerer
::StartRegistration()
{

}



//
//
//
void
Registerer
::ApplyParameters( std::vector< ImageType::Pointer > images ) const
{

  std::cout << "Applying transformation: " << m_Transform << std::endl;

  // Loop over all images
  for ( std::vector< ImageType::Pointer >::iterator  it = images.begin();
        it != images.end(); ++it )
  {
    // Get the image
    ImageType::Pointer  image = *it;

    // Get the original image-to-world tranform of the image
    typedef itk::AffineTransform< double, 3 >  AffineTransformType;
    AffineTransformType::MatrixType  scale;
    AffineTransformType::OffsetType  offset;
    for ( int i = 0; i < 3; i++ )
    {
      scale[ i ][ i ] = image->GetSpacing()[ i ];
      offset[ i ] = image->GetOrigin()[ i ];
    }
    AffineTransformType::Pointer  original = AffineTransformType::New();
    original->SetMatrix( image->GetDirection() * scale );
    original->SetOffset( offset );
    std::cout << "Original image-to-world: " << original << std::endl;

    // Pre-multiply it with the calculated transformation
    AffineTransformType::Pointer updated = AffineTransformType::New();
    updated->SetMatrix( m_Transform->GetMatrix() );
    updated->SetOffset( m_Transform->GetOffset() );
    updated->Compose( original, true );
    std::cout << "Updated image-to-world: " << updated << std::endl;

    // Now overwrite image information to get exactly that effect
    ImageType::PointType   newOrigin;
    ImageType::SpacingType  newSpacing;
    ImageType::DirectionType  newDirection;
    for ( int i = 0; i < 3; i++ )
    {
      // Offset part
      newOrigin[ i ] = updated->GetOffset()[ i ];

      // For every column, determine norm (which will be voxel spacing), and normalize direction
      double  normOfColumn = 0.0;
      for ( int j = 0; j < 3; j++ )
      {
        normOfColumn += pow( updated->GetMatrix()[ j ][ i ], 2 );
      }
      normOfColumn = sqrt( normOfColumn );
      newSpacing[ i ] = normOfColumn;
      for ( int j = 0; j < 3; j++ )
      {
        newDirection[ j ][ i ] = updated->GetMatrix()[ j ][ i ] / normOfColumn;
      }
    }
    image->SetOrigin( newOrigin );
    image->SetSpacing( newSpacing );
    image->SetDirection( newDirection );

  } // end loop over all images


}



} // end namespace kvl
