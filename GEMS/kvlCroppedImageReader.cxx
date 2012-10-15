/**
 * @file  kvlCroppedImageReader.cxx
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
#include "kvlCroppedImageReader.h"


#include "vnl/vnl_matlab_read.h"
#include "vnl/vnl_matrix.h"
#include "itkImageFileReader.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkShrinkImageFilter.h"
#include "itkMGHImageIO.h"


namespace kvl
{

//
//
//
CroppedImageReader
::CroppedImageReader()
{
  m_Image = 0;
  m_Transform = TransformType::New();
  m_WorldToImageTransform = TransformType::New();

  m_ExtraFraction = 0.0;
  m_DownSamplingFactor = 1;

  for ( int i = 0; i < 3; i++ )
  {
    m_BoundingBoxSize[ i ] = 0;
  }

}




//
//
//
CroppedImageReader
::~CroppedImageReader()
{

}


//
//
//
CroppedImageReader::TransformType::Pointer
CroppedImageReader
::GetTransformOfFileName( const std::string& filename )
{
  // Construct the default transform
  TransformType::Pointer  transform = TransformType::New();


  // Try first reading as MGH. If that fails, try reading .mat file
  itk::MGHImageIO::Pointer  io = itk::MGHImageIO::New();
  if ( io->CanReadFile( filename.c_str() ) )
  {
    std::cout << "Constructing image-to-world transform from MGH header information ("
              << filename << ")" << std::endl;

    io->SetFileName( filename.c_str() );
    io->ReadImageInformation();

    // Reconstruct the image-to-world tranform of the image
    TransformType::MatrixType  scale;
    TransformType::OffsetType  offset;
    itk::Matrix< float, 3, 3 >  direction;
    for ( int i = 0; i < 3; i++ )
    {
      scale[ i ][ i ] = io->GetSpacing( i );
      offset[ i ] = io->GetOrigin( i );

      const std::vector< double >  axis = io->GetDirection( i );
      for ( int j = 0; j< 3; j++ )
      {
        direction[ j ][ i ] = axis[ j ];
      }
    }

    transform->SetMatrix( direction * scale );
    transform->SetOffset( offset );
  }
  else
  {
    // Construct .mat file name
    const std::string::size_type  it = filename.find_last_of( "." );
    std::string  matFileName( filename, 0, it );
    matFileName += ".mat";


    // Try to read the contents of the .mat filename
    vcl_ifstream in( matFileName.c_str() );
    vnl_matlab_readhdr Mh( in );
    if ( Mh )
    {
      std::cout << "Slurping up contents of .mat filename " << matFileName << std::endl;

      if ( ( Mh.rows() == 4 ) && ( Mh.cols() == 4 ) )
      {
        vnl_matrix<double>  matlabMatrix( 4, 4 );
        if ( Mh.read_data( matlabMatrix.data_array() ) )
        {
          std::cout << "Using contents of .mat file!" << std::endl;

          TransformType::ParametersType  parameters( 12 );
          for ( unsigned int row = 0; row < 3; row++ )
          {
            for ( unsigned int col = 0; col < 3; col++ )
            {
              parameters[ row * 3 + col ] = matlabMatrix( row, col );
            }
            parameters[ 9 + row ] = matlabMatrix( row, 3 );
          }

          transform->SetParameters( parameters );

        }
      }
    }
    else
    {
      std::cout << "Couldn't reconstruct image-to-world transform for file " << filename << std::endl;
    }

  }


  std::cout << "Returning the following transform: " << std::endl;
  transform->Print( std::cout );

  return transform;
}



//
//
//
void
CroppedImageReader
::Read( const char* fileName, const char* boundingFileName )
{

  // Read the first image
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  ReaderType::Pointer  reader = ReaderType::New();
  reader->SetFileName( fileName );
  reader->Update();

  m_OriginalImageOriginalRegion = reader->GetOutput()->GetBufferedRegion();

#if 0
  // Make sure to unset spacing and origin as we don't look at it,
  // but VTK sure does!
  const float spacing[] = { 1.0f, 1.0f, 1.0f };
  const float origin[] = { 0.0f, 0.0f, 0.0f };
  reader->GetOutput()->SetSpacing( spacing );
  reader->GetOutput()->SetOrigin( origin );
#endif

  // If bounding file name is given, obtain a cropped image, and
  // construct the transform between the bounding image and the
  // cropped imagae
  if ( !boundingFileName )
  {
    m_Image = reader->GetOutput();
  }
  else
  {


    // The transformation is obtained from the image-to-world transforms of
    // the two images
    TransformType::Pointer  transform = TransformType::New();
    this->GetTransformOfFileName( fileName )->GetInverse( transform );
    m_Transform = this->GetTransformOfFileName( boundingFileName );
    m_Transform->Compose( transform );
    std::cout << "Transform before cropping: " << std::endl;
    m_Transform->Print( std::cout );

#if 0
    // Shift the inverse transform above by one voxel to map SPM's (1,1,1)-origin
    // coordinate system into (0,0,0)-based image grid
    m_WorldToImageTransform->Compose( transform );
    TransformType::OutputVectorType  originShift;
    originShift[ 0 ] = -1;
    originShift[ 1 ] = -1;
    originShift[ 2 ] = -1;
    m_WorldToImageTransform->Translate( originShift );
    std::cout << "WorldToImageTransform before cropping: " << std::endl;
    m_WorldToImageTransform->Print( std::cout );
#else
    m_WorldToImageTransform->Compose( transform );
    std::cout << "WorldToImageTransform before cropping: " << std::endl;
    m_WorldToImageTransform->Print( std::cout );
#endif

    // Obtain the bounding box size
    itk::ImageIOBase::Pointer   imageIO =
      itk::ImageIOFactory::CreateImageIO( boundingFileName, itk::ImageIOFactory::ReadMode );
    if ( imageIO.IsNull() )
    {
      itkExceptionMacro( << "Can't find a suitable filter to import image " << boundingFileName );
    }
    imageIO->SetFileName( boundingFileName );
    imageIO->ReadImageInformation();
    for ( int i = 0; i < 3; i++ )
    {
      m_BoundingBoxSize[ i ] = imageIO->GetDimensions( i );
    }
    std::cout << " m_BoundingBoxSize: [" << m_BoundingBoxSize[ 0 ] << ", "
              << m_BoundingBoxSize[ 1 ] << ","
              << m_BoundingBoxSize[ 2 ] << "]" << std::endl;



    // Map each of the corners of the bounding box, and record minima and maxima
    float  minimalMappedCoordinate[ 3 ];
    float  maximalMappedCoordinate[ 3 ];
    for ( int i = 0; i < 3; i++ )
    {
      minimalMappedCoordinate[ i ] = itk::NumericTraits< float >::max();
      maximalMappedCoordinate[ i ] = itk::NumericTraits< float >::min();
    }
    std::vector< TransformType::InputPointType >  meshCornerPoints;
    TransformType::InputPointType  corner;
    corner[ 2 ] = 0;
    corner[ 0 ] = 0;
    corner[ 1 ] = 0;
    meshCornerPoints.push_back( corner );
    corner[ 0 ] = m_BoundingBoxSize[ 0 ]-1;
    corner[ 1 ] = 0;
    meshCornerPoints.push_back( corner );
    corner[ 0 ] = m_BoundingBoxSize[ 0 ]-1;
    corner[ 1 ] = m_BoundingBoxSize[ 1 ]-1;
    meshCornerPoints.push_back( corner );
    corner[ 0 ] = 0;
    corner[ 1 ] = m_BoundingBoxSize[ 1 ]-1;
    meshCornerPoints.push_back( corner );
    corner[ 2 ] = m_BoundingBoxSize[ 2 ]-1;
    corner[ 0 ] = 0;
    corner[ 1 ] = 0;
    meshCornerPoints.push_back( corner );
    corner[ 0 ] = m_BoundingBoxSize[ 0 ]-1;
    corner[ 1 ] = 0;
    meshCornerPoints.push_back( corner );
    corner[ 0 ] = m_BoundingBoxSize[ 0 ]-1;
    corner[ 1 ] = m_BoundingBoxSize[ 1 ]-1;
    meshCornerPoints.push_back( corner );
    corner[ 0 ] = 0;
    corner[ 1 ] = m_BoundingBoxSize[ 1 ]-1;
    meshCornerPoints.push_back( corner );

    for ( std::vector< TransformType::InputPointType >::const_iterator  it = meshCornerPoints.begin();
          it != meshCornerPoints.end(); ++it )
    {
      TransformType::OutputPointType  mappedCorner =  m_Transform->TransformPoint( *it );
      for ( int i = 0; i < 3; i++ )
      {
        if ( mappedCorner[ i ] < minimalMappedCoordinate[ i ] )
        {
          minimalMappedCoordinate[ i ] = mappedCorner[ i ];
        }

        if ( mappedCorner[ i ] > maximalMappedCoordinate[ i ] )
        {
          maximalMappedCoordinate[ i ] = mappedCorner[ i ];
        }

      }
    }


    std::cout << " minimalMappedCoordinate: [" << minimalMappedCoordinate[ 0 ] << ", "
              << minimalMappedCoordinate[ 1 ] << ","
              << minimalMappedCoordinate[ 2 ] << "]" << std::endl;

    std::cout << " maximalMappedCoordinate: [" << maximalMappedCoordinate[ 0 ] << ", "
              << maximalMappedCoordinate[ 1 ] << ","
              << maximalMappedCoordinate[ 2 ] << "]" << std::endl;


    // Given the boundaries, crop the original image obtained from the reader
    // Loop over the image, and find the bounding box with content
    int  min[ 3 ];
    for ( int i = 0; i < 3; i++ )
    {
      min[ i ] = static_cast< int >( minimalMappedCoordinate[ i ] );
    }
    int  max[ 3 ];
    for ( int i = 0; i < 3; i++ )
    {
      max[ i ] = static_cast< int >( maximalMappedCoordinate[ i ] + 1 );
    }


    // Add some margin as a fraction of the found bounding box
    for ( int i = 0; i < 3; i++ )
    {
      const int  margin = static_cast< int >( m_ExtraFraction * ( max[ i ] - min[ i ] ) );

      min[ i ] -= margin;
      //if ( min[ i ] < 0 )
      //  min[ i ] = 0;

      max[ i ] += margin;
      //if ( max[ i ] > ( reader->GetOutput()->GetBufferedRegion().GetSize()[i] - 1 ) )
      //  max[ i ] = ( reader->GetOutput()->GetBufferedRegion().GetSize()[i] - 1 );
    }

    std::cout << "Cropping with min [" << min[ 0 ] << "  " << min[ 1 ] << "  " << min[ 2 ] << "]" << std::endl;
    std::cout << "          and max [" << max[ 0 ] << "  " << max[ 1 ] << "  " << max[ 2 ] << "]" << std::endl;


    // Gets the image regions
    ImageType::IndexType  originalImageIndex = {{ min[ 0 ], min[ 1 ], min[ 2 ] }};
    ImageType::SizeType  originalImageSize;
    for ( int i = 0; i < 3; i++ )
    {
      originalImageSize[ i ] = max[ i ] - min[ i ] + 1;
    }
    m_OriginalImageRegion = ImageType::RegionType( originalImageIndex, originalImageSize );

    ImageType::IndexType  croppedImageIndex = {{ 0, 0, 0 }};
    ImageType::SizeType  croppedImageSize = originalImageSize;
    m_CroppedImageRegion = ImageType::RegionType( croppedImageIndex, croppedImageSize );

    // Create an empty image, and fill it up with zeroes
    m_Image = ImageType::New();
    m_Image->SetRegions( m_CroppedImageRegion );
    m_Image->Allocate();
    m_Image->FillBuffer( 0 );

    // Copy the intensities by looping over the voxels of interest. Make sure you don't visit non-existing
    // voxels if the bounding box is bigger than the actual image...
    std::cout << "m_OriginalImageRegion before checking margins: " << m_OriginalImageRegion << std::endl;
    std::cout << "m_CroppedImageRegion before checking margins: " << m_CroppedImageRegion << std::endl;
    for ( int i = 0; i < 3; i++ )
    {
      const int  endMargin = ( reader->GetOutput()->GetBufferedRegion().GetSize()[i] - 1 ) -
                             ( m_OriginalImageRegion.GetIndex( i ) + m_OriginalImageRegion.GetSize( i ) - 1 );
      if ( endMargin < 0 )
      {
        m_OriginalImageRegion.SetSize( i, m_OriginalImageRegion.GetSize( i ) + endMargin );
        m_CroppedImageRegion.SetSize( i, m_CroppedImageRegion.GetSize( i ) + endMargin );
      }

      const int  startMargin = m_OriginalImageRegion.GetIndex( i );
      if ( startMargin < 0 )
      {
        m_OriginalImageRegion.SetIndex( i, m_OriginalImageRegion.GetIndex( i ) - startMargin );
        m_OriginalImageRegion.SetSize( i, m_OriginalImageRegion.GetSize( i ) + startMargin );

        m_CroppedImageRegion.SetIndex( i, m_CroppedImageRegion.GetIndex( i ) - startMargin );
        m_CroppedImageRegion.SetSize( i, m_CroppedImageRegion.GetSize( i ) + startMargin );
      }

    }
    std::cout << std::endl;
    std::cout << "m_OriginalImageRegion after checking margins: " << m_OriginalImageRegion << std::endl;
    std::cout << "m_CroppedImageRegion after checking margins: " << m_CroppedImageRegion << std::endl;



    itk::ImageRegionConstIterator< ImageType >  originalImageIterator( reader->GetOutput(), m_OriginalImageRegion );
    itk::ImageRegionIterator< ImageType >  croppedImageIterator( m_Image, m_CroppedImageRegion );
    for ( ; !originalImageIterator.IsAtEnd(); ++originalImageIterator, ++croppedImageIterator )
    {
      croppedImageIterator.Value() = originalImageIterator.Value();
    }


    // Finally, also adjust the transform to reflect the cropping
    TransformType::OutputVectorType  translation;
    translation[ 0 ] = -min[ 0 ];
    translation[ 1 ] = -min[ 1 ];
    translation[ 2 ] = -min[ 2 ];
    m_Transform->Translate( translation );
    std::cout << "Transform after cropping: " << std::endl;
    m_Transform->Print( std::cout );

    m_WorldToImageTransform->Translate( translation );
    std::cout << "WorldToImageTransform after cropping: " << std::endl;
    m_WorldToImageTransform->Print( std::cout );
  } // End test if bounding file name is given



  // Downsample image if necessary
  if ( m_DownSamplingFactor != 1 )
  {
    // Downsample image
    typedef itk::ShrinkImageFilter< ImageType, ImageType >  ShrinkerType;
    ShrinkerType::Pointer  shrinker = ShrinkerType::New();
    shrinker->SetInput( m_Image );
    shrinker->SetShrinkFactors( m_DownSamplingFactor );
    shrinker->Update();
    m_Image = shrinker->GetOutput();

    // Make sure to unset spacing and origin as we don't look at it,
    // but VTK sure does!
    const float spacing[] = { 1.0f, 1.0f, 1.0f };
    const float origin[] = { 0.0f, 0.0f, 0.0f };
    m_Image->SetSpacing( spacing );
    m_Image->SetOrigin( origin );

    // Adjust transformation
    m_Transform->Scale( 1 / static_cast< float >( m_DownSamplingFactor ) );
    std::cout << "Transform after scaling: " << std::endl;
    m_Transform->Print( std::cout );

    m_WorldToImageTransform->Scale( 1 / static_cast< float >( m_DownSamplingFactor ) );
    std::cout << "WorldToImageTransform after scaling: " << std::endl;
    m_WorldToImageTransform->Print( std::cout );
  }

}



} // end namespace kvl
