/**
 * @file  kvlGenerateToyExampleData.cxx
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
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileWriter.h"
#include "kvlAtlasParameterEstimator.h"
#include "kvlAtlasMeshSmoother.h"
#include "vnl/vnl_sample.h"



int main( int argc,  char* argv[] )
{

  // Input check
  if ( argc < 2 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " imageSize [ partialVolumingFactorX partialVolumingFactorY partialVolumingFactorZ ]" << std::endl;
    exit( - 1 );
  }

  // Retrieve the size
  std::istringstream  imageSizeStream( argv[ 1 ] );
  int  imageSize;
  imageSizeStream >> imageSize;

  // Retrieve the partial volume model specifications
  int  partialVolumingFactors[ 3 ] = { 1, 1, 1 };
  if ( argc > 2 )
  {
    for ( int i = 0; i < 3; i++ )
    {
      std::istringstream  inputStream( argv[ 2 + i ] );
      inputStream >> partialVolumingFactors[ i ];
    }
  }


  // Create an empty image
  typedef itk::Image< unsigned short, 3 >  ImageType;
  ImageType::SizeType  size;
  for ( int i = 0; i < 3; i++ )
  {
    size[ i ] = imageSize * partialVolumingFactors[ i ];
  }
  ImageType::Pointer  image = ImageType::New();
  image->SetRegions( size );
  image->Allocate();

  // Draw a binary object in the image, where each of the two
  // labels is obscured by Gaussian noise with a certain mean and variance
  const float  insideMean = 110.0f;
  const float  outsideMean = 130.0f;
  const float  insideSigma = 1.0f;
  const float  outsideSigma = 1.0f;

  float  center[ 3 ];
  for ( int i = 0; i < 3; i++ )
  {
    center[ i ] = ( size[ i ] - 1 ) / 2.0f;
  }
  float  radius[ 3 ];
  radius[ 0 ] = size[ 0 ] / 4.5f;
  radius[ 1 ] = size[ 1 ] / 3.5f;
  radius[ 2 ] = size[ 2 ] / 5.5f;
  for ( itk::ImageRegionIteratorWithIndex< ImageType >  it( image, image->GetBufferedRegion() );
        !it.IsAtEnd(); ++it )
  {
    float  MahalanobisDistance = 0;
    for ( int i = 0; i < 3; i++ )
    {
      MahalanobisDistance += pow( ( it.GetIndex()[ i ] - center[ i ] ) / radius[ i ], 2 );
    }
    MahalanobisDistance = sqrt( MahalanobisDistance );

    if ( MahalanobisDistance < 1 )
    {
      it.Value() = static_cast< ImageType::PixelType >( vnl_sample_normal( insideMean, insideSigma ) );
    }
    else
    {
      it.Value() = static_cast< ImageType::PixelType >( vnl_sample_normal( outsideMean, outsideSigma ) );
    }

  }

  typedef itk::ImageFileWriter< ImageType >  WriterType;
  WriterType::Pointer  writer = WriterType::New();
  writer->SetFileName( "toyImage.mhd" );
  writer->SetInput( image );
  writer->Update();


  // Now generate a mesh representation of an atlas that we could use to play with
  ImageType::SizeType  reducedSize;
  for ( int i = 0; i < 3; i++ )
  {
    reducedSize[ i ] = static_cast< int >( size[ i ] * 0.8 );
  }
  for ( int i = 0; i < 3; i++ )
  {
    center[ i ] = ( reducedSize[ i ] - 1 ) / 2.0f;
  }
  typedef itk::Image< unsigned char, 3 >  LabelImageType;
  LabelImageType::Pointer  labelImage = LabelImageType::New();
  labelImage->SetRegions( reducedSize );
  labelImage->Allocate();
  labelImage->FillBuffer( 0 );

  const float  fixedRadius = reducedSize[ 0 ] / 3.5f;
  for ( itk::ImageRegionIteratorWithIndex< LabelImageType >  it( labelImage, labelImage->GetBufferedRegion() );
        !it.IsAtEnd(); ++it )
  {
    float  distance = 0;
    for ( int i = 0; i < 3; i++ )
    {
      distance += pow( it.GetIndex()[ i ] - center[ i ], 2 );
    }
    distance = sqrt( distance );

    if ( distance < fixedRadius )
    {
      it.Value() = 1;
    }

  }

  typedef itk::ImageFileWriter< LabelImageType >  LabelWriterType;
  LabelWriterType::Pointer  labelWriter = LabelWriterType::New();
  labelWriter->SetFileName( "inputImage.mhd" );
  labelWriter->SetInput( labelImage );
  labelWriter->Update();


  // Construct a mesh collection, so that each node coincides exactly with each pixel
  // in the image
  kvl::AtlasMeshCollection::Pointer  collection =  kvl::AtlasMeshCollection::New();
  unsigned int  meshSize[ 3 ];
  unsigned int  domainSize[ 3 ];
  for ( int i = 0; i < 3; i++ )
  {
    meshSize[ i ] = reducedSize[ i ];
    domainSize[ i ] = reducedSize[ i ];
  }
  const unsigned int  numberOfClasses = 2;
  const unsigned int  numberOfMeshes = 1;
  collection->Construct( meshSize, domainSize, 1000,
                         numberOfClasses, numberOfMeshes );
  collection->Write( "initialMesh.txt" );


  // Estimate the model parameters for this mesh collection. This should yield probabilies
  // of either 0 or 1 in each of the nodes
  kvl::AtlasParameterEstimator::Pointer  estimator = kvl::AtlasParameterEstimator::New();
  std::vector< LabelImageType::ConstPointer >  labelImages;
  labelImages.push_back( labelImage.GetPointer() );
  estimator->SetLabelImages( labelImages );
  estimator->SetInitialMeshCollection( collection );
  estimator->Estimate();
  //collection = estimator->GetCurrentMeshCollection();
  collection->SetK( 0.01 );
  collection->Write( "estimatedMesh.txt" );


  // Smooth a little
  const float  sigma = reducedSize[ 0 ] / 5.0f;
  kvl::AtlasMeshSmoother::Pointer  smoother = kvl::AtlasMeshSmoother::New();
  smoother->SetMeshCollection( collection );
  smoother->SetSigma( sigma );
  kvl::AtlasMeshCollection::Pointer  smoothedCollection = smoother->GetSmoothedMeshCollection();
  smoothedCollection->Write( "smoothed.txt" );


  // Now translate a little to the center, so that the mesh is centered in the originally sized data
  typedef kvl::AtlasMeshCollection::TransformType  TransformType;
  TransformType::Pointer  transform = TransformType::New();
  TransformType::OutputVectorType  translation;
  for ( int i = 0; i < 3; i++ )
  {
    translation[ i ] = ( size[ i ] - reducedSize[ i ] ) / 2.0f;
  }
  transform->Translate( translation );
  //transform->Print( std::cout );
  for ( int meshNumber = -1; meshNumber < static_cast< int >( numberOfMeshes ); meshNumber++ )
  {
    smoothedCollection->Transform( meshNumber, transform );
  }
  smoothedCollection->Write( "toyMesh.txt" );





  if ( argc > 2 )
  {
    // Create an empty partial volumed image
    ImageType::SizeType  partialVolumedSize;
    int  numberOfSubvoxels = 1;
    for ( int i = 0; i < 3; i++ )
    {
      partialVolumedSize[ i ] = imageSize;
      numberOfSubvoxels *= partialVolumingFactors[ i ];
    }
    ImageType::Pointer  partialVolumedImage = ImageType::New();
    partialVolumedImage->SetRegions( partialVolumedSize );
    partialVolumedImage->Allocate();
    partialVolumedImage->FillBuffer( 0 );

    // Fill in the partial volumed data
    for ( itk::ImageRegionIteratorWithIndex< ImageType >  it( partialVolumedImage, partialVolumedImage->GetBufferedRegion() );
          !it.IsAtEnd(); ++it )
    {
      std::cout << "Filling in partial volume voxel with index " << it.GetIndex() << std::endl;

      ImageType::IndexType  index;
      for ( int xStep = 0; xStep < partialVolumingFactors[ 0 ]; xStep++ )
      {
        index[ 0 ] = it.GetIndex()[ 0 ] * partialVolumingFactors[ 0 ] + xStep;
        for ( int yStep = 0; yStep < partialVolumingFactors[ 1 ]; yStep++ )
        {
          index[ 1 ] = it.GetIndex()[ 1 ] * partialVolumingFactors[ 1 ] + yStep;
          for ( int zStep = 0; zStep < partialVolumingFactors[ 2 ]; zStep++ )
          {
            index[ 2 ] = it.GetIndex()[ 2 ] * partialVolumingFactors[ 2 ] + zStep;

            std::cout << "         Adding contribution of high-res voxel with index " << index << std::endl;
            it.Value() += image->GetPixel( index );
          }
        }
      }

      it.Value() = static_cast< unsigned int >( it.Value() / static_cast< float >( numberOfSubvoxels ) );
    }

    // Write out the partial volumed image
    writer->SetFileName( "toyPartialVolumedImage.mhd" );
    writer->SetInput( partialVolumedImage );
    writer->Update();



    // If we're using a partial volume model, the mesh is in the image grid space of a
    // higher-resolution image than the one we're working with. Re-adjust the position
    // accordingly. Beside an obvious scaling, there is also a small translation involved:
    // the voxel at the origin of the low-res image grid is supposed to be an average of
    // the voxels close to the origin, but with positive coordinates, of the high-res grid.
    // In one dimension, if we're averaging the N voxels with high-res positions 0, 1, 2, ... N-1,
    // the center of the "averaging" voxel in the low-res grid is ( 0 + 1 + 2 + ... + ( N-1) ) / N,
    // which is [ N*(N-1)/2 ] / N (cf. primary school), which is (N-1)/2

    // Build the transform from low-res into high-res coordinates
    transform->SetIdentity();
    TransformType::OutputVectorType  translation;
    TransformType::OutputVectorType  scaling;
    for ( int i = 0; i < 3; i++ )
    {
      translation[ i ] = ( partialVolumingFactors[ i ] - 1 ) / 2.0f;
      scaling[ i ] = partialVolumingFactors[ i ];
    }
    transform->Scale( scaling );
    transform->Translate( translation );
    transform->Print( std::cout );

    // Apply the inverse of this transform to the mesh
    TransformType::Pointer  inverseTransform = TransformType::New();
    transform->GetInverse( inverseTransform );
    inverseTransform->Print( std::cout );
    for ( int meshNumber = -1; meshNumber < static_cast< int >( numberOfMeshes ); meshNumber++ )
    {
      smoothedCollection->Transform( meshNumber, inverseTransform );
    }
    smoothedCollection->Write( "toyPartialVolumedMesh.txt" );


  }


  return 0;
};

