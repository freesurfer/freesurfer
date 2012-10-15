/**
 * @file  kvlPartialVolumeSegmentationTest.cxx
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
#include "itkImage.h"
#include "itkAutomaticTopologyMeshSource.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkImageFileWriter.h"
#include "kvlAtlasMesh.h"
#include "kvlAtlasMeshAlphaDrawer.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "kvlAtlasMeshCollection.h"
#include "kvlAtlasMeshVisitCounter.h"



int main( int argc, char* argv[] )
{

#if 0
  for ( kvl::AtlasMesh::PointsContainer::Iterator  pointIt = meshSource->GetOutput()->GetPoints()->Begin();
        pointIt != meshSource->GetOutput()->GetPoints()->End();
        ++pointIt )
  {
    kvl::AtlasMesh::PointType  p = pointIt.Value();

    kvl::AtlasAlphasType   alphas( 2 );
    float  squareDistanceToCenter =
      pow( p[ 0 ] - static_cast<float>( size[0] ), 2 ) +
      pow( p[ 1 ] - static_cast<float>( size[1] ), 2 );
    pow( p[ 2 ] - static_cast<float>( size[2] ), 2 );
    float  sigma = static_cast< float >( size[ 0 ] ) / 4.0f;
    alphas[ 0 ] = exp( -squareDistanceToCenter / 2 / sigma / sigma );
    alphas[ 1 ] = 1 - alphas[ 0 ];

    //std::cout << "Setting alpha of vertex with position " << p << " to " << alphas[ 0 ] << std::endl;
    kvl::AtlasMesh::PixelType  pointParameters;
    pointParameters.m_Alphas = alphas;
    meshSource->GetOutput()->SetPointData( pointIt.Index(), pointParameters );
  }


#endif



  // Create an empty image
  typedef itk::Image< unsigned char, 3 >  ImageType;
  ImageType::SizeType  size;
  size[ 0 ] = 10;
  size[ 1 ] = 10;
  size[ 2 ] = 10;
  ImageType::Pointer  image = ImageType::New();
  image->SetRegions( size );
  image->Allocate();

  // Create a mesh, in which there is exactly one vertex coinciding exactly
  // with each voxel in the image
  kvl::AtlasMeshCollection::Pointer  collection =  kvl::AtlasMeshCollection::New();
  unsigned int  meshSize[ 3 ];
  unsigned int  domainSize[ 3 ];
  for ( int i = 0; i < 3; i++ )
  {
    meshSize[ i ] = size[ i ];
    domainSize[ i ] = size[ i ];
  }
  const unsigned int  numberOfClasses = 2;
  const unsigned int  numberOfMeshes = 1;
  collection->Construct( meshSize, domainSize, 0.05,
                         numberOfClasses, numberOfMeshes );


  // Create a spherical "object" in the mesh
  const float  radius = static_cast< float >( size[ 0 ] ) / 4.0f;
  for ( kvl::AtlasMesh::PointsContainer::ConstIterator  pointIt = collection->GetReferencePosition()->Begin(),
        kvl::AtlasMesh::PointDataContainer::Iterator  paramIt = collection->GetPointParameters()->Begin();
        pointIt != collection->GetReferencePosition()->End();
        ++pointIt, ++paramIt )
  {
    kvl::AtlasMesh::PointType  p = pointIt.Value();

    kvl::AtlasAlphasType   alphas( 2 );
    const float  distanceToCenter = sqrt(
                                      pow( p[ 0 ] - static_cast<float>( size[0] ), 2 ) +
                                      pow( p[ 1 ] - static_cast<float>( size[1] ), 2 );
                                      pow( p[ 2 ] - static_cast<float>( size[2] ), 2 ) );
    if ( distanceToCenter < radius )
    {
      alphas[ 0 ] = 1;
    }
    else
    {
      alphas[ 0 ] = 0;
    }
    alphas[ 1 ] = 1 - alphas[ 0 ];

    //std::cout << "Setting alpha of vertex with position " << p << " to " << alphas[ 0 ] << std::endl;
    paramIt.Value().m_Alphas = alphas;
  }
  collection->Write( "meshToRasterize.txt" );



  // Rasterize mesh
  kvl::AtlasMeshAlphaDrawer::Pointer  alphaDrawer = kvl::AtlasMeshAlphaDrawer::New();
  alphaDrawer->SetLabelImage( image );
  alphaDrawer->SetLabelNumber( 1 );
  alphaDrawer->Rasterize( collection->GetReferenceMesh() );

  // Write out
  typedef itk::ImageFileWriter< kvl::AtlasMeshAlphaDrawer::AlphaImageType >  AlphaWriterType;
  AlphaWriterType::Pointer  alphaWriter = AlphaWriterType::New();
  alphaWriter->SetFileName( "testAlpha.mhd" );
  alphaWriter->SetInput( alphaDrawer->GetAlphaImage() );
  alphaWriter->Write();



#if 0
  // Now rasterize mesh again with VisitCounter
  kvl::AtlasMeshVisitCounter::Pointer  visitCounter = kvl::AtlasMeshVisitCounter::New();
  visitCounter->SetLabelImage( image );
  visitCounter->Rasterize( collection->GetReferenceMesh() );


  // Report any errors you find: voxels that have not been visited at all, and voxels that have been
  // visit more than once
  itk::ImageRegionConstIteratorWithIndex< kvl::AtlasMeshVisitCounter::CountImageType >  it( visitCounter->GetCountImage(),
      visitCounter->GetCountImage()->GetLargestPossibleRegion() );
  for( ; !it.IsAtEnd(); ++it )
  {
    if ( it.Value() != 1 )
    {
      // Don't warn on the sides. We know there is a problem there.
      if ( ( it.GetIndex()[ 0 ] != visitCounter->GetCountImage()->GetLargestPossibleRegion().GetSize()[ 0 ]-1 ) &&
           ( it.GetIndex()[ 1 ] != visitCounter->GetCountImage()->GetLargestPossibleRegion().GetSize()[ 1 ]-1 ) &&
           ( it.GetIndex()[ 2 ] != visitCounter->GetCountImage()->GetLargestPossibleRegion().GetSize()[ 2 ]-1 ) )
      {
        std::cout << "Got " << static_cast< int >( it.Value() ) << " visits in voxel with index " << it.GetIndex() << std::endl;
      }
    }

  }



  // Write out
  typedef itk::ImageFileWriter< kvl::AtlasMeshVisitCounter::CountImageType >  CountWriterType;
  CountWriterType::Pointer  countWriter = CountWriterType::New();
  countWriter->SetFileName( "testCount.mhd" );
  countWriter->SetInput( visitCounter->GetCountImage() );
  countWriter->Write();


  // Print out the images the old-fashed way
  std::cout << "Visting count: " << std::endl;
  ImageType::IndexType  index;
  for ( int z = 0; z < size[ 2 ]; z++ )
  {
    index[ 2 ] = z;
    for ( int y = 0; y < size[ 1 ]; y++ )
    {
      index[ 1 ] = y;
      for ( int x = 0; x < size[ 0 ]; x++ )
      {
        index[ 0 ] = x;
        std::cout << static_cast< int >( visitCounter->GetCountImage()->GetPixel( index ) ) << "  ";
      }
      std::cout << "\n" << std::endl;
    }

    std::cout << "\n\n\n" << std::endl;
  }


  std::cout << "\n\n\nAlpha image: " << std::endl;
  for ( int z = 0; z < size[ 2 ]; z++ )
  {
    index[ 2 ] = z;
    for ( int y = 0; y < size[ 1 ]; y++ )
    {
      index[ 1 ] = y;
      for ( int x = 0; x < size[ 0 ]; x++ )
      {
        index[ 0 ] = x;
        std::cout << alphaDrawer->GetAlphaImage()->GetPixel( index ) << "  ";
      }
      std::cout << "\n" << std::endl;
    }

    std::cout << "\n\n\n" << std::endl;
  }

#endif




  return 0;
};

