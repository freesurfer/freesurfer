/**
 * @file  kvlAtlasMeshRasterizorTest.cxx
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

  // Create an empty image to serve as a template for the alpha drawer
  typedef itk::Image< unsigned char, 3 >  ImageType;
  ImageType::SizeType  size;
  size[ 0 ] = 100;
  size[ 1 ] = 100;
  size[ 2 ] = 100;
  ImageType::Pointer  image = ImageType::New();
  image->SetRegions( size );
  image->Allocate();
  image->FillBuffer( 0 );

  // Use a mesh source to create the mesh
  typedef itk::AutomaticTopologyMeshSource< kvl::AtlasMesh >  MeshSourceType;
  MeshSourceType::Pointer  meshSource = MeshSourceType::New();
#if 1
  ImageType::SizeType  meshSize;
  meshSize[ 0 ] = 13;
  meshSize[ 1 ] = 7;
  meshSize[ 2 ] = 4;
  for ( unsigned int  x = 0; x < meshSize[ 0 ]; x++ )
  {
    for ( unsigned int  y = 0; y < meshSize[ 1 ]; y++ )
    {
      for ( unsigned int  z = 0; z < meshSize[ 2 ]; z++ )
      {
        float  x1 = static_cast< float >( x ) / static_cast< float >( meshSize[ 0 ] );
        float  y1 = static_cast< float >( y ) / static_cast< float >( meshSize[ 1 ] );
        float  z1 = static_cast< float >( z ) / static_cast< float >( meshSize[ 2 ] );

        float  x2 = static_cast< float >( x+1 ) / static_cast< float >( meshSize[ 0 ] );
        float  y2 = static_cast< float >( y+1 ) / static_cast< float >( meshSize[ 1 ] );
        float  z2 = static_cast< float >( z+1 ) / static_cast< float >( meshSize[ 2 ] );

        x1 = static_cast< unsigned int >( ( size[ 0 ] - 1 ) * x1 );
        y1 = static_cast< unsigned int >( ( size[ 1 ] - 1 ) * y1 );
        z1 = static_cast< unsigned int >( ( size[ 2 ] - 1 ) * z1 );

        x2 = static_cast< unsigned int >( ( size[ 0 ] - 1 ) * x2 );
        y2 = static_cast< unsigned int >( ( size[ 1 ] - 1 ) * y2 );
        z2 = static_cast< unsigned int >( ( size[ 2 ] - 1 ) * z2 );

        const float  p0[] = { x1, y1, z1 };
        const float  p1[] = { x2, y1, z1 };
        const float  p2[] = { x1, y2, z1 };
        const float  p3[] = { x2, y2, z1 };
        const float  p4[] = { x1, y1, z2 };
        const float  p5[] = { x2, y1, z2 };
        const float  p6[] = { x1, y2, z2 };
        const float  p7[] = { x2, y2, z2 };

        // Each cube will be filled by 5 tetrahedra. There are two different configurations however that should alternate in a checker-board pattern.
        bool flippedConfiguration = true;
        if ( x % 2 )
        {
          flippedConfiguration = !flippedConfiguration;
        }
        if ( y % 2 )
        {
          flippedConfiguration = !flippedConfiguration;
        }
        if ( z % 2 )
        {
          flippedConfiguration = !flippedConfiguration;
        }

        if ( flippedConfiguration )
        {
          meshSource->AddTetrahedron( meshSource->AddPoint( p0 ),
                                      meshSource->AddPoint( p1 ),
                                      meshSource->AddPoint( p4 ),
                                      meshSource->AddPoint( p2 ) );
          meshSource->AddTetrahedron( meshSource->AddPoint( p1 ),
                                      meshSource->AddPoint( p4 ),
                                      meshSource->AddPoint( p7 ),
                                      meshSource->AddPoint( p5 ) );
          meshSource->AddTetrahedron( meshSource->AddPoint( p1 ),
                                      meshSource->AddPoint( p4 ),
                                      meshSource->AddPoint( p2 ),
                                      meshSource->AddPoint( p7 ) );
          meshSource->AddTetrahedron( meshSource->AddPoint( p1 ),
                                      meshSource->AddPoint( p2 ),
                                      meshSource->AddPoint( p3 ),
                                      meshSource->AddPoint( p7 ) );
          meshSource->AddTetrahedron( meshSource->AddPoint( p2 ),
                                      meshSource->AddPoint( p7 ),
                                      meshSource->AddPoint( p4 ),
                                      meshSource->AddPoint( p6 ) );

        }
        else
        {
          meshSource->AddTetrahedron( meshSource->AddPoint( p3 ),
                                      meshSource->AddPoint( p1 ),
                                      meshSource->AddPoint( p0 ),
                                      meshSource->AddPoint( p5 ) );
          meshSource->AddTetrahedron( meshSource->AddPoint( p0 ),
                                      meshSource->AddPoint( p3 ),
                                      meshSource->AddPoint( p6 ),
                                      meshSource->AddPoint( p2 ) );
          meshSource->AddTetrahedron( meshSource->AddPoint( p3 ),
                                      meshSource->AddPoint( p5 ),
                                      meshSource->AddPoint( p6 ),
                                      meshSource->AddPoint( p7 ) );
          meshSource->AddTetrahedron( meshSource->AddPoint( p0 ),
                                      meshSource->AddPoint( p6 ),
                                      meshSource->AddPoint( p5 ),
                                      meshSource->AddPoint( p4 ) );
          meshSource->AddTetrahedron( meshSource->AddPoint( p0 ),
                                      meshSource->AddPoint( p3 ),
                                      meshSource->AddPoint( p5 ),
                                      meshSource->AddPoint( p6 ) );
        }

      }
    }
  }


  // Assign alphas reflecting distance from center
  for ( kvl::AtlasMesh::PointsContainer::Iterator  pointIt = meshSource->GetOutput()->GetPoints()->Begin();
        pointIt != meshSource->GetOutput()->GetPoints()->End();
        ++pointIt )
  {
    kvl::AtlasMesh::PointType  p = pointIt.Value();

    kvl::AtlasAlphasType   alphas( 2 );
    float  squareDistanceToCenter =
      pow( p[ 0 ] - static_cast<float>( size[0] ) / 2.0f, 2 ) +
      pow( p[ 1 ] - static_cast<float>( size[1] ) / 2.0f, 2 ) +
      pow( p[ 2 ] - static_cast<float>( size[2] ) / 2.0f, 2 );
    float  sigma = static_cast< float >( size[ 0 ] ) / 4.0f;
    alphas[ 0 ] = exp( -squareDistanceToCenter / 2 / sigma / sigma );
    alphas[ 1 ] = 1 - alphas[ 0 ];

    //std::cout << "Setting alpha of vertex with position " << p << " to " << alphas[ 0 ] << std::endl;
    kvl::AtlasMesh::PixelType  pointParameters;
    pointParameters.m_Alphas = alphas;
    meshSource->GetOutput()->SetPointData( pointIt.Index(), pointParameters );
  }
#else

  const float  p0[] = { 50, 50, 70 };
  const float  p1[] = { 30, 50, 30 };
  const float  p2[] = { 70, 30, 34 };
  const float  p3[] = { 70, 70, 38 };
  meshSource->AddTetrahedron( meshSource->AddPoint( p0 ),
                              meshSource->AddPoint( p1 ),
                              meshSource->AddPoint( p2 ),
                              meshSource->AddPoint( p3 ) );

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

  // Rasterize mesh
  kvl::AtlasMeshAlphaDrawer::Pointer  alphaDrawer = kvl::AtlasMeshAlphaDrawer::New();
  alphaDrawer->SetLabelImage( image );
  alphaDrawer->SetLabelNumber( 0 );
  alphaDrawer->Rasterize( meshSource->GetOutput() );

  // Write out
  typedef itk::ImageFileWriter< kvl::AtlasMeshAlphaDrawer::AlphaImageType >  AlphaWriterType;
  AlphaWriterType::Pointer  alphaWriter = AlphaWriterType::New();
  alphaWriter->SetFileName( "testAlpha.mhd" );
  alphaWriter->SetInput( alphaDrawer->GetAlphaImage() );
  alphaWriter->Write();


  // Now rasterize mesh again with VisitCounter
  kvl::AtlasMeshVisitCounter::Pointer  visitCounter = kvl::AtlasMeshVisitCounter::New();
  visitCounter->SetLabelImage( image );
  visitCounter->Rasterize( meshSource->GetOutput() );

#endif


  // Create an empty image
  typedef itk::Image< unsigned char, 3 >  ImageType;
  ImageType::SizeType  size;
  size[ 0 ] = 3;
  size[ 1 ] = 3;
  size[ 2 ] = 3;
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
  collection->Construct( meshSize, domainSize, 1000,
                         numberOfClasses, numberOfMeshes );
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
      if ( ( it.GetIndex()[ 0 ] !=
             static_cast< long >( visitCounter->GetCountImage()->GetLargestPossibleRegion().GetSize()[ 0 ]-1 ) ) &&
           ( it.GetIndex()[ 1 ] !=
             static_cast< long >( visitCounter->GetCountImage()->GetLargestPossibleRegion().GetSize()[ 1 ]-1 ) ) &&
           ( it.GetIndex()[ 2 ] !=
             static_cast< long >( visitCounter->GetCountImage()->GetLargestPossibleRegion().GetSize()[ 2 ]-1 ) ) )
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
  for ( unsigned int z = 0; z < size[ 2 ]; z++ )
  {
    index[ 2 ] = z;
    for ( unsigned int y = 0; y < size[ 1 ]; y++ )
    {
      index[ 1 ] = y;
      for ( unsigned int x = 0; x < size[ 0 ]; x++ )
      {
        index[ 0 ] = x;
        std::cout << static_cast< int >( visitCounter->GetCountImage()->GetPixel( index ) ) << "  ";
      }
      std::cout << "\n" << std::endl;
    }

    std::cout << "\n\n\n" << std::endl;
  }


  std::cout << "\n\n\nAlpha image: " << std::endl;
  for ( unsigned int z = 0; z < size[ 2 ]; z++ )
  {
    index[ 2 ] = z;
    for ( unsigned int y = 0; y < size[ 1 ]; y++ )
    {
      index[ 1 ] = y;
      for ( unsigned int x = 0; x < size[ 0 ]; x++ )
      {
        index[ 0 ] = x;
        std::cout << alphaDrawer->GetAlphaImage()->GetPixel( index ) << "  ";
      }
      std::cout << "\n" << std::endl;
    }

    std::cout << "\n\n\n" << std::endl;
  }

  return 0;
};

