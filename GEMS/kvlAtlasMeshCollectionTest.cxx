/**
 * @file  kvlAtlasMeshCollectionTest.cxx
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:38 $
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
#include "kvlAtlasMeshCollection.h"
#include "kvlAtlasMeshCollectionValidator.h"

#if 0
#include "kvlImageViewer.h"
#endif


int main( int, char** )
{

  // Create
  kvl::AtlasMeshCollection::Pointer  collection =  kvl::AtlasMeshCollection::New();
  const unsigned int  meshSize[] = { 10, 10, 10 };
  const unsigned int  domainSize[] = {100, 100, 100 };
  const float  initialStiffness = 5.0f;
  const unsigned int  numberOfClasses = 2;
  const unsigned int  numberOfMeshes = 4;
  collection->Construct( meshSize, domainSize, initialStiffness,
                         numberOfClasses, numberOfMeshes );

  // Write the mesh collection
  collection->Write( "collection.txt" );

  // Reconstruct the collection from the file
  kvl::AtlasMeshCollection::Pointer  collection2 =  kvl::AtlasMeshCollection::New();
  collection2->Read( "collection.txt" );

  // Write the reconstructed mesh collection again
  collection2->Write( "collection2.txt" );

  // Region-grow
  kvl::AtlasMeshCollection::Pointer  regionGrown = collection->GetRegionGrown( 730 /*34*/, 1 );
  regionGrown->Write( "regionGrown1.txt" );
  regionGrown = collection->GetRegionGrown( 730 /*34*/, 2 );
  regionGrown->Write( "regionGrown2.txt" );


  // Check validaty of the region-grown mesh collection
  kvl::AtlasMeshCollectionValidator::Pointer  validator = kvl::AtlasMeshCollectionValidator::New();
  std::cout << "Valdating regionGrown" << std::endl;
  validator->Validate( regionGrown );



  typedef kvl::AtlasMeshCollection::TransformType  TransformType;
  TransformType::Pointer  transform = TransformType::New();
  TransformType::OutputVectorType  translation;
  translation[ 0 ] = 10;
  translation[ 1 ] = 10;
  translation[ 2 ] = 10;
  transform->Translate( translation );
  transform->Scale( 0.7 );

  kvl::AtlasMeshCollection::Pointer  regionGrown2 =  kvl::AtlasMeshCollection::New();
  regionGrown2->Read( "regionGrown2.txt" );
  for ( int i = -1; i < static_cast< int >( numberOfMeshes ); i++ )
  {
    regionGrown2->Transform( i, transform );
  }
  regionGrown2->Write( "regionGrown2_readAndWritten.txt" );


  // Collapse edge
  kvl::AtlasMeshCollection::Pointer  collapsed;
  std::set< kvl::AtlasMesh::CellIdentifier >  disappearingCells;
  kvl::AtlasMesh::CellIdentifier  unifiedVertexId;
  //collection->GetCollapsed( 730, collapsed, disappearingCells, unifiedVertexId );
  regionGrown->GetCollapsed( 730, collapsed, disappearingCells, unifiedVertexId );
  collapsed->Write( "collapsed.txt" );

  // Validate collapsed edge mesh collection
  std::cout << "Validating collapsed" << std::endl;
  validator->Validate( collapsed );



  /*  // Display the last mesh of the original collection
    typedef kvl::ImageViewer::ImageType  ImageType;
    ImageType::SizeType  size;
    size[ 0 ] = domainSize[ 0 ];
    size[ 1 ] = domainSize[ 1 ];
    ImageType::Pointer  image = ImageType::New();
    image->SetRegions( size );
    image->Allocate();
    image->FillBuffer( 0 );

    const unsigned long  triangleIdToHighlight = 86;
    //const unsigned long edgeIdToHighlight = 89;
    //const unsigned long edgeIdToHighlight = 82;
    //const unsigned long edgeIdToHighlight = 50;
    const unsigned long edgeIdToHighlight = 9;

    kvl::ImageViewer  viewer( 100, 100, 256, 256 );
    viewer.SetImage( image );
    viewer.SetMesh( collection->GetMesh( numberOfMeshes-1 ) );
    viewer.SetTriangleIdToHighlight( triangleIdToHighlight );
    viewer.SetEdgeIdToHighlight( edgeIdToHighlight );
    viewer.show();

    // Display the last mesh of the reconstructed collection
    kvl::ImageViewer  viewer2( 400, 100, 256, 256 );
    viewer2.SetImage( image );
    viewer2.SetMesh( collection2->GetMesh( numberOfMeshes-1 ) );
    viewer2.SetTriangleIdToHighlight( triangleIdToHighlight );
    viewer2.SetEdgeIdToHighlight( edgeIdToHighlight );
    viewer2.show();*/
  /*
    // Collapse edge
    kvl::AtlasMeshCollection::Pointer  collapsed;
    std::set< kvl::AtlasMesh::CellIdentifier >  disappearingCells;
    kvl::AtlasMesh::CellIdentifier  unifiedVertexId;
    collection->GetCollapsed( edgeIdToHighlight, collapsed, disappearingCells, unifiedVertexId );
    collapsed->Write( "collapsed.txt" );

    // Display the last collapsed mesh
    kvl::ImageViewer  viewer3( 700, 100, 256, 256 );
    viewer3.SetImage( image );
    viewer3.SetMesh( collapsed->GetMesh( numberOfMeshes-1 ) );
    viewer3.SetTriangleIdToHighlight( triangleIdToHighlight );
    viewer3.SetEdgeIdToHighlight( edgeIdToHighlight );
    viewer3.show();

    // Region grow a vertex
    kvl::AtlasMeshCollection::Pointer  regionGrown = collapsed->GetRegionGrown( 43, 2 );
    regionGrown->Write( "regionGrown.txt" );

    // Display the last region grown mesh
    kvl::ImageViewer  viewer4( 100, 400, 256, 256 );
    viewer4.SetImage( image );
    viewer4.SetMesh( regionGrown->GetMesh( numberOfMeshes-1 ) );
    viewer4.SetTriangleIdToHighlight( triangleIdToHighlight );
    viewer4.SetEdgeIdToHighlight( edgeIdToHighlight );
    viewer4.show();

    // Upsample the original mesh
    kvl::AtlasMeshCollection::Pointer  upsampled = collection->GetUpsampled( domainSize, initialStiffness );
    upsampled->Write( "upsampled.txt" );

    // Display the last upsampled mesh
    kvl::ImageViewer  viewer5( 400, 400, 256, 256 );
    viewer5.SetImage( image );
    viewer5.SetMesh( upsampled->GetMesh( numberOfMeshes-1 ) );
    viewer5.show();

    // Edge-split the original mesh
    kvl::AtlasMesh::CellIdentifier  newVertexId = 10000;
    kvl::AtlasMesh::PointIdentifier  newPointId = 10000;
    kvl::AtlasMeshCollection::Pointer  splitted = collection->GetEdgeSplitted( edgeIdToHighlight, newVertexId, newPointId );
    splitted->Write( "splitted.txt" );

    // Display the edge-split mesh
    kvl::ImageViewer  viewer6( 700, 400, 256, 256 );
    viewer6.SetImage( image );
    viewer6.SetMesh( splitted->GetMesh( numberOfMeshes-1 ) );
    viewer6.show();

    // Edge-swap the original mesh
    kvl::AtlasMeshCollection::Pointer  swapped = collection->GetEdgeSwapped( edgeIdToHighlight );
    swapped->Write( "swapped.txt" );

    // Display the edge-swapped mesh
    kvl::ImageViewer  viewer7( 100, 700, 256, 256 );
    viewer7.SetImage( image );
    viewer7.SetMesh( swapped->GetMesh( numberOfMeshes-1 ) );
    viewer7.Show();*/


  return 0;
};


