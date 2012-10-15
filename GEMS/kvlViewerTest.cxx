/**
 * @file  kvlViewerTest.cxx
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


#include "itkVTKImageExport.h"
#include "vtkImageImport.h"
#include "vtkSmartPointer.h"
#include "vtkImageData.h"
#include "vtkCellArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkMetaImageWriter.h"
#include "vtkDataSetWriter.h"
#include "itkGaussianImageSource.h"
#include "itkImageFileWriter.h"

#include "kvlImageViewer.h"



/**
 * This function will connect the given itk::VTKImageExport filter to
 * the given vtkImageImport filter.
 */
template <typename ITK_Exporter, typename VTK_Importer>
void ConnectPipelines(ITK_Exporter exporter, VTK_Importer importer)
{
  importer->SetUpdateInformationCallback(exporter->GetUpdateInformationCallback());
  importer->SetPipelineModifiedCallback(exporter->GetPipelineModifiedCallback());
  importer->SetWholeExtentCallback(exporter->GetWholeExtentCallback());
  importer->SetSpacingCallback(exporter->GetSpacingCallback());
  importer->SetOriginCallback(exporter->GetOriginCallback());
  importer->SetScalarTypeCallback(exporter->GetScalarTypeCallback());
  importer->SetNumberOfComponentsCallback(exporter->GetNumberOfComponentsCallback());
  importer->SetPropagateUpdateExtentCallback(exporter->GetPropagateUpdateExtentCallback());
  importer->SetUpdateDataCallback(exporter->GetUpdateDataCallback());
  importer->SetDataExtentCallback(exporter->GetDataExtentCallback());
  importer->SetBufferPointerCallback(exporter->GetBufferPointerCallback());
  importer->SetCallbackUserData(exporter->GetCallbackUserData());
}



vtkSmartPointer< vtkUnstructuredGrid >  ConvertToVTKUnstructedGrid( const kvl::AtlasMesh* mesh )
{
  // Construct a VTK points container
  vtkSmartPointer< vtkPoints >  vPoints = vtkSmartPointer< vtkPoints >::New();
  vPoints->SetNumberOfPoints( mesh->GetNumberOfPoints() );
  for( kvl::AtlasMesh::PointsContainer::ConstIterator it = mesh->GetPoints()->Begin();
       it != mesh->GetPoints()->End(); ++it )
  {
    vPoints->SetPoint( it->Index(), it->Value().GetDataPointer() );
  }

  std::cout << "Converted points" << std::endl;

  // Construct a VTK type array and cell array
  int*  vTypes = new int[ mesh->GetNumberOfCells() ];
  vtkSmartPointer< vtkCellArray >  vCells = vtkSmartPointer< vtkCellArray >::New();
  vCells->EstimateSize( mesh->GetNumberOfCells(), 3 );
  int vCellCount = 0;
  for ( kvl::AtlasMesh::CellsContainer::ConstIterator  it = mesh->GetCells()->Begin();
        it != mesh->GetCells()->End(); ++it )
  {
    kvl::AtlasMesh::CellType*  cell = it.Value();
    if( cell->GetType() != kvl::AtlasMesh::CellType::TRIANGLE_CELL )
    {
      continue;
    }

    // Insert a triangle cell
    //vCells->InsertNextCell( 3,  static_cast< vtkIdType* >( cell->PointIdsBegin() ) );
    vCells->InsertNextCell( 3,  (vtkIdType*)cell->PointIdsBegin() );
    vTypes[ vCellCount ] = VTK_TRIANGLE;
    vCellCount++;

  } // End loop over all cells in the mesh

  std::cout << "Converted cells and types" << std::endl;


  // Put everything together into a VTK unstructed grid
  vtkSmartPointer< vtkUnstructuredGrid >  vGrid = vtkSmartPointer< vtkUnstructuredGrid >::New();
  vGrid->SetPoints( vPoints );
  vGrid->SetCells( vTypes, vCells );

  std::cout << "Putted everything" << std::endl;

  return vGrid;
}





int main( int argc, char* argv[] )
{


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
  std::cout << "Generating the mesh" << std::endl;
  typedef itk::AutomaticTopologyMeshSource< kvl::AtlasMesh >  MeshSourceType;
  MeshSourceType::Pointer  meshSource = MeshSourceType::New();
#if 1
  ImageType::SizeType  meshSize;
  meshSize[ 0 ] = 30;
  meshSize[ 1 ] = 30;
  meshSize[ 2 ] = 30;
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
  std::cout << "Rasterizing the mesh" << std::endl;
  kvl::AtlasMeshAlphaDrawer::Pointer  alphaDrawer = kvl::AtlasMeshAlphaDrawer::New();
  alphaDrawer->SetLabelImage( image );
  alphaDrawer->SetLabelNumber( 0 );
  alphaDrawer->Rasterize( meshSource->GetOutput() );

  // Make the image available as a VTK object
  std::cout << "Makding and writing VTK image" << std::endl;
  typedef itk::VTKImageExport< kvl::AtlasMeshAlphaDrawer::AlphaImageType > ExportFilterType;
  ExportFilterType::Pointer  exporter = ExportFilterType::New();
  exporter->SetInput( alphaDrawer->GetAlphaImage() );
  vtkSmartPointer< vtkImageImport >  vtkImporter = vtkSmartPointer< vtkImageImport >::New();
  ConnectPipelines( exporter, vtkImporter );

  // Write the VTK image
  vtkSmartPointer< vtkMetaImageWriter >  vtkImageWriter = vtkSmartPointer< vtkMetaImageWriter >::New();
  //vtkImporter->GetOutput()->Print( std::cout );
  vtkImageWriter->SetInput( vtkImporter->GetOutput() );
  vtkImageWriter->SetFileName( "alphaImage.mhd" );
  vtkImageWriter->Write();



  // Make the mesh's vertices, edges, and triangles available as a VTK object
  std::cout << "Makding and writing VTK grid" << std::endl;
  vtkSmartPointer< vtkUnstructuredGrid >  vGrid = ConvertToVTKUnstructedGrid( meshSource->GetOutput() );
  std::cout << "Still here" << std::endl;
  vGrid->Print( std::cout );
  std::cout << "Being out" << std::endl;

  // Write the VTK unstructured grid
  vtkSmartPointer< vtkDataSetWriter >  gridWriter = vtkSmartPointer< vtkDataSetWriter >::New();
  gridWriter->SetInput( vGrid );
  gridWriter->SetFileName( "grid.vtk" );
  gridWriter->Write();
  std::cout << "now here" << std::endl;

  // Create some image to serve as an overlay and write it out
  typedef itk::GaussianImageSource< ImageType >  GaussianSourceType;
  GaussianSourceType::ArrayType  mean;
  mean[ 0 ] = size[ 0 ] / 4;
  mean[ 1 ] = size[ 1 ] / 4;
  mean[ 2 ] = size[ 2 ] / 4;
  GaussianSourceType::ArrayType  sigma;
  sigma[ 0 ] = size[ 0 ] / 4;
  sigma[ 1 ] = size[ 1 ] / 4;
  sigma[ 2 ] = size[ 2 ] / 4;
  GaussianSourceType::Pointer  source = GaussianSourceType::New();
  source->SetSize( size );
  source->SetMean( mean );
  source->SetSigma( sigma );
  source->SetScale( 255 );

  typedef itk::ImageFileWriter< ImageType >  WriterType;
  WriterType::Pointer  writer = WriterType::New();
  writer->SetInput( source->GetOutput() );
  writer->SetFileName( "overlay.mhd" );
  writer->Write();







  // Now try the viewer
  typedef itk::IntensityWindowingImageFilter< kvl::AtlasMeshAlphaDrawer::AlphaImageType,
          kvl::ImageViewer::ImageType >   WindowerType;
  WindowerType::Pointer  windower = WindowerType::New();
  windower->SetInput( alphaDrawer->GetAlphaImage() );
  windower->SetWindowMinimum( 0 );
  windower->SetWindowMaximum( 1.0f );
  windower->SetOutputMinimum( 0 );
  windower->SetOutputMaximum( 255 );
  windower->Update();

  kvl::ImageViewer  viewer( 100, 100, 1500, 750 );
  viewer.SetImage( windower->GetOutput() );
  viewer.SetOverlayImage( source->GetOutput() );
  viewer.SetMesh( meshSource->GetOutput() );
  viewer.show();

  Fl::run();



  return 0;
};

