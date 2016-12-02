#include "kvlImageViewer.h"

#include "itkMinimumMaximumImageCalculator.h"
#include "vtkImageImport.h"
#include "vtkImageMapToColors.h"
#include "vtkImageData.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderer.h"
#include "vtkProperty.h"
#include "vtkClipPolyData.h"
#include "vtkTubeFilter.h"
#include "vtkRenderWindow.h"
#include "vtkInteractorStyleSwitch.h"
#include "vtkCellArray.h"
#include "vtkWindowToImageFilter.h"
#include "vtkPNGWriter.h"


namespace kvl
{



//
//
//
template <typename ITK_Exporter, typename VTK_Importer>
void ConnectPipelines( ITK_Exporter exporter, VTK_Importer importer )
{
  importer->SetUpdateInformationCallback(exporter->GetUpdateInformationCallback());
  importer->SetPipelineModifiedCallback(exporter->GetPipelineModifiedCallback());
  importer->SetWholeExtentCallback(exporter->GetWholeExtentCallback());
  //importer->SetSpacingCallback(exporter->GetSpacingCallback());
  //importer->SetOriginCallback(exporter->GetOriginCallback());
  importer->SetScalarTypeCallback(exporter->GetScalarTypeCallback());
  importer->SetNumberOfComponentsCallback(exporter->GetNumberOfComponentsCallback());
  importer->SetPropagateUpdateExtentCallback(exporter->GetPropagateUpdateExtentCallback());
  importer->SetUpdateDataCallback(exporter->GetUpdateDataCallback());
  importer->SetDataExtentCallback(exporter->GetDataExtentCallback());
  importer->SetBufferPointerCallback(exporter->GetBufferPointerCallback());
  importer->SetCallbackUserData(exporter->GetCallbackUserData());
}



//
//
//
ImageViewer
::ImageViewer( int x, int y, int w, int h, const char *l ) : vtkFlRenderWindowInteractor(x,y,w,h,l)
{ 
  m_Image = 0;
  m_OverlayImage = 0;
  m_Mesh = 0;

  m_OverlayColor[ 0 ] = 0.3;
  m_OverlayColor[ 1 ] = 0.3;
  m_OverlayColor[ 2 ] = 1.0;
  m_OverlayAlpha = 1.0;

  for ( int i = 0; i < 3; i++ )
    {
    m_MaximumImageIndex[ i ] = 0;
    }

  m_SagittalSliceNumber = 0;
  m_AxialSliceNumber = 0;
  m_CoronalSliceNumber = 0;

  m_ImageExporter = 0;
  m_RGBAImageExporter = 0;
  m_OverlayImageExporter = 0;
  m_RGBAOverlayImageExporter = 0;

  // Create points of connection to VTK pipelines
  m_OutlineFilter = vtkSmartPointer< vtkOutlineFilter >::New();
  m_EdgeExtracter = vtkSmartPointer< vtkExtractEdges >::New();


  // Create lookup tables
  m_ImageLookupTable = vtkSmartPointer< vtkLookupTable >::New();
  m_ImageLookupTable->SetTableRange( 0, 255 );
  m_ImageLookupTable->SetSaturationRange( 0, 0 );
  m_ImageLookupTable->SetHueRange( 0, 0 );
  m_ImageLookupTable->SetValueRange( 0, 1 );
  m_ImageLookupTable->Build();
  //m_ImageLookupTable->SetTableValue( 0, 0.0, 0.0, 0.0, 0.0 ); // Make black voxels transparent

  m_OverlayImageLookupTable = vtkSmartPointer< vtkLookupTable >::New();
  m_OverlayImageLookupTable->SetTableRange( 0, 255 );
  m_OverlayImageLookupTable->SetHueRange( 0.75, 0.75 );  // The actual color
  m_OverlayImageLookupTable->SetSaturationRange( 1, 1 );
  m_OverlayImageLookupTable->SetValueRange( 1, 1 );
  m_OverlayImageLookupTable->SetAlphaRange( 0, 1 );
  m_OverlayImageLookupTable->Build();


  // Set up implicit functions used for cutting and clipping
  m_SagittalPlane = vtkSmartPointer< vtkPlane >::New();
  m_SagittalPlane->SetNormal( -1, 0, 0 );

  m_CoronalPlane = vtkSmartPointer< vtkPlane >::New();
  m_CoronalPlane->SetNormal( 0, -1, 0 );

  m_AxialPlane = vtkSmartPointer< vtkPlane >::New();
  m_AxialPlane->SetNormal( 0, 0, -1 );



  // Set up VTK pipeline for the m_Image
  m_Importer = vtkSmartPointer< vtkImageImport >::New();

  m_SagittalColors = vtkSmartPointer< vtkImageMapToColors >::New();
  m_SagittalColors->SetInput( m_Importer->GetOutput() );
  m_SagittalColors->SetLookupTable( m_ImageLookupTable );
  m_SagittalColors->PassAlphaToOutputOn();

  m_CoronalColors = vtkSmartPointer< vtkImageMapToColors >::New();
  m_CoronalColors->SetInput( m_Importer->GetOutput() );
  m_CoronalColors->SetLookupTable( m_ImageLookupTable );
  m_CoronalColors->PassAlphaToOutputOn();

  m_AxialColors = vtkSmartPointer< vtkImageMapToColors >::New();
  m_AxialColors->SetInput( m_Importer->GetOutput() );
  m_AxialColors->SetLookupTable( m_ImageLookupTable );
  m_AxialColors->PassAlphaToOutputOn();



  // Set up VTK pipeline for the m_ImageOverlay
  m_OverlayImporter = vtkSmartPointer< vtkImageImport >::New();

  m_SagittalOverlayColors = vtkSmartPointer< vtkImageMapToColors >::New();
  m_SagittalOverlayColors->SetInput( m_OverlayImporter->GetOutput() );
  //m_SagittalOverlayColors->SetLookupTable( m_OverlayImageLookupTable );
  m_SagittalOverlayColors->PassAlphaToOutputOn();

  m_CoronalOverlayColors = vtkSmartPointer< vtkImageMapToColors >::New();
  m_CoronalOverlayColors->SetInput( m_OverlayImporter->GetOutput() );
  //m_CoronalOverlayColors->SetLookupTable( m_OverlayImageLookupTable );
  m_CoronalOverlayColors->PassAlphaToOutputOn();

  m_AxialOverlayColors = vtkSmartPointer< vtkImageMapToColors >::New();
  m_AxialOverlayColors->SetInput( m_OverlayImporter->GetOutput() );
  //m_AxialOverlayColors->SetLookupTable( m_OverlayImageLookupTable );
  m_AxialOverlayColors->PassAlphaToOutputOn();

  this->SetOverlayImageLookupTable( m_OverlayImageLookupTable );

  // Set up VTK pipeline for the blended m_Image and m_ImageOverlay
  m_SagittalBlender = vtkSmartPointer< vtkImageBlend >::New();
  m_SagittalActor = vtkSmartPointer< vtkImageActor >::New();
  m_SagittalActor->SetInput( m_SagittalBlender->GetOutput() );
  m_SagittalActor->InterpolateOff();

  m_CoronalBlender = vtkSmartPointer< vtkImageBlend >::New();
  m_CoronalActor = vtkSmartPointer< vtkImageActor >::New();
  m_CoronalActor->SetInput( m_CoronalBlender->GetOutput() );
  m_CoronalActor->InterpolateOff();

  m_AxialBlender = vtkSmartPointer< vtkImageBlend >::New();
  m_AxialActor = vtkSmartPointer< vtkImageActor >::New();
  m_AxialActor->SetInput( m_AxialBlender->GetOutput() );
  m_AxialActor->InterpolateOff();


  // Set up VTK pipeline for the m_Mesh: outline
  vtkSmartPointer< vtkPolyDataMapper >  outlineMapper = vtkSmartPointer< vtkPolyDataMapper >::New();
  outlineMapper->SetInput( m_OutlineFilter->GetOutput() );

  m_OutlineActor = vtkSmartPointer< vtkActor >::New();
  m_OutlineActor->SetMapper( outlineMapper );
  m_OutlineActor->GetProperty()->SetColor( 0, 0, 0 );


  // Set up VTK pipeline for the m_Mesh: extract edges to display
  vtkSmartPointer< vtkClipPolyData >  axialEdgeClipper = vtkSmartPointer< vtkClipPolyData >::New();
  axialEdgeClipper->SetInput( m_EdgeExtracter->GetOutput() );
  axialEdgeClipper->SetClipFunction( m_AxialPlane );
  axialEdgeClipper->SetValue( 0.0f );

  vtkSmartPointer< vtkTubeFilter >  axialEdgeTuber = vtkSmartPointer< vtkTubeFilter >::New();
  axialEdgeTuber->SetInput( axialEdgeClipper->GetOutput() );
  axialEdgeTuber->SetRadius( 0.1f );
  axialEdgeTuber->SetNumberOfSides( 6 );

  vtkSmartPointer< vtkPolyDataMapper >  axialEdgeMapper = vtkSmartPointer< vtkPolyDataMapper >::New();
  axialEdgeMapper->SetInput( axialEdgeTuber->GetOutput() );
  m_AxialEdgeActor = vtkSmartPointer< vtkActor >::New();
  m_AxialEdgeActor->SetMapper( axialEdgeMapper );
  m_AxialEdgeActor->GetProperty()->SetColor( 1, 0, 0 );
  m_AxialEdgeActor->GetProperty()->SetSpecularColor( 1, 1, 1 );
  m_AxialEdgeActor->GetProperty()->SetSpecular( 0.3 );
  m_AxialEdgeActor->GetProperty()->SetSpecularPower( 20 );
  m_AxialEdgeActor->GetProperty()->SetAmbient( 0.2 );
  m_AxialEdgeActor->GetProperty()->SetDiffuse( 0.8 );


  // Set up VTK pipeline for the m_Mesh: extract the 2-D cut
  m_SagittalCutter = vtkSmartPointer< vtkCutter >::New();
  m_SagittalCutter->SetCutFunction( m_SagittalPlane );

  vtkSmartPointer< vtkPolyDataMapper >  sagittalCutMapper = vtkSmartPointer< vtkPolyDataMapper >::New();
  sagittalCutMapper->SetInput( m_SagittalCutter->GetOutput() );

  m_SagittalCutActor = vtkSmartPointer< vtkActor >::New();
  m_SagittalCutActor->SetMapper( sagittalCutMapper );
  m_SagittalCutActor->GetProperty()->SetColor( 0, 1, 0 );
  m_SagittalCutActor->GetProperty()->SetSpecularColor( 1, 1, 1 );
  m_SagittalCutActor->GetProperty()->SetSpecular( 0.3 );
  m_SagittalCutActor->GetProperty()->SetSpecularPower( 20 );
  m_SagittalCutActor->GetProperty()->SetAmbient( 0.2 );
  m_SagittalCutActor->GetProperty()->SetDiffuse( 0.8 );


  m_CoronalCutter = vtkSmartPointer< vtkCutter >::New();
  m_CoronalCutter->SetCutFunction( m_CoronalPlane );

  vtkSmartPointer< vtkPolyDataMapper >  coronalCutMapper = vtkSmartPointer< vtkPolyDataMapper >::New();
  coronalCutMapper->SetInput( m_CoronalCutter->GetOutput() );

  m_CoronalCutActor = vtkSmartPointer< vtkActor >::New();
  m_CoronalCutActor->SetMapper( coronalCutMapper );
  m_CoronalCutActor->GetProperty()->SetColor( 0, 1, 0 );
  m_CoronalCutActor->GetProperty()->SetSpecularColor( 1, 1, 1 );
  m_CoronalCutActor->GetProperty()->SetSpecular( 0.3 );
  m_CoronalCutActor->GetProperty()->SetSpecularPower( 20 );
  m_CoronalCutActor->GetProperty()->SetAmbient( 0.2 );
  m_CoronalCutActor->GetProperty()->SetDiffuse( 0.8 );


  m_AxialCutter = vtkSmartPointer< vtkCutter >::New();
  m_AxialCutter->SetCutFunction( m_AxialPlane );

  vtkSmartPointer< vtkPolyDataMapper >  axialCutMapper = vtkSmartPointer< vtkPolyDataMapper >::New();
  axialCutMapper->SetInput( m_AxialCutter->GetOutput() );

  m_AxialCutActor = vtkSmartPointer< vtkActor >::New();
  m_AxialCutActor->SetMapper( axialCutMapper );
  m_AxialCutActor->GetProperty()->SetColor( 0, 1, 0 );
  m_AxialCutActor->GetProperty()->SetSpecularColor( 1, 1, 1 );
  m_AxialCutActor->GetProperty()->SetSpecular( 0.3 );
  m_AxialCutActor->GetProperty()->SetSpecularPower( 20 );
  m_AxialCutActor->GetProperty()->SetAmbient( 0.2 );
  m_AxialCutActor->GetProperty()->SetDiffuse( 0.8 );


  // Set up renderers
  static double  backgroundColor[] = /* { 0.8, 0.8, 1 } */ { 0.2, 0.2, 0.25 };

  m_ThreeDRenderer = vtkSmartPointer< vtkRenderer >::New();
  m_ThreeDRenderer->SetActiveCamera( vtkSmartPointer< vtkCamera >::New() );
  m_ThreeDRenderer->SetBackground( backgroundColor );

  m_SagittalCamera = vtkSmartPointer< vtkCamera >::New();
  m_SagittalCamera->SetViewUp( 0, 0, 1 );

  m_CoronalCamera = vtkSmartPointer< vtkCamera >::New();
  m_CoronalCamera->SetViewUp( 0, 0, 1 );

  m_AxialCamera = vtkSmartPointer< vtkCamera >::New();
  m_AxialCamera->SetViewUp( 0, -1, 0 );

  m_SagittalRenderer = vtkSmartPointer< vtkRenderer >::New();
  m_SagittalRenderer->SetActiveCamera( m_SagittalCamera );
  m_SagittalRenderer->SetLayer( 0 );
  m_SagittalRenderer->SetBackground( backgroundColor );

  m_SagittalRenderer2 = vtkSmartPointer< vtkRenderer >::New();
  m_SagittalRenderer2->SetActiveCamera( m_SagittalCamera );
  m_SagittalRenderer2->SetLayer( 1 );
  m_SagittalRenderer2->SetBackground( backgroundColor );

  m_CoronalRenderer = vtkSmartPointer< vtkRenderer >::New();
  m_CoronalRenderer->SetActiveCamera( m_CoronalCamera );
  m_CoronalRenderer->SetLayer( 0 );
  m_CoronalRenderer->SetBackground( backgroundColor );

  m_CoronalRenderer2 = vtkSmartPointer< vtkRenderer >::New();
  m_CoronalRenderer2->SetActiveCamera( m_CoronalCamera );
  m_CoronalRenderer2->SetLayer( 1 );
  m_CoronalRenderer2->SetBackground( backgroundColor );

  m_AxialRenderer = vtkSmartPointer< vtkRenderer >::New();
  m_AxialRenderer->SetActiveCamera( m_AxialCamera );
  m_AxialRenderer->SetLayer( 0 );
  m_AxialRenderer->SetBackground( backgroundColor );

  m_AxialRenderer2 = vtkSmartPointer< vtkRenderer >::New();
  m_AxialRenderer2->SetActiveCamera( m_AxialCamera );
  m_AxialRenderer2->SetLayer( 1 );
  m_AxialRenderer2->SetBackground( backgroundColor );


  // Set up render window
  vtkSmartPointer< vtkRenderWindow >  renderWindow = vtkSmartPointer< vtkRenderWindow >::New();
  renderWindow->SetSize( w, h );
  renderWindow->AddRenderer( m_ThreeDRenderer );
  renderWindow->SetNumberOfLayers( 2 );
  renderWindow->AddRenderer( m_SagittalRenderer );
  renderWindow->AddRenderer( m_SagittalRenderer2 );
  renderWindow->AddRenderer( m_CoronalRenderer );
  renderWindow->AddRenderer( m_CoronalRenderer2 );
  renderWindow->AddRenderer( m_AxialRenderer );
  renderWindow->AddRenderer( m_AxialRenderer2 );

  this->LookAt( 5 );


  // Set up interactor (which is ourselves)
  vtkSmartPointer< vtkInteractorStyleSwitch >  interactorStyleSwitch = vtkSmartPointer< vtkInteractorStyleSwitch >::New();
  this->SetInteractorStyle( interactorStyleSwitch );
  interactorStyleSwitch->SetCurrentStyleToTrackballCamera();
  this->SetRenderWindow( renderWindow );
  this->Initialize();



  // Set to the correct initial location, set the correct overlay value etc
  //this->SetSliceLocation( m_SagittalSliceNumber, m_CoronalSliceNumber, m_AxialSliceNumber );
  this->SetOverlayAlpha( m_OverlayAlpha );

}



//
//
//
ImageViewer
::~ImageViewer()
{

}


//
//
//
void
ImageViewer
::SetImage( const ImageBaseType* image )
{
//   // Remember if this is the first time an image is set. If so,
//   // we'll automatically set the 3D camera for the user
//   bool  reset3DCamera = false;
//   if ( !m_Image )
//     {
//     reset3DCamera = true;
//     }

  m_Image = image;

  if ( !m_Image )
    {
    // Remove the actors and return
    m_ThreeDRenderer->RemoveActor( m_SagittalActor );
    m_ThreeDRenderer->RemoveActor( m_AxialActor );
    m_ThreeDRenderer->RemoveActor( m_CoronalActor );

    m_SagittalRenderer->RemoveActor( m_SagittalActor );
    m_CoronalRenderer->RemoveActor( m_CoronalActor );
    m_AxialRenderer->RemoveActor( m_AxialActor );

    return;
    }

  // Connect the VTK pipeline bits and pieces
  if ( dynamic_cast< const ImageType* >( m_Image.GetPointer() ) )
    {
    ImageType::ConstPointer castImage = static_cast< const ImageType* >( m_Image.GetPointer() );
    m_ImageExporter = ExporterType::New();
    m_ImageExporter->SetInput( castImage );
    ConnectPipelines( m_ImageExporter, m_Importer );

    m_SagittalColors->SetLookupTable( m_ImageLookupTable );
    m_CoronalColors->SetLookupTable( m_ImageLookupTable );
    m_AxialColors->SetLookupTable( m_ImageLookupTable );
    }
  else if ( dynamic_cast< const RGBAImageType* >( m_Image.GetPointer() ) )
    {
    RGBAImageType::ConstPointer castImage = static_cast< const RGBAImageType* >( m_Image.GetPointer() );
    m_RGBAImageExporter = RGBAExporterType::New();
    m_RGBAImageExporter->SetInput( castImage );
    ConnectPipelines( m_RGBAImageExporter, m_Importer );

    m_SagittalColors->SetLookupTable( 0 );
    m_CoronalColors->SetLookupTable( 0 );
    m_AxialColors->SetLookupTable( 0 );
    }
  else
    {
    std::cerr << "Unsupported image type" << std::endl;
    return;
    }




  m_SagittalBlender->SetInput( 0, m_SagittalColors->GetOutput() );
  m_CoronalBlender->SetInput( 0, m_CoronalColors->GetOutput() );
  m_AxialBlender->SetInput( 0, m_AxialColors->GetOutput() );

  m_ThreeDRenderer->AddActor( m_SagittalActor );
  m_ThreeDRenderer->AddActor( m_AxialActor );
  m_ThreeDRenderer->AddActor( m_CoronalActor );

  m_SagittalRenderer->AddActor( m_SagittalActor );
  m_CoronalRenderer->AddActor( m_CoronalActor );
  m_AxialRenderer->AddActor( m_AxialActor );


  // Calculate default slice location and set the viewer to it
  for ( int i = 0; i < 3; i++ )
    {
    m_MaximumImageIndex[ i ] = m_Image->GetLargestPossibleRegion().GetSize()[ i ] - 1;
    }

  unsigned int  sagittalSliceNumber = static_cast< unsigned int >( ( m_Image->GetLargestPossibleRegion().GetSize()[ 0 ] - 1 ) / 2.0f );
  unsigned int  coronalSliceNumber = static_cast< unsigned int >( ( m_Image->GetLargestPossibleRegion().GetSize()[ 1 ] - 1 ) / 2.0f );
  unsigned int  axialSliceNumber = static_cast< unsigned int >( ( m_Image->GetLargestPossibleRegion().GetSize()[ 2 ] - 1 ) / 2.0f );
  this->SetSliceLocation( sagittalSliceNumber, coronalSliceNumber, axialSliceNumber );

  // Set up the 3D view camera based on the dimensions of the dataset
  m_ThreeDRenderer->GetActiveCamera()->SetFocalPoint( ( m_Image->GetLargestPossibleRegion().GetSize()[ 0 ] - 1 ) / 2.0f,
                                 ( m_Image->GetLargestPossibleRegion().GetSize()[ 1 ] - 1 ) / 2.0f,
                                 ( m_Image->GetLargestPossibleRegion().GetSize()[ 2 ] - 1 ) / 2.0f );
  m_ThreeDRenderer->GetActiveCamera()->SetPosition( ( m_Image->GetLargestPossibleRegion().GetSize()[ 0 ] - 1 ) * 10.0f,
                               ( m_Image->GetLargestPossibleRegion().GetSize()[ 1 ] - 1 ) * 10.0f,
                               ( m_Image->GetLargestPossibleRegion().GetSize()[ 2 ] - 1 ) * 10.0f );
  m_ThreeDRenderer->GetActiveCamera()->SetViewUp( 0, 0, 1 );

/*  if ( reset3DCamera )
    {*/
    // Reset the camera for the 3D renderer: this automatically sets up the camera based on the visible actors.
    m_ThreeDRenderer->ResetCamera();
//     }

}



//
void
ImageViewer
::SetOverlayImage( const ImageBaseType* overlayImage )
{
  // Can't have an overlay image without an image
  if ( !m_Image )
    {
    return;
    }

  m_OverlayImage = overlayImage;

  // Connect the VTK pipeline bits and pieces
  if ( dynamic_cast< const ImageType* >( m_OverlayImage.GetPointer() ) )
    {
    ImageType::ConstPointer castOverlayImage = static_cast< const ImageType* >( m_OverlayImage.GetPointer() );
    m_OverlayImageExporter = ExporterType::New();
    m_OverlayImageExporter->SetInput( castOverlayImage );
    ConnectPipelines( m_OverlayImageExporter, m_OverlayImporter );

    this->SetOverlayImageLookupTable( m_OverlayImageLookupTable );

    }
  else if ( dynamic_cast< const RGBAImageType* >( m_OverlayImage.GetPointer() ) )
    {
    RGBAImageType::ConstPointer castOverlayImage = static_cast< const RGBAImageType* >( m_OverlayImage.GetPointer() );
    m_RGBAOverlayImageExporter = RGBAExporterType::New();
    m_RGBAOverlayImageExporter->SetInput( castOverlayImage );
    ConnectPipelines( m_RGBAOverlayImageExporter, m_OverlayImporter );

    this->SetOverlayImageLookupTable( 0 );
    }
  else
    {
    std::cerr << "Unsupported image type" << std::endl;
    return;
    }

  m_SagittalBlender->SetInput( 1, m_SagittalOverlayColors->GetOutput() );
  m_CoronalBlender->SetInput( 1, m_CoronalOverlayColors->GetOutput() );
  m_AxialBlender->SetInput( 1, m_AxialOverlayColors->GetOutput() );


}



//
//
//
void
ImageViewer
::SetMesh( const AtlasMesh* mesh )
{
  m_Mesh = mesh;

  if ( !m_Mesh )
    {
    m_ThreeDRenderer->RemoveActor( m_OutlineActor );

    m_ThreeDRenderer->RemoveActor( m_SagittalCutActor );
    m_ThreeDRenderer->RemoveActor( m_CoronalCutActor );
    m_ThreeDRenderer->RemoveActor( m_AxialCutActor );

    m_ThreeDRenderer->RemoveActor( m_AxialEdgeActor );

    m_SagittalRenderer2->RemoveActor( m_SagittalCutActor );
    m_CoronalRenderer2->RemoveActor( m_CoronalCutActor );
    m_AxialRenderer2->RemoveActor( m_AxialCutActor );
    
    return;
    }

  // Convert into VTK object
  vtkSmartPointer< vtkUnstructuredGrid >  vGrid = this->GetVTKUnstructedGrid( m_Mesh );

  // Connect the VTK pipeline bits and pieces
  m_OutlineFilter->SetInput( vGrid );
  m_EdgeExtracter->SetInput( vGrid );

  m_SagittalCutter->SetInput( vGrid );
  m_CoronalCutter->SetInput( vGrid );
  m_AxialCutter->SetInput( vGrid );

  //m_ThreeDRenderer->AddActor( m_OutlineActor );
  m_ThreeDRenderer->AddActor( m_AxialCutActor );
  m_ThreeDRenderer->AddActor( m_AxialEdgeActor );

  m_SagittalRenderer2->AddActor( m_SagittalCutActor );
  m_CoronalRenderer2->AddActor( m_CoronalCutActor );
  m_AxialRenderer2->AddActor( m_AxialCutActor );

}



//
//
//
void
ImageViewer
::SetOverlayAlpha( float overlayAlpha )
{
  m_OverlayAlpha = overlayAlpha;

  m_SagittalBlender->SetOpacity( 1, m_OverlayAlpha );
  m_CoronalBlender->SetOpacity( 1, m_OverlayAlpha );
  m_AxialBlender->SetOpacity( 1, m_OverlayAlpha );

}


//
//
//
void
ImageViewer
::SetSliceLocation( unsigned int sagittalSliceNumber, unsigned int coronalSliceNumber, unsigned int axialSliceNumber )
{
  m_SagittalSliceNumber = sagittalSliceNumber;
  m_CoronalSliceNumber = coronalSliceNumber;
  m_AxialSliceNumber = axialSliceNumber;

  if ( !m_Image )
    {
    return;
    }

  m_SagittalActor->SetDisplayExtent( m_SagittalSliceNumber, m_SagittalSliceNumber,
                                     0, m_MaximumImageIndex[ 1 ],
                                     0, m_MaximumImageIndex[ 2 ] );


  m_CoronalActor->SetDisplayExtent( 0, m_MaximumImageIndex[ 0 ],
                                    m_CoronalSliceNumber, m_CoronalSliceNumber,
                                    0, m_MaximumImageIndex[ 2 ] );

  m_AxialActor->SetDisplayExtent( 0, m_MaximumImageIndex[ 0 ],
                                  0, m_MaximumImageIndex[ 1 ],
                                  m_AxialSliceNumber, m_AxialSliceNumber );

  m_SagittalPlane->SetOrigin( sagittalSliceNumber * m_Image->GetSpacing()[ 0 ], 0, 0 );
  m_CoronalPlane->SetOrigin( 0, coronalSliceNumber * m_Image->GetSpacing()[ 1 ], 0 );
  //m_AxialPlane->SetOrigin( 0, 0, axialSliceNumber * m_Image->GetSpacing()[ 2 ] );
  m_AxialPlane->SetOrigin( 0, 0, ( axialSliceNumber + 0.1 ) * m_Image->GetSpacing()[ 2 ] );

  const float  scale = vnl_math_max( vnl_math_max( m_MaximumImageIndex[ 0 ] / 2.0, m_MaximumImageIndex[ 1 ] / 2.0 ), m_MaximumImageIndex[ 2 ] / 2.0 );
  m_SagittalCamera->SetFocalPoint( m_SagittalSliceNumber,
                                   m_MaximumImageIndex[ 1 ] / 2.0, m_MaximumImageIndex[ 2 ] / 2.0 );
  m_SagittalCamera->ParallelProjectionOn();
  //m_SagittalCamera->SetParallelScale( vnl_math_max( m_MaximumImageIndex[ 1 ] / 2.0, m_MaximumImageIndex[ 2 ] / 2.0 ) );
  m_SagittalCamera->SetParallelScale( scale );
  //m_SagittalCamera->SetPosition( 0, m_MaximumImageIndex[ 1 ] / 2.0, m_MaximumImageIndex[ 2 ] / 2.0 );
  m_SagittalCamera->SetPosition( m_MaximumImageIndex[ 0 ], m_MaximumImageIndex[ 1 ] / 2.0, m_MaximumImageIndex[ 2 ] / 2.0 );


  m_CoronalCamera->SetFocalPoint( m_MaximumImageIndex[ 0 ] / 2.0, m_CoronalSliceNumber,
                                  m_MaximumImageIndex[ 2 ] / 2.0 );
  m_CoronalCamera->ParallelProjectionOn();
  //m_CoronalCamera->SetParallelScale( vnl_math_max( m_MaximumImageIndex[ 0 ] / 2.0, m_MaximumImageIndex[ 2 ] / 2.0 ) );
  m_CoronalCamera->SetParallelScale( scale );
  m_CoronalCamera->SetPosition( m_MaximumImageIndex[ 0 ] / 2.0, 0, m_MaximumImageIndex[ 2 ] / 2.0 );


  m_AxialCamera->SetFocalPoint( m_MaximumImageIndex[ 0 ] / 2.0, m_MaximumImageIndex[ 1 ] / 2.0,
                                m_AxialSliceNumber );
  m_AxialCamera->ParallelProjectionOn();
  //m_AxialCamera->SetParallelScale( vnl_math_max( m_MaximumImageIndex[ 0 ] / 2.0, m_MaximumImageIndex[ 1 ] / 2.0 ) );
  m_AxialCamera->SetParallelScale( scale );
  m_AxialCamera->SetPosition( m_MaximumImageIndex[ 0 ] / 2.0, m_MaximumImageIndex[ 1 ] / 2.0, 0 );
  //m_AxialCamera->SetPosition( m_MaximumImageIndex[ 0 ] / 2.0, m_MaximumImageIndex[ 1 ] / 2.0, m_MaximumImageIndex[ 2 ] + 1 );

}



//
//
//
void
ImageViewer
::SetScale( float scale )
{
  //std::cout << "Setting imageLookupTable table range: 0 -> " << 255.0f / scale << std::endl;
  m_ImageLookupTable->SetTableRange( 0, 255.0f / scale );

}

//
//
//
void
ImageViewer
::SetOverlayScale( float overlayScale )
{

  m_OverlayImageLookupTable->SetTableRange( 0, 255.0f / overlayScale );

}






//
//
//
void 
ImageViewer
::SetScaleToFillRange()
{

  if ( !m_Image )
    return;

  if ( dynamic_cast< const ImageType* >( m_Image.GetPointer() ) )
    {
    ImageType::ConstPointer castImage = static_cast< const ImageType* >( m_Image.GetPointer() );
    this->SetScale( 255.0f / this->CalculateMaximum( castImage ) );
    }

}


//
//
//
void 
ImageViewer
::SetOverlayScaleToFillRange()
{

  if ( !m_OverlayImage )
    return;

  if ( dynamic_cast< const ImageType* >( m_OverlayImage.GetPointer() ) )
    {
    ImageType::ConstPointer castOverlayImage = static_cast< const ImageType* >( m_OverlayImage.GetPointer() );

    this->SetOverlayScale( 255.0f / this->CalculateMaximum( castOverlayImage ) );
    }

}

//
//
//
float 
ImageViewer
::CalculateMaximum( const ImageType* image ) const
{
  typedef itk::MinimumMaximumImageCalculator< ImageType >  RangeCalculatorType;
  RangeCalculatorType::Pointer  rangeCalculator = RangeCalculatorType::New();
  rangeCalculator->SetImage( image );
  rangeCalculator->ComputeMaximum();

  return static_cast< float >( rangeCalculator->GetMaximum() );

}






//
//
//
vtkSmartPointer< vtkUnstructuredGrid >
ImageViewer
::GetVTKUnstructedGrid( const kvl::AtlasMesh* mesh ) const
{
  // Construct a VTK points container
  // VTK is not as smart as ITK. As a result, we have to allocate a point container
  // as big as the largest index in ITK, even if some (or even most) of the intermediate
  // indices don't exist because of our mesh operations...
  vtkSmartPointer< vtkPoints >  vPoints = vtkSmartPointer< vtkPoints >::New();
  //vPoints->SetNumberOfPoints( mesh->GetNumberOfPoints() );
  kvl::AtlasMesh::PointsContainer::ConstIterator  endIt = mesh->GetPoints()->End();
  --endIt;
  //std::cout << "Highest index in point container: " << endIt->Index() << std::endl;
  vPoints->SetNumberOfPoints( endIt->Index() + 1 );
  for( kvl::AtlasMesh::PointsContainer::ConstIterator it = mesh->GetPoints()->Begin();
       it != mesh->GetPoints()->End(); ++it )
    {
    vPoints->SetPoint( it->Index(), it->Value().GetDataPointer() );
    }

  //std::cout << "Converted points" << std::endl;

  //std::cout << "mesh->GetNumberOfCells(): " << mesh->GetNumberOfCells() << std::endl;

  // Construct a VTK type array and cell array
  int*  vTypes = new int[ mesh->GetNumberOfCells() ];
  vtkSmartPointer< vtkCellArray >  vCells = vtkSmartPointer< vtkCellArray >::New();
#if 1
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
#else
  vCells->EstimateSize( mesh->GetNumberOfCells(), 2 );
  int vCellCount = 0;
  for ( kvl::AtlasMesh::CellsContainer::ConstIterator  it = mesh->GetCells()->Begin(); 
        it != mesh->GetCells()->End(); ++it )
    {
    kvl::AtlasMesh::CellType*  cell = it.Value();
    if( cell->GetType() != kvl::AtlasMesh::CellType::LINE_CELL )
      {
      continue;
      }

    // Insert a line cell
    //vCells->InsertNextCell( 3,  static_cast< vtkIdType* >( cell->PointIdsBegin() ) );
    vCells->InsertNextCell( 2,  (vtkIdType*)cell->PointIdsBegin() );
    vTypes[ vCellCount ] = VTK_LINE;
    vCellCount++;

    } // End loop over all cells in the mesh
#endif

  //std::cout << "Converted cells and types" << std::endl;


  // Put everything together into a VTK unstructed grid
  vtkSmartPointer< vtkUnstructuredGrid >  vGrid = vtkSmartPointer< vtkUnstructuredGrid >::New();
  vGrid->SetPoints( vPoints );
  vGrid->SetCells( vTypes, vCells );

  //std::cout << "Put everything" << std::endl;
  delete[] vTypes;

  return vGrid;
}



//
//
//
void
ImageViewer
::LookAt( int quadrantNumber )
{

  switch ( quadrantNumber )
    {
    case 1:
      m_SagittalRenderer->SetViewport( 0, 0, 0, 0 );
      m_SagittalRenderer2->SetViewport( 0, 0, 0, 0 );
      m_CoronalRenderer->SetViewport( 0, 0, 1, 1 );
      m_CoronalRenderer2->SetViewport( 0, 0, 1, 1 );
      m_AxialRenderer->SetViewport( 0, 0, 0, 0 );
      m_AxialRenderer2->SetViewport( 0, 0, 0, 0 );
      m_ThreeDRenderer->SetViewport( 0, 0, 0, 0 );
      break;
    case 2:
      m_SagittalRenderer->SetViewport( 0, 0, 1, 1 );
      m_SagittalRenderer2->SetViewport( 0, 0, 1, 1 );
      m_CoronalRenderer->SetViewport( 0, 0, 0, 0 );
      m_CoronalRenderer2->SetViewport( 0, 0, 0, 0 );
      m_AxialRenderer->SetViewport( 0, 0, 0, 0 );
      m_AxialRenderer2->SetViewport( 0, 0, 0, 0 );
      m_ThreeDRenderer->SetViewport( 0, 0, 0, 0 );
      break;
    case 3:
      m_SagittalRenderer->SetViewport( 0, 0, 0, 0 );
      m_SagittalRenderer2->SetViewport( 0, 0, 0, 0 );
      m_CoronalRenderer->SetViewport( 0, 0, 0, 0 );
      m_CoronalRenderer2->SetViewport( 0, 0, 0, 0 );
      m_AxialRenderer->SetViewport( 0, 0, 1, 1 );
      m_AxialRenderer2->SetViewport( 0, 0, 1, 1 );
      m_ThreeDRenderer->SetViewport( 0, 0, 0, 0 );
      break;
    case 4:
      m_SagittalRenderer->SetViewport( 0, 0, 0, 0 );
      m_SagittalRenderer2->SetViewport( 0, 0, 0, 0 );
      m_CoronalRenderer->SetViewport( 0, 0, 0, 0 );
      m_CoronalRenderer2->SetViewport( 0, 0, 0, 0 );
      m_AxialRenderer->SetViewport( 0, 0, 0, 0 );
      m_AxialRenderer2->SetViewport( 0, 0, 0, 0 );
      m_ThreeDRenderer->SetViewport( 0, 0, 1, 1 );
      break;
    default:
      m_SagittalRenderer->SetViewport( 0.5, 0.5, 1, 1 );
      m_SagittalRenderer2->SetViewport( 0.5, 0.5, 1, 1 );
      m_CoronalRenderer->SetViewport( 0, 0.5, 0.5, 1 );
      m_CoronalRenderer2->SetViewport( 0, 0.5, 0.5, 1 );
      m_AxialRenderer->SetViewport( 0, 0, 0.5, 0.5 );
      m_AxialRenderer2->SetViewport( 0, 0, 0.5, 0.5 );
      m_ThreeDRenderer->SetViewport( 0.5, 0, 1, 0.5 );
    }


}



//
//
//
void
ImageViewer
::SetOverlayImageLookupTable( vtkLookupTable*  overlayImageLookupTable )
{
  m_OverlayImageLookupTable = overlayImageLookupTable;
  
  m_SagittalOverlayColors->SetLookupTable( m_OverlayImageLookupTable );
  m_CoronalOverlayColors->SetLookupTable( m_OverlayImageLookupTable );
  m_AxialOverlayColors->SetLookupTable( m_OverlayImageLookupTable );

}



//
//
//
void
ImageViewer
::WriteScreenShot( const std::string&  fileName )
{

  vtkSmartPointer< vtkWindowToImageFilter >  windowToImage = vtkWindowToImageFilter::New();
  windowToImage->SetInput( this->GetRenderWindow() );

  vtkSmartPointer< vtkPNGWriter >  writer = vtkPNGWriter::New();
  writer->SetInput( windowToImage->GetOutput() );
  writer->SetFileName( fileName.c_str() );
  writer->Write();

  std::cout << "Wrote screenshot to file " << fileName << std::endl;

}


} // End namespace kvl

