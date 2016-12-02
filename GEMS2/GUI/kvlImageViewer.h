#ifndef __kvlImageViewer_h
#define __kvlImageViewer_h


#include "vtkFlRenderWindowInteractor.h"
#include "itkImage.h"
#include "itkVTKImageExport.h"
#include "kvlAtlasMesh.h"
#include "vtkSmartPointer.h"
#include "vtkOutlineFilter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkOutlineFilter.h"
#include "vtkExtractEdges.h"
#include "vtkCutter.h"
#include "vtkLookupTable.h"
#include "vtkImageBlend.h"
#include "vtkImageActor.h"
#include "vtkActor.h"
#include "vtkPlane.h"
#include "vtkCamera.h"
#include "vtkImageMapToColors.h"
#include "itkRGBAPixel.h"


class vtkImageImport;

namespace kvl
{



class ImageViewer : public vtkFlRenderWindowInteractor
{
public:

  typedef itk::ImageBase< 3 >  ImageBaseType;
  typedef itk::Image< unsigned char, 3 >  ImageType;
  typedef itk::Image< itk::RGBAPixel< unsigned char >, 3 >  RGBAImageType;

  ImageViewer( int x, int y, int w, int h, const char *l=0 );
  ~ImageViewer();

  void SetImage( const ImageBaseType* image );

  const ImageBaseType*  GetImage() const
    { return m_Image; }

  void SetScale( float scale );

  void SetScaleToFillRange();

  void SetOverlayImage( const ImageBaseType* overlayImage );

  const ImageBaseType*  GetOverlayImage() const
    { return m_OverlayImage; }

  void SetOverlayScale( float overlayScale );

  void SetOverlayScaleToFillRange();

  void SetOverlayAlpha( float overlayAlpha );

  float  GetOverlayAlpha() const
    { return m_OverlayAlpha; }

  void SetMesh( const AtlasMesh* mesh );

  const AtlasMesh* GetMesh() const
    { return m_Mesh; }

  void  SetSliceLocation( unsigned int sagittalSliceNumber, unsigned int coronalSliceNumber, unsigned int axialSliceNumber );

  void  GetSliceLocation( unsigned int& sagittalSliceNumber, unsigned int& coronalSliceNumber, unsigned int& axialSliceNumber )
    {
    sagittalSliceNumber = m_SagittalSliceNumber;
    coronalSliceNumber = m_CoronalSliceNumber;
    axialSliceNumber = m_AxialSliceNumber;
    }

  const int*  GetMaximumImageIndex() const
    {
    return m_MaximumImageIndex;
    }

  void SetPositionGradient( const AtlasPositionGradientContainerType* positionGradient ) {}

  vtkLookupTable*  GetOverlayImageLookupTable()
    { return m_OverlayImageLookupTable; }

  void  SetOverlayImageLookupTable( vtkLookupTable*  overlayImageLookupTable );

  void  LookAt( int quadrantNumber = 5 );

  void  WriteScreenShot( const std::string&  fileName );

protected :

  float CalculateMaximum( const ImageType* image ) const;

  //
  typedef itk::VTKImageExport< ImageType > ExporterType;
  typedef itk::VTKImageExport< RGBAImageType >  RGBAExporterType;

  //
  vtkSmartPointer< vtkUnstructuredGrid >  GetVTKUnstructedGrid( const kvl::AtlasMesh* mesh ) const;


private :

  ImageBaseType::ConstPointer  m_Image;
  ImageBaseType::ConstPointer  m_OverlayImage;
  AtlasMesh::ConstPointer  m_Mesh;


  float  m_OverlayColor[3];
  float  m_OverlayAlpha;

  int  m_MaximumImageIndex[ 3 ];

  unsigned int  m_SagittalSliceNumber;
  unsigned int  m_AxialSliceNumber;
  unsigned int  m_CoronalSliceNumber;

  ExporterType::Pointer  m_ImageExporter;
  RGBAExporterType::Pointer  m_RGBAImageExporter;
  ExporterType::Pointer  m_OverlayImageExporter;
  RGBAExporterType::Pointer  m_RGBAOverlayImageExporter;

  vtkSmartPointer< vtkImageImport >  m_Importer;
  vtkSmartPointer< vtkImageImport >  m_OverlayImporter;

  vtkSmartPointer< vtkOutlineFilter >  m_OutlineFilter;
  vtkSmartPointer< vtkExtractEdges >  m_EdgeExtracter;

  vtkSmartPointer< vtkCutter >  m_SagittalCutter;
  vtkSmartPointer< vtkCutter >  m_CoronalCutter;
  vtkSmartPointer< vtkCutter >  m_AxialCutter;

  vtkSmartPointer< vtkLookupTable >  m_ImageLookupTable;
  vtkSmartPointer< vtkLookupTable >  m_OverlayImageLookupTable;

  vtkSmartPointer< vtkImageMapToColors >  m_SagittalColors;
  vtkSmartPointer< vtkImageMapToColors >  m_CoronalColors;
  vtkSmartPointer< vtkImageMapToColors >  m_AxialColors;

  vtkSmartPointer< vtkImageMapToColors >  m_SagittalOverlayColors;
  vtkSmartPointer< vtkImageMapToColors >  m_CoronalOverlayColors;
  vtkSmartPointer< vtkImageMapToColors >  m_AxialOverlayColors;

  vtkSmartPointer< vtkImageBlend >  m_SagittalBlender;
  vtkSmartPointer< vtkImageBlend >  m_CoronalBlender;
  vtkSmartPointer< vtkImageBlend >  m_AxialBlender;

  vtkSmartPointer< vtkImageActor >  m_SagittalActor;
  vtkSmartPointer< vtkImageActor >  m_CoronalActor;
  vtkSmartPointer< vtkImageActor >  m_AxialActor;

  vtkSmartPointer< vtkPlane >  m_SagittalPlane;
  vtkSmartPointer< vtkPlane >  m_CoronalPlane;
  vtkSmartPointer< vtkPlane >  m_AxialPlane;

  vtkSmartPointer< vtkActor >  m_OutlineActor;

  vtkSmartPointer< vtkActor >  m_SagittalCutActor;
  vtkSmartPointer< vtkActor >  m_CoronalCutActor;
  vtkSmartPointer< vtkActor >  m_AxialCutActor;

  vtkSmartPointer< vtkActor >  m_AxialEdgeActor;

  vtkSmartPointer< vtkCamera >  m_SagittalCamera;
  vtkSmartPointer< vtkCamera >  m_CoronalCamera;
  vtkSmartPointer< vtkCamera >  m_AxialCamera;

  vtkSmartPointer< vtkRenderer >  m_ThreeDRenderer;

  vtkSmartPointer< vtkRenderer >  m_SagittalRenderer;
  vtkSmartPointer< vtkRenderer >  m_CoronalRenderer;
  vtkSmartPointer< vtkRenderer >  m_AxialRenderer;

  vtkSmartPointer< vtkRenderer >  m_SagittalRenderer2;
  vtkSmartPointer< vtkRenderer >  m_CoronalRenderer2;
  vtkSmartPointer< vtkRenderer >  m_AxialRenderer2;


};


} // End namespace kvl

#endif


