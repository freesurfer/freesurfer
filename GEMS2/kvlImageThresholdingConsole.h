#ifndef __kvlImageThresholdingConsole_h
#define __kvlImageThresholdingConsole_h

#include "kvlImageThresholdingConsoleGUI.h"
#include "itkBinaryThresholdImageFilter.h"


namespace kvl
{

class ImageThresholdingConsole : public kvlImageThresholdingConsoleGUI
{

public: 

  //
  ImageThresholdingConsole();

  //
  virtual ~ImageThresholdingConsole();

  //
  void LoadImages( const std::vector< std::string >& fileNames );

  //
  void SetImageToThreshold( int imageToThreshold );

  //
  void SetImageToMask( int imageToMask );

  //
  void Show();

  //
  void ShowSelectedView();

  //
  void SetThresholds( float lowerThreshold, float upperThreshold );

  //
  void MaskImage();

  //
  void WriteMask();

  //
  void ShowInverseMask( bool showInverseMask );

protected:

  //
  void Draw();

  //
  void SetOverlayOpacity( float overlayOpacity );

  //
  void SetSliceLocation( unsigned int  sagittalSliceNumber,
                         unsigned int  coronalSliceNumber,
                         unsigned int  axialSliceNumber );


private:

  typedef itk::Image< short, 3 >  ImageType;
  std::vector< ImageType::ConstPointer >  m_Images;

  typedef itk::BinaryThresholdImageFilter< ImageType, ImageViewer::ImageType >  ThresholderType;
  ThresholderType::Pointer  m_Thresholder;

};



} // end namespace kvl

#endif
