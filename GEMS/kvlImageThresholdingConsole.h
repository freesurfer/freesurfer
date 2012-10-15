/**
 * @file  kvlImageThresholdingConsole.h
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
