/**
 * @file  kvlAtlasMeshSegmenterConsole.h
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
#ifndef __kvlAtlasMeshSegmenterConsole_h
#define __kvlAtlasMeshSegmenterConsole_h

#include "kvlAtlasMeshSegmenterConsoleGUI.h"
#include <vector>
#include <string>
#include "kvlAtlasMeshSegmentationDriver.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"


namespace kvl
{

class AtlasMeshSegmenterConsole : public kvlAtlasMeshSegmenterConsoleGUI
{

public:

  //
  AtlasMeshSegmenterConsole();

  //
  virtual ~AtlasMeshSegmenterConsole();

  //
  void  SetUp( const std::string& setUpFileName );

  //
  void Show();

  //
  void Segment( bool  useAffine );

  //
  void HandleEMSegmenterEvent( itk::Object* object, const itk::EventObject & event );

  //
  void HandleSegmenterEvent( itk::Object* object, const itk::EventObject & event );

  //
  void HandleDriverEvent( itk::Object* object, const itk::EventObject & event );

  //
  void SetGradientDescentStepSize( float stepSize );

  //
  void SetPositionUpdatingStopCriterion( float stopCriterion );

  //
  void ShowSelectedView();

  //
  void SetSliceLocation( unsigned int  sagittalSliceNumber,
                         unsigned int  coronalSliceNumber,
                         unsigned int  axialSliceNumber );

protected:

  //
  void DisplayLabel( unsigned int labelNumber );

  //
  void Step();

  //
  void Interrupt();

  //
  void Continue();

  //
  void  ApplyAffineParameters();

  //
  void  ResetAffineParameters();

  //
  typedef itk::Image< float, 3 >  ExponentiatedImageType;
  template< class TPixel >
  static ExponentiatedImageType::Pointer GetExponentiatedImage( const itk::Image< TPixel, 3 >* image )
  {
    ExponentiatedImageType::Pointer  exponentiatedImage = ExponentiatedImageType::New();
    exponentiatedImage->SetRegions( image->GetBufferedRegion() );
    exponentiatedImage->Allocate();

    itk::ImageRegionConstIterator< itk::Image< TPixel, 3 > >  imageIt( image,
        image->GetLargestPossibleRegion() );
    itk::ImageRegionIterator< ExponentiatedImageType >  expIt( exponentiatedImage,
        exponentiatedImage->GetLargestPossibleRegion() );
    for ( ; !imageIt.IsAtEnd(); ++imageIt, ++expIt )
    {
      if ( imageIt.Value() > 0 )
      {
        std::cout << "Exponentiating " << static_cast< ExponentiatedImageType::PixelType >( imageIt.Value() )
                  << " to " <<  exp( static_cast< ExponentiatedImageType::PixelType >( imageIt.Value() ) ) << std::endl;
        expIt.Value() = exp( static_cast< ExponentiatedImageType::PixelType >( imageIt.Value() ) );
      }
      else
      {
        expIt.Value() = 0;
      }
    }

    return exponentiatedImage;
  }


private:

  AtlasMeshSegmentationDriver::Pointer  m_Driver;

  std::string  m_ProgressLabel;

  bool  m_Interrupted;
  bool  m_Stepping;


};



} // end namespace kvl

#endif
