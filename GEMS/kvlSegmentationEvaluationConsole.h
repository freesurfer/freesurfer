/**
 * @file  kvlSegmentationEvaluationConsole.h
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
#ifndef __kvlSegmentationEvaluationConsole_h
#define __kvlSegmentationEvaluationConsole_h

#include "kvlSegmentationEvaluationConsoleGUI.h"
#include "kvlAtlasMeshCollection.h"


namespace kvl
{

class SegmentationEvaluationConsole : public kvlSegmentationEvaluationConsoleGUI
{

public:

  //
  SegmentationEvaluationConsole();

  //
  virtual ~SegmentationEvaluationConsole();

  //
  void LoadImage( const char* fileName );

  //
  void LoadOverlayImage1( const char* fileName );

  //
  void LoadOverlayImage2( const char* fileName );

  //
  void LoadOverlayLookupTable( const char* fileName );

  //
  void Show();

  //
  void ShowSelectedView();

  //
  void Swap();

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

};



} // end namespace kvl

#endif
