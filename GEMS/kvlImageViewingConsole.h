/**
 * @file  kvlImageViewingConsole.h
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
#ifndef __kvlImageViewingConsole_h
#define __kvlImageViewingConsole_h

#include "kvlImageViewingConsoleGUI.h"
#include "kvlAtlasMeshCollection.h"


namespace kvl
{

class ImageViewingConsole : public kvlImageViewingConsoleGUI
{

public:

  //
  ImageViewingConsole();

  //
  virtual ~ImageViewingConsole();

  //
  void LoadImage( const char* fileName );

  //
  void LoadOverlayImage( const char* fileName );

  //
  void LoadOverlayLookupTable( const char* fileName );

  //
  void LoadMesh( const char* fileName, int meshNumber=-1 );

  //
  void Show();

  //
  void ShowSelectedView();

  //
  void GetScreenShot();

  //
  void GetScreenShotSeries( int directionNumber );

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
