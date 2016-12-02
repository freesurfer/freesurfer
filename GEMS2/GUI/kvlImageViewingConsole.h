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
