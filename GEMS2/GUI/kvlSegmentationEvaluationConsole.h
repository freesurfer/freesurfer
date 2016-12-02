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
