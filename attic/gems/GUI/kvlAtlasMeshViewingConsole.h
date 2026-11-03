#ifndef __kvlAtlasMeshViewingConsole_h
#define __kvlAtlasMeshViewingConsole_h

#include "kvlAtlasMeshViewingConsoleGUI.h"
#include "kvlAtlasMeshCollection.h"
#include "kvlCompressionLookupTable.h"


namespace kvl
{

class AtlasMeshViewingConsole : public kvlAtlasMeshViewingConsoleGUI
{

public: 

  //
  AtlasMeshViewingConsole();
  
  //
  virtual ~AtlasMeshViewingConsole();

  //
  void LoadMeshCollection( const char* fileName,
                           const int* templateSize = 0,
                           const std::string& backgroundImageFileName = std::string() );
  
  // 
  void Show();
  
  //
  void ShowSelectedView();

  //
  void  GetScreenShot();

  //
  void  GetScreenShotSeries();

  //
  void  DumpImage();

  //
  void  GetScreenShotSeries( int directionNumber );

protected:

  //
  void Draw();
  
  //
  void SetSliceLocation( unsigned int  sagittalSliceNumber,
                         unsigned int  coronalSliceNumber,
                         unsigned int  axialSliceNumber );

  //
  typedef kvl::ImageViewer::ImageType  ImageType;
  static ImageType::Pointer  ReadBackgroundImage( const std::string&  backgroundImageFileName );
  
private:  
  
  AtlasMeshCollection::Pointer  m_MeshCollection;
  //AtlasMesh::CellIdentifier  m_EdgeIdToHighlight;
  ImageType::Pointer  m_BackgroundImage;

  CompressionLookupTable::Pointer  m_Compressor;
};



} // end namespace kvl

#endif
