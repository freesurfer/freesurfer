/**
 * @file  kvlAtlasMeshViewingConsole.h
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
  void SelectTriangleContainingPoint( float x, float y );

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
