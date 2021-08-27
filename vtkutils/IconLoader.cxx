/**
 * @brief Loads TIF and GIF icons and sets KW Menus and Pushbuttons
 *
 * Allows you to load icons from a text file list or at run time from
 * TIF and GIF files, and lets you set menu icons and button icons.
 * 
 */
/*
 * Original Author: Kevin Teich
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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


#include <sstream>
#include <string>
#include <stdexcept>
#include "tiffio.h"
#include "IconLoader.h"
#include "vtkKWApplication.h"
#include "vtkKWCheckButton.h"
#include "vtkKWMenu.h"
#include "vtkKWPushButton.h"

using namespace std;

vtkKWApplication* IconLoader::mApp;
map<string,bool> IconLoader::mabTkIconLoaded;
map<string,unsigned char*> IconLoader::maTIFFData;
map<string,int> IconLoader::maWidth;
map<string,int> IconLoader::maHeight;

void
IconLoader::Initialize ( vtkKWApplication* iApp ) {

  mApp = iApp;
  mabTkIconLoaded.clear();
  maTIFFData.clear();
  maWidth.clear();
  maHeight.clear();
}

void
IconLoader::ShutDown () {

  map<string,unsigned char*>::iterator tTIFFData;
  for ( tTIFFData = maTIFFData.begin();
	tTIFFData != maTIFFData.end();
	++tTIFFData )
    if ( NULL != tTIFFData->second )
      _TIFFfree( tTIFFData->second );
  
}

void
IconLoader::SetPushButtonIcon ( const char* isKey,
				vtkKWPushButton* iButton ) {

  map<string,unsigned char*>::iterator tTIFFData = maTIFFData.find( isKey );
  if( maTIFFData.end() == tTIFFData )
    throw runtime_error( string("Icon key ") + 
			 string(isKey) + 
			 string(" not found in data map.") );

  map<string,int>::iterator tWidth = maWidth.find( isKey );
  if( maWidth.end() == tWidth )
    throw runtime_error( string("Icon key ") + 
			 string(isKey) + 
			 string(" not found in width map.") );

  map<string,int>::iterator tHeight = maHeight.find( isKey );
  if( maWidth.end() == tHeight )
    throw runtime_error( string("Icon key ") + 
			 string(isKey) + 
			 string(" not found in height map.") );

  if ( NULL == tTIFFData->second )
    throw runtime_error( string("Icon with key ") + 
			 string(isKey) + 
			 string(" had NULL data.") );

  iButton->SetImageToPixels( tTIFFData->second, 
			     tWidth->second, tHeight->second, 4, 0 );
}

void
IconLoader::SetCheckButtonIcon ( const char* isKey,
				 vtkKWCheckButton* iButton ) {

  map<string,unsigned char*>::iterator tTIFFData = maTIFFData.find( isKey );
  if( maTIFFData.end() == tTIFFData )
    throw runtime_error( string("Icon key ") + 
			 string(isKey) + 
			 string(" not found in data map.") );

  map<string,int>::iterator tWidth = maWidth.find( isKey );
  if( maWidth.end() == tWidth )
    throw runtime_error( string("Icon key ") + 
			 string(isKey) + 
			 string(" not found in width map.") );

  map<string,int>::iterator tHeight = maHeight.find( isKey );
  if( maWidth.end() == tHeight )
    throw runtime_error( string("Icon key ") + 
			 string(isKey) + 
			 string(" not found in height map.") );

  if ( NULL == tTIFFData->second )
    throw runtime_error( string("Icon with key ") + 
			 string(isKey) + 
			 string(" had NULL data.") );

  iButton->SetImageToPixels( tTIFFData->second, 
			     tWidth->second, tHeight->second, 4, 0 );
}

void
IconLoader:: SetMenuItemIcon ( const char* isKey,
			       vtkKWMenu* iMenu, int inItem ) {

  map<string,bool>::iterator tbLoaded = mabTkIconLoaded.find( isKey );
  if( mabTkIconLoaded.end() == tbLoaded )
    throw runtime_error( string("Icon key ") + 
			 string(isKey) + 
			 string(" not found in loaded map.") );
  
  if ( !tbLoaded->second )
    throw runtime_error( string("Icon with key ") + 
			 string(isKey) + 
			 string(" is not loaded.") );
  
  iMenu->SetItemImage( inItem, isKey );
}
 
int
IconLoader::LoadIconsFromFile ( const char* ifn ) {

  int errs=0;

  // Try to open the file.
  ifstream fIcons( ifn, ios::in );
  if( !fIcons || fIcons.bad() ) {
    throw runtime_error( string("Couldn't read icons file ") + ifn );
  }

  // Go through the file.
  int nLine = 1;
  string sLine;
  while( !fIcons.eof() ) {

    // Get a line.
    getline( fIcons, sLine );

    // If that read gave us a EOL, we're done.
    if( fIcons.eof() )
      break;

    // Try to get three strings out of this line.
    stringstream ssLine( sLine );
    string sKey, sfnTIFF, sfnGIF;
    if( ssLine >> sKey >> sfnTIFF >> sfnGIF ) {

      // Try to load the icon with the strings that we got.
      try { 
        //cout << "LoadIcon: " << sKey.c_str() << ", " 
        //   << sfnTIFF.c_str() << ", " 
        //   << sfnGIF.c_str() << endl;
        errs += LoadIcon( sKey.c_str(), sfnTIFF.c_str(), sfnGIF.c_str() );
      } 
      catch( exception& e ) {
        cerr << "Error reading " << ifn << ", line " << nLine << ": "
             << e.what() << endl;
        errs++;
      }

    } else {
      
      cerr << "Error reading " << ifn << ", line " << nLine 
           << ": Malformed line, unexpected number of entries." << endl;
      errs++;
    }

    nLine++;
  }

  fIcons.close();

  return errs;
}

int
IconLoader::LoadIcon ( const char* isKey,
		       const char* ifnTIFF,
		       const char* ifnGIF ) {
  
  int errs=0;

  // We'll substitute any IMAGEDIR in the file name with these
  // strings.
  string fnImageDir;
  char* FREESURFER_HOME = getenv( "FREESURFER_HOME" );
  if (NULL == FREESURFER_HOME) fnImageDir= "../images";
  else fnImageDir = string(FREESURFER_HOME) + "/lib/images";

  string fnTIFF = ifnTIFF;
  string fnGIF = ifnGIF;

  size_t nFound = fnTIFF.find("IMAGEDIR");
  if( nFound != string::npos )
    fnTIFF.replace( nFound, strlen("IMAGEDIR"), fnImageDir );

  nFound = fnGIF.find("IMAGEDIR");
  if( nFound != string::npos )
    fnGIF.replace( nFound, strlen("IMAGEDIR"), fnImageDir );

  // Load the data icon.
  TIFF* tiff = TIFFOpen( fnTIFF.c_str(), "r" );
  if ( tiff ) {

    uint32 zImageWidth, zImageHeight;
    size_t cPixels;
    uint32* tempData;
    uint32* data;

    TIFFGetField( tiff, TIFFTAG_IMAGEWIDTH, &zImageWidth );
    TIFFGetField( tiff, TIFFTAG_IMAGELENGTH, &zImageHeight );
    cPixels = zImageWidth * zImageHeight;

    // OK, originally I used TIFFReadRGBAImageOriented to read the
    // image so it culd flip it on the way in, but libtiff on rh9
    // doesn't have that function, so we have to use the older
    // TIFFReadRGBAImage and flip it manually.
    tempData = (uint32*)_TIFFmalloc( cPixels * sizeof(uint32) );
    data = (uint32*)_TIFFmalloc( cPixels * sizeof(uint32) );
    if ( data ) {
//       if( TIFFReadRGBAImageOriented( tiff, zImageWidth, zImageHeight,
//          data, ORIENTATION_TOPLEFT, 1 ) ) {
      if ( TIFFReadRGBAImage( tiff, zImageWidth, zImageHeight,
                              tempData, 1 ) ) {

        // We read it into tempData, now copy it in reverse row order
        // into data.
        for ( int nRow = zImageHeight-1; nRow >= 0; nRow-- ) {
          memmove( (void*)&data[(zImageHeight-1-nRow) * zImageWidth],
                   (void*)&tempData[nRow * zImageWidth],
                   sizeof(uint32) * zImageWidth );
        } 

        maTIFFData[isKey] = (unsigned char*) data;
        maWidth[isKey] = zImageWidth;
        maHeight[isKey] = zImageHeight;

      } else {
        errs++;
        _TIFFfree( data );
        _TIFFfree( tempData );
      }
    } else errs++;

    _TIFFfree( tempData );
    TIFFClose( tiff );

  } else {
    errs++;
    string sError;
    sError = string("Couldn't open ") + fnTIFF + ", missing?";
    throw runtime_error( sError );
  }

  //  string sCmd = string("image create photo ") + sTkIconName + " -file [file join " + sFreeSurferHome + " lib images " + fnIconTk + "]";
  //  mApp->Script( sCmd.c_str() );

  mApp->Script( "image create photo %s -file %s", isKey, fnGIF.c_str());

  mabTkIconLoaded[isKey] = true;

  return errs;
}
