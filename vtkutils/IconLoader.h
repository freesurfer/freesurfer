/**
 * @file  IconLoader.h
 * @brief Loads TIF and GIF icons and sets KW Menus and Pushbuttons
 *
 * Allows you to load icons from a text file list or at run time from
 * TIF and GIF files, and lets you set menu icons and button icons.
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/03/01 19:56:45 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


// .NAME IconLoader - static interface for loading and setting icons
// .SECTION Description
// This class provides an enumeration of custom icons and can set the
// icons in KW classes (current vtkKWPushButtons and vtkKWMenu items)
// to those icons. Internally it loads the icons from .gif and .tif
// data. The locations and file names for those icons are hard coded
// within.

#ifndef IconLoader_h
#define IconLoader_h

#include <string>
#include <map>
class vtkKWApplication;
class vtkKWCheckButton;
class vtkKWMenu;
class vtkKWPushButton;

class IconLoader {

public:

  // Description:
  // Initialize the class with a pointer to an application so it can
  // use its Sript() function to set icons.
  static void Initialize ( vtkKWApplication* iApp );
  static void ShutDown ();

  // Description:
  // Load the icons from a file. The file should be a list of three
  // text strings: the key, the TIFF filename, and the GIF
  // filename. It will substitute IMAGEDIR in the filename with the
  // FREESURFER_HOME/lib/images, or a static file name can be used.
  // Ex:
  // home IMAGEDIR/home.tif IMAGEDIR/home.gif
  // back /usr/share/tif/back.tif /usr/share/gif/back.gif
  // forwards ./forward.tif ./forward.gif
  static void LoadIconsFromFile ( const char* ifn );
  
  // Description:
  // Call this to load icons at run time. The key is a string key (no
  // spaces) that can be used to refer to the icon later. A TIFF file
  // is used to load image data, and a GIF file will be loaded into
  // Tk. It will substitute IMAGEDIR/ in the filename with the
  // FREESURFER_HOME/lib/images/, or a static file name can be used.
  static void LoadIcon ( const char* isKey,
			 const char* ifnTIFF,
			 const char* ifnGIF );

  // Description:
  // Set the icon for a vtkKWPushButton. Internally, it will load the
  // icon if not loaded from a .tif resource, extract the pixel data,
  // and set the vtkKWPushButton's image data. Can throw an exception
  // if inIcon is invalid.
  static void SetPushButtonIcon ( const char* isKey,
				  vtkKWPushButton* iButton );

  // Description:
  // Set the icon for a vtkKWCheckButton. Internally, it will load the
  // icon if not loaded from a .tif resource, extract the pixel data,
  // and set the vtkKWPushButton's image data. Can throw an exception
  // if inIcon is invalid.
  static void SetCheckButtonIcon ( const char* isKey,
				   vtkKWCheckButton* iButton );

  // Description:
  // Set the icon for a menu item. Interally, it will load the icon as
  // a Tk image, and set vtkKWMenu item's icon using that Tk image
  // name. Can throw an exception if inIcon is invalid.
  static void SetMenuItemIcon ( const char* isKey,
				vtkKWMenu* iMenu, int inItem );

protected:

  static vtkKWApplication* mApp;

  static std::map<std::string,bool> mabTkIconLoaded;
  static std::map<std::string,unsigned char*> maTIFFData;
  static std::map<std::string,int> maWidth;
  static std::map<std::string,int> maHeight;
};


#endif
