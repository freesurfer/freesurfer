/**
 * @brief Command line parsing and registry
 *
 * Application code that parses the command line options and passes
 * them to the main window. Also gets registry options pertaining to
 * the window size.
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


#ifndef __vtkKWQdecApp_h
#define __vtkKWQdecApp_h

#include <string>
#include "vtkKWApplication.h"
#include "vtkSmartPointer.h"

class vtkKWDialog;
class vtkKWQdecWindow;
class vtkKWTopLevel;

class vtkKWQdecApp : public vtkKWApplication {

public:

  static vtkKWQdecApp* New ();
  vtkTypeRevisionMacro ( vtkKWQdecApp, vtkKWApplication );

  // Override to show our window.
  virtual void Start ( int argc, char* argv[] );
  
  // Override to clean up.
  virtual int Exit ();

  // Override to add text to 'Help About' box
  virtual void AddAboutText ( ostream &os );

  // Thse are just passed right to the window.
  void LoadDataTable ( const char* ifnDataTable );
  void LoadProjectFile ( const char* ifnProject );
  void LoadSurface ( const char* ifnSurface );
  void LoadGDFFile ( const char* ifnGDF );
  void LoadSurfaceScalars ( const char* ifnScalars );
  void LoadSurfaceCurvatureScalars ( const char* ifnScalars );
  void LoadAnnotation ( const char* ifnAnnotation );
  void LoadSurfaceOverlayScalars ( const char* ifnScalars,
				   const char* ifnColors );
  void LoadLabel ( const char* ifnLabel );
  void SetAverageSubject ( const char* isAvgSubj );

  // Display an error message dialog and log the message.
  virtual void ErrorMessage ( const char* isMessage );
  
  // Show our help dialog.
  void DisplayHelpDialog ( vtkKWTopLevel* iTop );

protected:

  vtkKWQdecApp ();
  ~vtkKWQdecApp ();

  //BTX

  // The main application window.
  vtkSmartPointer<vtkKWQdecWindow> mWindow;

  // Help dialog.
  vtkSmartPointer<vtkKWDialog> mDlogHelp;

  //ETX
};

#endif
