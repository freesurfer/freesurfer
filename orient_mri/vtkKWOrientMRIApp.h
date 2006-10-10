#ifndef __vtkKWOrientMRIApp_h
#define __vtkKWOrientMRIApp_h

#include <string>
#include "vtkKWApplication.h"

class vtkKWOrientMRIWindow;

class vtkKWOrientMRIApp : public vtkKWApplication {

 public:

  static vtkKWOrientMRIApp* New ();
  vtkTypeRevisionMacro ( vtkKWOrientMRIApp, vtkKWApplication );
  
  // Override to show our window.
  virtual void Start ( int argc, char* argv[] );

  // Override to save/load our prefs.
  virtual void RestoreApplicationSettingsFromRegistry ();
  virtual void SaveApplicationSettingsToRegistry ();

  // Registry key constants.
  static const char* sMainWindowWidthRegKey;
  static const char* sMainWindowHeightRegKey;

  void LoadVolume ( const char* ifnVolume );

 protected:
  
  vtkKWOrientMRIApp ();
  ~vtkKWOrientMRIApp ();

  vtkKWOrientMRIWindow* mWindow;

  int mMainWindowWidth;
  int mMainWindowHeight;
};

#endif
