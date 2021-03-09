/**
 * @brief A KWWidgets progress dialog
 *
 * A simple KWWidgets dialog box with a progress bar tht listens to
 * update commands.
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


#ifndef vtkKWProgressDialog_h
#define vtkKWProgressDialog_h

#include "vtkSetGet.h"
#include "vtkCommand.h"

class vtkObject;
class vtkKWDialog;
class vtkKWProgressGauge;
class vtkKWApplication;

class vtkKWProgressDialog : public vtkCommand {

public:

  static vtkKWProgressDialog* New ();
  vtkTypeRevisionMacro( vtkKWProgressDialog, vtkCommand );

  void Execute ( vtkObject* iCaller, unsigned long iEvent, void* iCallData );

  void SetApplication ( vtkKWApplication* iApplication );
  void SetWindowTitle ( const char* isTitle );

protected:

  vtkKWProgressDialog ();
  virtual ~vtkKWProgressDialog ();

  vtkKWApplication* mApplication;
  vtkKWDialog* mDialog;
  vtkKWProgressGauge* mProgressGauge;

  char msWindowTitle[1024];
};

#endif
