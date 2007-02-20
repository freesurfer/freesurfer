/**
 * @file  vtkKWProgressDialog.h
 * @brief A KWWidgets progress dialog
 *
 * A simple KWWidgets dialog box with a progress bar tht listens to
 * update commands.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/02/20 22:18:59 $
 *    $Revision: 1.2 $
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
