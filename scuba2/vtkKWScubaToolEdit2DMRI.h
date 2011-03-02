/**
 * @file  vtkKWScubaToolEdit2DMRI.h
 * @brief A tool for editing MRI volumes.
 *
 * Tool implementation that works on MRI volumes (2D and 3D).
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.2 $
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


#ifndef vtkKWScubaToolEdit2DMRI_h
#define vtkKWScubaToolEdit2DMRI_h

#include "vtkKWScubaTool.h"

class vtkKWScaleWithEntry;

class vtkKWScubaToolEdit2DMRI : public vtkKWScubaTool {

public:

  static vtkKWScubaToolEdit2DMRI* New ();
  vtkTypeRevisionMacro( vtkKWScubaToolEdit2DMRI, vtkKWScubaTool );

  // Description:
  // We need picks during our mouse drags so we can edit.
  virtual bool SuspendPickEvents () { return false; }

  // Description:
  // The new value to use in editing.
  int GetNewValue ();
  void SetNewValue ( int iNewValue );

protected:

  vtkKWScubaToolEdit2DMRI ();
  ~vtkKWScubaToolEdit2DMRI ();

  // Description:
  // We add and remove a slider for our new value.
  void AddControls ( vtkKWWidget* iPanel );
  void RemoveControls ();

  // Description:
  virtual void DoMouseDrag  ( vtkKWScubaWindow* iWindow, vtkKWScubaView* iView,
                              vtkKWScubaLayer* iLayer, float iRAS[3] );

  int mNewValue;
  vtkKWScaleWithEntry* mScaleNewValue;

};

#endif
