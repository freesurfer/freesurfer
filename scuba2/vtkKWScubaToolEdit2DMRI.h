/**
 * @file  vtkKWScubaToolEdit2DMRI.h
 * @brief A tool for editing MRI volumes.
 *
 * Tool implementation that works on MRI volumes (2D and 3D).
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/04/06 22:23:06 $
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
