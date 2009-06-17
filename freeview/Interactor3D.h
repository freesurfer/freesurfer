/**
 * @file  Interactor3D.h
 * @brief Base Interactor class for 3D render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/06/17 20:41:17 $
 *    $Revision: 1.5 $
 *
 * Copyright (C) 2008-2009,
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

#ifndef Interactor3D_h
#define Interactor3D_h

#include "Interactor.h"

class Interactor3D : public Interactor
{
public:
  Interactor3D();
  virtual ~Interactor3D();

  // return true if to have parent Interactor3D continue processing the event
  // return false to stop event from further processing
  virtual bool ProcessMouseDownEvent( wxMouseEvent& event, RenderView* view );
  virtual bool ProcessMouseUpEvent( wxMouseEvent& event, RenderView* view );
  virtual bool ProcessMouseMoveEvent( wxMouseEvent& event, RenderView* view );
  virtual bool ProcessKeyDownEvent( wxKeyEvent& event, RenderView* view );

protected:
  int  m_nMousePosX;
  int  m_nMousePosY;

  bool m_bWindowLevel;
};

#endif


