/**
 * @file  Interactor2D.h
 * @brief Interactor2D to manage mouse and key input on render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/06/04 20:43:24 $
 *    $Revision: 1.2.2.1 $
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
 
#ifndef Interactor2D_h
#define Interactor2D_h

#include "Interactor.h"

class Interactor2D : public Interactor
{
public:
	Interactor2D();
	virtual ~Interactor2D();
		
	// return true if to have parent Interactor2D continue processing the event
	// return false to stop event from further processing
	virtual bool ProcessMouseDownEvent( wxMouseEvent& event, RenderView* view );
	virtual bool ProcessMouseUpEvent( wxMouseEvent& event, RenderView* view );
	virtual bool ProcessMouseMoveEvent( wxMouseEvent& event, RenderView* view );
	virtual bool ProcessKeyDownEvent( wxKeyEvent& event, RenderView* view );
		
	virtual void ProcessPostMouseWheelEvent( wxMouseEvent& event, RenderView* view );
	virtual void ProcessPostMouseMoveEvent( wxMouseEvent& event, RenderView* view );
	
protected:	
	int		m_nMousePosX;
	int		m_nMousePosY;
	
	bool	m_bWindowLevel;
	bool 	m_bChangeSlice;
};

#endif 


