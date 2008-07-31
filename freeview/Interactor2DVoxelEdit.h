/**
 * @file  Interactor2DVoxelEdit.h
 * @brief Interactor2DVoxelEdit to manage mouse and key input on render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/07/31 22:51:45 $
 *    $Revision: 1.5 $
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
 
#ifndef Interactor2DVoxelEdit_h
#define Interactor2DVoxelEdit_h

#include "Interactor2D.h"
#include <vector>

class wxWindow;

class Interactor2DVoxelEdit : public Interactor2D
{
	public:
		Interactor2DVoxelEdit();
		virtual ~Interactor2DVoxelEdit();
	
		enum EditMode { EM_Freehand = 0, EM_Fill, EM_Polyline };
	
	// return true if to have parent Interactor2DROIEdit continue processing the event
	// return false to stop event from further processing
		virtual bool ProcessMouseDownEvent( wxMouseEvent& event, RenderView* view );
		virtual bool ProcessMouseUpEvent( wxMouseEvent& event, RenderView* view );
		virtual bool ProcessMouseMoveEvent( wxMouseEvent& event, RenderView* view );
		virtual bool ProcessKeyDownEvent( wxKeyEvent& event, RenderView* view );
		virtual bool ProcessKeyUpEvent( wxKeyEvent& event, RenderView* view );
		
	protected:	
		void UpdateCursor( wxEvent& event, wxWindow* wnd );
		
		bool	m_bEditing;
	
		std::vector<double>		m_dPolylinePoints;
};

#endif 


