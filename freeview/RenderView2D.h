/**
 * @file  RenderView2D.h
 * @brief View class for 2D image rendering.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/03/27 20:39:00 $
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
 
#ifndef RenderView2D_h
#define RenderView2D_h

#include "RenderView.h"

class Annotation2D;
class Cursor2D;

class VTK_RENDERING_EXPORT RenderView2D : public RenderView
{
	friend class Interactor2D;
	
  DECLARE_DYNAMIC_CLASS(RenderView2D)

public:
	RenderView2D();
    RenderView2D(wxWindow *parent, int id);
    virtual ~RenderView2D();	
    
	enum InteractionMode { IM_Navigate = 0, IM_VoxelEdit, IM_ROIEdit }; 
	
    static RenderView2D * New();
    void PrintSelf(ostream& os, vtkIndent indent);

	virtual void RefreshAllActors();
	virtual void DoListenToMessage ( std::string const iMessage, void* const iData );
	
	virtual void TriggerContextMenu( const wxPoint& pos );
	
	void SetViewPlane( int nPlane );	
	int GetViewPlane();
	
	void UpdateViewByWorldCoordinate();
	
	void UpdateMouseRASPosition( int nX, int nY );
	void UpdateCursorRASPosition( int nX, int nY, bool bConnectPrevious = false );
	
	void SetInteractionMode( int nMode );
	
	void MousePositionToRAS( int nX, int nY, double* ras );
	
	Cursor2D* GetCursor()
		{ return m_cursor2D; }
	
	void UpdateCursor();
	
	void MoveLeft();
	void MoveRight();
	void MoveUp();
	void MoveDown();
    
protected:
	void Initialize2D();
	void UpdateAnnotation();
	
	int		m_nViewPlane;
	
	int		m_nCursorPosX;
	int		m_nCursorPosY;
	
	Annotation2D*		m_annotation2D;
	Cursor2D*			m_cursor2D;
	
    // any class wishing to process wxWindows events must use this macro
    DECLARE_EVENT_TABLE()

};

#endif 


