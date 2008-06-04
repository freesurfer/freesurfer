/**
 * @file  RenderView3D.h
 * @brief View class for 2D image rendering.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/06/04 20:43:25 $
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
 
#ifndef RenderView3D_h
#define RenderView3D_h

#include "RenderView.h"

class VTK_RENDERING_EXPORT RenderView3D : public RenderView
{
  DECLARE_DYNAMIC_CLASS(RenderView3D)

public:
	RenderView3D();
    RenderView3D(wxWindow *parent, int id);
    virtual ~RenderView3D();	
    
    static RenderView3D * New();
    void PrintSelf(ostream& os, vtkIndent indent);
	
	virtual void RefreshAllActors();
	
	void UpdateViewByWorldCoordinate();
	
private:

    // any class wishing to process wxWindows events must use this macro
    DECLARE_EVENT_TABLE()

};

#endif 


