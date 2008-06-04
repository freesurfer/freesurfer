/**
 * @file  Annotation2D.h
 * @brief Annotation class for 2D view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/06/04 20:43:23 $
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
 
#ifndef Annotation2D_h
#define Annotation2D_h

#include "RenderView.h"
#include "vtkSmartPointer.h"

class vtkTextActor;
class vtkRenderer;

class Annotation2D
{
public:
	Annotation2D();
    virtual ~Annotation2D();	

	void Update( vtkRenderer* renderer, int nPlane );   
	
	void AppendAnnotations( vtkRenderer* renderer );
	
private:
	vtkSmartPointer<vtkTextActor>	m_actorCoordinates[5];
};

#endif 


