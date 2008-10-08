/**
 * @file  FSLabel.h
 * @brief Interactor to manage mouse and key input on render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/10/08 19:14:35 $
 *    $Revision: 1.3 $
 *
 * Copyright (C) 2002-2009,
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
 
#ifndef FSLabel_h
#define FSLabel_h

#include "vtkImageData.h"
#include "vtkMatrix4x4.h"

extern "C" {
#include "label.h"
}

class FSVolume;
class wxWindow;
class wxCommandEvent;

class FSLabel 
{
public:
	FSLabel();
	virtual ~FSLabel();
	
	bool LabelRead( const char* filename );	
	bool LabelWrite( const char* filename );
		
	void UpdateLabelFromImage( vtkImageData* rasImage_in, FSVolume* ref_vol, wxWindow* wnd, wxCommandEvent& event );
	void UpdateRASImage( vtkImageData* rasImage_out, FSVolume* ref_vol );
	
protected:	
	LABEL*			m_label;

};

#endif 


