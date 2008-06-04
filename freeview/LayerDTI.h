/**
 * @file  LayerDTI.h
 * @brief Layer data object for MRI volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/06/04 20:43:24 $
 *    $Revision: 1.3.2.1 $
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
 
#ifndef LayerDTI_h
#define LayerDTI_h

#include "LayerMRI.h"
#include "vtkSmartPointer.h"
#include <string>


class FSVectorVolume;
class wxWindow;
class wxCommandEvent;
class LayerPropertiesDTI;

class LayerDTI : public LayerMRI
{
	public:
		LayerDTI();
		virtual ~LayerDTI();
					
		bool LoadDTIFromFile( wxWindow* wnd, wxCommandEvent& event );
		
		void SetVectorFileName( std::string filename )
			{ m_sVectorFileName = filename; }
		
		const char* GetVectorFileName()
			{ return m_sVectorFileName.c_str(); }
		
		LayerPropertiesDTI*	GetProperties();

	protected:
		void UpdateColorMap();
		void InitializeDTIColorMap( wxWindow* wnd, wxCommandEvent& event );
		
		FSVolume*		m_vectorSource;		
		std::string		m_sVectorFileName;
};

#endif 


