/**
 * @file  FSSurface.h
 * @brief Interactor to manage mouse and key input on render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/05/20 16:28:31 $
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
 
#ifndef FSSurface_h
#define FSSurface_h

#include "vtkSmartPointer.h"
#include "vtkImageData.h"
#include "vtkMatrix4x4.h"

extern "C" {
#include "mrisurf.h"
}

class wxWindow;
class wxCommandEvent;
class vtkTransform;

class FSSurface 
{
public:
	FSSurface();
	virtual ~FSSurface();
	
	bool MRISRead( const char* filename, wxWindow* wnd, wxCommandEvent& event );	

	void GetBounds ( float oRASBounds[6] );
	
	void ConvertSurfaceToRAS ( float iSurfX, float iSurfY, float iSurfZ,
							   float& oRASX, float& oRASY, float& oRASZ ) const;
	void ConvertSurfaceToRAS ( double iSurfX, double iSurfY, double iSurfZ,
							   double& oRASX, double& oRASY, double& oRASZ ) const;
	void ConvertRASToSurface ( float iRASX, float iRASY, float iRASZ,
							   float& oSurfX, float& oSurfY, float& oSurfZ) const;
	void ConvertRASToSurface ( double iRASX, double iRASY, double iRASZ,
							   double& oSurfX, double& oSurfY, double& oSurfZ) const;
	void ConvertSurfaceToRAS ( float const iSurf[3], float oRAS[3] ) const;
	void ConvertSurfaceToRAS ( double const iSurf[3], double oRAS[3] ) const;
	void ConvertRASToSurface ( float const iRAS[3], float oSurf[3] ) const;
	void ConvertRASToSurface ( double const iRAS[3], double oSurf[3] ) const;

	
	int GetNumberOfVertices () const;
	
	MRIS* GetMRIS() { return m_MRIS; }
	
protected:	
	MRIS*			m_MRIS;
	
	double			m_SurfaceToRASMatrix[16];
	vtkSmartPointer<vtkTransform> m_SurfaceToRASTransform;

	
	// RAS bounds.
	bool  	m_bBoundsCacheDirty;
	float 	m_RASBounds[6];
	float	m_RASCenter[3];
	

	// Hash table so we can look up vertices. Uses v->x,y,z.
	MRIS_HASH_TABLE* m_HashTable;
};

#endif 


