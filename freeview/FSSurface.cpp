/**
 * @file  FSSurface.h
 * @brief Interactor to manage mouse and key input on render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/06/04 20:43:24 $
 *    $Revision: 1.1.2.1 $
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

#include <wx/wx.h>
#include "FSSurface.h"
#include <stdexcept>
#include "vtkShortArray.h"
#include "vtkLongArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkIntArray.h"
#include "vtkSmartPointer.h"
#include "vtkImageReslice.h"
#include "vtkMatrix4x4.h"
#include "vtkTransform.h"
#include "vtkImageChangeInformation.h"

using namespace std;

FSSurface::FSSurface() :
	m_MRIS( NULL ),
	m_bBoundsCacheDirty( true ),
	m_HashTable( NULL )
{
}
	
FSSurface::~FSSurface()
{
	if ( m_MRIS )
		::MRISfree( &m_MRIS );
	
	if ( m_HashTable )
		MHTfree( &m_HashTable );
}
		
bool FSSurface::MRISRead( const char* filename, wxWindow* wnd, wxCommandEvent& event  )
{
	if ( m_MRIS )
		::MRISfree( &m_MRIS );
	
	event.SetInt( event.GetInt() + 1 );
	wxPostEvent( wnd, event );
	
	char* fn = strdup( filename );
	m_MRIS = ::MRISread( fn );
	free( fn );	
	
	if ( m_MRIS == NULL ) {
		cerr << "MRISread failed";
		return false;
	}	
	cout << "MRISread finished" << endl;
	
  // Get some info from the MRIS. This can either come from the volume
  // geometry data embedded in the surface; this is done for newer
  // surfaces. Or it can come from the source information in the
  // transform. We use it to get the RAS center offset for the
  // surface->RAS transform.
	if ( m_MRIS->vg.valid ) 
	{
		m_SurfaceToRASMatrix[0] = 1;
		m_SurfaceToRASMatrix[1] = 0;
		m_SurfaceToRASMatrix[2] = 0;
		m_SurfaceToRASMatrix[3] = m_MRIS->vg.c_r;
		m_SurfaceToRASMatrix[4] = 0;
		m_SurfaceToRASMatrix[5] = 1;
		m_SurfaceToRASMatrix[6] = 0;
		m_SurfaceToRASMatrix[7] = m_MRIS->vg.c_a;
		m_SurfaceToRASMatrix[8] = 0;
		m_SurfaceToRASMatrix[9] = 0;
		m_SurfaceToRASMatrix[10] = 1;
		m_SurfaceToRASMatrix[11] = m_MRIS->vg.c_s;
		m_SurfaceToRASMatrix[12] = 0;
		m_SurfaceToRASMatrix[13] = 0;
		m_SurfaceToRASMatrix[14] = 0;
		m_SurfaceToRASMatrix[15] = 1;
    
	} 
	else if ( m_MRIS->lta ) 
	{
		m_SurfaceToRASMatrix[0] = 1;
		m_SurfaceToRASMatrix[1] = 0;
		m_SurfaceToRASMatrix[2] = 0;
		m_SurfaceToRASMatrix[3] = -m_MRIS->lta->xforms[0].src.c_r;
		m_SurfaceToRASMatrix[4] = 0;
		m_SurfaceToRASMatrix[5] = 1;
		m_SurfaceToRASMatrix[6] = 0;
		m_SurfaceToRASMatrix[7] = -m_MRIS->lta->xforms[0].src.c_a;
		m_SurfaceToRASMatrix[8] = 0;
		m_SurfaceToRASMatrix[9] = 0;
		m_SurfaceToRASMatrix[10] = 1;
		m_SurfaceToRASMatrix[11] = -m_MRIS->lta->xforms[0].src.c_s;
		m_SurfaceToRASMatrix[12] = 0;
		m_SurfaceToRASMatrix[13] = 0;
		m_SurfaceToRASMatrix[14] = 0;
		m_SurfaceToRASMatrix[15] = 1;
	}

  // Make our transform object and set the matrix.
	m_SurfaceToRASTransform = vtkSmartPointer<vtkTransform>::New();
	m_SurfaceToRASTransform->SetMatrix( m_SurfaceToRASMatrix );
  
  // Make the hash table. This makes it with v->x,y,z.
	if ( m_HashTable )
		MHTfree( &m_HashTable );
	m_HashTable = MHTfillVertexTableRes( m_MRIS, NULL, CURRENT_VERTICES, 2.0 );
		
	return true;
}


void FSSurface::ConvertSurfaceToRAS ( float iX, float iY, float iZ,
		float& oX, float& oY, float& oZ ) const 
{
	float surface[3];
	float ras[3];
	
	surface[0] = iX;
	surface[1] = iY;
	surface[2] = iZ;

	this->ConvertSurfaceToRAS( surface, ras );

	oX = ras[0];
	oY = ras[1];
	oZ = ras[2];
}

void FSSurface::ConvertSurfaceToRAS ( double iX, double iY, double iZ,
				double& oX, double& oY, double& oZ ) const 
{
	double surface[3];
	double ras[3];

	surface[0] = iX;
	surface[1] = iY;
	surface[2] = iZ;

	this->ConvertSurfaceToRAS( surface, ras );

	oX = ras[0];
	oY = ras[1];
	oZ = ras[2];
}

void FSSurface::ConvertRASToSurface ( float iX, float iY, float iZ,
						float& oX, float& oY, float& oZ ) const 
{
	float ras[3];
	float surface[3];
  
	ras[0] = iX;
	ras[1] = iY;
	ras[2] = iZ;
  
	this->ConvertRASToSurface( ras, surface );

	oX = surface[0];
	oY = surface[1];
	oZ = surface[2];
}

void FSSurface::ConvertRASToSurface ( double iX, double iY, double iZ,
								double& oX, double& oY, double& oZ ) const 
{
	double ras[3];
	double surface[3];

	ras[0] = iX;
	ras[1] = iY;
	ras[2] = iZ;
  
	this->ConvertRASToSurface( ras, surface );

	oX = surface[0];
	oY = surface[1];
	oZ = surface[2];
}

void FSSurface::ConvertSurfaceToRAS ( float const iSurf[3], float oRAS[3] ) const 
{
	m_SurfaceToRASTransform->TransformPoint( iSurf, oRAS );
}

void FSSurface::ConvertSurfaceToRAS ( double const iSurf[3], double oRAS[3] ) const 
{
	m_SurfaceToRASTransform->TransformPoint( iSurf, oRAS );
}

void FSSurface::ConvertRASToSurface ( float const iRAS[3], float oSurf[3] ) const 
{
	m_SurfaceToRASTransform->GetInverse()->TransformPoint( iRAS, oSurf );
}

void FSSurface::ConvertRASToSurface ( double const iRAS[3], double oSurf[3] ) const 
{
	m_SurfaceToRASTransform->GetInverse()->TransformPoint( iRAS, oSurf );
}

void FSSurface::GetBounds ( float oRASBounds[6] ) 
{
	if ( NULL == m_MRIS ) 
	{
		oRASBounds[0] = oRASBounds[1] = oRASBounds[2] =
				oRASBounds[3] = oRASBounds[4] = oRASBounds[5] = 0;
		return;
	}

	if ( m_bBoundsCacheDirty ) 
	{
		m_RASBounds[0] = m_RASBounds[2] = m_RASBounds[4] = 999999;
		m_RASBounds[1] = m_RASBounds[3] = m_RASBounds[5] = -999999;

    // Find the bounds.
		for ( int vno = 0; vno < m_MRIS->nvertices; vno++ ) {

      // Translate to actual RAS coords.
			float rasX, rasY, rasZ;
			this->ConvertSurfaceToRAS( m_MRIS->vertices[vno].x,
									   m_MRIS->vertices[vno].y,
									   m_MRIS->vertices[vno].z,
									   rasX, rasY, rasZ );

			if ( rasX < m_RASBounds[0] ) m_RASBounds[0] = rasX;
			if ( rasX > m_RASBounds[1] ) m_RASBounds[1] = rasX;
			if ( rasY < m_RASBounds[2] ) m_RASBounds[2] = rasY;
			if ( rasY > m_RASBounds[3] ) m_RASBounds[3] = rasY;
			if ( rasZ < m_RASBounds[4] ) m_RASBounds[4] = rasZ;
			if ( rasZ > m_RASBounds[5] ) m_RASBounds[5] = rasZ;

		}

		m_bBoundsCacheDirty = false;
	}

	oRASBounds[0] = m_RASBounds[0];
	oRASBounds[1] = m_RASBounds[1];
	oRASBounds[2] = m_RASBounds[2];
	oRASBounds[3] = m_RASBounds[3];
	oRASBounds[4] = m_RASBounds[4];
	oRASBounds[5] = m_RASBounds[5];
}


int FSSurface::GetNumberOfVertices () const 
{
	if ( m_MRIS )
		return m_MRIS->nvertices;
	else
		return 0;
}
