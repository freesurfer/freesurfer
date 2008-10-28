/**
 * @file  FSSurface.h
 * @brief Interactor to manage mouse and key input on render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/10/28 17:26:46 $
 *    $Revision: 1.9 $
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

#include <wx/wx.h>
#include <wx/filename.h>
#include "FSSurface.h"
#include <stdexcept>
#include "vtkShortArray.h"
#include "vtkLongArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkCellArray.h"
#include "vtkIntArray.h"
#include "vtkSmartPointer.h"
#include "vtkImageReslice.h"
#include "vtkMatrix4x4.h"
#include "vtkTransform.h"
#include "vtkImageChangeInformation.h"
#include "vtkPolyData.h"
#include "FSVolume.h"

using namespace std;

FSSurface::FSSurface( FSVolume* ref ) :
	m_MRIS( NULL ),
	m_bBoundsCacheDirty( true ),
	m_bCurvatureLoaded( false ),
	m_nActiveSurface( SurfaceMain ), 
	m_volumeRef( ref )
{
	m_polydata = vtkPolyData::New();
	
	for ( int i = 0; i < NUM_OF_VSETS; i++ )
	{
		m_fVertexSets[i] = NULL;
		m_fNormalSets[i] = NULL;
		m_bSurfaceLoaded[i] = false;
		m_HashTable[i] = NULL;
	}
}
	
FSSurface::~FSSurface()
{
	if ( m_MRIS )
		::MRISfree( &m_MRIS );
	
	for ( int i = 0; i < NUM_OF_VSETS; i++ )
	{
		if ( m_fNormalSets[i] )
			delete[] m_fNormalSets[i];
		if ( m_fVertexSets[i] )
			delete[] m_fVertexSets[i];
		if ( m_HashTable[i] )
			MHTfree( &m_HashTable[i] );
	}
		
	m_polydata->Delete();
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
		cerr << "MRISread failed" << endl;
		return false;
	}	
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
	if ( m_HashTable[0] )
		MHTfree( &m_HashTable[0] );
	m_HashTable[0] = MHTfillVertexTableRes( m_MRIS, NULL, CURRENT_VERTICES, 2.0 );
		
	UpdatePolyData();
	
	// Save main vertices and normals
	SaveVertices( m_MRIS, SurfaceMain );
	SaveNormals	( m_MRIS, SurfaceMain );
	m_bSurfaceLoaded[SurfaceMain] = true;
			
	LoadSurface	( "inflated", 	SurfaceInflated );
	LoadSurface	( "white", 		SurfaceWhite );
	LoadSurface	( "pial", 		SurfacePial );
	LoadSurface	( "orig", 		SurfaceOriginal );
	LoadCurvature();
	
//	cout << "MRISread finished" << endl;
		
	return true;
}

bool FSSurface::LoadSurface( const char* filename, int nSet )
{
	if ( ::MRISreadVertexPositions( m_MRIS, (char*)filename ) != 0 )
	{
		cerr << "could not load surface from " << filename << endl;
		m_bSurfaceLoaded[nSet] = false;
		return false;
	}
	else
	{
		if ( m_HashTable[nSet] )
				MHTfree( &m_HashTable[nSet] );
		m_HashTable[nSet] = MHTfillVertexTableRes( m_MRIS, NULL, CURRENT_VERTICES, 2.0 );
		ComputeNormals();		
		SaveVertices( m_MRIS, nSet );
		SaveNormals	( m_MRIS, nSet );
		m_bSurfaceLoaded[nSet] = true;
		return true;
	}	
}

bool FSSurface::LoadCurvature( const char* filename )
{
	if ( ::MRISreadCurvatureFile( m_MRIS, (char*)(filename ? filename : "curv" ) ) != 0 )
	{
		cerr << "could not read curvature from " << (filename ? filename : "curv") << endl;
		m_bCurvatureLoaded = false;
		return false;
	}
	else
	{
		int cVertices = m_MRIS->nvertices;

		vtkSmartPointer<vtkFloatArray> curvs =
				vtkSmartPointer<vtkFloatArray>::New();
		curvs->Allocate( cVertices );
		curvs->SetNumberOfComponents( 1 );
		curvs->SetName( "Curvature" );

 		for ( int vno = 0; vno < cVertices; vno++ ) 
		{
			curvs->InsertNextValue( m_MRIS->vertices[vno].curv );
		}
		m_polydata->GetPointData()->SetScalars( curvs );
		m_bCurvatureLoaded = true;
		return true;
	}
}

void FSSurface::SaveVertices( MRIS* mris, int nSet )
{
	if ( !mris || nSet >= NUM_OF_VSETS )
		return;
	
	int nvertices = mris->nvertices;
	VERTEX *v;
	
	if ( m_fVertexSets[nSet] == NULL )
	{
		m_fVertexSets[nSet] = new VertexItem[nvertices];
		if ( !m_fVertexSets[nSet] )
		{
			cerr << "Can not allocate memory for vertex sets." << endl;
			return;
		}
	}
	for ( int vno = 0; vno < nvertices; vno++ )
	{
		v = &mris->vertices[vno];
		m_fVertexSets[nSet][vno].x = v->x;
		m_fVertexSets[nSet][vno].y = v->y;
		m_fVertexSets[nSet][vno].z = v->z;
	}
}
	
void FSSurface::RestoreVertices( MRIS* mris, int nSet )
{
	if ( !mris || nSet >= NUM_OF_VSETS || m_fVertexSets[nSet] == NULL)
		return;
	
	int nvertices = mris->nvertices;
	VERTEX *v;

	for ( int vno = 0; vno < nvertices; vno++ )
	{
		v = &mris->vertices[vno];
		v->x = m_fVertexSets[nSet][vno].x;
		v->y = m_fVertexSets[nSet][vno].y;
		v->z = m_fVertexSets[nSet][vno].z;
	}
}


void FSSurface::SaveNormals( MRIS* mris, int nSet )
{
	if ( !mris || nSet >= NUM_OF_VSETS )
		return;
	
	int nvertices = mris->nvertices;
	VERTEX *v;
	
	if ( m_fNormalSets[nSet] == NULL )
	{
		m_fNormalSets[nSet] = new VertexItem[nvertices];
		if ( !m_fNormalSets[nSet] )
		{
			cerr << "Can not allocate memory for normal sets." << endl;
			return;
		}
	}
	for ( int vno = 0; vno < nvertices; vno++ )
	{
		v = &mris->vertices[vno];
		m_fNormalSets[nSet][vno].x = v->nx;
		m_fNormalSets[nSet][vno].y = v->ny;
		m_fNormalSets[nSet][vno].z = v->nz;
	}
}
	
void FSSurface::RestoreNormals( MRIS* mris, int nSet )
{
	if ( !mris || nSet >= NUM_OF_VSETS || m_fNormalSets[nSet] == NULL)
		return;
	
	int nvertices = mris->nvertices;
	VERTEX *v;

	for ( int vno = 0; vno < nvertices; vno++ )
	{
		v = &mris->vertices[vno];
		v->nx = m_fNormalSets[nSet][vno].x;
		v->ny = m_fNormalSets[nSet][vno].y;
		v->nz = m_fNormalSets[nSet][vno].z;
	}
}


void FSSurface::UpdatePolyData()
{
  // Allocate all our arrays.
	int cVertices = m_MRIS->nvertices;
	int cFaces = m_MRIS->nfaces;

	vtkSmartPointer<vtkPoints> newPoints = 
			vtkSmartPointer<vtkPoints>::New();
	newPoints->Allocate( cVertices );

	vtkSmartPointer<vtkCellArray> newPolys = 
			vtkSmartPointer<vtkCellArray>::New();
	newPolys->Allocate( newPolys->EstimateSize(cFaces,VERTICES_PER_FACE) );

	vtkSmartPointer<vtkFloatArray> newNormals =
			vtkSmartPointer<vtkFloatArray>::New();
	newNormals->Allocate( cVertices );
	newNormals->SetNumberOfComponents( 3 );
	newNormals->SetName( "Normals" );

  // Go through the surface and copy the vertex and normal for each
  // vertex. We have to transform them from surface RAS into normal
  // RAS.
	float point[3], normal[3], surfaceRAS[3];
	for ( int vno = 0; vno < cVertices; vno++ ) {

		surfaceRAS[0] = m_MRIS->vertices[vno].x;
		surfaceRAS[1] = m_MRIS->vertices[vno].y;
		surfaceRAS[2] = m_MRIS->vertices[vno].z;
		this->ConvertSurfaceToRAS( surfaceRAS, point );
		if ( m_volumeRef )
			m_volumeRef->RASToTarget( point, point );
		newPoints->InsertNextPoint( point );

		normal[0] = m_MRIS->vertices[vno].nx;
		normal[1] = m_MRIS->vertices[vno].ny;
		normal[2] = m_MRIS->vertices[vno].nz;
		newNormals->InsertNextTuple( normal );
	}

  // Go through and add the face indices.
	vtkIdType face[VERTICES_PER_FACE];
	for ( int fno = 0; fno < cFaces; fno++ ) {

		face[0] = m_MRIS->faces[fno].v[0];
		face[1] = m_MRIS->faces[fno].v[1];
		face[2] = m_MRIS->faces[fno].v[2];
		newPolys->InsertNextCell( 3, face );
	}

	m_polydata->SetPoints( newPoints );
	m_polydata->GetPointData()->SetNormals( newNormals );
	newPolys->Squeeze(); // since we've estimated size; reclaim some space
	m_polydata->SetPolys( newPolys );	
//	m_polydata->Update();
}

void FSSurface::UpdateVerticesAndNormals()
{
  // Allocate all our arrays.
	int cVertices = m_MRIS->nvertices;

	vtkSmartPointer<vtkPoints> newPoints = 
			vtkSmartPointer<vtkPoints>::New();
	newPoints->Allocate( cVertices );

	vtkSmartPointer<vtkFloatArray> newNormals =
			vtkSmartPointer<vtkFloatArray>::New();
	newNormals->Allocate( cVertices );
	newNormals->SetNumberOfComponents( 3 );
	newNormals->SetName( "Normals" );

  // Go through the surface and copy the vertex and normal for each
  // vertex. We have to transform them from surface RAS into normal
  // RAS.
	float point[3], normal[3], surfaceRAS[3];
	for ( int vno = 0; vno < cVertices; vno++ ) {

		surfaceRAS[0] = m_MRIS->vertices[vno].x;
		surfaceRAS[1] = m_MRIS->vertices[vno].y;
		surfaceRAS[2] = m_MRIS->vertices[vno].z;
		this->ConvertSurfaceToRAS( surfaceRAS, point );
		if ( m_volumeRef )
			m_volumeRef->RASToTarget( point, point );
		newPoints->InsertNextPoint( point );

		normal[0] = m_MRIS->vertices[vno].nx;
		normal[1] = m_MRIS->vertices[vno].ny;
		normal[2] = m_MRIS->vertices[vno].nz;
		newNormals->InsertNextTuple( normal );
	}

	m_polydata->SetPoints( newPoints );
	m_polydata->GetPointData()->SetNormals( newNormals );
//	m_polydata->Update();
}

bool FSSurface::SetActiveSurface( int nIndex )
{
	if ( nIndex == m_nActiveSurface )
		return false;
	
	m_nActiveSurface = nIndex;
	
	RestoreVertices( m_MRIS, nIndex );
	
	if ( m_fNormalSets[nIndex] == NULL )
	{
		ComputeNormals();
		SaveNormals( m_MRIS, nIndex );
	}
	else
	{
		RestoreNormals( m_MRIS, nIndex );
	}
	
	UpdateVerticesAndNormals();
	
	return true;
}

void FSSurface::NormalFace(int fac, int n, float *norm)
{
	int n0,n1;
	FACE *f;
	MRIS* mris = m_MRIS;
	float v0[3],v1[3];

	n0 = (n==0)                   ? VERTICES_PER_FACE-1 : n-1;
	n1 = (n==VERTICES_PER_FACE-1) ? 0                   : n+1;
	f = &mris->faces[fac];
	v0[0] = mris->vertices[f->v[n]].x - mris->vertices[f->v[n0]].x;
	v0[1] = mris->vertices[f->v[n]].y - mris->vertices[f->v[n0]].y;
	v0[2] = mris->vertices[f->v[n]].z - mris->vertices[f->v[n0]].z;
	v1[0] = mris->vertices[f->v[n1]].x - mris->vertices[f->v[n]].x;
	v1[1] = mris->vertices[f->v[n1]].y - mris->vertices[f->v[n]].y;
	v1[2] = mris->vertices[f->v[n1]].z - mris->vertices[f->v[n]].z;
	Normalize(v0);
	Normalize(v1);
	norm[0] = -v1[1]*v0[2] + v0[1]*v1[2];
	norm[1] = v1[0]*v0[2] - v0[0]*v1[2];
	norm[2] = -v1[0]*v0[1] + v0[0]*v1[1];
  /*
	printf("[%5.2f,%5.2f,%5.2f] x [%5.2f,%5.2f,%5.2f] = [%5.2f,%5.2f,%5.2f]\n",
	v0[0],v0[1],v0[2],v1[0],v1[1],v1[2],norm[0],norm[1],norm[2]);
  */
}

float FSSurface::TriangleArea(int fac, int n)
{
	int n0,n1;
	FACE *f;
	MRIS* mris = m_MRIS;
	float v0[3],v1[3],d1,d2,d3;

	n0 = (n==0)                   ? VERTICES_PER_FACE-1 : n-1;
	n1 = (n==VERTICES_PER_FACE-1) ? 0                   : n+1;
	f = &mris->faces[fac];
	v0[0] = mris->vertices[f->v[n]].x - mris->vertices[f->v[n0]].x;
	v0[1] = mris->vertices[f->v[n]].y - mris->vertices[f->v[n0]].y;
	v0[2] = mris->vertices[f->v[n]].z - mris->vertices[f->v[n0]].z;
	v1[0] = mris->vertices[f->v[n1]].x - mris->vertices[f->v[n]].x;
	v1[1] = mris->vertices[f->v[n1]].y - mris->vertices[f->v[n]].y;
	v1[2] = mris->vertices[f->v[n1]].z - mris->vertices[f->v[n]].z;
	d1 = -v1[1]*v0[2] + v0[1]*v1[2];
	d2 = v1[0]*v0[2] - v0[0]*v1[2];
	d3 = -v1[0]*v0[1] + v0[0]*v1[1];
	return sqrt(d1*d1+d2*d2+d3*d3)/2;
}

void FSSurface::Normalize( float v[3] )
{
	float d;

	d = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	if (d>0)
	{
		v[0] /= d;
		v[1] /= d;
		v[2] /= d;
	}
}

void FSSurface::ComputeNormals()  
{
	if ( !m_MRIS )
		return;
	
	MRIS* mris = m_MRIS;
	int k,n;
	VERTEX *v;
	FACE *f;
	float norm[3],snorm[3];

	for (k=0;k<mris->nfaces;k++)
		if (mris->faces[k].ripflag)
	{
		f = &mris->faces[k];
		for (n=0;n<VERTICES_PER_FACE;n++)
			mris->vertices[f->v[n]].border = TRUE;
	}
	for (k=0;k<mris->nvertices;k++)
		if (!mris->vertices[k].ripflag)
	{
		v = &mris->vertices[k];
		snorm[0]=snorm[1]=snorm[2]=0;
		v->area = 0;
		for (n=0;n<v->num;n++)
			if (!mris->faces[v->f[n]].ripflag)
		{
			NormalFace(v->f[n],v->n[n],norm);
			snorm[0] += norm[0];
			snorm[1] += norm[1];
			snorm[2] += norm[2];
			v->area += TriangleArea(v->f[n],v->n[n]);
			/* Note: overest. area by 2! */
		}
		Normalize( snorm );

		if (v->origarea<0)
			v->origarea = v->area;

		v->nx = snorm[0];
		v->ny = snorm[1];
		v->nz = snorm[2];
	}
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


int FSSurface::FindVertexAtRAS ( float const iRAS[3], float* oDistance ) 
{
	float surf[3];
	this->ConvertRASToSurface( iRAS, surf );

	return this->FindVertexAtSurfaceRAS( surf, oDistance );
}

int	FSSurface::FindVertexAtRAS ( double const iRAS[3], double* oDistance ) 
{
	double surf[3];
	this->ConvertRASToSurface( iRAS, surf );

	return this->FindVertexAtSurfaceRAS( surf, oDistance );
}

int	FSSurface::FindVertexAtSurfaceRAS ( float const iSurfaceRAS[3],
						float* oDistance ) 
{
	VERTEX v;
	v.x = iSurfaceRAS[0];
	v.y = iSurfaceRAS[1];
	v.z = iSurfaceRAS[2];
	float distance;
	int nClosestVertex =
	MHTfindClosestVertexNo( m_HashTable[m_nActiveSurface], m_MRIS, &v, &distance );

	if ( -1 == nClosestVertex ) 
	{
	//	cerr << "No vertices found." << endl;
		return -1;
	}

	if ( NULL != oDistance ) 
	{
		*oDistance = distance;
	}

	return nClosestVertex;
}

int	FSSurface::FindVertexAtSurfaceRAS ( double const iSurfaceRAS[3],
								double* oDistance ) 
{
	VERTEX v;
	v.x = static_cast<float>(iSurfaceRAS[0]);
	v.y = static_cast<float>(iSurfaceRAS[1]);
	v.z = static_cast<float>(iSurfaceRAS[2]);
	float distance;
	int nClosestVertex = MHTfindClosestVertexNo( m_HashTable[m_nActiveSurface], m_MRIS, &v, &distance );
	if ( -1 == nClosestVertex ) 
	{
	//	cerr << "No vertices found.";
	}

	if ( NULL != oDistance ) 
	{
		*oDistance = static_cast<double>(distance);
	}

	return nClosestVertex;
}

bool FSSurface::GetRASAtVertex ( int inVertex, float ioRAS[3] ) 
{
	float surfaceRAS[3];
	if ( this->GetSurfaceRASAtVertex( inVertex, surfaceRAS ) )
	{
		this->ConvertSurfaceToRAS( surfaceRAS, ioRAS );
		return true;
	}
	else
		return false;
}

bool FSSurface::GetRASAtVertex ( int inVertex, double ioRAS[3] ) 
{
	double surfaceRAS[3];
	if ( this->GetSurfaceRASAtVertex( inVertex, surfaceRAS ) )
	{
		this->ConvertSurfaceToRAS( surfaceRAS, ioRAS );
		return true;
	}
	else
		return false;
}

bool FSSurface::GetSurfaceRASAtVertex ( int inVertex, float ioRAS[3] ) 
{
	if( m_MRIS == NULL )
//		throw runtime_error( "GetRASAtVertex: m_MRIS was NULL" );
		return false;

	if( inVertex < 0 || inVertex >= m_MRIS->nvertices )
//		throw runtime_error( "GetRASAtVertex: inVertex was invalid" );
		return false;

	if ( m_nActiveSurface >= 0 && m_fVertexSets[m_nActiveSurface] != NULL )
	{
		ioRAS[0] = m_fVertexSets[m_nActiveSurface][inVertex].x;
		ioRAS[1] = m_fVertexSets[m_nActiveSurface][inVertex].y;
		ioRAS[2] = m_fVertexSets[m_nActiveSurface][inVertex].z;
	}
	else
	{
		ioRAS[0] = m_MRIS->vertices[inVertex].x;
		ioRAS[1] = m_MRIS->vertices[inVertex].y;
		ioRAS[2] = m_MRIS->vertices[inVertex].z;
	}
	
	return true;
}

bool FSSurface::GetSurfaceRASAtVertex ( int inVertex, double ioRAS[3] ) 
{
	if( m_MRIS == NULL )
//		throw runtime_error( "GetRASAtVertex: m_MRIS was NULL" );
		return false;
	
	if( inVertex < 0 || inVertex >= m_MRIS->nvertices )
//		throw runtime_error( "GetRASAtVertex: inVertex was invalid" );
		return false;

	if ( m_nActiveSurface >= 0 && m_fVertexSets[m_nActiveSurface] != NULL )
	{
		ioRAS[0] = m_fVertexSets[m_nActiveSurface][inVertex].x;
		ioRAS[1] = m_fVertexSets[m_nActiveSurface][inVertex].y;
		ioRAS[2] = m_fVertexSets[m_nActiveSurface][inVertex].z;
	}
	else
	{	ioRAS[0] = static_cast<double>(m_MRIS->vertices[inVertex].x);
		ioRAS[1] = static_cast<double>(m_MRIS->vertices[inVertex].y);
		ioRAS[2] = static_cast<double>(m_MRIS->vertices[inVertex].z);
	}
	
	return true;
}
