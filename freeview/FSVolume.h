/**
 * @file  FSVolume.h
 * @brief Interactor to manage mouse and key input on render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/03/27 18:12:14 $
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
 
#ifndef FSVolume_h
#define FSVolume_h

#include "vtkSmartPointer.h"
#include "vtkImageData.h"
#include "vtkMatrix4x4.h"

extern "C" {
#include "mri.h"
}

class wxWindow;
class wxCommandEvent;

class FSVolume 
{
public:
	FSVolume();
	virtual ~FSVolume();
		
	void Create( FSVolume* src, bool bCopyVoxelData );
	
	bool MRIRead( const char* filename, wxWindow* wnd, wxCommandEvent& event );	
	bool MRIWrite( const char* filename );
	bool MRIWrite();
		
	int OriginalIndexToRAS( float iIdxX, float iIdxY, float iIdxZ,
					float& oRASX, float& oRASY, float& oRASZ );
	int RASToOriginalIndex ( float iRASX, float iRASY, float iRASZ,
					float& oIdxX, float& oIdxY, float& oIdxZ );
	int RASToOriginalIndex ( float iRASX, float iRASY, float iRASZ,
					int& oIdxX, int& oIdxY, int& oIdxZ );
	
	void UpdateMRIFromImage( vtkImageData* rasImage, wxWindow* wnd, wxCommandEvent& event );
	
	vtkImageData* GetImageOutput();

	void GetBounds ( float oRASBounds[6] );
	
	void GetPixelSize( double* pixelSize );
	
	int GetNumberOfFrames();
	
	double GetMinValue ()
		{ return m_fMinValue; }
	double GetMaxValue ()
		{ return m_fMaxValue; }
	
	double* GetRASToVoxelMatrix () 
		{ return m_RASToVoxelMatrix; }
	
	double* GetVoxelToRASMatrix()
		{ return m_VoxelToRASMatrix; }
	
	bool GetResampleToRAS()
		{ return m_bResampleToRAS; }
	
	void SetResampleToRAS( bool bRemap )
		{ m_bResampleToRAS = bRemap; }
	
	void RemapPositionToRealRAS( const double* pos_in, double* pos_out );
	void RemapPositionToRealRAS( double x_in, double y_in, double z_in, 
								 double& x_out, double& y_out, double& z_out );
	
	void RASToOutputIndex( const double* pos_in, int* index_out );
	
protected:	
	void MapMRIToImage( wxWindow* wnd, wxCommandEvent& event );
	void CopyMatricesFromMRI();	
	void SetOriginalOrigin( double* origin );
	
	vtkSmartPointer<vtkImageData>	m_imageData;
			
	MRI*			m_MRI;
	
	double			m_RASToVoxelMatrix[16];
	double			m_VoxelToRASMatrix[16];
	
	float			m_fMinValue;
	float			m_fMaxValue;
	
	double			m_fOriginalOrigin[3];
	
	bool			m_bResampleToRAS;	
	double			m_MRIToImageMatrix[16];		
	
	// RAS bounds.
	bool  	m_bBoundsCacheDirty;
	float 	m_RASBounds[6];
};

#endif 


