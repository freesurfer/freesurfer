/**
 * @file  FSVolume.h
 * @brief Interactor to manage mouse and key input on render view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2009/04/29 22:53:50 $
 *    $Revision: 1.2.2.3 $
 *
 * Copyright (C) 2008-2009,
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
#include "CommonDataStruct.h"

extern "C"
{
#include "mri.h"
}

class wxWindow;
class wxCommandEvent;

class FSVolume
{
public:
  FSVolume( FSVolume* ref );
  virtual ~FSVolume();

  void Create( FSVolume* src, bool bCopyVoxelData );

  bool MRIRead( const char* filename, const char* reg_filename, wxWindow* wnd, wxCommandEvent& event );
  bool MRIWrite( const char* filename );
  bool MRIWrite();

  int OriginalIndexToRAS( float iIdxX, float iIdxY, float iIdxZ,
                          float& oRASX, float& oRASY, float& oRASZ );
  int RASToOriginalIndex ( float iRASX, float iRASY, float iRASZ,
                           float& oIdxX, float& oIdxY, float& oIdxZ );
  int RASToOriginalIndex ( float iRASX, float iRASY, float iRASZ,
                           int& oIdxX, int& oIdxY, int& oIdxZ );
  
  double GetVoxelValue( int i, int j, int k, int frame );
  
  void UpdateMRIFromImage( vtkImageData* rasImage, wxWindow* wnd, wxCommandEvent& event );

  vtkImageData* GetImageOutput();

  void GetBounds ( float oRASBounds[6] );

  void GetBounds( MRI* mri, float oRASBounds[6] );

  void GetPixelSize( double* pixelSize );

  int GetNumberOfFrames();

  double GetMinValue ()
  {
    return m_fMinValue;
  }
  double GetMaxValue ()
  {
    return m_fMaxValue;
  }

  double* GetRASToVoxelMatrix ()
  {
    return m_RASToVoxelMatrix;
  }

  double* GetVoxelToRASMatrix()
  {
    return m_VoxelToRASMatrix;
  }

  bool GetResampleToRAS()
  {
    return m_bResampleToRAS;
  }

  void SetResampleToRAS( bool bRemap )
  {
    m_bResampleToRAS = bRemap;
  }

  void TargetToRAS( const double* pos_in, double* pos_out );

  void TargetToRAS( const float* pos_in, float* pos_out );

  void TargetToRAS( double x_in, double y_in, double z_in,
                    double& x_out, double& y_out, double& z_out );

  void RASToTarget( const double* pos_in, double* pos_out );

  void RASToTarget( const float* pos_in, float* pos_out );

  void RASToTargetIndex( const double* pos_in, int* index_out );

  void RASToNativeRAS( const double* pos_in, double* pos_out ); // when there is registration/transformation involved,
                                                                // ras is not native ras!
  void NativeRASToRAS( const double* pos_in, double* pos_out );

  void TkRegToNativeRAS( const double* pos_in, double* pos_out );

  void NativeRASToTkReg( const double* pos_in, double* pos_out );

  void SetMRITarget( MRI* mri );

  void SetMRI( MRI*& mri_out, MRI* mri_in );

  MRI* GetMRITarget()
  {
    return m_MRITarget;
  }

  MRI* GetMRI()
  {
    return m_MRI;
  }

  MATRIX* GetRegMatrix()
  {
    return m_matReg;
  }

  bool HasOriginalTarget()
  {
    return m_MRIOrigTarget != NULL;
  }

  bool Rotate( std::vector<RotationElement>& rotations, wxWindow* wnd, wxCommandEvent& event );

protected:
  bool LoadRegistrationMatrix( const char* filename );
  void MapMRIToImage( wxWindow* wnd, wxCommandEvent& event );
  void CopyMRIDataToImage( MRI* mri, vtkImageData* image, wxWindow* wnd, wxCommandEvent& event );
  void CopyMatricesFromMRI();
  void CreateImage( MRI* mri, wxWindow* wnd, wxCommandEvent& event );
  void UpdateRASToRASMatrix();

  MATRIX* GetRotationMatrix( int nPlane, double angle, double* origin );

  vtkSmartPointer<vtkImageData> m_imageData;

  MRI*      m_MRI;
  MRI*      m_MRITarget;  // target space. header only
  MRI*      m_MRIRef;   // reference target space, can also serve as the registration target. header only
  MRI*      m_MRIOrigTarget; // orignal target space, header only
  MATRIX*   m_matReg;

  FSVolume* m_volumeRef;

  double    m_RASToVoxelMatrix[16];
  double    m_VoxelToRASMatrix[16];
  double    m_VoxelToVoxelMatrix[16]; // native to target
  double    m_RASToRASMatrix[16];  // native to target
  double    m_RASToTkRegMatrix[16];

  float     m_fMinValue;
  float     m_fMaxValue;

  bool      m_bResampleToRAS;
  double    m_MRIToImageMatrix[16];

  // RAS bounds.
  bool      m_bBoundsCacheDirty;
  float     m_RASBounds[6];
};

#endif


