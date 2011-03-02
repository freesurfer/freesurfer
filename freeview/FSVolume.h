/**
 * @file  FSVolume.h
 * @brief Base volume class that takes care of I/O and data conversion.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:02 $
 *    $Revision: 1.31 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
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
#include "colortab.h"
}

class wxWindow;
class wxCommandEvent;
class vtkTransform;

class FSVolume
{
public:
  FSVolume( FSVolume* ref );
  virtual ~FSVolume();

  bool Create( FSVolume* src, bool bCopyVoxelData, int data_type );

  bool MRIRead( const char* filename, const char* reg_filename, wxWindow* wnd, wxCommandEvent& event );
  bool MRIWrite( const char* filename, int nSampleMethod = SAMPLE_NEAREST );
  bool MRIWrite();
  bool SaveRegistration( const char* filename );
  bool Restore( const char* filename, const char* reg_filename, wxWindow* wnd, wxCommandEvent& event );

  int OriginalIndexToRAS( float iIdxX, float iIdxY, float iIdxZ,
                          float& oRASX, float& oRASY, float& oRASZ );
  int RASToOriginalIndex ( float iRASX, float iRASY, float iRASZ,
                           float& oIdxX, float& oIdxY, float& oIdxZ );
  int RASToOriginalIndex ( float iRASX, float iRASY, float iRASZ,
                           int& oIdxX, int& oIdxY, int& oIdxZ );
  
  double GetVoxelValue( int i, int j, int k, int frame );
  
  bool UpdateMRIFromImage( vtkImageData* rasImage, wxWindow* wnd, wxCommandEvent& event, 
                           bool resampleToOriginal = true );

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

  bool Rotate( std::vector<RotationElement>& rotations, wxWindow* wnd, wxCommandEvent& event, int nSampleMethod = -1 );
  
  int GetInterpolationMethod()
  {
    return m_nInterpolationMethod;
  }
  
  void SetInterpolationMethod( int nMethod );

  int GetDataType();
  
  MATRIX* GetTargetToRASMatrix();
  
  COLOR_TABLE*  GetEmbeddedColorTable()
  {
    return m_ctabEmbedded;
  }
  
  const char* GetOrientationString()
  {
    return m_strOrientation;
  }
  
  void SetCroppingBounds( double* bound );
  
  vtkTransform* GetTransform();
  
  void SetConform( bool bConform );
  
protected:
  bool LoadMRI( const char* filename, const char* reg_filename, wxWindow* wnd, wxCommandEvent& event );
  bool LoadRegistrationMatrix( const char* filename );
  bool MapMRIToImage( wxWindow* wnd, wxCommandEvent& event );
  void CopyMRIDataToImage( MRI* mri, vtkImageData* image, wxWindow* wnd, wxCommandEvent& event );
  void CopyMatricesFromMRI();
  bool CreateImage( MRI* mri, wxWindow* wnd, wxCommandEvent& event );
  bool ResizeRotatedImage( MRI* mri, MRI* refTarget, vtkImageData* refImageData, double* rasPoint,
                    wxWindow* wnd, wxCommandEvent& event );
  void UpdateRASToRASMatrix();
  MRI* CreateTargetMRI( MRI* src, MRI* refTarget, bool AllocatePixel = true, bool bConform = false );

  MATRIX* GetRotationMatrix( int nPlane, double angle, double* origin );

  vtkSmartPointer<vtkImageData> m_imageData;
  vtkSmartPointer<vtkTransform> m_transform;

  MRI*      m_MRI;
  MRI*      m_MRITarget;      // target space. header only
  MRI*      m_MRIRef;         // reference target space, can also serve as the registration target. header only
  MRI*      m_MRIOrigTarget;  // orignal target space, header only
  MRI*      m_MRITemp;        // temp mri for saving
  MATRIX*   m_matReg;   
  COLOR_TABLE*  m_ctabEmbedded;
  
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
  
  int       m_nInterpolationMethod;
  bool      m_bConform;
  char      m_strOrientation[4];
  
  double    m_dBounds[6];
  bool      m_bCrop;
};

#endif


