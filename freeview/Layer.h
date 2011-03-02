/**
 * @file  Layer.h
 * @brief Base Layer class. A layer is an independent data object with 2D and 3D graphical representations.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:02 $
 *    $Revision: 1.18 $
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

#ifndef Layer_h
#define Layer_h

#include "Listener.h"
#include "Broadcaster.h"
#include "CommonDataStruct.h"
#include <string>
#include <vector>

class vtkRenderer;
class vtkProp;
class wxWindow;
class wxCommandEvent;
class LayerProperties;

class Layer : public Listener, public Broadcaster
{
public:
  Layer();
  virtual ~Layer();

  const char* GetName()
  {
    return m_strName.c_str();
  }

  void SetName( const char* name );

  virtual void Append2DProps( vtkRenderer* renderer, int nPlane ) = 0;
  virtual void Append3DProps( vtkRenderer* renderer, bool* bPlaneVisibility = NULL ) = 0;

  virtual bool HasProp( vtkProp* prop ) = 0;

  virtual void SetVisible( bool bVisible = true ) = 0;
  virtual bool IsVisible() = 0;

  bool Rotate( std::vector<RotationElement>& rotations, wxWindow* wnd, wxCommandEvent& event );
  
  bool Translate( double x, double y, double z );  
  bool Translate( double* dPos );
  
  void Scale( double* scale, int nSampleMethod = 1 /* SAMPLE_TRILINEAR */ );
  
  void Restore();
  
  void ResetTranslate()
  {
    m_dTranslate[0] = 0;
    m_dTranslate[1] = 0;
    m_dTranslate[2] = 0;
  }
  
  void GetTranslate( double* pos )
  {
    pos[0] = m_dTranslate[0];
    pos[1] = m_dTranslate[1];
    pos[2] = m_dTranslate[2];
  }

  void ResetScale()
  {
    m_dScale[0] = 1;
    m_dScale[1] = 1;
    m_dScale[2] = 1;
  }
  
  void GetScale( double* scale )
  {
    scale[0] = m_dScale[0];
    scale[1] = m_dScale[1];
    scale[2] = m_dScale[2];
  }
  
  double* GetWorldOrigin();
  void GetWorldOrigin( double* origin );
  void SetWorldOrigin( double* origin );

  double* GetWorldVoxelSize();
  void GetWorldVoxelSize( double* voxelsize );
  void SetWorldVoxelSize( double* voxelsize );

  double* GetWorldSize();
  void GetWorldSize( double* size );
  void SetWorldSize( double* size );

  double* GetSlicePosition();
  void GetSlicePosition( double* slicePos );
  void SetSlicePosition( double* slicePos );
  void SetSlicePosition( int nPlane, double slicePos );

  void RASToVoxel( const double* pos, int* n );
  void VoxelToRAS( const int* n, double* pos );

  virtual void OnSlicePositionChanged( int nPlane ) = 0;

  bool IsTypeOf( std::string tname );

  std::string GetErrorString()
  {
    return m_strError;
  }

  void SetErrorString( std::string msg )
  {
    m_strError = msg;
  }

  bool IsLocked()
  {
    return m_bLocked;
  }

  void Lock( bool bLock );

  std::string GetEndType();
  
  inline LayerProperties* GetProperties()
  {
    return mProperties;
  }
  
  virtual void GetBounds( double* bounds );
  
  virtual void GetDisplayBounds( double* bounds );

  void Show()
  {
    SetVisible( true );
  }
  
  void Hide() 
  {
    SetVisible( false );
  }
  
protected:
  virtual void DoListenToMessage( std::string const iMessage, void* iData, void* sender );
  virtual bool DoRotate( std::vector<RotationElement>& rotations, wxWindow* wnd, wxCommandEvent& event ) 
  { 
    return true; 
  }
  virtual void DoRestore() {}
  
  virtual void DoTranslate( double* offset ) {}
  virtual void DoScale( double* scale, int nSampleMethod ) {}
  
  std::string  m_strName;
  double    m_dSlicePosition[3];
  double    m_dWorldOrigin[3];
  double    m_dWorldVoxelSize[3];
  double    m_dWorldSize[3];

  // translate and scale are for volume transformation
  double    m_dTranslate[3];
  double    m_dScale[3];
  
  bool   m_bLocked;

  LayerProperties*  mProperties;  
      
  std::string  m_strError;
  std::vector<std::string> m_strTypeNames;
};

#endif


