/**
 * @file  Layer.h
 * @brief Base Layer class. A layer is an independent data object with 2D and 3D graphical representations.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/12 00:28:48 $
 *    $Revision: 1.19 $
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

#ifndef Layer_h
#define Layer_h

#include <QObject>
#include <QString>
#include <QStringList>
#include <vector>
#include "CommonDataStruct.h"

class vtkRenderer;
class vtkProp;
class LayerProperty;

class Layer : public QObject
{
    Q_OBJECT
public:
  Layer( QObject* parent = 0 );
  virtual ~Layer();

  QString GetName() const
  {
    return m_strName;
  }

  void SetName( const QString& name );

  virtual void Append2DProps( vtkRenderer* renderer, int nPlane ) = 0;
  virtual void Append3DProps( vtkRenderer* renderer, bool* bPlaneVisibility = NULL ) = 0;

  virtual bool HasProp( vtkProp* prop ) = 0;

  virtual void SetVisible( bool bVisible = true )
  {
      emit VisibilityChanged(bVisible);
  }

  virtual bool IsVisible() = 0;

  bool Rotate( std::vector<RotationElement>& rotations );
  
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


  virtual void OnSlicePositionChanged( int nPlane ) = 0;

  bool IsTypeOf( const QString& tname );

  bool IsLocked()
  {
    return m_bLocked;
  }

  void Lock( bool bLock );

  QString GetEndType() const;
  
  inline LayerProperty* GetProperty()
  {
    return mProperty;
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
  
  QString GetFileName()
  {
    return m_sFilename;
  }

  void SetFileName( const QString& fn )
  {
    m_sFilename = fn;
  }

  void SendActorUpdated()
  {
      emit ActorUpdated();
  }

Q_SIGNALS:
  void NameChanged( const QString& name );
  void Transformed();
  void Locked( bool bLock );
  void ActorUpdated();
  void ActorChanged();
  void VisibilityChanged(bool bVisible);

protected:
  virtual bool DoRotate( std::vector<RotationElement>& rotations )
  { 
    return true; 
  }

  virtual void DoRestore() {}
  
  virtual void DoTranslate( double* offset ) {}
  virtual void DoScale( double* scale, int nSampleMethod ) {}
  
  QString   m_strName;
  double    m_dSlicePosition[3];
  double    m_dWorldOrigin[3];
  double    m_dWorldVoxelSize[3];
  double    m_dWorldSize[3];

  // translate and scale are for volume transformation
  double    m_dTranslate[3];
  double    m_dScale[3];
  
  bool      m_bLocked;

  LayerProperty*  mProperty;
      
  QString   m_sFilename;
  QStringList m_strTypeNames;

  int       m_nID;
  static int m_nLastID;
};

#endif


