/**
 * @file  Layer.h
 * @brief Base Layer class. A layer is an independent data object with 2D and 3D graphical representations.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2012/05/03 19:50:01 $
 *    $Revision: 1.24 $
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

  bool Transform (double* mat, int sample_method);
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

  void ResetRotate()
  {
    m_dRotate[0] = 0;
    m_dRotate[1] = 0;
    m_dRotate[2] = 0;
  }

  void GetRotate( double* rotate)
  {
    rotate[0] = m_dRotate[0];
    rotate[1] = m_dRotate[1];
    rotate[2] = m_dRotate[2];
  }

  void SetRotate(double* rotate, bool bAroundCenter = true);
  void SetTranslate(double* offset);
  void SetScale(double* scale);
  void UpdateTransform(int sample_method = 0);

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
  virtual void DoTransform( double* mat, int sample_method ) {}
  virtual void DoRotate( double* rotate, double* pos) {}

  // new transform scheme
  virtual void DoTransform(int sample_method) {}

  QString   m_strName;
  double    m_dSlicePosition[3];
  double    m_dWorldOrigin[3];
  double    m_dWorldVoxelSize[3];
  double    m_dWorldSize[3];

  // translate and scale are for volume transformation
  double    m_dTranslate[3];
  double    m_dScale[3];
  double    m_dRotate[3];
  bool      m_bRotateAroundCenter;

  bool      m_bLocked;

  LayerProperty*  mProperty;

  QString   m_sFilename;
  QStringList m_strTypeNames;

  int       m_nID;
  static int m_nLastID;
};

#endif


