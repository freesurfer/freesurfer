/**
 * @brief Base Layer class. A layer is an independent data object with 2D and 3D graphical representations.
 *
 */
/*
 * Original Author: Ruopeng Wang
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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
#include "vtkSmartPointer.h"

class vtkRenderer;
class vtkProp;
class LayerProperty;
class vtkTransform;

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

  void GetFlip( bool* flip)
  {
    flip[0] = m_bFlip[0];
    flip[1] = m_bFlip[1];
    flip[2] = m_bFlip[2];
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

  void GetRotationCenter( double* c_pos)
  {
    c_pos[0] = m_dRotationCenter[0];
    c_pos[1] = m_dRotationCenter[1];
    c_pos[2] = m_dRotationCenter[2];
  }

  void SetRotate(double* rotate, bool bAroundCenter = true);
  void SetTranslate(double* offset);
  void SetTranslateByCenterPosition(double* c_pos);
  void SetScale(double* scale);
  void SetFlip(bool* flip);
  void SetRotationCenter(double* c_pos);
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

  QString GetPrimaryType() const;

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

  virtual void SetFileName( const QString& fn )
  {
    m_sFilename = fn;
  }

  void SendActorUpdated()
  {
    emit ActorUpdated();
  }

  QString GetSubjectName()
  {
    return m_sSubjectName;
  }

  void ParseSubjectName(const QString& file_path);

  void SetLayerIndex(int n)
  {
    m_nLayerIndex = n;
  }

  int GetLayerIndex()
  {
    return m_nLayerIndex;
  }

  int GetID()
  {
    return m_nID;
  }

  void SetID(int n)
  {
    m_nID = n;
    if (m_nLastID <= m_nID)
      m_nLastID = m_nID+1;
  }

  bool GetAboutToDelete()
  {
    return m_bAboutToDelete;
  }

  void MarkAboutToDelete()
  {
    m_bAboutToDelete = true;
  }

  void CopyTransformation(Layer* layer);

  void AppendTransform(vtkTransform* t);

  void UpdateLastTransform(vtkTransform* t);

  void AppendIdentityTransform();

  void UpdateTotalTransform();

  bool HasTransformUndo()
  {
    return !m_listTransform.isEmpty();
  }

Q_SIGNALS:
  void NameChanged( const QString& name );
  void Transformed();
  void Locked( bool bLock );
  void ActorUpdated();
  void ActorChanged();
  void VisibilityChanged(bool bVisible);

public slots:
  void ResetTransform();
  void UndoLastTransform();

protected:
  virtual bool DoRotate( std::vector<RotationElement>& rotations )
  {
    Q_UNUSED( rotations );
    return true;
  }

  virtual void DoRestore() {}

  virtual void DoTranslate( double* offset ) { Q_UNUSED(offset); }
  virtual void DoScale( double* scale, int nSampleMethod ) { Q_UNUSED(scale); Q_UNUSED(nSampleMethod); }
  virtual void DoTransform( double* mat, int sample_method ) { Q_UNUSED(mat); Q_UNUSED(sample_method); }
  virtual void DoRotate( double* rotate, double* pos) { Q_UNUSED(rotate); Q_UNUSED(pos); }

  // new transform scheme
  virtual void DoTransform(int sample_method) { Q_UNUSED(sample_method); }

  QString   m_strName;
  double    m_dSlicePosition[3];
  double    m_dWorldOrigin[3];
  double    m_dWorldVoxelSize[3];
  double    m_dWorldSize[3];

  // translate and scale are for volume transformation
  double    m_dTranslate[3];
  double    m_dScale[3];
  double    m_dRotate[3];
  bool      m_bFlip[3];
  double    m_dRotationCenter[3];
  bool      m_bRotateAroundCenter;
  bool      m_bUseRotationCenter;

  bool      m_bLocked;
  bool      m_bAboutToDelete;

  LayerProperty*  mProperty;

  QString   m_sFilename;
  QString   m_sSubjectName;
  QStringList m_strTypeNames;
  QString   m_sPrimaryType;
  QList<vtkTransform*> m_listTransform;
  vtkSmartPointer<vtkTransform> m_transform;

  int       m_nID;
  static int m_nLastID;
  int       m_nLayerIndex;
};

#endif


