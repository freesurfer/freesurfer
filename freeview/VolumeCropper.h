/**
 * @brief Class to crop volume.
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

#ifndef VolumeCropper_h
#define VolumeCropper_h

#include "vtkSmartPointer.h"
#include <QObject>

class vtkBox;
class vtkCubeSource;
class vtkSphereSource;
class vtkPlaneSource;
class vtkActor;
class vtkProp;
class vtkRenderer;
class vtkClipPolyData;
class LayerMRI;
class RenderView;
class RenderView2D;

class VolumeCropper : public QObject
{
  Q_OBJECT
public:
  VolumeCropper( QObject* parent = NULL );
  virtual ~VolumeCropper();

  void SetVolume( LayerMRI* mri );

  LayerMRI* GetVolume()
  {
    return m_mri;
  }

  void SetEnabled( bool bEnable );

  bool GetEnabled()
  {
    return m_bEnabled;
  }

  void Show( bool bShow = true );

  vtkActor* GetProp();

  void Append3DProps( vtkRenderer* renderer );
  void Append2DProps( vtkRenderer* renderer, int n );

  bool IsShown();

  bool PickActiveBound( vtkProp* prop );
  void MoveActiveBound( RenderView* view, int nx, int ny );
  void ReleaseActiveBound();

  bool PickActiveBound2D( RenderView2D* view, int nX, int nY );
  void MoveActiveBound2D( RenderView2D* view, int x, int y );

  double* GetBounds()
  {
    return m_bounds;
  }

  int* GetExtent()
  {
    return m_extent;
  }

  void SetExtent( int nComp, int nValue );

  void UpdateProps();

  int GetActivePlane()
  {
    return m_nActivePlane;
  }

Q_SIGNALS:
  void CropBoundChanged( LayerMRI* );

public slots:
  void Reset();
  void Apply();

protected:
  void UpdateExtent();
  void UpdateSliceActorVisibility();
  void UpdateActivePlane();
  void ValidateActiveBound();

  vtkSmartPointer<vtkBox>         m_box;
  vtkSmartPointer<vtkCubeSource>  m_boxSource;
  vtkSmartPointer<vtkPlaneSource> m_planeSource;
  vtkSmartPointer<vtkActor>     m_actorBox;
  vtkSmartPointer<vtkActor>     m_actorFrame;
  vtkSmartPointer<vtkActor>     m_actorActivePlane;
  vtkSmartPointer<vtkActor>     m_actorSphere[6];
  vtkSmartPointer<vtkActor>     m_actorBox2D[3];
  vtkSmartPointer<vtkActor>     m_actorFrame2D[3];
  vtkSmartPointer<vtkActor>     m_actorActivePlane2D[3];
  vtkSmartPointer<vtkSphereSource>  m_sphereSource[6];
  vtkSmartPointer<vtkClipPolyData>  m_clipper;

  LayerMRI*         m_mri;
  double            m_bounds[6];
  int               m_extent[6];
  bool              m_bEnabled;

  int               m_nActivePlane;
};

#endif


