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

#include "VolumeCropper.h"
#include "LayerMRI.h"
#include "LayerPropertyMRI.h"
#include "RenderView2D.h"
#include <vtkBox.h>
#include <vtkCubeSource.h>
#include <vtkActor.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkClipPolyData.h>
#include <vtkImageData.h>
#include <vtkRenderer.h>
#include <vtkProperty.h>
#include <vtkSphereSource.h>
#include <vtkImageActor.h>
#include <vtkMath.h>
#include <vtkPlaneSource.h>

VolumeCropper::VolumeCropper( QObject* parent ) :
  QObject( parent ),
  m_mri( NULL ),
  m_bEnabled( false ),
  m_nActivePlane( -1 )
{
  m_actorBox = vtkSmartPointer<vtkActor>::New();
  m_actorBox->VisibilityOff();
  m_actorBox->GetProperty()->SetOpacity( 0.2 );
  m_actorBox->GetProperty()->SetColor( 1, 1, 0 );
  m_actorBox->PickableOff();

  double ratio = 1;
#if VTK_MAJOR_VERSION > 7
  QWidget* w = qobject_cast<QWidget*>(parent);
  if (w)
    ratio = w->devicePixelRatio();
#endif
  m_actorFrame = vtkSmartPointer<vtkActor>::New();
  m_actorFrame->VisibilityOff();
  m_actorFrame->GetProperty()->SetColor( 1, 1, 0 );
  m_actorFrame->GetProperty()->SetLineWidth( 2*ratio );
  m_actorFrame->GetProperty()->SetRepresentationToWireframe();
  m_actorFrame->GetProperty()->SetDiffuse( 0.0 );
  m_actorFrame->GetProperty()->SetAmbient( 1.0 );

  for ( int i = 0; i < 6; i++ )
  {
    m_actorSphere[i] = vtkSmartPointer<vtkActor>::New();
    m_actorSphere[i]->VisibilityOff();
    m_sphereSource[i] = vtkSmartPointer<vtkSphereSource>::New();
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection( m_sphereSource[i]->GetOutputPort() );
    m_actorSphere[i]->SetMapper( mapper );
    m_actorSphere[i]->GetProperty()->SetColor( 1, 1, 0 );
  }

  m_actorActivePlane = vtkSmartPointer<vtkActor>::New();
  m_actorActivePlane->VisibilityOff();
  m_actorActivePlane->GetProperty()->SetOpacity( 0.3 );
  m_actorActivePlane->GetProperty()->SetColor( 1, 0, 0 );
  m_planeSource = vtkSmartPointer<vtkPlaneSource>::New();
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection( m_planeSource->GetOutputPort() );
  m_actorActivePlane->SetMapper( mapper );

  for ( int i = 0; i < 3; i++ )
  {
    m_actorActivePlane2D[i] = vtkSmartPointer<vtkActor>::New();
    m_actorActivePlane2D[i]->VisibilityOff();
    m_actorActivePlane2D[i]->GetProperty()->SetColor( 1, 0, 0 );
    m_actorActivePlane2D[i]->GetProperty()->SetLineWidth( 2*ratio );
    m_actorActivePlane2D[i]->GetProperty()->SetRepresentationToWireframe();
    m_actorActivePlane2D[i]->GetProperty()->SetDiffuse( 0.0 );
    m_actorActivePlane2D[i]->GetProperty()->SetAmbient( 1.0 );
    mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection( m_planeSource->GetOutputPort() );
    m_actorActivePlane2D[i]->SetMapper( mapper );
  }

  m_boxSource = vtkSmartPointer<vtkCubeSource>::New();
  mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection( m_boxSource->GetOutputPort() );
  m_actorBox->SetMapper( mapper );
  mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection( m_boxSource->GetOutputPort() );
  m_actorFrame->SetMapper( mapper );

  m_box = vtkSmartPointer<vtkBox>::New();

  m_clipper = vtkSmartPointer<vtkClipPolyData>::New();
  m_clipper->SetClipFunction( m_box );
  m_clipper->InsideOutOn();

  for (int i = 0; i < 3; i++)
  {
    m_actorBox2D[i] = vtkSmartPointer<vtkActor>::New();
    m_actorBox2D[i]->VisibilityOff();
    m_actorBox2D[i]->SetProperty( m_actorBox->GetProperty() );
    m_actorFrame2D[i] = vtkSmartPointer<vtkActor>::New();
    m_actorFrame2D[i]->VisibilityOff();
    m_actorFrame2D[i]->SetProperty( m_actorFrame->GetProperty() );

    mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection( m_boxSource->GetOutputPort() );
    m_actorBox2D[i]->SetMapper( mapper );
    mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection( m_boxSource->GetOutputPort() );
    m_actorFrame2D[i]->SetMapper( mapper );
  }

#if VTK_MAJOR_VERSION > 5
  m_actorBox->ForceOpaqueOn();
  m_actorFrame->ForceOpaqueOn();
  m_actorActivePlane->ForceOpaqueOn();
  for (int i = 0; i < 3; i++)
  {
    m_actorActivePlane2D[i]->ForceOpaqueOn();
    m_actorFrame2D[i]->ForceOpaqueOn();
    m_actorBox2D[i]->ForceOpaqueOn();
  }
  for (int i = 0; i < 6; i++)
  {
    m_actorSphere[i]->ForceOpaqueOn();
  }
#endif
}

VolumeCropper::~VolumeCropper()
{}

vtkActor* VolumeCropper::GetProp()
{
  return m_actorBox;
}

void VolumeCropper::Append3DProps( vtkRenderer* renderer )
{
  renderer->AddViewProp( m_actorFrame );
  for ( int i = 0; i < 6; i++ )
  {
    renderer->AddViewProp( m_actorSphere[i] );
  }
  renderer->AddViewProp( m_actorActivePlane );
  renderer->AddViewProp( m_actorBox );
}

void VolumeCropper::Append2DProps( vtkRenderer* renderer, int n )
{
  renderer->AddViewProp( m_actorBox2D[n] );
  renderer->AddViewProp( m_actorFrame2D[n] );
  renderer->AddViewProp( m_actorActivePlane2D[n] );
}

void VolumeCropper::SetVolume( LayerMRI* mri)
{
  if ( m_mri != mri )
  {
    if (m_mri)
    {
      disconnect(m_mri, SIGNAL(IsoSurfaceUpdated()), this, SLOT(Apply()));
    }
    m_mri = mri;
    connect(m_mri, SIGNAL(IsoSurfaceUpdated()), this, SLOT(Apply()), Qt::UniqueConnection);
    Reset();
  }
}

void VolumeCropper::Reset()
{
  m_mri->GetDisplayBounds( m_bounds );
  m_boxSource->SetBounds( m_bounds );
  double dMax = 0;
  for ( int i = 0; i < 6; i+=2 )
  {
    if ( dMax < (m_bounds[i+1] - m_bounds[i]) )
    {
      dMax = m_bounds[i+1] - m_bounds[i];
    }
  }
  for ( int i = 0; i < 6; i++ )
  {
    m_sphereSource[i]->SetRadius( dMax/75.0 );
  }

  UpdateExtent();
  m_box->SetBounds( m_bounds );
  UpdateProps();
}

// apply contour clipper
void VolumeCropper::Apply()
{
  m_box->SetBounds( m_bounds );
}

void VolumeCropper::UpdateProps()
{
  UpdateSliceActorVisibility();

  double c[3];
  for ( int i = 0; i < 3; i++ )
  {
    c[i] = ( m_bounds[i*2] + m_bounds[i*2+1] ) / 2;
  }

  double offset = m_sphereSource[0]->GetRadius()/2;
  m_sphereSource[0]->SetCenter( m_bounds[0]-offset, c[1], c[2] );
  m_sphereSource[1]->SetCenter( m_bounds[1]+offset, c[1], c[2] );
  m_sphereSource[2]->SetCenter( c[0], m_bounds[2]-offset, c[2] );
  m_sphereSource[3]->SetCenter( c[0], m_bounds[3]+offset, c[2] );
  m_sphereSource[4]->SetCenter( c[0], c[1], m_bounds[4]-offset );
  m_sphereSource[5]->SetCenter( c[0], c[1], m_bounds[5]+offset );

  m_boxSource->SetBounds( m_bounds );

  int nScale = 1;
  if ( m_mri->GetProperty()->GetShowLabelOutline() )
  {
    nScale = 4;
  }
  else if ( m_mri->GetProperty()->GetUpSampleMethod() )
  {
    nScale = 2;
  }
  int ext[6];
  for ( int i = 0; i < 6; i++ )
  {
    ext[i] = nScale * m_extent[i];
  }
  m_mri->m_sliceActor2D[0]->SetDisplayExtent( ext[2], ext[3], ext[4], ext[5], 0, 0 );
  m_mri->m_sliceActor2D[1]->SetDisplayExtent( ext[0], ext[1], ext[4], ext[5], 0, 0 );
  m_mri->m_sliceActor2D[2]->SetDisplayExtent( ext[0], ext[1], ext[2], ext[3], 0, 0 );
  for ( int i = 0; i < 3; i++ )
  {
    m_mri->m_sliceActor3D[i]->SetDisplayExtent( m_mri->m_sliceActor2D[i]->GetDisplayExtent() );
  }

  vtkPolyDataMapper* mapper = vtkPolyDataMapper::SafeDownCast( m_mri->m_actorContour->GetMapper() );
  if ( mapper && mapper->GetInput() && m_clipper->GetOutput() != mapper->GetInput() )
  {
#if VTK_MAJOR_VERSION > 5
    m_clipper->SetInputData( mapper->GetInput() );
#else
    m_clipper->SetInput( mapper->GetInput() );
#endif
    mapper->SetInputConnection( m_clipper->GetOutputPort() );
  }

  emit CropBoundChanged( m_mri );
}

void VolumeCropper::UpdateExtent()
{
  vtkImageData* image = m_mri->GetImageData();
  double* orig = image->GetOrigin();
  double* voxelsize = image->GetSpacing();
  for ( int i = 0; i < 6; i++ )
  {
    m_extent[i] = (int)( (m_bounds[i] - orig[i/2])/voxelsize[i/2] + 0.5 );
  }

  m_mri->SetCroppingBounds( m_bounds );
}

void VolumeCropper::SetEnabled( bool bEnable )
{
  m_bEnabled = bEnable;
  if ( !bEnable )
  {
    Show( false );
  }
}

void VolumeCropper::Show( bool bShow )
{
  m_actorBox->SetVisibility( bShow );
  m_actorFrame->SetVisibility( bShow );
  for ( int i = 0; i < 6; i++ )
  {
    m_actorSphere[i]->SetVisibility( bShow );
  }
  for ( int i = 0; i < 3; i++ )
  {
    m_actorBox2D[i]->SetVisibility( bShow );
    m_actorFrame2D[i]->SetVisibility( bShow );
  }
  if ( !bShow )
  {
    m_actorActivePlane->SetVisibility( false );
    for ( int i = 0; i < 3; i++ )
    {
      m_actorActivePlane2D[i]->SetVisibility( false );
    }
  }
}

bool VolumeCropper::IsShown()
{
  return m_actorBox->GetVisibility();
}

bool VolumeCropper::PickActiveBound( vtkProp* prop )
{
  m_nActivePlane = -1;
  for ( int i = 0; i < 6; i++ )
  {
    if ( m_actorSphere[i].GetPointer() == prop )
    {
      m_nActivePlane = i;
      m_actorSphere[i]->GetProperty()->SetColor( 1, 0, 0 );
      UpdateActivePlane();
      m_actorActivePlane->VisibilityOn();
      return true;
    }
  }
  return false;
}

bool VolumeCropper::PickActiveBound2D( RenderView2D* view, int x, int y )
{
  int nPlane = view->GetViewPlane();
  if ( !IsShown() )
  {
    m_actorActivePlane2D[nPlane]->VisibilityOff();
    return false;
  }

  double pt[3];
  view->ScreenToWorld( x, y, 0, pt[0], pt[1], pt[2] );
  pt[nPlane] = ( m_bounds[nPlane*2] + m_bounds[nPlane*2+1] ) / 2;
  double dist = 1e10;
  int n = -1;
  for ( int i = 0; i < 6; i++ )
  {
    if ( i/2 != nPlane && fabs( m_bounds[i] - pt[i/2] ) < dist )
    {
      int m = 0;
      for ( m = 0; m < 3; m++ )
      {
        if ( m != nPlane && m != i/2 )
        {
          break;
        }
      }
      if ( pt[m] <= m_bounds[m*2+1] && pt[m] >= m_bounds[m*2] )
      {
        dist = fabs( m_bounds[i] - pt[i/2] );
        n = i;
      }
    }
  }

  int nOldActivePlane = m_nActivePlane;
  int nTolerance = 5; // in pixels
  if ( n >= 0 )
  {
    for ( int i = 0; i < 3; i++ )
    {
      pt[i] = m_bounds[i*2];
    }
    pt[n/2] = m_bounds[n];
    int x2, y2;
    view->WorldToScreen( pt[0], pt[1], pt[2], x2, y2 );
    if ( fabs( x - x2 ) <= nTolerance || fabs( y - y2 ) <= nTolerance )
    {
      m_nActivePlane = n;
      UpdateActivePlane();
      m_actorActivePlane2D[nPlane]->VisibilityOn();
      view->RequestRedraw();
      return true;
    }
  }
  m_nActivePlane = -1;
  m_actorActivePlane2D[nPlane]->VisibilityOff();
  if ( m_nActivePlane != nOldActivePlane )
  {
    view->RequestRedraw();
  }
  return false;
}

void VolumeCropper::ReleaseActiveBound()
{
  if ( m_nActivePlane >= 0 )
  {
    m_actorSphere[m_nActivePlane]->GetProperty()->SetColor( 1, 1, 0 );
  }

  m_nActivePlane = -1;
  // update clipper box
  m_box->SetBounds( m_bounds );
  m_actorActivePlane->VisibilityOff();
}

void VolumeCropper::MoveActiveBound( RenderView* view, int nx, int ny )
{
  if ( !m_mri )
  {
    return;
  }

  double dMaxBounds[6];
  m_mri->GetDisplayBounds( dMaxBounds );
  int x1, y1, x2, y2;
  double pt1[3], pt2[3];
  for ( int i = 0; i < 3; i++ )
  {
    pt1[i] = ( dMaxBounds[i*2] + dMaxBounds[i*2+1] ) / 2;
    pt2[i] = pt1[i];
  }
  int n = m_nActivePlane/2;
  pt1[n] = dMaxBounds[m_nActivePlane];
  if ( m_nActivePlane%2 == 0 )
  {
    pt2[n] = dMaxBounds[m_nActivePlane+1];
  }
  else
  {
    pt2[n] = dMaxBounds[m_nActivePlane-1];
  }
  view->WorldToScreen( pt1[0], pt1[1], pt1[2], x1, y1 );
  view->WorldToScreen( pt2[0], pt2[1], pt2[2], x2, y2 );
  double ratio = 0;
  if ( fabs( x1-x2 ) > fabs( y1-y2 ) )
  {
    ratio = (x2==x1) ? 0 : ( nx / sqrt( (x2-x1)*(x2-x1) ) * (x2>x1?1:-1) );
  }
  else
  {
    ratio = (y2==y1) ? 0 : ( ny / sqrt( (y2-y1)*(y2-y1) ) * (y2>y1?1:-1) );
  }

  m_bounds[m_nActivePlane] += ratio * (pt2[n] - pt1[n]);

  ValidateActiveBound();
  UpdateActivePlane();
  UpdateExtent();
  UpdateProps();
}

void VolumeCropper::MoveActiveBound2D( RenderView2D* view, int x, int y )
{
  if ( m_nActivePlane < 0 )
  {
    return;
  }

  double pos[3];
  view->ScreenToWorld( x, y, 0, pos[0], pos[1], pos[2] );
  m_bounds[m_nActivePlane] = pos[m_nActivePlane/2];

  ValidateActiveBound();
  UpdateActivePlane();
  UpdateExtent();
  UpdateProps();
}

void VolumeCropper::ValidateActiveBound()
{
  if ( m_nActivePlane < 0 || !m_mri )
  {
    return;
  }

  double dMaxBounds[6];
  m_mri->GetDisplayBounds( dMaxBounds );

  int n = m_nActivePlane/2;
  if ( m_bounds[m_nActivePlane] < dMaxBounds[n*2] )
  {
    m_bounds[m_nActivePlane] = dMaxBounds[n*2];
  }
  else if ( m_bounds[m_nActivePlane] > dMaxBounds[n*2+1] )
  {
    m_bounds[m_nActivePlane] = dMaxBounds[n*2+1];
  }

  double* voxelsize = m_mri->GetImageData()->GetSpacing();
  if ( m_nActivePlane%2 == 0 && m_bounds[m_nActivePlane] >= m_bounds[m_nActivePlane+1] )
  {
    m_bounds[m_nActivePlane] = m_bounds[m_nActivePlane+1] - voxelsize[n]/2;
  }
  else if ( m_nActivePlane%2 == 1 && m_bounds[m_nActivePlane] < m_bounds[m_nActivePlane-1] )
  {
    m_bounds[m_nActivePlane] = m_bounds[m_nActivePlane-1] + voxelsize[n]/2;
  }
}

void VolumeCropper::SetExtent( int nComp, int nValue )
{
  m_extent[nComp] = nValue;
  if ( nComp%2 == 0 && m_extent[nComp+1] < nValue )
  {
    m_extent[nComp] = m_extent[nComp+1];
  }
  else if ( nComp%2 == 1 && m_extent[nComp-1] > nValue )
  {
    m_extent[nComp] = m_extent[nComp-1];
  }

  vtkImageData* image = m_mri->GetImageData();
  double* orig = image->GetOrigin();
  double* voxelsize = image->GetSpacing();
  m_bounds[nComp] = m_extent[nComp]*voxelsize[nComp/2] - orig[nComp/2];

  m_mri->SetCroppingBounds( m_bounds );
  UpdateProps();
}

void VolumeCropper::UpdateSliceActorVisibility()
{
  if ( !m_mri || !IsShown() )
  {
    return;
  }

  double* pos = m_mri->GetSlicePosition();
  for ( int i = 0; i < 3; i++ )
  {
    if ( pos[i] < m_bounds[i*2] || pos[i] > m_bounds[i*2+1] )
    {
      m_mri->m_sliceActor2D[i]->SetVisibility( false );
      m_mri->m_sliceActor3D[i]->SetVisibility( false );
    }
    else
    {
      m_mri->m_sliceActor2D[i]->SetVisibility( true );
      m_mri->m_sliceActor3D[i]->SetVisibility( true );
    }
  }
}

void VolumeCropper::UpdateActivePlane()
{
  if ( m_nActivePlane >= 0 )
  {
    double pt[3], pt1[3], pt2[3];
    for ( int i = 0; i < 3; i++ )
    {
      pt[i] = m_bounds[i*2];
    }
    pt[m_nActivePlane/2] = m_bounds[m_nActivePlane];
    for ( int i = 0; i < 3; i++ )
    {
      pt1[i] = pt2[i] = pt[i];
    }
    switch ( m_nActivePlane/2 )
    {
    case 0:
      pt1[1] = m_bounds[3];
      pt2[2] = m_bounds[5];
      break;
    case 1:
      pt1[0] = m_bounds[1];
      pt2[2] = m_bounds[5];
      break;
    case 2:
      pt1[0] = m_bounds[1];
      pt2[1] = m_bounds[3];
      break;
    }
    m_planeSource->SetOrigin( 1000, 200, 50 );   // random set origin/points first to work around a dumb bug in vtk
    m_planeSource->SetPoint1( 3, 3, 0 );
    m_planeSource->SetPoint2( 200, 45, 10000 );
    m_planeSource->SetPoint1( pt1 );
    m_planeSource->SetPoint2( pt2 );
    m_planeSource->SetOrigin( pt );
  }
}


