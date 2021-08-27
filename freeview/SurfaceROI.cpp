/**
 * @brief Surface region from a surface selection in 3D view.
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

#include "SurfaceROI.h"
#include "vtkRenderer.h"
#include "vtkActor2D.h"
#include "vtkProperty.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"
#include "vtkXMLPolyDataWriter.h"
#include "MainWindow.h"
#include "RenderView3D.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkPoints.h"
#include "vtkSelectPolyData.h"
#include "vtkProperty.h"
#include "vtkClipPolyData.h"
#include "vtkBox.h"
#include "vtkMath.h"
#include "MyUtils.h"
#include "LayerSurface.h"
#include "LayerPropertySurface.h"
#include "vtkCleanPolyData.h"
#include "vtkAppendPolyData.h"
#include <QFile>

SurfaceROI::SurfaceROI( LayerSurface* owner ) :
  QObject( owner )
{
  m_actorOutline = vtkSmartPointer<vtkActor>::New();
  m_actorOutline->GetProperty()->SetColor( 0, 1, 0 );
  double ratio = 1;
#if VTK_MAJOR_VERSION > 7
  ratio = MainWindow::GetMainWindow()->devicePixelRatio();
#endif
  m_actorOutline->GetProperty()->SetLineWidth(4*ratio);

  m_points = vtkSmartPointer<vtkPoints>::New();
  m_mris = owner;
  SetColor(Qt::yellow);
}

SurfaceROI::~SurfaceROI()
{}


void SurfaceROI::InitializeOutline(double* pos)
{
  m_actorOutline->VisibilityOn();
  m_points->Reset();
  double* offset = m_mris->GetProperty()->GetPosition();
  double pt[3];
  for (int i = 0; i < 3; i++)
    pt[i] = pos[i] - offset[i];
  m_points->InsertNextPoint(pt);
}

void SurfaceROI::RebuildOutline( bool bClose )
{
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
  lines->InsertNextCell( m_points->GetNumberOfPoints() + (bClose?1:0) );
  for ( int i = 0; i < m_points->GetNumberOfPoints(); i++ )
  {
    lines->InsertCellPoint( i );
  }
  if (bClose)
    lines->InsertCellPoint(0);
  polydata->SetPoints( m_points );
  polydata->SetLines( lines );
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION > 5
  mapper->SetInputData( polydata );
#else
  mapper->SetInput( polydata );
#endif
  m_actorOutline->SetMapper( mapper );
}

void SurfaceROI::AddPoint( double* pos )
{
  double* offset = m_mris->GetProperty()->GetPosition();
  double pt[3];
  for (int i = 0; i < 3; i++)
    pt[i] = pos[i] - offset[i];
  m_points->InsertNextPoint( pt );
  RebuildOutline( false );
}

void SurfaceROI::Close()
{
  RebuildOutline( true );
}

void SurfaceROI::Update()
{}

void SurfaceROI::AppendProps( vtkRenderer* renderer )
{
  renderer->AddViewProp( m_actorOutline );
}

void SurfaceROI::SetColor( const QColor& color )
{
  m_color = color;
  m_actorOutline->GetProperty()->SetColor( color.redF(), color.greenF(), color.blueF() );
  emit ColorChanged( color );
}

QColor SurfaceROI::GetColor()
{
  return m_color;
}

void SurfaceROI::Show( bool bShow )
{
  m_actorOutline->SetVisibility( bShow?1:0 );
}

vtkActor* SurfaceROI::GetActor()
{
  return m_actorOutline;
}
