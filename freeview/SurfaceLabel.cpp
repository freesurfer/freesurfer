/**
 * @file  SurfaceLabel.cxx
 * @brief Implementation for surface label.
 *
 * In 2D, the MRI is viewed as a single slice, and controls are
 * provided to change the color table and other viewing options. In
 * 3D, the MRI is viewed in three planes in 3D space, with controls to
 * move each plane axially.
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2012/03/13 21:32:06 $
 *    $Revision: 1.10 $
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


#include "SurfaceLabel.h"
#include "LayerSurface.h"
#include "LayerPropertySurface.h"
#include "FSSurface.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkAppendPolyData.h"
#include "vtkCellArray.h"
#include "vtkMath.h"
#include "MyUtils.h"
#include <QFile>
#include <QDebug>

SurfaceLabel::SurfaceLabel ( LayerSurface* surf ) :
  QObject( surf ),
  m_label( NULL ),
  m_surface( surf )
{
  m_rgbColor[0] = 1.0;
  m_rgbColor[1] = 1.0;
  m_rgbColor[2] = 0.0;

  m_actorOutline = vtkSmartPointer<vtkActor>::New();
  m_actorOutline->SetPosition(surf->GetProperty()->GetPosition());
  m_actorOutline->GetProperty()->SetLineWidth(4);

  SetShowOutline(false);
  SetColor(1.0, 1.0, 0);
}

SurfaceLabel::~SurfaceLabel ()
{
  if ( m_label )
  {
    ::LabelFree( &m_label );
  }
}


QString SurfaceLabel::GetName()
{
  return m_strName;
}

void SurfaceLabel::SetName( const QString& name )
{
  m_strName = name;
}

bool SurfaceLabel::LoadLabel( const QString& filename )
{
  if ( m_label )
  {
    ::LabelFree( &m_label );
  }

  m_label = ::LabelRead( NULL, filename.toAscii().data() );

  if ( m_label == NULL )
  {
    cerr << "LabelRead failed";
    return false;
  }

  QFile file( filename );
  if ( !file.open( QIODevice::ReadOnly | QIODevice::Text ) )
  {
    return false;
  }

  QString strg;
  while (!file.atEnd())
  {
    strg += file.readLine();
  }

  if ( strg.contains( "vox2ras=", Qt::CaseInsensitive ) &&
       !strg.contains( "vox2ras=TkReg", Qt::CaseInsensitive ) )
  {
    m_bTkReg = false;
  }

  // create outline
  if (m_surface)
  {
    MRIS* mris = m_surface->GetSourceSurface()->GetMRIS();
    VERTEX *v;

    MRISclearMarks(mris);
    LabelMarkSurface(m_label, mris);
    QList<int> vertices;
    for (int n = 0 ; n < m_label->n_points ; n++)
    {
      if (m_label->lv[n].vno >= 0)
      {
        v = &mris->vertices[m_label->lv[n].vno] ;
        if (v->ripflag)
          continue ;

        for (int m = 0 ; m < v->vnum ; m++)
        {
          if (mris->vertices[v->v[m]].marked == 0)
          {
            if (!vertices.contains(m_label->lv[n].vno))
              vertices << m_label->lv[n].vno;
            break;
          }
        }
      }
    }
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    QList<int> neighbor_count;
    for (int i = 0; i < vertices.size(); i++)
    {
      v = &mris->vertices[vertices[i]];
      double pt[3] = { v->x, v->y, v->z };
      m_surface->GetTargetAtSurfaceRAS(pt, pt);
      points->InsertNextPoint(pt);
      int n = 0;
      for (int m = 0; m < v->vnum; m++)
      {
        if (vertices.contains(v->v[m]))
          n++;
      }
      neighbor_count << n;
    }
    QList< QPair<int, int> > pairs;
    for (int i = 0; i < vertices.size(); i++)
    {
      v = &mris->vertices[vertices[i]];
      for (int m = 0; m < v->vnum; m++)
      {
        int n = vertices.indexOf(v->v[m]);
        if (n >= 0)
        {
        //  if (neighbor_count[i] < 3 || neighbor_count[n] < 3)
          QPair<int, int> p(i, n);
          if (i > n)
            p = QPair<int, int>(n, i);
          if (!pairs.contains(p))
          {
            lines->InsertNextCell(2);
            lines->InsertCellPoint(i);
            lines->InsertCellPoint(n);
            pairs << p;
          }
        }
      }
    }

    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);
    polydata->SetLines(lines);
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInput( polydata);
    m_actorOutline->SetMapper( mapper );
  }

  return true;
}

vtkPolyData* SurfaceLabel::MakeEdgePolyData(const QList<int> &indices_in, const QList<int> &vertices)
{
  MRIS* mris = m_surface->GetSourceSurface()->GetMRIS();
  VERTEX *v;
  QList<int> indices = indices_in;
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkPolyData* polydata = vtkPolyData::New();
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
  lines->InsertNextCell( indices.size() + 1 );
  for ( int i = 0; i < indices.size(); i++ )
  {
    lines->InsertCellPoint( i );
    v = &mris->vertices[indices[i]];
    double pt[3] = { v->x, v->y, v->z };
    m_surface->GetTargetAtSurfaceRAS(pt, pt);
    points->InsertNextPoint(pt);
  }
  lines->InsertCellPoint(0);

  polydata->SetPoints( points );
  polydata->SetLines( lines );

  // Has a hole in the label, recursively add it to the overlay polydata
  if (indices.size() < vertices.size())
  {
    QList<int> sub_vertices = vertices;
    for (int i = 0; i < indices.size(); i++)
    {
      if (sub_vertices.contains(indices[i]))
        sub_vertices.removeOne(indices[i]);
    }

    indices.clear();
    indices << sub_vertices[0];
    indices = DoConnectEdgeVertices(indices, sub_vertices);
    vtkPolyData* polydata2 = MakeEdgePolyData(indices, sub_vertices);
    vtkSmartPointer<vtkAppendPolyData> append = vtkSmartPointer<vtkAppendPolyData>::New();
    append->AddInput(polydata);
    append->AddInput(polydata2);
    append->Update();
    vtkPolyData* polydata_all = vtkPolyData::New();
    polydata_all->DeepCopy(append->GetOutput());
    return polydata_all;
  }
  else
    return polydata;
}

QList<int> SurfaceLabel::DoConnectEdgeVertices(const QList<int> &prev_indices, const QList<int> &vertices)
{
  MRIS* mris = m_surface->GetSourceSurface()->GetMRIS();
  VERTEX *v;

//  qDebug() << vertices.size() << prev_indices.size();
  QList<int> indices = prev_indices;
  while (indices.size() < vertices.size())
  {
    v = &mris->vertices[indices.last()];
    QList<int> connects;
    for (int m = 0 ; m < v->vnum ; m++)
    {
      if (!indices.contains(v->v[m]) && vertices.contains(v->v[m]))
      {
        connects << v->v[m];
      }
    }
    if (!connects.isEmpty())
    {
      if (connects.size() == 1)
        indices << connects[0];
      else
      {
        int ncount = 0;
        int n = -1;
        QList<int> solutions[10];
        for (int i = 0; i < connects.size(); i++)
        {
          QList<int> sol_indices = indices;
          sol_indices << connects[i];
          solutions[i] = DoConnectEdgeVertices(sol_indices, vertices);
          if (solutions[i].size() > ncount)
          {
            ncount = solutions[i].size();
            n = i;
          }
        }
        if (n == -1)
        {
          qDebug() << n;  // should not happen
          n = 0;
        }
        return solutions[n];
      }
    }
    else
    {
      return indices;
    }
  }
  return indices;
}

void SurfaceLabel::SetColor( double r, double g, double b )
{
  m_rgbColor[0] = r;
  m_rgbColor[1] = g;
  m_rgbColor[2] = b;

  m_actorOutline->GetProperty()->SetColor(r, g, b);

  emit SurfaceLabelChanged();
}

void SurfaceLabel::MapLabel( unsigned char* colordata, int nVertexCount )
{
  if ( !m_label || m_bShowOutline)
  {
    return;
  }

  for ( int i = 0; i < m_label->n_points; i++ )
  {
    int vno = m_label->lv[i].vno;
    if ( vno < nVertexCount )
    {
      double opacity = m_label->lv[i].stat;
      if ( opacity > 1 )
      {
        opacity = 1;
      }
      else if ( opacity < 0 )
      {
        opacity = 0;
      }
      opacity = 1;    // ignore opacity for now
      colordata[vno*4]    = ( int )( colordata[vno*4]   * ( 1 - opacity ) + m_rgbColor[0] * 255 * opacity );
      colordata[vno*4+1]  = ( int )( colordata[vno*4+1] * ( 1 - opacity ) + m_rgbColor[1] * 255 * opacity );
      colordata[vno*4+2]  = ( int )( colordata[vno*4+2] * ( 1 - opacity ) + m_rgbColor[2] * 255 * opacity );
    }
  }
}

void SurfaceLabel::SetShowOutline(bool bOutline)
{
  this->m_bShowOutline = bOutline;
  m_actorOutline->SetVisibility(bOutline?1:0);
}

vtkActor* SurfaceLabel::GetOutlineActor()
{
  return m_actorOutline;
}
