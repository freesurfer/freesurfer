/**
 * @file  SurfaceAnnotation.cxx
 * @brief Implementation for surface annotation.
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
 *    $Date: 2011/08/02 15:58:25 $
 *    $Revision: 1.15 $
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


#include "SurfaceAnnotation.h"
#include "vtkLookupTable.h"
#include "vtkRGBAColorTransferFunction.h"
#include "vtkMath.h"
#include "LayerSurface.h"
#include "FSSurface.h"
#include <QFileInfo>
#include <vtkActor.h>
#include <vtkAppendPolyData.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <QDebug>

SurfaceAnnotation::SurfaceAnnotation ( LayerSurface* surf ) :
  QObject( surf ),
  m_nIndices( NULL ),
  m_nCenterVertices( NULL ),
  m_lut( NULL ),
  m_surface( surf )
{
  m_actorOutline = vtkSmartPointer<vtkActor>::New();
}

SurfaceAnnotation::~SurfaceAnnotation ()
{
  Reset();
}

void SurfaceAnnotation::Reset()
{
  if ( m_nIndices )
  {
    delete[] m_nIndices;
  }

  if ( m_nCenterVertices )
  {
    delete[] m_nCenterVertices;
  }

  m_nIndices = NULL;
  m_lut = NULL;
  m_nCenterVertices = NULL;
}

/*
void SurfaceAnnotation::DoListenToMessage ( std::string const iMessage, void* iData, void* sender )
{
  if ( iMessage == "ColorMapChanged" )
  {
    this->SendBroadcast( "AnnotationChanged", this );
  }
}
*/

bool SurfaceAnnotation::LoadAnnotation( const QString& fn )
{
  if ( m_surface )
  {
    MRIS* mris = m_surface->GetSourceSurface()->GetMRIS();
    m_nIndices = NULL;

    if ( MRISreadAnnotation( mris, fn.toAscii().data() ) != 0 )
    {
      cerr << "Could not load annotation from file " << qPrintable(fn) << ".\n";
      return false;
    }
    else
    {
      Reset();
      m_lut = mris->ct;
      m_nIndexSize = mris->nvertices;
      m_nIndices = new int[m_nIndexSize];
      CTABgetNumberOfValidEntries( m_lut, &m_nAnnotations );
      m_nCenterVertices = new int[m_nAnnotations];
      for ( int i = 0; i < m_nAnnotations; i++ )
      {
        m_nCenterVertices[i] = -1;
      }

      // convert annotations to lookup table indices
      // also find the "center vertex"
      double** pts = new double*[m_nAnnotations];
      int* vcount = new int[m_nAnnotations];
      memset( vcount, 0, sizeof( int )*m_nAnnotations );
      for ( int i = 0; i < m_nAnnotations; i++ )
      {
        pts[i] = new double[3];
        memset( pts[i], 0, sizeof( double ) * 3 );
      }

      int n;
      for ( int i = 0; i < m_nIndexSize; i++ )
      {
        if ( CTABfindAnnotation( m_lut, mris->vertices[i].annotation, &n ) == 0 )
        {
          m_nIndices[i] = n;

          if ( mris->vertices[i].ripflag == 0 &&
               m_nIndices[i] != -1 )
          {
            vcount[n]++;
            pts[n][0] += mris->vertices[i].x;
            pts[n][1] += mris->vertices[i].y;
            pts[n][2] += mris->vertices[i].z;
            m_nCenterVertices[n] = i;
          }
        }
        else
        {
          m_nIndices[i] = -1;
        }
      }

      for ( int i = 0; i < m_nAnnotations; i++ )
      {
        if ( vcount[i] > 0 )
        {
          pts[i][0] /= vcount[i];
          pts[i][1] /= vcount[i];
          pts[i][2] /= vcount[i];

          int nVertex = m_surface->GetSourceSurface()->FindVertexAtSurfaceRAS( pts[i], NULL );
          if ( nVertex >=0 && m_nIndices[nVertex] == i )
          {
            m_nCenterVertices[i] = nVertex;
          }
        }
      }

      delete[] vcount;
      for ( int i = 0; i < m_nAnnotations; i++ )
      {
        delete[] pts[i];
      }
      delete[] pts;

      // build outline actor, not work yet
      /*
      vtkSmartPointer<vtkAppendPolyData> append = vtkSmartPointer<vtkAppendPolyData>::New();
      for (int i = 1; i < m_nAnnotations; i++)
      {
        if (true)
        {
          VERTEX *v;
          MRISclearMarks(mris);
          LABEL* label = MRISannotation_to_label(mris, i);
          if (label)
          {

            LabelMarkSurface(label, mris);
            QList<int> vertices;
            for (int n = 0 ; n < label->n_points ; n++)
            {
              if (label->lv[n].vno >= 0)
              {
                v = &mris->vertices[label->lv[n].vno] ;
                if (v->ripflag)
                  continue ;

                for (int m = 0 ; m < v->vnum ; m++)
                {
                  if (mris->vertices[v->v[m]].marked == 0)
                  {
                    vertices << label->lv[n].vno;
                    break;
                  }
                }
              }
            }

            QList<int> indices;
            indices << vertices[0];
            indices = this->DoConnectEdgeVertices(indices, vertices);
            qDebug() << i;
            vtkSmartPointer<vtkPolyData> polydata = MakeEdgePolyData(indices, vertices);
            append->AddInput(polydata);
            LabelFree(&label);
          }
        }
      }

      vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
      mapper->SetInput( append->GetOutput() );
      m_actorOutline->SetMapper( mapper );
      m_actorOutline->GetProperty()->SetColor(1,1,1);*/
      return true;
    }
  }
  return false;
}

void SurfaceAnnotation::GetAnnotationPoint( int nIndex, double* pt_out )
{
  m_surface->GetTargetAtVertex( m_nCenterVertices[nIndex], pt_out );
}

QString SurfaceAnnotation::GetName()
{
  return m_strName;
}

void SurfaceAnnotation::SetName( const QString& name )
{
  m_strName = name;
}

int SurfaceAnnotation::GetIndexAtVertex( int nVertex )
{
  return m_nIndices[nVertex];
}

QString SurfaceAnnotation::GetAnnotationNameAtVertex( int nVertex )
{
  return GetAnnotationNameAtIndex( m_nIndices[nVertex] );
}

QString SurfaceAnnotation::GetAnnotationNameAtIndex( int nIndex )
{
  char name[128];
  int nValid = 0;
  int nTotalCount = 0;
  CTABgetNumberOfTotalEntries( m_lut, &nTotalCount );
  if ( nIndex < nTotalCount )
  {
    CTABisEntryValid( m_lut, nIndex, &nValid );
  }
  if ( nValid && CTABcopyName( m_lut, nIndex, name, 128 ) == 0 )
  {
    return name;
  }

  return "";
}

void SurfaceAnnotation::GetAnnotationColorAtIndex( int nIndex, int* rgb )
{
  int nValid = 0;
  int nTotalCount = 0;
  CTABgetNumberOfTotalEntries( m_lut, &nTotalCount );
  if ( nIndex < nTotalCount )
  {
    CTABisEntryValid( m_lut, nIndex, &nValid );
  }
  if ( nValid )
  {
    CTABrgbAtIndexi( m_lut, nIndex, rgb, rgb+1, rgb+2 );
  }
}

void SurfaceAnnotation::SetShowOutline(bool bOutline)
{
  this->m_bShowOutline = bOutline;
  m_actorOutline->SetVisibility(bOutline?1:0);
}

vtkActor* SurfaceAnnotation::GetOutlineActor()
{
  return m_actorOutline;
}


vtkPolyData* SurfaceAnnotation::MakeEdgePolyData(const QList<int> &indices_in, const QList<int> &vertices)
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

QList<int> SurfaceAnnotation::DoConnectEdgeVertices(const QList<int> &prev_indices, const QList<int> &vertices)
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
