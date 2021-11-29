#include "SurfacePath.h"
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
#include "FSSurface.h"
#include <QDebug>

SurfacePath::SurfacePath( LayerSurface* owner ) :
  QObject( owner ), m_bPathMade(false), m_bClosed(false),
  m_bCutLineMade(false)
{
  m_actorOutline = vtkSmartPointer<vtkActor>::New();
  m_actorOutline->GetProperty()->SetColor( 0, 1, 0 );
  double ratio = 1;
#if VTK_MAJOR_VERSION > 7
  ratio = MainWindow::GetMainWindow()->devicePixelRatio();
#endif
  m_actorOutline->GetProperty()->SetLineWidth(4*ratio);
  m_actorOutline->GetProperty()->SetPointSize(8*ratio);
  m_points = vtkSmartPointer<vtkPoints>::New();
  m_mris = owner;
  m_actorOutline->SetPosition(m_mris->GetProperty()->GetPosition());
  SetColor(Qt::yellow);
}

SurfacePath::~SurfacePath()
{}


void SurfacePath::Reset()
{
  m_actorOutline->VisibilityOn();
  m_listVertices.clear();
  m_points->Reset();
}

void SurfacePath::RebuildActor()
{
  if (m_bCutLineMade)   // if surface is already cut, do not display anything
  {
    m_actorOutline->SetMapper( vtkSmartPointer<vtkPolyDataMapper>::New() );
    emit Updated();
    return;
  }

  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints( m_points );
  if (m_bPathMade)
  {
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    lines->InsertNextCell( m_points->GetNumberOfPoints() + (m_bClosed?1:0) );
    for ( int i = 0; i < m_points->GetNumberOfPoints(); i++ )
    {
      lines->InsertCellPoint( i );
    }
    if (m_bClosed)
      lines->InsertCellPoint( 0 );
    polydata->SetLines( lines );
  }
  else
  {
    vtkIdType n = 0;
    vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();
    for ( int i = 0; i < m_points->GetNumberOfPoints(); i++ )
    {
      verts->InsertNextCell(1, &n);
      n++;
    }
    polydata->SetVerts(verts);
  }
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION > 5
  mapper->SetInputData( polydata );
#else
  mapper->SetInput( polydata );
#endif
  m_actorOutline->SetMapper( mapper );
  emit Updated();
}

bool SurfacePath::AddPoint( double* pos )
{
  int nVert = m_mris->GetVertexIndexAtTarget(pos, NULL);
  return AddPoint(nVert);
}

bool SurfacePath::AddPoint(int nVert)
{
  if (nVert >= 0 && !m_listVertices.contains(nVert))
  {
    // first find the closest point
//    int n = m_listVertices.size();
//    double d = 1e10;
//    double pt0[3];
//    m_mris->GetSurfaceRASAtVertex(nVert, pt0);
//    for (int i = 0; i < m_listVertices.size(); i++)
//    {
//      double pt[3];
//      m_mris->GetSurfaceRASAtVertex(m_listVertices[n], pt);
//      double val = vtkMath::Distance2BetweenPoints(pt0, pt);
//      if (val < d)
//      {
//        d = val;
//        n = i;
//      }
//    }

//    // then determine where to insert the point
//    int n1, n2;
//    if ( n == 0 )
//    {
//      n1 = 0;
//      n2 = 1;
//    }
//    else
//    {
//      n1 = n-1;
//      n2 = n;
//    }

//    double pt1[3], pt2[3];
//    m_mris->GetSurfaceRASAtVertex(m_listVertices[n1], pt1);
//    m_mris->GetSurfaceRASAtVertex(m_listVertices[n2], pt2);
//    double d1 = vtkMath::Distance2BetweenPoints( pt0, pt1 );
//    double d2 = vtkMath::Distance2BetweenPoints( pt0, pt2 );
//    double d3 = vtkMath::Distance2BetweenPoints( pt1, pt2 );

//    if ( d3 >= d1 && d3 >= d2 )
//    {
//      n = n2;
//    }
//    else if ( d1 < d2 )
//    {
//      n = n1;
//    }
//    else
//    {
//      n = n2 + 1;
//    }

//    if ( n >= m_listVertices.size() )
//    {
//      m_listVertices << nVert;
//    }
//    else
//    {
//      m_listVertices.insert( n, nVert );
//    }

    m_listVertices << nVert;
    Update();
    return true;
  }
  return false;
}

void SurfacePath::Update()
{
  UpdatePoints();
  RebuildActor();
}

bool SurfacePath::RemovePoint(int nVert)
{
  if (nVert >= 0 && m_listVertices.contains(nVert))
  {
    int n = m_listVertices.indexOf(nVert);
    m_listVertices.remove(m_listVertices.indexOf(nVert));
    Update();
    return true;
  }
  return false;
}

void SurfacePath::RemoveLastPoint()
{
  if (!m_listVertices.isEmpty())
  {
    m_listVertices.removeLast();
    Update();
  }
}

bool SurfacePath::RemovePoint( double* pos )
{
  int nVert = m_mris->GetVertexIndexAtTarget(pos, NULL);
  return RemovePoint(nVert);
}

void SurfacePath::Clear()
{
  m_listVertices.clear();
  Update();
}

void SurfacePath::UpdatePoints()
{
  m_points->Reset();
  double* offset = m_mris->GetProperty()->GetPosition();
  for (int n = 0; n < m_listVertices.size(); n++)
  {
    double pt[3];
    m_mris->GetTargetAtVertex(m_listVertices[n], pt);
    for (int i = 0; i < 3; i++)
      pt[i] = pt[i] - offset[i];
    m_points->InsertNextPoint(pt);
  }
}

bool find_path ( MRIS* mris, int* vert_vno, int num_vno, int max_path_length,
                 int* path, int* path_length, bool flag2d = false)
{
  if (num_vno < 2)
    return false;

  int cur_vert_vno;
  int src_vno;
  int dest_vno;
  int vno;
  char* check;
  float* dist;
  int* pred;
  char done;
  VERTEX* v;
  VERTEX_TOPOLOGY* vt;
  VERTEX* u;
  VERTEX_TOPOLOGY* ut;
  float closest_dist;
  int closest_vno;
  int neighbor;
  int neighbor_vno;
  float dist_uv;
  int path_vno;
  int num_path = 0;
  int num_checked;
  float vu_x, vu_y, vu_z;

  dist = (float*) calloc (mris->nvertices, sizeof(float));
  pred = (int*) calloc (mris->nvertices, sizeof(int));
  check = (char*) calloc (mris->nvertices, sizeof(char));
  num_path = 0;
  num_checked = 0;
  (*path_length) = 0;

  for (cur_vert_vno = 0; cur_vert_vno < num_vno-1; cur_vert_vno++)
  {
    /* clear everything */
    for (vno = 0; vno < mris->nvertices; vno++)
    {
      dist[vno] = 999999;
      pred[vno] = -1;
      check[vno] = FALSE;
    }

    /* Set src and dest */
    src_vno = vert_vno[cur_vert_vno+1];
    dest_vno = vert_vno[cur_vert_vno];

    /* make sure both are in range. */
    if (src_vno < 0 || src_vno >= mris->nvertices ||
        dest_vno < 0 || dest_vno >= mris->nvertices)
      continue;

    if (src_vno == dest_vno)
      continue;

    /* pull the src vertex in. */
    dist[src_vno] = 0;
    pred[src_vno] = vno;
    check[src_vno] = TRUE;

    done = FALSE;
    while (!done)
    {
      /* find the vertex with the shortest edge. */
      closest_dist = 999999;
      closest_vno = -1;
      for (vno = 0; vno < mris->nvertices; vno++)
        if (check[vno])
          if (dist[vno] < closest_dist)
          {
            closest_dist = dist[vno];
            closest_vno = vno;
          }
      v = &(mris->vertices[closest_vno]);
      vt = &(mris->vertices_topology[closest_vno]);
      check[closest_vno] = FALSE;

      /* if this is the dest node, we're done. */
      if (closest_vno == dest_vno)
      {
        done = TRUE;
      }
      else
      {
        /* relax its neighbors. */
        for (neighbor = 0; neighbor < vt->vnum; neighbor++)
        {
          neighbor_vno = vt->v[neighbor];
          u = &(mris->vertices[neighbor_vno]);

          /* calc the vector from u to v. */
          vu_x = u->x - v->x;
          vu_y = u->y - v->y;
          if (flag2d)
            vu_z = 0;
          else
            vu_z = u->z - v->z;

          /* recalc the weight. */
          if (flag2d)
            dist_uv = sqrt(((v->x - u->x) * (v->x - u->x)) +
                           ((v->y - u->y) * (v->y - u->y)));
          else
            dist_uv = sqrt(((v->x - u->x) * (v->x - u->x)) +
                           ((v->y - u->y) * (v->y - u->y)) +
                           ((v->z - u->z) * (v->z - u->z)));

          /* if this is a new shortest path, update the predecessor,
             weight, and add it to the list of ones to check next. */
          if (dist_uv + dist[closest_vno] < dist[neighbor_vno])
          {
            pred[neighbor_vno] = closest_vno;
            dist[neighbor_vno] = dist_uv + dist[closest_vno];
            check[neighbor_vno] = TRUE;
          }
        }
      }
      num_checked++;
      if ((num_checked % 100) == 0)
      {

      }
    }

    /* add the predecessors from the dest to the src to the path. */
    path_vno = dest_vno;
    path[(*path_length)++] = dest_vno;
    while (pred[path_vno] != src_vno &&
           (*path_length) < max_path_length )
    {
      path[(*path_length)++] = pred[path_vno];
      path_vno = pred[path_vno];
    }
  }
  free (dist);
  free (pred);
  free (check);

  return true;
}

QVector<int> SurfacePath::DoMakePath(const QVector<int> &verts)
{
  QVector<int> verts_out;

  int nverts = verts.size();
  int* pverts = new int[nverts];
  for (int i = 0; i < nverts; i++)
    pverts[i] = verts[i];
  MRIS* mris = m_mris->GetSourceSurface()->GetMRIS();
  int* path = new int[mris->nvertices];
  int path_length = 0;
  if (find_path(mris, pverts, nverts, mris->nvertices, path, &path_length))
  {
    for (int i = 0; i < path_length; i++)
      verts_out << path[i];
  }
  return verts_out;
}

bool SurfacePath::MakePath(bool bClosed)
{
  if (m_listVertices.size() < 2)
    return false;

  QVector<int> verts = m_listVertices;
  if (bClosed)
    verts << verts.first();

  verts = DoMakePath(verts);
  if (!verts.isEmpty())
  {
    m_listVertices = verts;
    m_bClosed = bClosed;
    m_bPathMade = true;
    Update();
    emit PathMade();
    return true;
  }
  else
    return false;
}

bool SurfacePath::MakeCutLine(bool bClosed)
{
  if (m_listVertices.size() < 2)
    return false;

  QVector<int> verts = m_listVertices;
  if (bClosed)
    verts << verts.first();

  verts = DoMakePath(verts);
  if (!verts.isEmpty())
  {
    m_listVertices = verts;
    m_bClosed = bClosed;
    m_bPathMade = true;
    m_bCutLineMade = true;
    Update();
    emit CutLineMade();
    return true;
  }
  else
    return false;
}

void SurfacePath::AppendProps( vtkRenderer* renderer )
{
  renderer->AddViewProp( m_actorOutline );
}

void SurfacePath::SetColor( const QColor& color )
{
  m_color = color;
  m_actorOutline->GetProperty()->SetColor( color.redF(), color.greenF(), color.blueF() );
  emit ColorChanged( color );
}

QColor SurfacePath::GetColor()
{
  return m_color;
}

void SurfacePath::Show( bool bShow )
{
  m_actorOutline->SetVisibility( bShow?1:0 );
}

vtkActor* SurfacePath::GetActor()
{
  return m_actorOutline;
}

bool SurfacePath::Contains(int nvo)
{
  return m_listVertices.contains(nvo);
}

double SurfacePath::GetLength()
{
    if (m_bPathMade)
    {
        double dist = 0;
        for (int i = 0; i < m_listVertices.size()-1; i++)
        {
            double pt0[3], pt1[3];
            m_mris->GetTargetAtVertex(m_listVertices[i], pt0);
            m_mris->GetTargetAtVertex(m_listVertices[i+1], pt1);
            dist += sqrt(vtkMath::Distance2BetweenPoints(pt0, pt1));
        }
        return dist;
    }
    else
        return 0;
}

bool SurfacePath::SaveAsControlPoints(const QString &filename)
{
  QFile file( filename );
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  {
    QString strg = file.errorString();
    if (strg.isEmpty())
      cerr << "Can not open file for writing\n";
    else
      cerr << qPrintable(strg) << "\n";
    return false;
  }

  QTextStream out(&file);
  for (int n = 0; n < m_listVertices.size(); n++)
  {
    double pt[3];
    m_mris->GetRASAtVertex(m_listVertices[n], pt, m_mris->IsInflated()?FSSurface::SurfaceWhite:-1);
    out << QString("%1 %2 %3\n").arg(pt[0]).arg(pt[1]).arg(pt[2]);
  }
  out << QString("info\nnumpoints %1\nuseRealRAS 1\n").arg( m_listVertices.size() );
  return true;
}

int SurfacePath::GetNumberOfPoints()
{
  return m_listVertices.size();
}
