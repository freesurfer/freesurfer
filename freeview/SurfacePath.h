#ifndef SURFACEPATH_H
#define SURFACEPATH_H

#include "vtkSmartPointer.h"
#include <QObject>
#include <QColor>
#include <QList>
#include <QVector>

class vtkRenderer;
class vtkActor;
class vtkPolyData;
class vtkPoints;
class vtkSelectPolyData;
class vtkBox;
class vtkProp;
class RenderView3D;
class LayerSurface;

class SurfacePath : public QObject
{
  Q_OBJECT
public:
  SurfacePath( LayerSurface* owner );
  virtual ~SurfacePath();

  bool AddPoint( double* pt );

  bool AddPoint( int nvo );

  bool RemovePoint( double* pt);

  bool RemovePoint( int nvo );

  void RemoveLastPoint();

  void Clear();

  QColor GetColor();
  void SetColor( const QColor& color );

  void Update();

  void AppendProps( vtkRenderer* renderer );

  void Show( bool bShow = true );

  vtkActor* GetActor();

  LayerSurface* GetSurface()
  {
    return m_mris;
  }

  bool IsPathMade()
  {
    return m_bPathMade;
  }

  bool IsCutLineMade()
  {
    return m_bCutLineMade;
  }

  bool Contains(int nvo);

  QVector<int> GetPathVerts()
  {
    return m_listVertices;
  }

  void SetUndoVerts(const QVector<int>& verts)
  {
    m_undoVertices = verts;
  }

  QVector<int> GetUndoVerts()
  {
    return m_undoVertices;
  }

  double GetLength();

  int GetNumberOfPoints();

  bool SaveAsControlPoints(const QString& filename);

signals:
  void ColorChanged( const QColor& );
  void Progress(int n);
  void Updated();
  void CutLineMade();
  void PathMade();

public slots:
  bool MakePath(bool bClosed);
  bool MakeCutLine(bool bClosed);

  void Reset();

private:
  void RebuildActor();
  QVector<int> DoMakePath(const QVector<int>& verts);
  void UpdatePoints();

  vtkSmartPointer<vtkActor>   m_actorOutline;
  vtkSmartPointer<vtkPoints>  m_points;
  QVector<int>        m_listVertices;
  QVector<int>        m_undoVertices;
  bool  m_bPathMade;
  bool  m_bCutLineMade;
  bool  m_bClosed;

  LayerSurface*   m_mris;
  QColor      m_color;
};

#endif // SURFACEPATH_H
