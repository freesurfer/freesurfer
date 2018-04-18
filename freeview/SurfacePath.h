#ifndef SURFACEPATH_H
#define SURFACEPATH_H

#include "vtkSmartPointer.h"
#include <QObject>
#include <QColor>
#include <QList>

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

  QList<int> GetPathVerts()
  {
    return m_listVertices;
  }

  void SetUndoVerts(const QList<int>& verts)
  {
    m_undoVertices = verts;
  }

  QList<int> GetUndoVerts()
  {
    return m_undoVertices;
  }

signals:
  void ColorChanged( const QColor& );
  void Progress(int n);
  void Updated();
  void CutLineMade();

public slots:
  bool MakePath(bool bClosed);
  bool MakeCutLine(bool bClosed);

  void Reset();

private:
  void RebuildActor();
  QList<int> DoMakePath(const QList<int>& verts);
  void UpdatePoints();

  vtkSmartPointer<vtkActor>   m_actorOutline;
  vtkSmartPointer<vtkPoints>  m_points;
  QList<int>        m_listVertices;
  QList<int>        m_undoVertices;
  bool  m_bPathMade;
  bool  m_bCutLineMade;
  bool  m_bClosed;

  LayerSurface*   m_mris;
  QColor      m_color;
};

#endif // SURFACEPATH_H
