#ifndef REGION3D_H
#define REGION3D_H

#include <QObject>
#include <QColor>
#include "vtkSmartPointer.h"

class vtkRenderer;
class vtkActor;
class vtkPolyData;
class vtkPoints;
class vtkSelectPolyData;
class vtkBox;
class vtkProp;
class vtkClipPolyData;
class vtkCleanPolyData;
class RenderView3D;
class LayerMRI;

class Region3D : public QObject
{
  Q_OBJECT
public:
  explicit Region3D(LayerMRI* owner);
  virtual ~Region3D();

  void AddPoint( double* pt );

  bool Close();

  void ResetOutline();

  QColor GetColor();
  void SetColor( const QColor& color );

  void Update();

  void AppendProps( vtkRenderer* renderer );

  void Show( bool bShow = true );

  LayerMRI* GetMRI()
  {
    return m_mri;
  }

  void Highlight( bool bHighlight );

  bool HasPoint(double *pt, double dist = 0.5);

  static bool WriteHeader( FILE* fp, LayerMRI* mri, int nNum );

  bool WriteBody( FILE* fp );

  bool Load(FILE* fp);

signals:
  void ColorChanged( const QColor& );

private:
  void RebuildOutline(bool bInterpolate);

  vtkSmartPointer<vtkActor>   m_actor;
  vtkSmartPointer<vtkPoints>  m_points;

  LayerMRI*   m_mri;
  QColor      m_color;
};

#endif // REGION3D_H
