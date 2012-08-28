#ifndef SURFACESPLINE_H
#define SURFACESPLINE_H

#include <QObject>
#include <QColor>
#include "vtkSmartPointer.h"

extern "C"
{
#include "mri.h"
#include "colortab.h"
}

class vtkActor;
class LayerSurface;
class vtkRenderer;
class vtkPoints;

class SurfaceSpline : public QObject
{
    Q_OBJECT
public:
    explicit SurfaceSpline(LayerSurface *parent = 0);
    virtual ~SurfaceSpline();

    bool Load(const QString& filename);
    bool IsVisible();

    bool IsValid()
    {
      return m_mri;
    }

    bool GetProjection()
    {
      return m_bProjection;
    }

    QColor GetColor()
    {
      return m_color;
    }

    QString GetName()
    {
      return m_strName;
    }

    void AppendProp3D(vtkRenderer* ren);

    void AppendProp2D(vtkRenderer* ren, int nPlane);

signals:
    void SplineChanged();

public slots:
    void SetActiveVertex(int n);
    void SetVisible(bool visible);
    void RebuildActors();
    void SetColor(const QColor& c);
    void SetProjection(bool bProjection );

private:
    void BuildSphereActor(vtkActor* actor, vtkPoints* pts);
    MRI*    m_mri;
    MRI*    m_mriSurf;
    vtkSmartPointer<vtkActor> m_actor;
    vtkSmartPointer<vtkActor> m_actorSpheres;
    vtkSmartPointer<vtkActor> m_actor2D[3];
    vtkSmartPointer<vtkActor> m_actor2DSpheres[3];
    QColor  m_color;
    int     m_nActiveVertex;
    bool    m_bProjection;
    QString m_strName;
};

#endif // SURFACESPLINE_H
