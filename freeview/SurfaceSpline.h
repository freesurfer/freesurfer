#ifndef SURFACESPLINE_H
#define SURFACESPLINE_H

#include <QObject>
#include "vtkSmartPointer.h"

extern "C"
{
#include "mri.h"
#include "colortab.h"
}

class vtkActor;
class LayerSurface;

class SurfaceSpline : public QObject
{
    Q_OBJECT
public:
    explicit SurfaceSpline(LayerSurface *parent = 0);
    virtual ~SurfaceSpline();

    vtkActor* GetActor();
    vtkActor* GetActor2D(int nPlane);

    bool Load(const QString& filename);

signals:

public slots:
    void SetActiveVertex(int n);
    void SetVisible(bool visible);
    void RebuildActors();

private:
    MRI*    m_mri;
    MRI*    m_mriSurf;
    vtkSmartPointer<vtkActor> m_actor;
    vtkSmartPointer<vtkActor> m_actor2D[3];
    int     m_nActiveVertex;
};

#endif // SURFACESPLINE_H
