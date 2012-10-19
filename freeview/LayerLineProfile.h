#ifndef LAYERLINEPROFILE_H
#define LAYERLINEPROFILE_H

#include "Layer.h"
#include "vtkSmartPointer.h"
#include "LineProf.h"

class vtkActor;
class vtkPolyDataMapper;
class vtkCellArray;
class vtkPoints;
class LayerPointSet;
class LayerMRI;

class LayerLineProfile : public Layer
{
  Q_OBJECT
public:
    LayerLineProfile(int nPlane, QObject* parent = NULL, LayerPointSet* line1 = NULL, LayerPointSet* line2 = NULL);
    ~LayerLineProfile();

    void Append2DProps( vtkRenderer* renderer, int nPlane );
    void Append3DProps( vtkRenderer* renderer, bool* bPlaneVisibility = NULL ) {}

    bool HasProp( vtkProp* prop );

    bool IsVisible();

    void SetVisible( bool bVisible = true );

    LayerPointSet*  GetLine1()
    {
      return m_line1;
    }

    LayerPointSet*  GetLine2()
    {
      return m_line2;
    }

    void SetSourceLayers(LayerPointSet* line1, LayerPointSet* line2);

    int GetPlane()
    {
      return m_nPlane;
    }

    bool Solve(double resolution);

    bool Export(const QString &filename, LayerMRI *mri, int nSample);

public slots:
    void UpdateActors();

protected slots:
    void OnSourceLineDestroyed();

protected:
    void OnSlicePositionChanged(int nPlane);

private:
    std::vector < std::vector < double > > Points3DToSpline2D(std::vector<double> pts3d, double distance);
    std::vector < std::vector < double > > Points2DToSpline3D(std::vector < std::vector<double> > pts2d, int nSample);
    void MakeFlatTube(vtkPoints* points, vtkCellArray* lines, vtkActor* actor_in);

    LayerPointSet*  m_line1;
    LayerPointSet*  m_line2;
    vtkSmartPointer<vtkActor> m_endLines;
    vtkSmartPointer<vtkActor> m_profileLines;
    int         m_nPlane;
    std::vector < std::vector < std::vector < double > > > m_ptsProfile;
    double      m_dSliceLocation;
};

#endif // LAYERLINEPROFILE_H
