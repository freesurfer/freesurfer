#ifndef LAYERLINEPROFILE_H
#define LAYERLINEPROFILE_H

#include "Layer.h"
#include "vtkSmartPointer.h"

#if !defined(ARM64) && !defined(DISABLE_LINEPROF)
#include "LineProf.h"
#endif 

class vtkActor;
class vtkPolyDataMapper;
class vtkCellArray;
class vtkPoints;
class LayerPointSet;
class LayerMRI;
class LayerPropertyLineProfile;

class LayerLineProfile : public Layer
{
  Q_OBJECT
public:
  LayerLineProfile(int nPlane, QObject* parent = NULL, LayerPointSet* line1 = NULL, LayerPointSet* line2 = NULL);
  ~LayerLineProfile();

  inline LayerPropertyLineProfile* GetProperty()
  {
    return (LayerPropertyLineProfile*)mProperty;
  }

  void Append2DProps( vtkRenderer* renderer, int nPlane );
  void Append3DProps( vtkRenderer* renderer, bool* bPlaneVisibility = NULL )
  { Q_UNUSED(renderer); Q_UNUSED(bPlaneVisibility); }

  bool HasProp( vtkProp* prop );

  bool IsVisible();

  void SetVisible( bool bVisible = true );

  void SetSourceLayers(LayerPointSet* line1, LayerPointSet* line2);

  int GetPlane()
  {
    return m_nPlane;
  }

  bool Solve(double profileSpacing, double referenceSize, double laplaceResolution, double offset);

  bool Export(const QString &filename, LayerMRI *mri, int nSample);

  bool ExportThickness(const QString& filename, const QList<LayerPointSet*>& layers, int nSample = 100);

  bool Save(const QString& filename);

  static LayerLineProfile* Load(const QString& filename, LayerMRI* ref);

  double GetResultion()
  {
    return m_dResolution;
  }

  double GetSpacing()
  {
    return m_dSpacing;
  }

  int GetNumberOfSamples()
  {
    return m_nSamples;
  }

  double GetOffset()
  {
    return m_dOffset;
  }

  LayerPointSet* GetSpline0()
  {
    return m_spline0;
  }

  LayerPointSet* GetSpline1()
  {
    return m_spline1;
  }

  vtkActor* GetLineProfileActor();

public slots:
  void UpdateActors();
  void SetActiveLineId(int nId);

protected slots:
  void OnSourceLineDestroyed();

  void OnSlicePositionChanged(int nPlane);
  void UpdateOpacity();
  void UpdateColor();
  void UpdateActiveLine();

private:
  std::vector < std::vector < double > > Points3DToSpline2D(std::vector<double> pts3d, double distance);
  std::vector < std::vector < double > > Points2DToSpline3D(std::vector < std::vector<double> > pts2d, int nSample);
  void MakeFlatTube(vtkPoints* points, vtkCellArray* lines, vtkActor* actor_in, double radius);
  QList<double> ComputeThicknessAlongLineProfile(std::vector < std::vector < double > >& line2d,
                                                 const QList<LayerPointSet*>& splines, int nSample);

  LayerPointSet*  m_spline0;
  LayerPointSet*  m_spline1;
  vtkSmartPointer<vtkActor> m_endLines;
  vtkSmartPointer<vtkActor> m_profileLines;
  vtkSmartPointer<vtkActor> m_activeLine;
  int         m_nPlane;
  std::vector < std::vector < std::vector < double > > > m_ptsProfile;
  double      m_dSliceLocation;
  double      m_dResolution;
  double      m_dSpacing;
  double      m_dOffset;
  int         m_nSamples;
  double      m_dReferenceSize;

  int         m_nActiveLineId;
};

#endif // LAYERLINEPROFILE_H
