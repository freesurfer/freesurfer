#ifndef LAYERODF_H
#define LAYERODF_H

#include "LayerMRI.h"
#include "vtkSmartPointer.h"

class LayerPropertyODF;
class vtkActor;
class vtkPolyData;

class LayerODF : public LayerMRI
{
  Q_OBJECT
public:
  LayerODF(LayerMRI* layerMRI, QObject* parent = NULL );
  virtual ~LayerODF();

  bool Load(const QString& fn, const QString& vertex_fn = "", const QString& face_fn = "");

  inline LayerPropertyODF* GetProperty()
  {
    return (LayerPropertyODF*)mProperty;
  }

  void Append2DProps( vtkRenderer* renderer, int nPlane );
  void Append3DProps( vtkRenderer* renderer, bool* bSliceVisibility = NULL );

  void OnSlicePositionChanged( int nPlane );
  void SetVisible( bool bVisible = true );
  bool IsVisible();

  LayerMRI* GetOdfMask()
  {
    return m_mask;
  }

  void SetOdfMask(LayerMRI* mri, bool bRefresh = true);
  void GetOdfMaskThreshold(double* th);
  void SetOdfMaskThreshold(double* th);

signals:
  void UpdateActorRequested(int n = -1);

protected slots:
  void UpdateActors(int n = -1);
  void OnColorCodeChanged();

protected:
  void BuildSlice(int nPlane = -1);

  vtkSmartPointer<vtkActor> m_actor3D;
  vtkSmartPointer<vtkPolyData> m_polydata[3];
  double m_dScalarRange[2];

  float   m_odfVector[4096][3];
  int     m_odfMesh[7200][3];

  LayerMRI*  m_mask;
  double  m_odfMaskThreshold[2];

  bool    m_bDtkFormat;
  int     m_nVectors;
  int     m_nMesh;
};

#endif // LAYERODF_H
