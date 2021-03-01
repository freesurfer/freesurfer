#ifndef LAYERODF_H
#define LAYERODF_H

#include "LayerMRI.h"
#include "vtkSmartPointer.h"

class LayerPropertyODF;
class vtkActor;

class LayerODF : public LayerMRI
{
  Q_OBJECT
public:
  LayerODF(LayerMRI* layerMRI, QObject* parent = NULL );
  virtual ~LayerODF();

  bool Load(const QString& fn);

  inline LayerPropertyODF* GetProperty()
  {
    return (LayerPropertyODF*)mProperty;
  }

  void Append2DProps( vtkRenderer* renderer, int nPlane );
  void Append3DProps( vtkRenderer* renderer, bool* bSliceVisibility = NULL );

  void OnSlicePositionChanged( int nPlane );
  void SetVisible( bool bVisible = true );
  bool IsVisible();

  void SetMask(LayerMRI* mri);

protected:
  void BuildSlice(int nPlane);
  void UpdateActors();

  vtkSmartPointer<vtkActor> m_actor3D;
  double m_dScalarRange[2];

  float   m_odfVector[181][3];
  int     m_odfMesh[720][3];
  bool    m_bInvert[3];

  LayerMRI*  m_mask;
};

#endif // LAYERODF_H
