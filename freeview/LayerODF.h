#ifndef LAYERODF_H
#define LAYERODF_H

#include "LayerMRI.h"
#include "vtkSmartPointer.h"
#include <QList>

//#define USE_ACTOR_LIST

class LayerPropertyODF;
class vtkActor;
class vtkPolyData;
class vtkPoints;
class vtkCellArray;
class vtkFloatArray;

class LayerODF : public LayerMRI
{
  Q_OBJECT
public:
  LayerODF(LayerMRI* layerMRI, QObject* parent = NULL, int nMainview = 0 );
  virtual ~LayerODF();

  bool Load(const QString& fn, const QString& vertex_fn = "", const QString& face_fn = "", bool bPermute = false, bool bHemisphere = false);

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

public slots:
  void UpdateActors(int n = -1);
  void OnColorCodeChanged();
  void OnMainViewChanged(int nView);
  void OnShowInAllChanged();

protected:
  void BuildSlice(int nPlane = -1);
  void ClearSliceActors(int nPlane = -1);
  vtkPolyData* BuildActor(vtkActor* actor2D, vtkActor* actor3D, vtkPoints* pts, vtkCellArray* polys, vtkUnsignedCharArray* scalars, vtkFloatArray* scalars_2);

  vtkSmartPointer<vtkActor> m_actor3D;
  vtkSmartPointer<vtkPolyData> m_polydata[3];
  double m_dScalarRange[2];

  float   m_odfVector[4096][3];
  int     m_odfMesh[7200][3];

  LayerMRI*  m_mask;
  double  m_odfMaskThreshold[2];

  bool    m_bDtkFormat;
  bool    m_bHemisphere;
  int     m_nVectors;
  int     m_nMesh;

  int     m_nMainView;

#ifdef USE_ACTOR_LIST
  QList<vtkActor*> m_listActor2D[3];
  QList<vtkActor*> m_listActor3D[3];
#endif
};

#endif // LAYERODF_H
