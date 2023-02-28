#ifndef LAYERFCD_H
#define LAYERFCD_H

#include "LayerVolumeBase.h"
#include "vtkSmartPointer.h"
#include <cstddef>
#include "colortab.h"
#include "fcd.h"

class LayerPropertyFCD;
class QThread;
class LayerMRI;
class LayerSurface;
class LayerFCDWorkerThread;
class vtkImageReslice;
class vtkImageMapToColors;
class vtkImageActor;
class vtkImageData;
class vtkProp;

class LayerFCD : public LayerVolumeBase
{
  friend class LayerFCDWorkerThread;

  Q_OBJECT
public:
  LayerFCD(LayerMRI* mri, QObject* parent = NULL);
  ~LayerFCD();

  bool LoadFromFile();
  bool Load(const QString& subdir, const QString& subject, const QString& suffix = "");

  virtual void Append2DProps( vtkRenderer* renderer, int nPlane );
  virtual void Append3DProps( vtkRenderer* renderer, bool* bSliceVisibility = NULL );

  bool HasProp( vtkProp* prop );

  void SetVisible( bool bVisible = true );
  bool IsVisible();

  inline LayerPropertyFCD* GetProperty()
  {
    return (LayerPropertyFCD*)mProperty;
  }

  FCD_DATA* GetFCDData()
  {
    return m_fcd;
  }

  void GetLabelCentroidPosition(int nLabelIndex, double* pos_out);

  QList<bool> GetLabelVisibility()
  {
    return m_labelVisibility;
  }

  void SetLabelVisible(int n, bool visible);

  QList<LayerMRI*> GetMRILayers();
  QList<LayerSurface*> GetSurfaceLayers();

  void SetMRILayerCTAB(COLOR_TABLE* ctab);

  bool IsBusy();

  QThread* GetWorkerThread();

  bool GoToContralateralPoint(double* pos, double* pos_out);

signals:
  void LabelsChanged();
  void StatusChanged();
  void LayerMRICreated(LayerMRI* mri);

public slots:
  void Recompute();
  void SaveFCDLabels(const QString &dir);
  void SetDisplayInNeurologicalView(bool b);

protected slots:
  void UpdateOpacity();
  void UpdateColorMap();
  void OnLayerDestroyed();

protected:
  void DoCompute(bool resetProgress = true);
  void InitializeData();
  void InitializeActors();
  void UpdateRASImage(vtkImageData* rasImage);
  void MakeAllLayers();

  virtual void OnSlicePositionChanged( int nPlane );

  FCD_DATA* m_fcd;
  vtkSmartPointer<vtkImageReslice>   mReslice[3];
  vtkSmartPointer<vtkImageMapToColors>  mColorMap[3];
  QList<bool> m_labelVisibility;

  LayerMRI*   m_layerSource;
  LayerMRI*   m_mri_norm;
  LayerMRI*   m_mri_flair;
  LayerMRI*   m_mri_t2;
  LayerMRI*   m_mri_aseg;
  LayerMRI*   m_mri_difference;
  LayerSurface* m_surf_lh;
  LayerSurface* m_surf_rh;
  LayerSurface* m_surf_lh_pial;
  LayerSurface* m_surf_rh_pial;
  LayerSurface* m_surf_lh_sphere_d1;
  LayerSurface* m_surf_rh_sphere_d1;

  vtkSmartPointer<vtkImageActor>  m_sliceActor2D[3];
  vtkSmartPointer<vtkImageActor>  m_sliceActor3D[3];

  QString   m_sSubjectDir;
  QString   m_sSubject;
  QString   m_sSuffix;

  LayerFCDWorkerThread* m_worker;
};

#endif // LAYERFCD_H
