#ifndef LAYERFCD_H
#define LAYERFCD_H

#include "LayerVolumeBase.h"
#include "vtkSmartPointer.h"

extern "C"
{
#include "fcd.h"
}

class LayerPropertyFCD;
class LayerMRI;
class vtkImageReslice;
class vtkImageMapToColors;
class vtkImageActor;
class vtkImageData;
class vtkProp;

class LayerFCD : public LayerVolumeBase
{
  Q_OBJECT
public:
  LayerFCD(LayerMRI* mri, QObject* parent = NULL);
  ~LayerFCD();

  bool LoadFromFile();
  bool Load(const QString& subdir, const QString& subject);

  virtual void Append2DProps( vtkRenderer* renderer, int nPlane );
  virtual void Append3DProps( vtkRenderer* renderer, bool* bSliceVisibility = NULL );

  bool HasProp( vtkProp* prop );

  void SetVisible( bool bVisible = true );
  bool IsVisible();

  inline LayerPropertyFCD* GetProperty()
  {
    return (LayerPropertyFCD*)mProperty;
  }

protected slots:
  void UpdateOpacity();
  void UpdateColorMap();
  void Recompute();

protected:
  void InitializeActors();
  void UpdateRASImage(vtkImageData* rasImage);

  virtual void OnSlicePositionChanged( int nPlane );

  FCD_DATA* m_fcd;
  vtkSmartPointer<vtkImageReslice>   mReslice[3];
  vtkSmartPointer<vtkImageMapToColors>  mColorMap[3];

  LayerMRI*  m_layerSource;

  vtkSmartPointer<vtkImageActor>  m_sliceActor2D[3];
  vtkSmartPointer<vtkImageActor>  m_sliceActor3D[3];

  QString   m_sSubjectDir;
  QString   m_sSubject;
};

#endif // LAYERFCD_H
