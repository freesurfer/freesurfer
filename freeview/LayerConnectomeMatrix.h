#ifndef LAYERCONNECTOMEMATRIX_H
#define LAYERCONNECTOMEMATRIX_H

#include "Layer.h"
extern "C"
{
#include "cmat.h"
#include "colortab.h"
}

class LayerMRI;
class LayerPropertyConnectomeMatrix;

class LayerConnectomeMatrix : public Layer
{
  Q_OBJECT
public:
  explicit LayerConnectomeMatrix(LayerMRI* ref, QObject *parent = 0);
  virtual ~LayerConnectomeMatrix();

  inline LayerPropertyConnectomeMatrix* GetProperty()
  {
    return (LayerPropertyConnectomeMatrix*)mProperty;
  }

  bool LoadFromFile(const QString& fn_cmat, const QString& fn_parcel, const QString& fn_ctab);
  bool LoadFromFile();

  void Append2DProps(vtkRenderer *renderer, int nPlane);
  void Append3DProps(vtkRenderer *renderer, bool *bPlaneVisibility);
  bool HasProp(vtkProp *prop);
  bool IsVisible();
  void OnSlicePositionChanged(int nPlane);

signals:

public slots:
  void SetParcelFilename(const QString& filename)
  {
    m_sFilenameParcel = filename;
  }

  void SetColorTableFilename(const QString& fn)
  {
    m_sFilenameCTAB = fn;
  }

private:
  LayerMRI* m_mriRef;
  LayerMRI* m_mriParcel;
  CMAT* m_cmat;
  COLOR_TABLE* m_ctab;
  QString   m_sFilenameParcel;
  QString   m_sFilenameCTAB;
};

#endif // LAYERCONNECTOMEMATRIX_H
