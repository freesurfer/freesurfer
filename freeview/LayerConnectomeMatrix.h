#ifndef LAYERCONNECTOMEMATRIX_H
#define LAYERCONNECTOMEMATRIX_H

#include "Layer.h"
#include <QVariantMap>
#include <QList>
#include <vtkSmartPointer.h>

extern "C"
{
#include "cmat.h"
#include "colortab.h"
}

class LayerMRI;
class LayerPropertyConnectomeMatrix;
class vtkActor;

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

  bool LoadFromFile(const QString& fn_cmat, const QString& fn_parcel);
  bool LoadFromFile();

  void Append2DProps(vtkRenderer *renderer, int nPlane);
  void Append3DProps(vtkRenderer *renderer, bool *bPlaneVisibility);
  bool HasProp(vtkProp *prop);
  bool IsVisible();
  void OnSlicePositionChanged(int nPlane);

  QList<int> GetLabelList();

  QString GetLabelName(int n);
  QColor GetLabelColor(int n);

  bool HasConnection(int i, int j);

signals:

public slots:
  void SetParcelFilename(const QString& filename)
  {
    m_sFilenameParcel = filename;
  }

  void SetColorTable(COLOR_TABLE* ctab);

  void SetFromLabelIndex(int n);
  void SetToLabelIndex(int n);
  void SetToAllLabels(bool bAll);

  void RebuildSplineActors();

private:
  void BuildLabelActors();
  void UpdateLabelActors();

  LayerMRI* m_mriRef;
  LayerMRI* m_mriParcel;
  CMAT* m_cmat;
  COLOR_TABLE* m_ctab;
  QString   m_sFilenameParcel;
  QList<int> m_listLabels;
  int     m_nFromLabelIndex;
  int     m_nToLabelIndex;
  bool    m_bToAllLabels;

  vtkSmartPointer<vtkActor> m_actorSplines;
  QList< vtkSmartPointer<vtkActor> > m_actorLabels;
};

#endif // LAYERCONNECTOMEMATRIX_H
