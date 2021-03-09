/**
 * @brief Layer class for structural landmarks.
 *
 */
/*
 * Original Author: Ruopeng Wang
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 *
 */

#ifndef LAYERLANDMARKS_H
#define LAYERLANDMARKS_H

#include "LayerEditable.h"
#include <QColor>
#include <QList>
#include <QPointer>
#include "vtkSmartPointer.h"

class vtkActor;
class LayerMRI;

struct Landmark
{
  Landmark();

  double pos[3];
  QColor color;
  bool   valid;
  vtkSmartPointer<vtkActor> actorSphere;
  vtkSmartPointer<vtkActor> actorSlice[3];
};


class LayerLandmarks : public LayerEditable
{
  Q_OBJECT

public:
  LayerLandmarks(QObject* parent = NULL);
  ~LayerLandmarks();

  void Append2DProps( vtkRenderer* renderer, int nPlane );
  void Append3DProps( vtkRenderer* renderer, bool* bPlaneVisibility = NULL );

  bool HasProp( vtkProp* prop );

  bool IsVisible();

  void SetLandmarkPosition(int n, double x, double y, double z);
  void SetLandmarkPosition(int n, double* pos);

  void SetLandmarkColor(int n, const QColor& color);

  Landmark& GetLandmark(int n);

  void SetRadius(double dRadius);

  void SetVisible( bool bVisible = true );

signals:
  void LandmarkChanged(int n, const Landmark& lm);
  void LandmarkAdded();

public slots:
  void SetMRIRef(LayerMRI* mri);

protected:
  void OnSlicePositionChanged(int nPlane);
  void DoTransform(double *mat, int sample_method);
  void DoRestore();

private:
  bool MakeSureLandmarkExist(int n);
  void UpdateActors(bool bBuild3D = true);

  QList<Landmark> m_landmarks;
  QList<Landmark> m_landmarksOriginal;
  double        m_dRadius;
  QPointer<LayerMRI>     m_mriRef;
};

#endif // LAYERLANDMARKS_H
