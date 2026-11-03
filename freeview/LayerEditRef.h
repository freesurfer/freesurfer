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

#ifndef LayerEditRef_H
#define LayerEditRef_H

#include "LayerEditable.h"
#include <QColor>
#include <QList>
#include <QPointer>
#include <QVector>
#include "vtkSmartPointer.h"
#include "vtkPoints.h"

class vtkActor;
class LayerMRI;

class LayerEditRef : public LayerEditable
{
  Q_OBJECT

public:
  LayerEditRef(QObject* parent = NULL);
  ~LayerEditRef();

  void Append2DProps( vtkRenderer* renderer, int nPlane );
  void Append3DProps( vtkRenderer* renderer, bool* bPlaneVisibility = NULL );

  bool HasProp( vtkProp* prop );

  bool IsVisible();

  void SetStartPosition(double* pos);
  void SetEndPosition(double* pos);


  int GetNumberOfMarks();

  LayerMRI* GetMRIRef()
  {
    return m_mriRef;
  }


signals:

public slots:
  void SetColor(const QColor& color);
  void SetVisible( bool bVisible = true );
  void SetMRIRef(LayerMRI* mri);
  void Reset();
  void ApplyToMRI(LayerMRI* mri_in = NULL);

protected:
  void OnSlicePositionChanged(int nPlane);

private:
  void UpdateActors(bool bBuild3D = true);

  vtkSmartPointer<vtkPoints>  m_points;
  QPointer<LayerMRI>     m_mriRef;
  vtkSmartPointer<vtkActor> m_actorSlice[3];
  vtkSmartPointer<vtkActor> m_actor;

  struct POINT {
    int n[3];
    friend bool operator==(const POINT& lhs, const POINT& rhs)
    {
      return (lhs.n[0] == rhs.n[0] && lhs.n[1] == rhs.n[1] && lhs.n[2] == rhs.n[2]);
    }
  };

  QVector<POINT> m_voxels;
  QColor  m_color;
};

#endif // LayerEditRef_H
