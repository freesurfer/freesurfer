#ifndef __vtkInteractorStyleMyTrackballCamera_h
#define __vtkInteractorStyleMyTrackballCamera_h

#include "vtkInteractorStyleTrackballCamera.h"


class vtkInteractorStyleMyTrackballCamera : public vtkInteractorStyleTrackballCamera
{
public:
  static vtkInteractorStyleMyTrackballCamera *New();
  vtkTypeMacro(vtkInteractorStyleMyTrackballCamera,vtkInteractorStyleTrackballCamera);

  void SetRotateByPoint(bool b, double* dPos = NULL);

  void Rotate() override;

  bool GetRotateByPoint()
  {
    return m_bRotateAroundPoint;
  }

protected:
 vtkInteractorStyleMyTrackballCamera();
 virtual ~vtkInteractorStyleMyTrackballCamera();

private:
  vtkInteractorStyleMyTrackballCamera(const vtkInteractorStyleMyTrackballCamera&);
  void operator=(const vtkInteractorStyleMyTrackballCamera&);

  bool    m_bRotateAroundPoint;
  double  m_dCenterPoint[3];
};

#endif
