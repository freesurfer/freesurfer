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
 */
#ifndef _GenericRenderView_h
#define _GenericRenderView_h
#include "vtkVersion.h"

#if VTK_MAJOR_VERSION > 7
#include "QVTKOpenGLNativeWidget.h"
#else
#include "QVTKWidget.h"
#endif
#include <vtkSmartPointer.h>
#include <QPoint>
#include <QPair>
#include <QList>

class vtkGenericRenderView;
class QKeyEvent;
class QWheelEvent;
class QMouseEvent;
class QString;
class QColor;
class vtkCamera;
class vtkActor;
class vtkProp;
class vtkPropCollection;
class vtkLightKit;
class vtkRenderer;

typedef QPair<QString, double> CameraOperation;
typedef QList<CameraOperation> CameraOperations;

#if VTK_MAJOR_VERSION > 7
class GenericRenderView : public QVTKOpenGLNativeWidget
#else
class GenericRenderView : public QVTKWidget
#endif
{
  Q_OBJECT

public:
  // constructor & deconstructor
  GenericRenderView(QWidget* parent = NULL, Qt::WindowFlags f = 0);
  virtual ~GenericRenderView();

  // call render window to render
  virtual void Render();
  // overload key press event,
  virtual void keyPressEvent(QKeyEvent* event);

  // overloaded mouse press handler
  virtual void mousePressEvent(QMouseEvent* event);
  // overloaded mouse release handler
  virtual void mouseReleaseEvent(QMouseEvent* event);

  virtual void wheelEvent(QWheelEvent* event);

  void RenderSelf();

  QColor GetBackgroundColor();

  int GetAntialiasing();
  int GetStereoRender();

  int GetStereoPairAngle()
  {
    return m_nStereoPairAngle;
  }

  vtkCamera* GetCamera();
  void SetCamera(vtkCamera* camera);

  void Zoom( double dZoomFactor );

  inline vtkRenderer* GetRenderer()
  {
    return m_renderer;
  }

  bool SaveImage(const QString& filename, bool bAntiAliasing, int nMag = 1);

  inline bool IsDualRendererOn()
  {
    return m_renderer2 != NULL;
  }

  inline void EnableRender(bool b)
  {
    m_bEnableRender = b;
  }

  double GetKeyLightIntensity();
  double GetFillLightIntensity();
  double GetHeadLightIntensity();
  double GetBackLightIntensity();

  void GetVisibleProps(vtkPropCollection* propc);
  vtkProp* PickObject(const QPoint& point, vtkPropCollection* propc = NULL, double* pickpos = NULL);

  bool SetCameraOperations(CameraOperations ops);

signals:
  void RenderTriggeredByWheel();
  void MouseReleasedWithoutMove(QMouseEvent*);
  void BackgroundColorChanged(const QColor&);
  void ActorsUpdated();

  void KeyLightIntensityChanged(double);
  void FillLightIntensityChanged(double);
  void HeadLightIntensityChanged(double);
  void BackLightIntensityChanged(double);

public slots:
  void ResetCameraClippingRange();
  void SetBackgroundColor(const QColor& qc);
  void SetAntialiasing(int n, bool redraw = true);
  void SetAntialiasing(bool b, bool redraw = true)
  {
    SetAntialiasing(b?1:0, redraw);
  }
  void CopyToClipboard();
  void EnableInteractor(bool bEnable);

  void SetStereoRender(bool bOn);
  inline void StereoRenderOff()
  {
    SetStereoRender(false);
  }
  void SetStereoTypeToAnaglyph();
  void SetStereoTypeToRedBlue();
  void SetStereoTypeToInterlaced();
  void SetStereoTypeToDresden();
  void SetStereoTypeToCrystalEyes();
  void SetStereoTypeToLeftRight(bool b = true);
  void SetStereoPairAngle(int nAngle);

  void SetKeyLightIntensity(double d, bool redraw = true);
  void SetFillLightIntensity(double d, bool redraw = true);
  void SetHeadLightIntensity(double d, bool redraw = true);
  void SetBackLightIntensity(double d, bool redraw = true);
  void SetLightIntensity(double key, double head, double fill, double back);
  void SetLightToDefault();

  inline void EmitRenderTriggeredByInteractor()
  {
    //emit RenderTriggeredByInteractor();
  }

  inline void Refresh()
  {
    Render();
  }

  // UpdateProps
  virtual void RefreshAllActors(bool bForScreenShot = false);

private slots:
  void UpdateRenderer2();
  void UpdateCamera2();

protected:

  // main renderer
  vtkRenderer*    m_renderer;
  vtkRenderer*    m_renderer2;

private:
  vtkSmartPointer<vtkLightKit>  m_lightKit;
  int     m_nStereoPairAngle;
  QPoint    ptOld;
  bool    m_bEnableRender;
};


#endif
