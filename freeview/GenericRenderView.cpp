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
#include "GenericRenderView.h"
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkCamera.h>
#include <vtkLightKit.h>
#include <vtkImageWriter.h>
#include <vtkBMPWriter.h>
#include <vtkJPEGWriter.h>
#include <vtkPNGWriter.h>
#include <vtkTIFFWriter.h>
#include <vtkPostScriptWriter.h>
#include <vtkRenderLargeImage.h>
#include <vtkVRMLExporter.h>
#include <vtkRenderWindow.h>
#include <vtkPropCollection.h>
#include <vtkPropPicker.h>
#include <vtkAssemblyPath.h>
#include <vtkAssemblyNode.h>
#include <vtkCellPicker.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <QApplication>
#include <QClipboard>
#include <QImage>
#include <QMouseEvent>
#include <QString>
#include <QFileInfo>
#include <QDateTime>
#include <QDebug>

// fix for retina screens
#ifdef Q_OS_OSX
#include "MacRetina.h"
#endif

#define MAX_KEY_LIGHT   1.8
#define MIN_KEY_LIGHT   0.2
#define MIN_RATIO_LIGHT   1.0/100

#define DEFAULT_KEY_LIGHT 0.4
#define DEFAULT_HEAD_LIGHT  0.3
#define DEFAULT_FILL_LIGHT  0.25
#define DEFAULT_BACK_LIGHT  0.2

#if VTK_MAJOR_VERSION > 7
GenericRenderView::GenericRenderView(QWidget* parent, Qt::WindowFlags f) : QVTKOpenGLNativeWidget(parent, f)
#else
GenericRenderView::GenericRenderView(QWidget* parent, Qt::WindowFlags f) : QVTKWidget(parent, f)
#endif
{
#if VTK_MAJOR_VERSION > 7
  setAutoFillBackground(false);
  this->setAttribute(Qt::WA_NoBackground);
#endif

  m_renderer = vtkRenderer::New();
#if VTK_MAJOR_VERSION > 7
  vtkSmartPointer<vtkGenericOpenGLRenderWindow> renWin = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
  SetRenderWindow(renWin);
  setEnableHiDPI(true);
#else
  vtkRenderWindow* renWin = GetRenderWindow();
  // fix for retina screens
#ifdef Q_OS_OSX
  disableGLHiDPI(this->winId());
#endif
#endif
  renWin->AddRenderer(m_renderer);

  m_renderer2 = NULL;
  m_bEnableRender = true;
  m_nStereoPairAngle = -3;

  m_lightKit = vtkSmartPointer<vtkLightKit>::New();
  m_lightKit->AddLightsToRenderer(m_renderer);

  SetLightToDefault();
}

GenericRenderView::~GenericRenderView()
{
  if (m_renderer)
  {
    m_renderer->Delete();
  }
  if (m_renderer2)
  {
    m_renderer2->Delete();
  }
}

void GenericRenderView::SetLightIntensity(double key, double head, double fill, double back)
{
  SetKeyLightIntensity(key, false);
  SetHeadLightIntensity(head, false);
  SetFillLightIntensity(fill, false);
  SetBackLightIntensity(back, false);
  Render();
}

void GenericRenderView::SetLightToDefault()
{
  SetLightIntensity(DEFAULT_KEY_LIGHT, DEFAULT_HEAD_LIGHT, DEFAULT_FILL_LIGHT, DEFAULT_BACK_LIGHT);
}

double GenericRenderView::GetKeyLightIntensity()
{
  return (m_lightKit->GetKeyLightIntensity() - MIN_KEY_LIGHT) / (MAX_KEY_LIGHT - MIN_KEY_LIGHT);
}

void GenericRenderView::SetKeyLightIntensity(double d, bool redraw)
{
  m_lightKit->SetKeyLightIntensity(MIN_KEY_LIGHT + (MAX_KEY_LIGHT-MIN_KEY_LIGHT)*d);
  if (redraw)
  {
    Render();
  }
  emit KeyLightIntensityChanged(d);
}

double GenericRenderView::GetHeadLightIntensity()
{
  return (1.0/m_lightKit->GetKeyToHeadRatio() - MIN_RATIO_LIGHT) / (1 - MIN_RATIO_LIGHT);
}

void GenericRenderView::SetHeadLightIntensity(double d, bool redraw)
{
  m_lightKit->SetKeyToHeadRatio(1.0 / (MIN_RATIO_LIGHT + (1 - MIN_RATIO_LIGHT)*d));
  if (redraw)
  {
    Render();
  }
  emit HeadLightIntensityChanged(d);
}

double GenericRenderView::GetFillLightIntensity()
{
  return (1.0/m_lightKit->GetKeyToFillRatio() - MIN_RATIO_LIGHT) / (1 - MIN_RATIO_LIGHT);
}

void GenericRenderView::SetFillLightIntensity(double d, bool redraw)
{
  m_lightKit->SetKeyToFillRatio(1.0 / (MIN_RATIO_LIGHT + (1 - MIN_RATIO_LIGHT)*d));
  if (redraw)
  {
    Render();
  }
  emit FillLightIntensityChanged(d);
}

double GenericRenderView::GetBackLightIntensity()
{
  return (1.0/m_lightKit->GetKeyToBackRatio() - MIN_RATIO_LIGHT) / (1 - MIN_RATIO_LIGHT);
}

void GenericRenderView::SetBackLightIntensity(double d, bool redraw)
{
  m_lightKit->SetKeyToBackRatio(1.0 / (MIN_RATIO_LIGHT + (1 - MIN_RATIO_LIGHT)*d));
  if (redraw)
  {
    Render();
  }
  emit BackLightIntensityChanged(d);
}

// call render window to render immediately
void GenericRenderView::Render()
{
  if (!m_bEnableRender)
  {
    return;
  }

  //  QCursor old_cursor = cursor();
  //  setCursor(Qt::WaitCursor);

  //  if (GetInteractor()->GetEnabled())
  //    GetInteractor()->Render();
  //  else
  //  qDebug() << "render" << this << QDateTime::currentDateTime().toMSecsSinceEpoch();
  if (isVisible())
    GetRenderWindow()->Render();

  //  setCursor(old_cursor);
}

void GenericRenderView::RenderSelf()
{
  if (!m_bEnableRender)
  {
    return;
  }

  GetRenderWindow()->Render();
}

void GenericRenderView::RefreshAllActors(bool bForScreenshot)
{
  Q_UNUSED(bForScreenshot);
  emit ActorsUpdated();
}

// avoid sending key event to QVTKOpenGLNativeWidget because of a bug in QVTKInteractor
void GenericRenderView::keyPressEvent(QKeyEvent* event)
{
  QWidget::keyPressEvent(event);
  //  QVTKOpenGLNativeWidget::keyPressEvent(event);
}

vtkCamera* GenericRenderView::GetCamera()
{
  return m_renderer->GetActiveCamera();
}

void GenericRenderView::SetCamera(vtkCamera* camera)
{
  m_renderer->SetActiveCamera(camera);
}

QColor GenericRenderView::GetBackgroundColor()
{
  double c[3];
  m_renderer->GetBackground(c);
  QColor qc;
  qc.setRgbF(c[0], c[1], c[2]);
  return qc;
}

void GenericRenderView::SetBackgroundColor(const QColor& qc)
{
  m_renderer->SetBackground(qc.redF(), qc.greenF(), qc.blueF());
  if (m_renderer2)
  {
    m_renderer2->SetBackground(qc.redF(), qc.greenF(), qc.blueF());
  }
  emit BackgroundColorChanged(qc);
}

void GenericRenderView::wheelEvent(QWheelEvent* event)
{
  // remove horizontal scrolling

  QWheelEvent* e = event;
  if (qAbs(e->pixelDelta().x()) > 0)
  {
    QPoint pixelDelta = e->pixelDelta();
    pixelDelta.setX(0);
    QPoint angleDelta = e->angleDelta();
    angleDelta.setX(0);
    e = new QWheelEvent(e->posF(), e->globalPosF(),
                        pixelDelta, angleDelta, e->buttons(), e->modifiers(), e->phase(), e->inverted(), e->source());
  }

#if VTK_MAJOR_VERSION > 7
  QVTKOpenGLNativeWidget::wheelEvent(e);
#else
  QVTKWidget::wheelEvent(e);
#endif
  emit RenderTriggeredByWheel();
}

void GenericRenderView::mousePressEvent(QMouseEvent* event)
{
  ptOld = event->pos();
#if VTK_MAJOR_VERSION > 7
  QVTKOpenGLNativeWidget::mousePressEvent(event);
#else
  QVTKWidget::mousePressEvent(event);
#endif
}

void GenericRenderView::mouseReleaseEvent(QMouseEvent* event)
{
  if (ptOld == event->pos())
  {
    emit MouseReleasedWithoutMove(event);
  }
#if VTK_MAJOR_VERSION > 7
  QVTKOpenGLNativeWidget::mouseReleaseEvent(event);
#else
  QVTKWidget::mouseReleaseEvent(event);
#endif
}

bool GenericRenderView::SaveImage(const QString& filename, bool bAntiAliasing, int nMag)
{
  QFileInfo fi(filename);
  QString ext = fi.suffix().toLower();
  vtkImageWriter* writer = 0;
  QString fn = filename;
  if (ext == "wrl")
  {
    vtkVRMLExporter* exporter = vtkVRMLExporter::New();
    exporter->SetFileName(fn.toLatin1().data());
    exporter->SetRenderWindow(GetRenderWindow());
    exporter->Write();
    exporter->Delete();
  }
  else if (ext == "jpg"|| ext == "jpeg")
  {
    writer = vtkJPEGWriter::New();
  }
  else if (ext == "bmp")
  {
    writer = vtkBMPWriter::New();
  }
  else if (ext == "ps")
  {
    writer = vtkPostScriptWriter::New();
  }
  else if (ext == "tif" || ext == "tiff")
  {
    writer = vtkTIFFWriter::New();
  }
  else
  {
    writer = vtkPNGWriter::New();
    if (ext != "png")
    {
      fn += ".png";
    }
  }
  if (writer)
  {
    bool bCurrentAA = GetAntialiasing() > 0;
    SetAntialiasing( (bAntiAliasing?true:bCurrentAA), false);
    vtkRenderLargeImage* image = vtkRenderLargeImage::New();
    image->SetInput(m_renderer);
    image->SetMagnification(nMag);
#if VTK_MAJOR_VERSION > 5
    writer->SetInputConnection(image->GetOutputPort());
#else
    writer->SetInput(image->GetOutput());
#endif
    writer->SetFileName(fn.toLatin1().data());
    writer->Write();
    image->Delete();
    writer->Delete();
    SetAntialiasing(bCurrentAA, false);
  }
  return true;
}

int GenericRenderView::GetAntialiasing()
{
#if VTK_MAJOR_VERSION > 5
  return GetRenderWindow()->GetMultiSamples() > 0 ? 1: 0;
#else
  return GetRenderWindow()->GetAAFrames() > 0 ? 1: 0;
#endif
}

void GenericRenderView::SetAntialiasing(int bSet, bool redraw)
{
#if VTK_MAJOR_VERSION > 5
  GetRenderWindow()->SetMultiSamples(bSet? 8 : 0);
#else
  GetRenderWindow()->SetAAFrames(bSet > 0 ? 6 : 0);
#endif
  if (redraw)
  {
    Render();
  }
}

void GenericRenderView::CopyToClipboard()
{
  QClipboard* clipboard = QApplication::clipboard();
  unsigned char* p = this->GetRenderWindow()->GetRGBACharPixelData(0, 0, this->width()-1, this->height()-1, 0);
  QImage image(p, width(), height(), QImage::Format_RGB32);

  unsigned char ch[2] = {0, 1};
  unsigned short* a = (unsigned short*)ch;
  if (*a == 1)  // Big Endian
  {
    int nsize = width()*height();
    unsigned char* ptr = p;
    for (int i = 0; i < nsize; i++)
    {
      qSwap(*ptr, *(ptr+3));
      qSwap(*(ptr+1), *(ptr+2));
      ptr+=4;
    }
  }
  clipboard->setImage(image.mirrored().rgbSwapped(), QClipboard::Clipboard);
}

void GenericRenderView::EnableInteractor(bool b)
{
  vtkRenderWindowInteractor* iren = GetRenderWindow()->GetInteractor();
  if (iren)
  {
    if (b)
    {
      iren->Enable();
    }
    else
    {
      iren->Disable();
    }
  }
}

int GenericRenderView::GetStereoRender()
{
  return GetRenderWindow()->GetStereoRender();
}

void GenericRenderView::SetStereoRender(bool bOn)
{
  GetRenderWindow()->SetStereoRender(bOn?1:0);
  //#ifdef Q_OS_MAC
  Render();
  //#endif
}

void GenericRenderView::SetStereoTypeToAnaglyph()
{
  vtkRenderWindow* wnd = GetRenderWindow();
  wnd->SetStereoTypeToAnaglyph();
  wnd->SetAnaglyphColorSaturation(0.6);
  SetStereoRender(true);
}

void GenericRenderView::SetStereoTypeToRedBlue()
{
  vtkRenderWindow* wnd = GetRenderWindow();
  wnd->SetStereoTypeToRedBlue();
  SetStereoRender(true);
}

void GenericRenderView::SetStereoTypeToInterlaced()
{
  vtkRenderWindow* wnd = GetRenderWindow();
  wnd->SetStereoTypeToInterlaced();
  SetStereoRender(true);
}

void GenericRenderView::SetStereoTypeToDresden()
{
  vtkRenderWindow* wnd = GetRenderWindow();
  wnd->SetStereoTypeToDresden();
  SetStereoRender(true);
}

void GenericRenderView::SetStereoTypeToCrystalEyes()
{
  vtkRenderWindow* wnd = GetRenderWindow();
  wnd->SetStereoTypeToCrystalEyes();
  SetStereoRender(true);
}

void GenericRenderView::SetStereoTypeToLeftRight(bool b)
{
  if (b)
  {
    if (!m_renderer2)
    {
      m_renderer2 = vtkRenderer::New();
      m_renderer2->SetLayer(0);
      m_lightKit->AddLightsToRenderer(m_renderer2);
    }
    m_renderer->SetViewport(0, 0, 0.5, 1);
    m_renderer2->SetViewport(0.5, 0, 1, 1);
    m_renderer2->InteractiveOff();
    m_renderer2->SetBackground(m_renderer->GetBackground());
    GetRenderWindow()->AddRenderer(m_renderer2);
    UpdateRenderer2();
    GetRenderWindow()->SetStereoTypeToRight();
    SetStereoRender(true);
    connect(this, SIGNAL(PropsUpdated()), this, SLOT(UpdateRenderer2()));
    connect(this, SIGNAL(RenderTriggeredByInteractor()), this, SLOT(UpdateCamera2()));
  }
  else
  {
    disconnect(this, SIGNAL(PropsUpdated()), this, SLOT(UpdateRenderer2()));
    disconnect(this, SIGNAL(RenderTriggeredByInteractor()), this, SLOT(UpdateCamera2()));
    if (m_renderer2)
    {
      GetRenderWindow()->RemoveRenderer(m_renderer2);
      GetRenderWindow()->SetStereoTypeToInterlaced(); // call to avoid a vtk bug
      m_renderer2->Delete();
      m_renderer2 = NULL;
    }
    m_renderer->SetViewport(0, 0, 1, 1);
    Render();
  }
}

void GenericRenderView::UpdateRenderer2()
{
  if (!m_renderer2)
  {
    return;
  }

  vtkPropCollection* props = m_renderer->GetViewProps(),
      * props2 = m_renderer2->GetViewProps();
  props2->RemoveAllItems();

  props->InitTraversal();
  vtkProp* prop = props->GetNextProp();
  while (prop)
  {
    props2->AddItem(prop);
    prop = props->GetNextProp();
  }
  UpdateCamera2();
}

void GenericRenderView::UpdateCamera2()
{
  vtkCamera* cam = m_renderer2->GetActiveCamera();
  //  MyVTKUtils::CopyCamera(m_renderer->GetActiveCamera(), cam);
  cam->Azimuth(m_nStereoPairAngle);
  m_renderer2->ResetCameraClippingRange();
  m_renderer2->Render();
}

void GenericRenderView::SetStereoPairAngle(int nAngle)
{
  if (m_nStereoPairAngle == nAngle)
  {
    return;
  }

  m_nStereoPairAngle = nAngle;
  if (m_renderer2)
  {
    UpdateCamera2();
    Render();
  }
}

vtkProp* GenericRenderView::PickObject(const QPoint& point, vtkPropCollection* propc, double* pickpos)
{
  if (GetRenderWindow()->GetStereoRender())
  {
    return NULL;
  }

  vtkPropCollection* p = propc;
  if (!p)
  {
    p = m_renderer->GetViewProps();
  }
  vtkPropPicker* picker = vtkPropPicker::New();
  QRect rc = rect();
  vtkProp* actor = NULL;
  if (picker->PickProp(point.x(), rc.height() - point.y(), m_renderer, p))
  {
    actor = picker->GetViewProp();
    if (pickpos)
    {
      picker->GetPickPosition(pickpos);
    }
  }

  picker->Delete();
  return actor;
}

void GenericRenderView::GetVisibleProps(vtkPropCollection* propc)
{
  vtkPropCollection* props = m_renderer->GetViewProps();
  props->InitTraversal();
  vtkProp* prop = props->GetNextProp();
  while (prop)
  {
    if (prop->GetVisibility())
    {
      propc->AddItem(prop);
    }
    prop = props->GetNextProp();
  }
}

void GenericRenderView::Zoom( double dZoomFactor )
{
  vtkCamera* cam = GetCamera();
  cam->Dolly( dZoomFactor );
}

bool GenericRenderView::SetCameraOperations(CameraOperations ops)
{
  vtkCamera* cam = GetCamera();
  for (int i = 0; i < ops.size(); i++)
  {
    if (ops[i].first.toLower() == "azimuth")
    {
      cam->Azimuth(ops[i].second);
    }
    else if (ops[i].first.toLower() == "dolly" || ops[i].first.toLower() == "zoom")
    {
      cam->Dolly(ops[i].second);
    }
    else if (ops[i].first.toLower() == "elevation" || ops[i].first.toLower() == "elevate")
    {
      cam->Elevation(ops[i].second);
    }
    else if (ops[i].first.toLower() == "roll" || ops[i].first.toLower() == "elevate")
    {
      cam->Roll(ops[i].second);
    }
    else if (ops[i].first.toLower() == "yaw" )
    {
      cam->Yaw(ops[i].second);
    }
    else
    {
      cerr << "Unrecognized camera operation: " << qPrintable(ops[i].first) << "\n";
      return false;
    }

    cam->OrthogonalizeViewUp();
  }
  return true;
}

void GenericRenderView::ResetCameraClippingRange()
{
  GetRenderer()->ResetCameraClippingRange();
}
