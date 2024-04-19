/**
 * @brief Base Interactor class manage mouse and key input in render view.
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

#ifndef Interactor_h
#define Interactor_h

#include <QObject>
#include <QMouseEvent>
#include <QKeyEvent>

class vtkRenderer;
class RenderView;
class QEvent;
class QWheelEvent;
class QMouseEvent;
class QKeyEvent;

class Interactor : public QObject
{
  Q_OBJECT
public:
  Interactor( QObject* parent = NULL );
  virtual ~Interactor();

  enum MeasureMode
  {
    MM_Line = 0, MM_Polyline, MM_Spline, MM_Rectangle, MM_Label, MM_SurfaceRegion, MM_DrawOnSurface
  };

  enum EditMode
  {
    EM_Freehand = 0, EM_Fill, EM_Polyline, EM_Livewire, EM_ColorPicker, EM_Contour, EM_Shift, EM_Clone, EM_GeoSeg, EM_ScribblePrompt
  };

  int GetAction();
  void SetAction( int nAction );

  // return true if to have parent interactor continue processing the event
  // return false to stop event from further processing
  virtual bool ProcessMouseDownEvent    ( QMouseEvent* event, RenderView* view );
  virtual bool ProcessMouseUpEvent      ( QMouseEvent* event, RenderView* view );
  virtual bool ProcessMouseMoveEvent    ( QMouseEvent* event, RenderView* view );
  virtual bool ProcessMouseWheelEvent   ( QWheelEvent* event, RenderView* view );
  virtual bool ProcessMouseEnterEvent   ( QEvent* event, RenderView* view );
  virtual bool ProcessMouseLeaveEvent   ( QEvent* event, RenderView* view );
  virtual bool ProcessKeyDownEvent      ( QKeyEvent* event, RenderView* view );
  virtual bool ProcessKeyUpEvent        ( QKeyEvent* event, RenderView* view );
  virtual void ProcessPostMouseMoveEvent( QMouseEvent* event, RenderView* view )
  { Q_UNUSED(event); Q_UNUSED(view); }
  virtual void ProcessPostMouseWheelEvent( QWheelEvent* event, RenderView* view )
  { Q_UNUSED(event); Q_UNUSED(view); }

  static void SetUseCommandControl(bool b);

protected:
  virtual void UpdateCursor( QEvent* event, QWidget* wnd );

  int  m_nDownPosX;
  int  m_nDownPosY;

  int  m_nAction;
  static Qt::KeyboardModifier CONTROL_MODIFIER;
  static Qt::Key CONTROL_KEY;
};

#endif


