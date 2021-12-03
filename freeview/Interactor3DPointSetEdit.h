#ifndef INTERACTOR3DPOINTSETEDIT_H
#define INTERACTOR3DPOINTSETEDIT_H

#include "Interactor3D.h"

class Interactor3DPointSetEdit : public Interactor3D
{
public:
  Interactor3DPointSetEdit( QObject* parent );

  virtual bool ProcessMouseDownEvent( QMouseEvent* event, RenderView* view );
//  virtual bool ProcessMouseMoveEvent( QMouseEvent* event, RenderView* view );
  virtual bool ProcessMouseUpEvent( QMouseEvent* event, RenderView* view );
  virtual bool ProcessKeyDownEvent( QKeyEvent* event, RenderView* renderview );

protected:
  bool  m_bEditAttempted;
};
#endif // INTERACTOR3DPOINTSETEDIT_H
