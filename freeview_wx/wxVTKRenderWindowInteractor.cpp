/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: wxVTKRenderWindowInteractor.cpp,v $
  Language:  C++
  Date:      $Date: 2011/03/11 23:27:43 $
  Version:   $Revision: 1.1 $

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include <assert.h>

#include "wxVTKRenderWindowInteractor.h"

// DRG - added to workaround build errors on Cocoa (need to define __WXCOCOA__)
#ifdef __WXOSX_COCOA__
#ifndef __WXCOCOA__
#define __WXCOCOA__
#endif
#endif

//This is needed for vtk 3.1 :
#ifndef VTK_MAJOR_VERSION
#  include "vtkVersion.h"
#endif

#if VTK_MAJOR_VERSION > 4 || (VTK_MAJOR_VERSION == 4 && VTK_MINOR_VERSION > 0)
#  include "vtkCommand.h"
#else
#  include "vtkInteractorStyle.h"
#endif
#include "vtkDebugLeaks.h"

#ifdef __WXMAC__
#ifdef __WXCOCOA__
#include "vtkCocoaRenderWindow.h"
#else
#include "vtkCarbonRenderWindow.h"
#endif
#endif


// function to get VTK keysyms from ascii characters
static const char* ascii_to_key_sym(int);
// function to get VTK keysyms from Qt keys
static const char* wx_key_to_key_sym(int);



//Keep this for compatibilty reason, this was introduced in wxGTK 2.4.0
#if (!wxCHECK_VERSION(2, 4, 0))
wxWindow* wxGetTopLevelParent(wxWindow *win)
{
    while ( win && !win->IsTopLevel() )
         win = win->GetParent();
    return win;
}
#endif

// To access objc calls on cocoa
#ifdef __WXCOCOA__
#ifdef VTK_USE_COCOA
// DRG - Commented this out, it is not necessary on Cocoa (DRG - 9/14/10)
//#import <Cocoa/Cocoa.h>
// This trick is no longer need in VTK CVS, should get rid of that:
#define id Id
#else
#error Build mismatch you need both wxWidgets and VTK to be configure against Cocoa to work
#endif //VTK_USE_COCOA
#endif //__WXCOCOA__

#ifdef __WXGTK__
#    include <gdk/gdkx.h> // GDK_WINDOW_XWINDOW is found here in wxWidgets 2.8.0
#    include "gdk/gdkprivate.h"
#if wxCHECK_VERSION(2, 8, 0)
#ifdef __WXGTK20__
#if  wxCHECK_VERSION(2, 9, 0)
#include <wx/gtk/private/win_gtk.h>
#else
#include <wx/gtk/win_gtk.h>
#endif
#else
#include <wx/gtk1/win_gtk.h>
#endif
#else
#include <wx/gtk/win_gtk.h>
#endif
#if  wxCHECK_VERSION(2, 9, 0)
#define piz(wxwin) WX_PIZZA((wxwin)->m_wxwindow)
#define GetXWindow(wxwin) (wxwin)->m_wxwindow ? \
GDK_WINDOW_XWINDOW(((GtkWidget*)piz(wxwin))->window) : \
                        GDK_WINDOW_XWINDOW((wxwin)->m_widget->window)
#else
#define GetXWindow(wxwin) (wxwin)->m_wxwindow ? \
GDK_WINDOW_XWINDOW(GTK_PIZZA((wxwin)->m_wxwindow)->bin_window) : \
                         GDK_WINDOW_XWINDOW((wxwin)->m_widget->window)
#endif
#endif

#ifdef __WXX11__
#include "wx/x11/privx.h"
#define GetXWindow(wxwin)   ((Window)(wxwin)->GetHandle())
#endif


//For more info on this class please go to:
//http://wxvtk.sf.net
//This hack is for some buggy wxGTK version:
#if wxCHECK_VERSION(2, 3, 2) && !wxCHECK_VERSION(2, 4, 1) && defined(__WXGTK__)
#  define WX_USE_X_CAPTURE 0
#else
#  define WX_USE_X_CAPTURE 1
#endif

#define ID_wxVTKRenderWindowInteractor_TIMER 1001

#if defined(__WXGTK__) && defined(wxUSE_GLCANVAS)
IMPLEMENT_DYNAMIC_CLASS(wxVTKRenderWindowInteractor, wxGLCanvas)
#else
IMPLEMENT_DYNAMIC_CLASS(wxVTKRenderWindowInteractor, wxWindow)
#endif  //__WXGTK__

//---------------------------------------------------------------------------
#if defined(__WXGTK__) && defined(wxUSE_GLCANVAS)
BEGIN_EVENT_TABLE(wxVTKRenderWindowInteractor, wxGLCanvas)
#else
BEGIN_EVENT_TABLE(wxVTKRenderWindowInteractor, wxWindow)
#endif //__WXGTK__
  //refresh window by doing a Render
  EVT_PAINT       (wxVTKRenderWindowInteractor::OnPaint)
  EVT_ERASE_BACKGROUND(wxVTKRenderWindowInteractor::OnEraseBackground)
  EVT_MOTION      (wxVTKRenderWindowInteractor::OnMotion)

  //Bind the events to the event converters
  EVT_LEFT_DOWN   (wxVTKRenderWindowInteractor::OnButtonDown)
  EVT_LEFT_DCLICK (wxVTKRenderWindowInteractor::OnMouseDoubleCilck)
  EVT_MIDDLE_DOWN (wxVTKRenderWindowInteractor::OnButtonDown)
  EVT_RIGHT_DOWN  (wxVTKRenderWindowInteractor::OnButtonDown)
  EVT_LEFT_UP     (wxVTKRenderWindowInteractor::OnButtonUp)
  EVT_MIDDLE_UP   (wxVTKRenderWindowInteractor::OnButtonUp)
  EVT_RIGHT_UP    (wxVTKRenderWindowInteractor::OnButtonUp)
#if !(VTK_MAJOR_VERSION == 3 && VTK_MINOR_VERSION == 1)
  EVT_ENTER_WINDOW(wxVTKRenderWindowInteractor::OnEnter)
  EVT_LEAVE_WINDOW(wxVTKRenderWindowInteractor::OnLeave)
  EVT_MOUSEWHEEL  (wxVTKRenderWindowInteractor::OnMouseWheel)
#if wxCHECK_VERSION(2, 8, 0)
  EVT_MOUSE_CAPTURE_LOST(wxVTKRenderWindowInteractor::OnMouseCaptureLost)
#endif
  EVT_KEY_DOWN    (wxVTKRenderWindowInteractor::OnKeyDown)
  EVT_KEY_UP      (wxVTKRenderWindowInteractor::OnKeyUp)
  EVT_CHAR        (wxVTKRenderWindowInteractor::OnChar)
#endif
  EVT_TIMER       (ID_wxVTKRenderWindowInteractor_TIMER, wxVTKRenderWindowInteractor::OnTimer)
  EVT_SIZE        (wxVTKRenderWindowInteractor::OnSize)
END_EVENT_TABLE()

vtkCxxRevisionMacro(wxVTKRenderWindowInteractor, "$Revision: 1.1 $")
vtkInstantiatorNewMacro(wxVTKRenderWindowInteractor)

#if defined(__WXGTK__) && defined(wxUSE_GLCANVAS)
static int wxvtk_attributes[]={
  WX_GL_DOUBLEBUFFER,
  WX_GL_RGBA,
  WX_GL_DEPTH_SIZE,
  16,
  0
};
#endif

//---------------------------------------------------------------------------
wxVTKRenderWindowInteractor::wxVTKRenderWindowInteractor()
#if defined(__WXGTK__) && defined(wxUSE_GLCANVAS)
  #if  wxCHECK_VERSION(2, 9, 0)
    : wxGLCanvas(0, -1, wxvtk_attributes, wxDefaultPosition, wxDefaultSize, 0, wxT("wxVTKRenderWindowInteractor")), 
  #else
    : wxGLCanvas(0, -1, wxDefaultPosition, wxDefaultSize, 0, wxT("wxVTKRenderWindowInteractor"), wxvtk_attributes),
  #endif
      vtkRenderWindowInteractor()
#else
    : wxWindow(), vtkRenderWindowInteractor()
#endif //__WXGTK__
      , timer(this, ID_wxVTKRenderWindowInteractor_TIMER)
      , ActiveButton(wxEVT_NULL)
      , Stereo(0)
      , Handle(0)
      , Created(true)
      , RenderWhenDisabled(1)
      , UseCaptureMouse(0)
{
#ifdef VTK_DEBUG_LEAKS
  vtkDebugLeaks::ConstructClass("wxVTKRenderWindowInteractor");
#endif
#if defined(__WXMSW__) || defined(__WXMAC__)
  this->SetUseCaptureMouse(1);
#endif
#if defined(__WXGTK__) && defined(wxUSE_GLCANVAS) && wxCHECK_VERSION(2, 9, 0)
  this->GLContext = new wxGLContext (this);
#endif
  this->RenderWindow = NULL;
  this->SetRenderWindow(vtkRenderWindow::New());
  this->RenderWindow->Delete();
}

//---------------------------------------------------------------------------
wxVTKRenderWindowInteractor::wxVTKRenderWindowInteractor(wxWindow *parent,
                                                         wxWindowID id,
                                                         const wxPoint &pos,
                                                         const wxSize &size,
                                                         long style,
                                                         const wxString &name)
#if defined(__WXGTK__) && defined(wxUSE_GLCANVAS)
  #if  wxCHECK_VERSION(2, 9, 0)
    : wxGLCanvas(parent, id, wxvtk_attributes, pos, size, style, name),
  #else
    : wxGLCanvas(parent, id, pos, size, style, name, wxvtk_attributes),
  #endif
      vtkRenderWindowInteractor()
#else
      : wxWindow(parent, id, pos, size, style, name), vtkRenderWindowInteractor()
#endif //__WXGTK__
      , timer(this, ID_wxVTKRenderWindowInteractor_TIMER)
      , ActiveButton(wxEVT_NULL)
      , Stereo(0)
      , Handle(0)
      , Created(true)
      , RenderWhenDisabled(1)
      , UseCaptureMouse(0)
{
#ifdef VTK_DEBUG_LEAKS
  vtkDebugLeaks::ConstructClass("wxVTKRenderWindowInteractor");
#endif
#if defined(__WXMSW__) || defined(__WXMAC__)
  this->SetUseCaptureMouse(1);
#endif
#if defined(__WXGTK__) && defined(wxUSE_GLCANVAS) && wxCHECK_VERSION(2, 9, 0)
  this->GLContext = new wxGLContext (this);
#endif
  this->RenderWindow = NULL;
  this->SetRenderWindow(vtkRenderWindow::New());
  this->RenderWindow->Delete();
#ifdef __WXMAC__
  // On Mac (Carbon) we don't get notified of the initial window size with an EVT_SIZE event,
  // so we update the size information of the interactor/renderwindow here
  this->UpdateSize(size.x, size.y);
#endif
}
//---------------------------------------------------------------------------
wxVTKRenderWindowInteractor::~wxVTKRenderWindowInteractor()
{
#if defined(__WXGTK__) && defined(wxUSE_GLCANVAS) && wxCHECK_VERSION(2, 9, 0)
  delete this->GLContext;
#endif
  SetRenderWindow(NULL);
  SetInteractorStyle(NULL);
}
//---------------------------------------------------------------------------
wxVTKRenderWindowInteractor * wxVTKRenderWindowInteractor::New()
{
  // we don't make use of the objectfactory, because we're not registered
  return new wxVTKRenderWindowInteractor;
}
//---------------------------------------------------------------------------
void wxVTKRenderWindowInteractor::Initialize()
{
  int *size = RenderWindow->GetSize();
  // enable everything and start rendering
  Enable();
  //RenderWindow->Start();

  // set the size in the render window interactor
  Size[0] = size[0];
  Size[1] = size[1];

  // this is initialized
  Initialized = 1;
}
//---------------------------------------------------------------------------
void wxVTKRenderWindowInteractor::Enable()
{
  // if already enabled then done
  if (Enabled)
    return;

  // that's it
  Enabled = 1;
#if defined(__WXGTK__) && defined(wxUSE_GLCANVAS)
#if wxCHECK_VERSION(2, 9, 0)
  SetCurrent( *GLContext );
#else
  SetCurrent();
#endif
#endif
  Modified();
}
//---------------------------------------------------------------------------
bool wxVTKRenderWindowInteractor::Enable(bool enable)
{
#if defined(__WXGTK__) && defined(wxUSE_GLCANVAS)
  return wxGLCanvas::Enable(enable);
#else
  return wxWindow::Enable(enable);
#endif
}
//---------------------------------------------------------------------------
void wxVTKRenderWindowInteractor::Disable()
{
  // if already disabled then done
  if (!Enabled)
    return;

  // that's it (we can't remove the event handler like it should be...)
  Enabled = 0;
  Modified();
}
//---------------------------------------------------------------------------
void wxVTKRenderWindowInteractor::Start()
{
  // the interactor cannot control the event loop
  vtkErrorMacro( << "wxVTKRenderWindowInteractor::Start() "
    "interactor cannot control event loop.");
}
//---------------------------------------------------------------------------
void wxVTKRenderWindowInteractor::UpdateSize(int x, int y)
{
  if( RenderWindow )
  {
    // if the size changed tell render window
    if ( x != Size[0] || y != Size[1] )
    {
      // adjust our (vtkRenderWindowInteractor size)
      Size[0] = x;
      Size[1] = y;
      // and our RenderWindow's size
      RenderWindow->SetSize(x, y);
    }
  }
}
//---------------------------------------------------------------------------
int wxVTKRenderWindowInteractor::CreateTimer(int WXUNUSED(timertype))
{
  // it's a one shot timer
  if (!timer.Start(10, TRUE))
    return 0;

  return 1;
  
}
#if VTK_MAJOR_VERSION > 5 || (VTK_MAJOR_VERSION == 5 && VTK_MINOR_VERSION >= 2)
//------------------------------------------------------------------
int wxVTKRenderWindowInteractor::InternalCreateTimer(int timerId, int timerType,
                                                     unsigned long duration)
{
  if (!timer.Start(duration, timerType == OneShotTimer))
    return 0;
    
  return ID_wxVTKRenderWindowInteractor_TIMER;
}
//------------------------------------------------------------------
int wxVTKRenderWindowInteractor::InternalDestroyTimer(int platformTimerId)
{
  timer.Stop();
  return 1;
}
#endif
//---------------------------------------------------------------------------
int wxVTKRenderWindowInteractor::DestroyTimer()
{
  // do nothing
  return 1;
}
//---------------------------------------------------------------------------
void wxVTKRenderWindowInteractor::OnTimer(wxTimerEvent& WXUNUSED(event))
{
  if (!Enabled)
    return;
    
#if VTK_MAJOR_VERSION > 4 || (VTK_MAJOR_VERSION == 4 && VTK_MINOR_VERSION > 0)
  // new style
#if VTK_MAJOR_VERSION > 5 || (VTK_MAJOR_VERSION == 5 && VTK_MINOR_VERSION >= 2)
  // pass the right timer id
  int timerId = this->GetCurrentTimerId();
  this->InvokeEvent(vtkCommand::TimerEvent, &timerId);
#else
  this->InvokeEvent(vtkCommand::TimerEvent, NULL);
#endif
#else
  // old style
  InteractorStyle->OnTimer();
#endif
}

//---------------------------------------------------------------------------
// NOTE on implementation:
// Bad luck you ended up in the only tricky place of this code.
// A few note, wxWidgets still refuse to provide such convenient method
// so I have to maintain it myself, eventhough this is completely integrated
// in wxPython...
// Anyway if this happen to break for you then compare to a recent version of wxPython
// and look for the function long wxPyGetWinHandle(wxWindow* win)
// in wxPython/src/helpers.cpp
long wxVTKRenderWindowInteractor::GetHandleHack()
{
  //helper function to hide the MSW vs GTK stuff
  long handle_tmp = 0;

// __WXMSW__ is for Win32
//__WXMAC__ stands for using Carbon C-headers, using either the CarbonLib/CFM or the native Mach-O builds (which then also use the latest features available)
// __WXGTK__ is for both gtk 1.2.x and gtk 2.x
#if defined(__WXMSW__) || defined(__WXMAC__)
    handle_tmp = (long)this->GetHandle();
#endif //__WXMSW__

//__WXCOCOA__ stands for using the objective-c Cocoa API
#ifdef __WXCOCOA__
   // Here is how to find the NSWindow
   wxTopLevelWindow* toplevel = dynamic_cast<wxTopLevelWindow*>(
     wxGetTopLevelParent( this ) );
   if (toplevel != NULL )    
   {
     // DRG - commented this out for Mac OS X Cocoa build - this does not compile and
     // also is not needed (9/14/10)
     //     handle_tmp = (long)toplevel->GetNSWindow();
   }
   // The NSView will be deducted from 
   // [(NSWindow*)Handle contentView]
   // if only I knew how to write that in c++
#endif //__WXCOCOA__

    // Find and return the actual X-Window.
#if defined(__WXGTK__) || defined(__WXX11__)
    return (long)GetXWindow(this);
#endif

//#ifdef __WXMOTIF__
//    handle_tmp = (long)this->GetXWindow();
//#endif

  return handle_tmp;
}
//---------------------------------------------------------------------------
void wxVTKRenderWindowInteractor::OnPaint(wxPaintEvent& WXUNUSED(event))
{
  //must always be here
  wxPaintDC pDC(this);

  //do it here rather than in the cstor: this is safer.
  if(!Handle)
  {
    Handle = GetHandleHack();
    RenderWindow->SetWindowId(reinterpret_cast<void *>(Handle));
// Cocoa
// this->GetNSView() <-> DisplayId
// this->GetTopLevel()->GetNSWindow() <-> WindowId
#ifdef __WXMSW__
    RenderWindow->SetParentId(reinterpret_cast<void *>(this->GetParent()->GetHWND()));
#endif //__WXMSW__
      
    // This is another hack to prevent the VTK Render Window from closing the display.
    // If VTK closes the display, ~wxContext chashes while trying to destroy its
    // glContext (because the display is closed). The Get -> Set makes this VTK
    // object think someone else is responsible for the display. 
    this->RenderWindow->SetDisplayId(this->RenderWindow->GetGenericDisplayId());
  }
  // get vtk to render to the wxWindows
  Render();
#ifdef __WXMAC__
  // This solves a problem with repainting after a window resize
  // See also: http://sourceforge.net/mailarchive/forum.php?thread_id=31690967&forum_id=41789
#ifdef __WXCOCOA__
  vtkCocoaRenderWindow * rwin = vtkCocoaRenderWindow::SafeDownCast(RenderWindow);
  if( rwin )
  {
    rwin->UpdateContext();
  }
#else
  vtkCarbonRenderWindow* rwin = vtkCarbonRenderWindow::SafeDownCast(RenderWindow);
  if( rwin )
  {
#if VTK_MAJOR_VERSION > 4 || (VTK_MAJOR_VERSION == 4 && VTK_MINOR_VERSION > 4)
    // Must be somewhere after VTK 4.4
    rwin->UpdateGLRegion();
#endif
  }
#endif
#endif
}
//---------------------------------------------------------------------------
void wxVTKRenderWindowInteractor::OnEraseBackground(wxEraseEvent &event)
{
  //turn off background erase to reduce flickering on MSW
  event.Skip(false);
}
//---------------------------------------------------------------------------
void wxVTKRenderWindowInteractor::OnSize(wxSizeEvent& WXUNUSED(event))
{
  int w, h;
  GetClientSize(&w, &h);
  UpdateSize(w, h);

  if (!Enabled) 
    {
    return;
    }

#if VTK_MAJOR_VERSION > 4 || (VTK_MAJOR_VERSION == 4 && VTK_MINOR_VERSION > 0)
  InvokeEvent(vtkCommand::ConfigureEvent, NULL);
#endif
  //this will check for Handle
  //Render();
}
//---------------------------------------------------------------------------
void wxVTKRenderWindowInteractor::OnMotion(wxMouseEvent &event)
{
 if (!Enabled) 
    {
    return;
    }
#if VTK_MAJOR_VERSION > 4 || (VTK_MAJOR_VERSION == 4 && VTK_MINOR_VERSION > 0)
  SetEventInformationFlipY(event.GetX(), event.GetY(), 
    event.CmdDown(), event.ShiftDown(), '\0', 0, NULL);

  InvokeEvent(vtkCommand::MouseMoveEvent, NULL);
#else
  InteractorStyle->OnMouseMove(event.CmdDown(), event.ShiftDown(),
    event.GetX(), Size[1] - event.GetY() - 1);
#endif
}
//---------------------------------------------------------------------------
#if !(VTK_MAJOR_VERSION == 3 && VTK_MINOR_VERSION == 1)
void wxVTKRenderWindowInteractor::OnEnter(wxMouseEvent &event)
{
  if (!Enabled) 
    {
    return;
    }

#if VTK_MAJOR_VERSION > 4 || (VTK_MAJOR_VERSION == 4 && VTK_MINOR_VERSION > 0)
    // new style
  SetEventInformationFlipY(event.GetX(), event.GetY(), 
      event.CmdDown(), event.ShiftDown(), '\0', 0, NULL);

  InvokeEvent(vtkCommand::EnterEvent, NULL);
#else
    // old style
  InteractorStyle->OnEnter(event.CmdDown(), event.ShiftDown(),
      event.GetX(), Size[1] - event.GetY() - 1);  
#endif
}
//---------------------------------------------------------------------------
void wxVTKRenderWindowInteractor::OnLeave(wxMouseEvent &event)
{
  if (!Enabled) 
    {
    return;
    }

#if VTK_MAJOR_VERSION > 4 || (VTK_MAJOR_VERSION == 4 && VTK_MINOR_VERSION > 0)
    // new style
  SetEventInformationFlipY(event.GetX(), event.GetY(), 
      event.CmdDown(), event.ShiftDown(), '\0', 0, NULL);

  InvokeEvent(vtkCommand::LeaveEvent, NULL);
#else
    // old style
  InteractorStyle->OnLeave(event.CmdDown(), event.ShiftDown(),
      event.GetX(), Size[1] - event.GetY() - 1);  
#endif
}
//---------------------------------------------------------------------------
void wxVTKRenderWindowInteractor::OnKeyDown(wxKeyEvent &event)
{
  if (!Enabled) 
    {
    return;
    }

#if VTK_MAJOR_VERSION > 4 || (VTK_MAJOR_VERSION == 4 && VTK_MINOR_VERSION > 0)
    // new style
  int keycode = event.GetKeyCode();
  char key = '\0';
  if (((unsigned int)keycode) < 256)
  {
    // TODO: Unicode in non-Unicode mode ??
    key = (char)keycode;
  }

  const char* keysym = ascii_to_key_sym( keycode );
  if(!keysym)
  {
    // get virtual keys
    keysym = wx_key_to_key_sym( keycode );
  }
  
  if(!keysym)
  {
    keysym = "None";
  }
  

  // we don't get a valid mouse position inside the key event on every platform
  // so we retrieve the mouse position explicitly and pass it along
  wxPoint mousePos = ScreenToClient(wxGetMousePosition());
  SetEventInformationFlipY(mousePos.x, mousePos.y, 
                           event.CmdDown(), event.ShiftDown(), key, 0, keysym);
  InvokeEvent(vtkCommand::KeyPressEvent, NULL);
#else
  InteractorStyle->OnKeyDown(event.CmdDown(), event.ShiftDown(), 
    event.GetKeyCode(), 1);
#endif
  event.Skip();
}
//---------------------------------------------------------------------------
void wxVTKRenderWindowInteractor::OnKeyUp(wxKeyEvent &event)
{
  if (!Enabled) 
    {
    return;
    }

#if VTK_MAJOR_VERSION > 4 || (VTK_MAJOR_VERSION == 4 && VTK_MINOR_VERSION > 0)
    // new style
  int keycode = event.GetKeyCode();
  char key = '\0';
  if (((unsigned int)keycode) < 256)
  {
    // TODO: Unicode in non-Unicode mode ??
    key = (char)keycode;
  }

  const char* keysym = ascii_to_key_sym( keycode );
  if(!keysym)
  {
    // get virtual keys
    keysym = wx_key_to_key_sym( keycode );
  }
  
  if(!keysym)
  {
    keysym = "None";
  }
  

  // we don't get a valid mouse position inside the key event on every platform
  // so we retrieve the mouse position explicitly and pass it along
  wxPoint mousePos = ScreenToClient(wxGetMousePosition());
  SetEventInformationFlipY(mousePos.x, mousePos.y, 
                           event.CmdDown(), event.ShiftDown(), key, 0, keysym);
  InvokeEvent(vtkCommand::KeyReleaseEvent, NULL);
#else
  InteractorStyle->OnKeyUp(event.CmdDown(), event.ShiftDown(), 
    event.GetKeyCode(), 1);
#endif
  event.Skip();
}
#endif //!(VTK_MAJOR_VERSION == 3 && VTK_MINOR_VERSION == 1)
 //---------------------------------------------------------------------------
void wxVTKRenderWindowInteractor::OnChar(wxKeyEvent &event)
{
  if (!Enabled) 
    {
    return;
    }
    
#if VTK_MAJOR_VERSION > 4 || (VTK_MAJOR_VERSION == 4 && VTK_MINOR_VERSION > 0)
  // new style
  int keycode = event.GetKeyCode();
  char key = '\0';
  if (((unsigned int)keycode) < 256)
    {
    // TODO: Unicode in non-Unicode mode ??
    key = (char)keycode;
    }

  const char* keysym = ascii_to_key_sym( keycode );
  if(!keysym)
  {
    // get virtual keys
    keysym = wx_key_to_key_sym( keycode );
  }
  
  if(!keysym)
  {
    keysym = "None";
  }
  
  // we don't get a valid mouse position inside the key event on every platform
  // so we retrieve the mouse position explicitly and pass it along
  wxPoint mousePos = ScreenToClient(wxGetMousePosition());
  SetEventInformationFlipY(mousePos.x, mousePos.y, 
                           event.CmdDown(), event.ShiftDown(), key, 0, keysym);
  InvokeEvent(vtkCommand::CharEvent, NULL);
#endif
  event.Skip();
}
//---------------------------------------------------------------------------
void wxVTKRenderWindowInteractor::OnButtonDown(wxMouseEvent &event)
{
  if (!Enabled || (ActiveButton != wxEVT_NULL))
    {
    return;
    }
  ActiveButton = event.GetEventType();

    // On Mac (Carbon) and Windows we don't automatically get the focus when
    // you click inside the window
    // we therefore set the focus explicitly
    // Apparently we need that on linux (GTK) too:
    this->SetFocus();

#if VTK_MAJOR_VERSION > 4 || (VTK_MAJOR_VERSION == 4 && VTK_MINOR_VERSION > 0)
  SetEventInformationFlipY(event.GetX(), event.GetY(), 
    event.CmdDown(), event.ShiftDown(), '\0', 0, NULL);
#endif

  if(event.RightDown())
  {
#if VTK_MAJOR_VERSION > 4 || (VTK_MAJOR_VERSION == 4 && VTK_MINOR_VERSION > 0)
    // new style
    InvokeEvent(vtkCommand::RightButtonPressEvent, NULL);
#else            
    // old style
    InteractorStyle->OnRightButtonDown(event.CmdDown(), event.ShiftDown(),
      event.GetX(), Size[1] - event.GetY() - 1);
#endif
  }
  else if(event.LeftDown())
  {
#if VTK_MAJOR_VERSION > 4 || (VTK_MAJOR_VERSION == 4 && VTK_MINOR_VERSION > 0)
    // new style
    InvokeEvent(vtkCommand::LeftButtonPressEvent, NULL);
#else            
    // old style
    InteractorStyle->OnLeftButtonDown(event.CmdDown(),  event.ShiftDown(),
      event.GetX(), Size[1] - event.GetY() - 1);
#endif
  }
  else if (event.LeftDClick())
  {
#if VTK_MAJOR_VERSION > 4 || (VTK_MAJOR_VERSION == 4 && VTK_MINOR_VERSION > 0)
    // new style
//     this->SetRepeatCount (1);
    InvokeEvent(vtkCommand::LeftButtonPressEvent, NULL);
//     this->SetRepeatCount (0);
#else            
    // old style
    InteractorStyle->OnLeftButtonDown(event.CmdDown(),  event.ShiftDown(),
      event.GetX(), Size[1] - event.GetY() - 1);
#endif
  }
  
  else if(event.MiddleDown())
  {
#if VTK_MAJOR_VERSION > 4 || (VTK_MAJOR_VERSION == 4 && VTK_MINOR_VERSION > 0)
    // new style
    InvokeEvent(vtkCommand::MiddleButtonPressEvent, NULL);
#else            
    // old style
    InteractorStyle->OnMiddleButtonDown(event.CmdDown(), event.ShiftDown(),
      event.GetX(), Size[1] - event.GetY() - 1);
#endif
  }
  //save the button and capture mouse until the button is released
  //Only if :
  //1. it is possible (WX_USE_X_CAPTURE)
  //2. user decided to.
  if ((ActiveButton != wxEVT_NULL) && WX_USE_X_CAPTURE && UseCaptureMouse)
  {
    CaptureMouse();
  }
}
//---------------------------------------------------------------------------
void wxVTKRenderWindowInteractor::OnButtonUp(wxMouseEvent &event)
{
  //EVT_xxx_DOWN == EVT_xxx_UP - 1
  //This is only needed if two mouse buttons are pressed at the same time.
  //In wxWindows 2.4 and later: better use of wxMOUSE_BTN_RIGHT or 
  //wxEVT_COMMAND_RIGHT_CLICK
  if (!Enabled || (ActiveButton != (event.GetEventType()-1))) 
    {
    return;
    }

  // See report by Shang Mu / Kerry Loux on wxVTK mailing list
    this->SetFocus();

#if VTK_MAJOR_VERSION > 4 || (VTK_MAJOR_VERSION == 4 && VTK_MINOR_VERSION > 0)
  SetEventInformationFlipY(event.GetX(), event.GetY(), 
    event.CmdDown(), event.ShiftDown(), '\0', 0, NULL);
#endif
  
  if(ActiveButton == wxEVT_RIGHT_DOWN)
  {
#if VTK_MAJOR_VERSION > 4 || (VTK_MAJOR_VERSION == 4 && VTK_MINOR_VERSION > 0)
    // new style
    InvokeEvent(vtkCommand::RightButtonReleaseEvent, NULL);
#else            
    // old style
    InteractorStyle->OnRightButtonUp(event.CmdDown(), event.ShiftDown(),
      event.GetX(), Size[1] - event.GetY() - 1);
#endif
  }
  else if(ActiveButton == wxEVT_LEFT_DOWN)
  {
#if VTK_MAJOR_VERSION > 4 || (VTK_MAJOR_VERSION == 4 && VTK_MINOR_VERSION > 0)
    // new style
    InvokeEvent(vtkCommand::LeftButtonReleaseEvent, NULL);
#else            
    // old style
    InteractorStyle->OnLeftButtonUp(event.CmdDown(), event.ShiftDown(),
      event.GetX(), Size[1] - event.GetY() - 1);
#endif
  }
  else if(ActiveButton == wxEVT_MIDDLE_DOWN)
  {
#if VTK_MAJOR_VERSION > 4 || (VTK_MAJOR_VERSION == 4 && VTK_MINOR_VERSION > 0)
    // new style
    InvokeEvent(vtkCommand::MiddleButtonReleaseEvent, NULL);
#else            
    // old style
    InteractorStyle->OnMiddleButtonUp(event.CmdDown(), event.ShiftDown(),
      event.GetX(), Size[1] - event.GetY() - 1);
#endif
  }
  //if the ActiveButton is realeased, then release mouse capture
  if ((ActiveButton != wxEVT_NULL) && WX_USE_X_CAPTURE && UseCaptureMouse)
  {
    ReleaseMouse();
  }
  ActiveButton = wxEVT_NULL;
}
//---------------------------------------------------------------------------
void wxVTKRenderWindowInteractor::OnMouseWheel(wxMouseEvent& event)
{
// Mouse wheel was only added after VTK 4.4 (I think...)
#if VTK_MAJOR_VERSION > 4 || (VTK_MAJOR_VERSION == 4 && VTK_MINOR_VERSION > 4)
  // new style
  //Set vtk event information ... The numebr of wheel rotations is stored in
  //the x varible.  y varible is zero
  SetEventInformationFlipY(event.GetX() , event.GetY(), 
                           event.CmdDown(), event.ShiftDown(), '\0', 0, NULL);
  if(event.GetWheelRotation() > 0)
    {
      //Send event to VTK
      InvokeEvent(vtkCommand::MouseWheelForwardEvent, NULL);
    }
  else
    {
      //Send event to VTK
      InvokeEvent(vtkCommand::MouseWheelBackwardEvent, NULL);
    }
#endif
    
}

//---------------------------------------------------------------------------
#if wxCHECK_VERSION(2, 8, 0)
void wxVTKRenderWindowInteractor::OnMouseCaptureLost(wxMouseCaptureLostEvent& event)
{
   if (ActiveButton != wxEVT_NULL)
   {
       //Maybe also invoke the button release event here
   }
   // Reset ActiveButton so that
   // 1. we do not process mouse button up events any more,
   // 2. the next button down event will be processed and call CaptureMouse().
   // Otherwise ReleaseMouse() will be called
   // without a previous CaptureMouse().
   ActiveButton = wxEVT_NULL;
}
#endif

//---------------------------------------------------------------------------
void wxVTKRenderWindowInteractor::Render()
{
#if wxCHECK_VERSION(2, 8, 0)
  int renderAllowed = !IsFrozen();
#else
  int renderAllowed = 1;
#endif
  if (renderAllowed && !RenderWhenDisabled)
    {
    //the user doesn't want us to render when the toplevel frame
    //is disabled - first find the top level parent
    wxWindow *topParent = wxGetTopLevelParent(this);
    if (topParent)
      {
      //if it exists, check whether it's enabled
      //if it's not enabeld, renderAllowed will be false
      renderAllowed = topParent->IsEnabled();
      }
    }

  if (renderAllowed)
    {
    if(Handle && (Handle == GetHandleHack()) )
      {
      RenderWindow->Render();
      }
#if VTK_MAJOR_VERSION > 4 || (VTK_MAJOR_VERSION == 4 && VTK_MINOR_VERSION > 2)
    else if(GetHandleHack())
      {
      //this means the user has reparented us; let's adapt to the
      //new situation by doing the WindowRemap dance
      //store the new situation
      Handle = GetHandleHack();
      RenderWindow->SetNextWindowId(reinterpret_cast<void *>(Handle));
      RenderWindow->WindowRemap();
      RenderWindow->Render();
      }
#endif
    }
}
//---------------------------------------------------------------------------
void wxVTKRenderWindowInteractor::SetRenderWhenDisabled(int newValue)
{
  //Change value of __RenderWhenDisabled ivar.
  //If __RenderWhenDisabled is false (the default), this widget will not
  //call Render() on the RenderWindow if the top level frame (i.e. the
  //containing frame) has been disabled.

  //This prevents recursive rendering during wxSafeYield() calls.
  //wxSafeYield() can be called during the ProgressMethod() callback of
  //a VTK object to have progress bars and other GUI elements updated -
  //it does this by disabling all windows (disallowing user-input to
  //prevent re-entrancy of code) and then handling all outstanding
  //GUI events.
        
  //However, this often triggers an OnPaint() method for wxVTKRWIs,
  //resulting in a Render(), resulting in Update() being called whilst
  //still in progress.

  RenderWhenDisabled = static_cast<bool>(newValue);
}
//---------------------------------------------------------------------------
//
// Set the variable that indicates that we want a stereo capable window
// be created. This method can only be called before a window is realized.
//
void wxVTKRenderWindowInteractor::SetStereo(int capable)
{
  if (Stereo != capable)
    {
    Stereo = capable;
    RenderWindow->StereoCapableWindowOn();
    RenderWindow->SetStereoTypeToCrystalEyes();
    Modified();
    }
}

//---------------------------------------------------------------------------
//
//
void wxVTKRenderWindowInteractor::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}



// ***** keysym stuff below  *****

static const char *AsciiToKeySymTable[] = {
  0, 0, 0, 0, 0, 0, 0, 0, 0, "Tab", 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  "space", "exclam", "quotedbl", "numbersign",
  "dollar", "percent", "ampersand", "quoteright",
  "parenleft", "parenright", "asterisk", "plus",
  "comma", "minus", "period", "slash",
  "0", "1", "2", "3", "4", "5", "6", "7",
  "8", "9", "colon", "semicolon", "less", "equal", "greater", "question",
  "at", "A", "B", "C", "D", "E", "F", "G",
  "H", "I", "J", "K", "L", "M", "N", "O",
  "P", "Q", "R", "S", "T", "U", "V", "W",
  "X", "Y", "Z", "bracketleft",
  "backslash", "bracketright", "asciicircum", "underscore",
  "quoteleft", "a", "b", "c", "d", "e", "f", "g",
  "h", "i", "j", "k", "l", "m", "n", "o",
  "p", "q", "r", "s", "t", "u", "v", "w",
  "x", "y", "z", "braceleft", "bar", "braceright", "asciitilde", "Delete",
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

const char* ascii_to_key_sym(int i)
{
  if(i >= 0 && i<256 )
  {
    return AsciiToKeySymTable[i];
  }
  return 0;
}


#define WXVTK_HANDLE(x,y) \
  case x : \
    ret = y; \
    break;

const char* wx_key_to_key_sym(int i)
{
  const char* ret = 0;
  switch(i)
  {
      WXVTK_HANDLE(WXK_CLEAR, "Clear")
      WXVTK_HANDLE(WXK_BACK, "BackSpace")
      WXVTK_HANDLE(WXK_TAB, "Tab")
      WXVTK_HANDLE(WXK_RETURN, "Return")
      WXVTK_HANDLE(WXK_SHIFT, "Shift_L")
      WXVTK_HANDLE(WXK_CONTROL, "Control_L")
      WXVTK_HANDLE(WXK_ALT, "Alt_L")
      WXVTK_HANDLE(WXK_PAUSE, "Pause")
      WXVTK_HANDLE(WXK_CAPITAL, "Caps_Lock")
      WXVTK_HANDLE(WXK_ESCAPE, "Escape")
      WXVTK_HANDLE(WXK_SPACE, "space")
      WXVTK_HANDLE(WXK_END, "End")
      WXVTK_HANDLE(WXK_HOME, "Home")
      WXVTK_HANDLE(WXK_LEFT, "Left")
      WXVTK_HANDLE(WXK_UP, "Up")
      WXVTK_HANDLE(WXK_RIGHT, "Right")
      WXVTK_HANDLE(WXK_DOWN, "Down")
      WXVTK_HANDLE(WXK_SNAPSHOT, "Snapshot")
      WXVTK_HANDLE(WXK_INSERT, "Insert")
      WXVTK_HANDLE(WXK_DELETE, "Delete")
      WXVTK_HANDLE(WXK_HELP, "Help")
      WXVTK_HANDLE(WXK_NUMPAD0, "0")
      WXVTK_HANDLE(WXK_NUMPAD1, "1")
      WXVTK_HANDLE(WXK_NUMPAD2, "2")
      WXVTK_HANDLE(WXK_NUMPAD3, "3")
      WXVTK_HANDLE(WXK_NUMPAD4, "4")
      WXVTK_HANDLE(WXK_NUMPAD5, "5")
      WXVTK_HANDLE(WXK_NUMPAD6, "6")
      WXVTK_HANDLE(WXK_NUMPAD7, "7")
      WXVTK_HANDLE(WXK_NUMPAD8, "8")
      WXVTK_HANDLE(WXK_NUMPAD9, "9")
      // KP_0 - KP_9
      WXVTK_HANDLE(WXK_MULTIPLY, "asterisk")
      WXVTK_HANDLE(WXK_ADD, "plus")
      // bar
      WXVTK_HANDLE(WXK_SUBTRACT, "minus")
      WXVTK_HANDLE(WXK_DECIMAL, "period")
      WXVTK_HANDLE(WXK_SEPARATOR, "slash")
      WXVTK_HANDLE(WXK_F1, "F1")
      WXVTK_HANDLE(WXK_F2, "F2")
      WXVTK_HANDLE(WXK_F3, "F3")
      WXVTK_HANDLE(WXK_F4, "F4")
      WXVTK_HANDLE(WXK_F5, "F5")
      WXVTK_HANDLE(WXK_F6, "F6")
      WXVTK_HANDLE(WXK_F7, "F7")
      WXVTK_HANDLE(WXK_F8, "F8")
      WXVTK_HANDLE(WXK_F9, "F9")
      WXVTK_HANDLE(WXK_F10, "F10")
      WXVTK_HANDLE(WXK_F11, "F11")
      WXVTK_HANDLE(WXK_F12, "F12")
      WXVTK_HANDLE(WXK_F13, "F13")
      WXVTK_HANDLE(WXK_F14, "F14")
      WXVTK_HANDLE(WXK_F15, "F15")
      WXVTK_HANDLE(WXK_F16, "F16")
      WXVTK_HANDLE(WXK_F17, "F17")
      WXVTK_HANDLE(WXK_F18, "F18")
      WXVTK_HANDLE(WXK_F19, "F19")
      WXVTK_HANDLE(WXK_F20, "F20")
      WXVTK_HANDLE(WXK_F21, "F21")
      WXVTK_HANDLE(WXK_F22, "F22")
      WXVTK_HANDLE(WXK_F23, "F23")
      WXVTK_HANDLE(WXK_F24, "F24")
      WXVTK_HANDLE(WXK_NUMLOCK, "Num_Lock")
      WXVTK_HANDLE(WXK_SCROLL, "Scroll_Lock")
      WXVTK_HANDLE(WXK_PAGEUP, "Page_Up")
      WXVTK_HANDLE(WXK_PAGEDOWN, "Page_Down")
#if WXWIN_COMPATIBILITY_2_6
      WXVTK_HANDLE(WXK_NUMPAD_PRIOR, "Page_Up")
      WXVTK_HANDLE(WXK_NUMPAD_NEXT, "Page_Down")
#endif	
      WXVTK_HANDLE(WXK_NUMPAD_END, "End")
	WXVTK_HANDLE(WXK_NUMPAD_BEGIN, "Begin")
	WXVTK_HANDLE(WXK_NUMPAD_INSERT, "Insert")
	WXVTK_HANDLE(WXK_NUMPAD_DELETE, "Delete")
	WXVTK_HANDLE(WXK_NUMPAD_EQUAL, "Equal")
	WXVTK_HANDLE(WXK_NUMPAD_MULTIPLY, "asterix")
	WXVTK_HANDLE(WXK_NUMPAD_ADD, "plus")
	WXVTK_HANDLE(WXK_NUMPAD_SEPARATOR, "slash")
	WXVTK_HANDLE(WXK_NUMPAD_SUBTRACT, "minus")
	WXVTK_HANDLE(WXK_NUMPAD_DECIMAL, "period")
	WXVTK_HANDLE(WXK_NUMPAD_DIVIDE, "slash")

	WXVTK_HANDLE(WXK_WINDOWS_LEFT, "Win_L")
	WXVTK_HANDLE(WXK_WINDOWS_RIGHT, "Win_R")
	WXVTK_HANDLE(WXK_WINDOWS_MENU, "Menu")
	WXVTK_HANDLE(WXK_COMMAND, "Command")


      default:
    break;
    }
  return ret;
}



//---------------------------------------------------------------------------
void wxVTKRenderWindowInteractor::OnMouseDoubleCilck(wxMouseEvent &event)
{
  if (!Enabled || (ActiveButton != wxEVT_NULL))
  {
    return;
  }
  ActiveButton = event.GetEventType();
  
  // On Mac (Carbon) and Windows we don't automatically get the focus when
  // you click inside the window
  // we therefore set the focus explicitly
  // Apparently we need that on linux (GTK) too:
  this->SetFocus();
  
#if VTK_MAJOR_VERSION > 4 || (VTK_MAJOR_VERSION == 4 && VTK_MINOR_VERSION > 0)
  SetEventInformationFlipY(event.GetX(), event.GetY(), 
			   event.CmdDown(), event.ShiftDown(), '\0', 1, NULL);
#endif
  
  if (event.LeftDClick())
  {
#if VTK_MAJOR_VERSION > 4 || (VTK_MAJOR_VERSION == 4 && VTK_MINOR_VERSION > 0)
    // new style
    //     this->SetRepeatCount (1);
    InvokeEvent(vtkCommand::LeftButtonPressEvent, NULL);
    InvokeEvent(vtkCommand::LeftButtonReleaseEvent, NULL);
    //     this->SetRepeatCount (0);
#else            
    // old style
    InteractorStyle->OnLeftButtonDown(event.CmdDown(),  event.ShiftDown(),
				      event.GetX(), Size[1] - event.GetY() - 1);
#endif
  }
  
  ActiveButton = wxEVT_NULL;
  
}

  
