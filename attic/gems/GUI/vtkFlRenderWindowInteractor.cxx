/*
 * vtkFlRenderWindowInteractor - class to enable VTK to render to and interact
 * with a FLTK window.
 * 
 * Copyright (c) 2002-2006 Charl P. Botha http://cpbotha.net/
 * Based on original code and concept copyright (c) 2000,2001 David Pont
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 * 
 */

/*  
 * You must not delete one of these classes.  Make use of the Delete()
 * method... this thing makes use of VTK reference counting.  Let me
 * repeat that: never "delete" an instance of this class, always use
 * ->Delete().
 */

#include "vtkFlRenderWindowInteractor.h"
// FLTK
#include <FL/x.H>
// vtk
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkInteractorStyle.h>
#include <vtkVersion.h>
#include <vtkCommand.h>

//---------------------------------------------------------------------------
vtkFlRenderWindowInteractor::vtkFlRenderWindowInteractor() : 
Fl_Gl_Window( 0, 0, 300, 300, "" ), vtkRenderWindowInteractor()
{
    // this is a subclass of Fl_Group, call end so children cant be added
    this->end();
}
//---------------------------------------------------------------------------
vtkFlRenderWindowInteractor::vtkFlRenderWindowInteractor( int x, int y, int w, int h, const char *l ) : 
Fl_Gl_Window( x, y, w, h, l ), vtkRenderWindowInteractor()
{
    // this is a subclass of Fl_Group, call end so children cant be added
    this->end();
}
//---------------------------------------------------------------------------
vtkFlRenderWindowInteractor::~vtkFlRenderWindowInteractor()
{
    // according to the fltk docs, destroying a widget does NOT remove it from
    // its parent, so we have to do that explicitly at destruction
    // (and remember, NEVER delete() an instance of this class)
    if (parent())
    {
        ((Fl_Group*)parent())->remove(*(Fl_Gl_Window*)this);
    }
}
//---------------------------------------------------------------------------
vtkFlRenderWindowInteractor * vtkFlRenderWindowInteractor::New()
{
    // we don't make use of the objectfactory, because we're not registered
    return new vtkFlRenderWindowInteractor;
}

//---------------------------------------------------------------------------
void vtkFlRenderWindowInteractor::Initialize()
{
    // if don't have render window then we can't do anything yet
    if (!RenderWindow)
    {
	vtkErrorMacro(<< "vtkFlRenderWindowInteractor::Initialize has no render window");
	return;
    }

    int *size = RenderWindow->GetSize();
    // enable everything and start rendering
    Enable();
    
    // We should NOT call ->Render yet, as it's entirely possible that
    // Initialize() is called before there's a valid Fl_Gl_Window!
    //RenderWindow->Render();

    // set the size in the render window interactor
    Size[0] = size[0];
    Size[1] = size[1];

    // this is initialized
    Initialized = 1;
}
//---------------------------------------------------------------------------
void vtkFlRenderWindowInteractor::Enable()
{
    // if already enabled then done
    if (Enabled)
      return;

    // that's it
    Enabled = 1;
    Modified();
}
//---------------------------------------------------------------------------
void vtkFlRenderWindowInteractor::Disable()
{
    // if already disabled then done
    if (!Enabled)
      return;

    // that's it (we can't remove the event handler like it should be...)
    Enabled = 0;
    Modified();
}
//---------------------------------------------------------------------------
void vtkFlRenderWindowInteractor::Start()
{
    // the interactor cannot control the event loop
    vtkErrorMacro(<<"vtkFlRenderWindowInteractor::Start() interactor cannot control event loop.");
}
//---------------------------------------------------------------------------
void vtkFlRenderWindowInteractor::SetRenderWindow(vtkRenderWindow *aren)
{
   vtkRenderWindowInteractor::SetRenderWindow(aren);
   // if a vtkFlRenderWindowInteractor has been shown already, and one
   // re-sets the RenderWindow, neither UpdateSize nor draw is called,
   // so we have to force the dimensions of the NEW RenderWindow to match
   // the our (vtkFlRWI) dimensions
   if (RenderWindow)
       RenderWindow->SetSize(this->w(), this->h());
}
//---------------------------------------------------------------------------
// this gets called during FLTK window draw()s and resize()s
void vtkFlRenderWindowInteractor::UpdateSize(int W, int H)
{
    if (RenderWindow != NULL)
    {
	// if the size changed tell render window
	if ( (W != Size[0]) || (H != Size[1]) )
	{
	    // adjust our (vtkRenderWindowInteractor size)
	    Size[0] = W;
	    Size[1] = H;
	    // and our RenderWindow's size
	    RenderWindow->SetSize(W, H);
	 
	    // FLTK can move widgets on resize; if that happened, make
	    // sure the RenderWindow position agrees with that of the
	    // Fl_Gl_Window
	    int *pos = RenderWindow->GetPosition();
	    if( pos[0] != x() || pos[1] != y() ) {
		RenderWindow->SetPosition( x(), y() );
	    }
	}
    }
}
//---------------------------------------------------------------------------
// FLTK needs global timer callbacks, but we set it up so that this global
// callback knows which instance OnTimer() to call
void OnTimerGlobal(void *p)
{
    if (p)
      ((vtkFlRenderWindowInteractor *)p)->OnTimer();
}

//---------------------------------------------------------------------------
int vtkFlRenderWindowInteractor::CreateTimer(int timertype)
{
    // to be called every 10 milliseconds, one shot timer
    // we pass "this" so that the correct OnTimer instance will be called
    if (timertype == VTKI_TIMER_FIRST)
      Fl::add_timeout(0.01, OnTimerGlobal, (void *)this);
    else
      Fl::repeat_timeout(0.01, OnTimerGlobal, (void *)this);
   
    return 1;
    // Fl::repeat_timer() is more correct, it doesn't measure the timeout
    // from now, but from when the system call that caused this timeout
    // elapsed.
}
//---------------------------------------------------------------------------
int vtkFlRenderWindowInteractor::DestroyTimer()
{
    // do nothing
    return 1;
}

void vtkFlRenderWindowInteractor::OnTimer(void)
{
    if (!Enabled)
      return;
    // this is all we need to do, InteractorStyle is stateful and will
    // continue with whatever it's busy
    
#if (VTK_MAJOR_VERSION >= 4)
    // new style
    this->InvokeEvent(vtkCommand::TimerEvent, NULL);
#else
    // old style
    InteractorStyle->OnTimer();
#endif
    
}


//---------------------------------------------------------------------------
void vtkFlRenderWindowInteractor::TerminateApp()
{
    ;
}

//---------------------------------------------------------------------------
// FLTK event handlers
//---------------------------------------------------------------------------
void vtkFlRenderWindowInteractor::flush(void)
{
    // err, we don't want to do any fansy pansy Fl_Gl_Window stuff, so we
    // bypass all of it (else we'll get our front and back buffers in all
    // kinds of tangles, and need extra glXSwapBuffers() calls and all that)
    draw();
}

//---------------------------------------------------------------------------
void vtkFlRenderWindowInteractor::draw(void){
    if (RenderWindow!=NULL)
    {
	// make sure the vtk part knows where and how large we are
	UpdateSize( this->w(), this->h() );
        
        // make sure the GL context exists and is current:
        // after a hide() and show() sequence e.g. there is no context yet
        // and the Render() will fail due to an invalid context.
        // see Fl_Gl_Window::show()
        make_current();
        
	RenderWindow->SetWindowId( (void *)fl_xid( this ) );
#if !defined(WIN32) && !defined(__APPLE__)
	RenderWindow->SetDisplayId( fl_display );
#endif
	// get vtk to render to the Fl_Gl_Window
	Render();
    }
}
//---------------------------------------------------------------------------
void vtkFlRenderWindowInteractor::resize( int x, int y, int w, int h ) {
    // make sure VTK knows about the new situation
    UpdateSize( w, h );
    // resize the FLTK window by calling ancestor method
    Fl_Gl_Window::resize( x, y, w, h ); 
}
//---------------------------------------------------------------------------
// main FLTK event handler
int vtkFlRenderWindowInteractor::handle( int event ) {
    if( !Enabled ) return 0;
    
#if (VTK_MAJOR_VERSION >= 4)
    // setup for new style
    // SEI(x, y, ctrl, shift, keycode, repeatcount, keysym)
    this->SetEventInformation(Fl::event_x(), this->h()-Fl::event_y()-1, 
                              Fl::event_state( FL_CTRL ), Fl::event_state( FL_SHIFT ),
                              Fl::event_key(), 1, NULL);
#endif    
    
    switch( event ) 
    {
      case FL_FOCUS:
      case FL_UNFOCUS:
	;   // Return 1 if you want keyboard events, 0 otherwise. Yes we do
	break;
      
      case FL_KEYBOARD:   // keypress
#if (VTK_MAJOR_VERSION >= 4)
        // new style
        this->InvokeEvent(vtkCommand::MouseMoveEvent, NULL);        
        this->InvokeEvent(vtkCommand::KeyPressEvent, NULL);
        this->InvokeEvent(vtkCommand::CharEvent, NULL);
#else
        // old style
	InteractorStyle->OnChar(Fl::event_state( FL_CTRL ), Fl::event_state( FL_SHIFT ), Fl::event_key(), 1);
#endif        
	// now for possible controversy: there is no way to find out if the InteractorStyle actually did
	// something with this event.  To play it safe (and have working hotkeys), we return "0", which indicates
	// to FLTK that we did NOTHING with this event.  FLTK will send this keyboard event to other children
	// in our group, meaning it should reach any FLTK keyboard callbacks (including hotkeys)
	return 0;
	break;
      
      case FL_PUSH: // mouse down
	this->take_focus();  // this allows key events to work
	switch( Fl::event_button() ) 
	{
	  case FL_LEFT_MOUSE:
#if (VTK_MAJOR_VERSION >= 4)
            // new style
            this->InvokeEvent(vtkCommand::LeftButtonPressEvent,NULL);
#else
            // old style
	    InteractorStyle->OnLeftButtonDown(Fl::event_state( FL_CTRL ), Fl::event_state( FL_SHIFT ), Fl::event_x(), this->h()-Fl::event_y()-1);
#endif            
	    break;
	  case FL_MIDDLE_MOUSE:
#if (VTK_MAJOR_VERSION >= 4)
            // new style
            this->InvokeEvent(vtkCommand::MiddleButtonPressEvent,NULL);
#else
            // old style
	    InteractorStyle->OnMiddleButtonDown(Fl::event_state( FL_CTRL ), Fl::event_state( FL_SHIFT ), Fl::event_x(), this->h()-Fl::event_y()-1);
#endif
	    break;
	  case FL_RIGHT_MOUSE:
#if (VTK_MAJOR_VERSION >= 4)
            // new style
            this->InvokeEvent(vtkCommand::RightButtonPressEvent,NULL);
#else
            // old style
	    InteractorStyle->OnRightButtonDown(Fl::event_state( FL_CTRL ), Fl::event_state( FL_SHIFT ), Fl::event_x(), this->h()-Fl::event_y()-1);
#endif            
            break;
	}
	break; // this break should be here, at least according to vtkXRenderWindowInteractor

	// we test for both of these, as fltk classifies mouse moves as with or
	// without button press whereas vtk wants all mouse movement (this bug took
	// a while to find :)
      case FL_DRAG:
      case FL_MOVE:
#if (VTK_MAJOR_VERSION >= 4)
        // new style
        this->InvokeEvent(vtkCommand::MouseMoveEvent, NULL);
#else        
	// old style
        InteractorStyle->OnMouseMove(Fl::event_state( FL_CTRL ), Fl::event_state( FL_SHIFT ), Fl::event_x(), this->h()-Fl::event_y()-1);
#endif        
	break;

      case FL_RELEASE:    // mouse up
	switch( Fl::event_button() ) {
	  case FL_LEFT_MOUSE:
#if (VTK_MAJOR_VERSION >= 4)
            // new style
            this->InvokeEvent(vtkCommand::LeftButtonReleaseEvent,NULL);
#else            
	    // old style
            InteractorStyle->OnLeftButtonUp(Fl::event_state( FL_CTRL ), Fl::event_state( FL_SHIFT ), Fl::event_x(), this->h()-Fl::event_y()-1);
#endif
	    break;
	  case FL_MIDDLE_MOUSE:
#if (VTK_MAJOR_VERSION >= 4)
            // new style
            this->InvokeEvent(vtkCommand::MiddleButtonReleaseEvent,NULL);
#else            
	    // old style
            InteractorStyle->OnMiddleButtonUp(Fl::event_state( FL_CTRL ), Fl::event_state( FL_SHIFT ), Fl::event_x(), this->h()-Fl::event_y()-1);
#endif
	    break;
	  case FL_RIGHT_MOUSE:
#if (VTK_MAJOR_VERSION >= 4)
            // new style
            this->InvokeEvent(vtkCommand::RightButtonReleaseEvent,NULL);
#else            
	    // old style
            InteractorStyle->OnRightButtonUp(Fl::event_state( FL_CTRL ), Fl::event_state( FL_SHIFT ), Fl::event_x(), this->h()-Fl::event_y()-1);
#endif            
	    break;
	}
	break;

      default:	// let the base class handle everything else 
	return Fl_Gl_Window::handle( event );

    } // switch(event)...

    return 1; // we handled the event if we didn't return earlier
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


const char *vtkFlRenderWindowInteractor_rcsid(void)
{
    return "vtkFlRenderWindowInteractor";
}
