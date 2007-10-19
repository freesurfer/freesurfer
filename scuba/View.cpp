/**
 * @file  View.cpp
 * @brief Basic View superclass, or a subpanel in a ScubaFrame.
 *
 * This is a base class used by ScubaFrame as an interface to views,
 * or the internal components of the Frame. Each Frame can have
 * multiple non-overlapping Views, and each View can have multiple
 * layers.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/19 22:31:57 $
 *    $Revision: 1.14 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#include "View.h"

using namespace std;

View::View() :
    Broadcaster( "View" ),
    Listener( "View" ),
    mWidth( 0 ),
    mHeight( 0 ),
    mbPostRedisplay( false ) {}

View::~View() {}

void
View::Draw () {
  this->DoDraw();
}

void
View::Reshape ( int iWidth, int iHeight ) {
  mWidth = iWidth;
  mHeight = iHeight;
  this->DoReshape( iWidth, iHeight );
}

void
View::Timer() {

  this->DoTimer();
}

void
View::MouseMoved( int iWindow[2], InputState& iInput, ScubaToolState& iTool ) {

  this->DoMouseMoved( iWindow, iInput, iTool );
}

void
View::MouseUp( int iWindow[2], InputState& iInput, ScubaToolState& iTool ) {

  this->DoMouseUp( iWindow, iInput, iTool );
}

void
View::MouseDown( int iWindow[2], InputState& iInput, ScubaToolState& iTool ) {

  this->DoMouseDown( iWindow, iInput, iTool );
}

void
View::KeyDown( int iWindow[2], InputState& iInput, ScubaToolState& iTool ) {

  this->DoKeyDown( iWindow, iInput, iTool );
}

void
View::KeyUp( int iWindow[2], InputState& iInput, ScubaToolState& iTool ) {

  this->DoKeyUp( iWindow, iInput, iTool );
}

TclCommandListener::TclCommandResult
View::DoListenToTclCommand ( char*, int, char** ) {
  return ok;
}

void
View::DoListenToMessage ( string, void* ) {}


void
View::DoDraw() {

  DebugOutput( << "View " << msLabel << ": DoDraw()" );
}

void
View::DoReshape( int, int ) {

  DebugOutput( << "View " << msLabel << ": DoReshape()" );
}

void
View::DoTimer() {

  DebugOutput( << "View " << msLabel << ": DoTimer()" );
}

void
View::DoMouseMoved( int[2], InputState&, ScubaToolState& ) {

  DebugOutput( << "View " << msLabel << ": DoMouseMoved()" );
}

void
View::DoMouseUp( int[2], InputState&, ScubaToolState& ) {

  DebugOutput( << "View " << msLabel << ": DoMouseUp()" );
}

void
View::DoMouseDown( int[2], InputState&, ScubaToolState& ) {

  DebugOutput( << "View " << msLabel << ": DoMouseDown()" );
}

void
View::DoKeyDown( int[2], InputState&, ScubaToolState& ) {

  DebugOutput( << "View " << msLabel << ": DoKeyDown()" );
}

void
View::DoKeyUp( int[2], InputState&, ScubaToolState& ) {

  DebugOutput( << "View " << msLabel << ": DoKeyUp()" );
}

