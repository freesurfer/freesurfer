/**
 * @file  vtkKWScubaApplicationSettingsInterface.cxx
 * @brief Preferences for Scuba
 *
 * Subclass of vtkKWApplicationSettingsInterface that adds our own
 * custom settings.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.3 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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

#include "vtkKWScubaApplicationSettingsInterface.h"

#include "vtkKWCheckButton.h"
#include "vtkObjectFactory.h"
#include "vtkKWFrameWithLabel.h"
#include "vtkKWScubaWindow.h"

vtkStandardNewMacro(vtkKWScubaApplicationSettingsInterface);
vtkCxxRevisionMacro(vtkKWScubaApplicationSettingsInterface, "$Revision: 1.3 $");

vtkKWScubaApplicationSettingsInterface::vtkKWScubaApplicationSettingsInterface() :
  mChkBtnAutoSizeInfoArea( NULL ) {
}

vtkKWScubaApplicationSettingsInterface::~vtkKWScubaApplicationSettingsInterface() {
  if( mChkBtnAutoSizeInfoArea )
    mChkBtnAutoSizeInfoArea->Delete();
}

void 
vtkKWScubaApplicationSettingsInterface::Create() {
  if( this->IsCreated() ) {
    vtkErrorMacro("The panel is already created.");
    return;
  }

  // Create the superclass instance (and set the application)
  this->Superclass::Create();

  // Put our checkbox in the interface settings frame.
  mChkBtnAutoSizeInfoArea = vtkKWCheckButton::New();
  mChkBtnAutoSizeInfoArea->SetParent( InterfaceSettingsFrame->GetFrame() );
  mChkBtnAutoSizeInfoArea->Create();
  mChkBtnAutoSizeInfoArea->SetText( "Automatically Resize Info Area" );
  mChkBtnAutoSizeInfoArea->
    SetBalloonHelpString( "Resize the info area according to how many rows "
			  "of information there are" );
  mChkBtnAutoSizeInfoArea->SetCommand( this, "AutoSizeInfoAreaCallback" );

  this->Script( "pack %s -side top -anchor w -expand no -fill none",
		mChkBtnAutoSizeInfoArea->GetWidgetName() );
  
  // Update
  this->Update();
}

void
vtkKWScubaApplicationSettingsInterface::SetScubaWindow ( vtkKWScubaWindow* iWindow ){

  if ( this->ScubaWindow == iWindow )
    return;

  ScubaWindow = iWindow;
  this->Modified();

  this->Update();
}

void
vtkKWScubaApplicationSettingsInterface::Update () {
  this->Superclass::Update();

  if( !this->IsCreated() || NULL == ScubaWindow )
    return;

  // Check our resize option.
  if( mChkBtnAutoSizeInfoArea )
    mChkBtnAutoSizeInfoArea->
      SetSelectedState( ScubaWindow->GetAutoSizeInfoArea() );
}

void
vtkKWScubaApplicationSettingsInterface::AutoSizeInfoAreaCallback ( int ibSize ) {

  if( NULL != ScubaWindow )
    ScubaWindow->SetAutoSizeInfoArea( ibSize );
}

void
vtkKWScubaApplicationSettingsInterface::UpdateEnableState () {

  this->Superclass::UpdateEnableState();

  if ( mChkBtnAutoSizeInfoArea )
    mChkBtnAutoSizeInfoArea->SetEnabled(this->GetEnabled());
}
