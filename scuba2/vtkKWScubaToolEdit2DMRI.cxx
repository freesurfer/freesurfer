/**
 * @file  vtkKWScubaToolEdit2DMRI.cxx
 * @brief A tool for editing MRI volumes.
 *
 * Tool implementation that works on MRI volumes (2D and 3D).
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.2 $
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


#include "vtkKWScubaToolEdit2DMRI.h"
#include "vtkObjectFactory.h"
#include "vtkKWScubaView.h"
#include "vtkKWScaleWithEntry.h"
#include "vtkKWScubaLayer2DMRI.h"
#include "vtkKWScubaLayer3DMRI.h"
#include "vtkFSVolumeSource.h"

using namespace std;

vtkStandardNewMacro( vtkKWScubaToolEdit2DMRI );
vtkCxxRevisionMacro( vtkKWScubaToolEdit2DMRI, "$Revision: 1.2 $" );

vtkKWScubaToolEdit2DMRI::vtkKWScubaToolEdit2DMRI () :
    mNewValue(0),
    mScaleNewValue(NULL) {

  msLabel = "Edit Voxels";
}

vtkKWScubaToolEdit2DMRI::~vtkKWScubaToolEdit2DMRI () {}

void
vtkKWScubaToolEdit2DMRI::AddControls ( vtkKWWidget* iPanel ) {

  mScaleNewValue = vtkKWScaleWithEntry::New();
  mScaleNewValue->SetParent( iPanel );
  mScaleNewValue->SetOrientationToHorizontal();
  mScaleNewValue->Create();
  mScaleNewValue->SetLabelText( "New Value: " );
  mScaleNewValue->SetRange( 0, 255 );
  mScaleNewValue->SetResolution( 1 );
  mScaleNewValue->SetEntryWidth( 3 );
  mScaleNewValue->SetCommand( this, "SetNewValue" );
  mScaleNewValue->SetValue( mNewValue );

  this->Script( "pack %s -side top -fill x -anchor nw",
                mScaleNewValue->GetWidgetName() );
}

void
vtkKWScubaToolEdit2DMRI::RemoveControls () {

  if( mScaleNewValue ) {
    mScaleNewValue->Delete();
    mScaleNewValue = NULL;
  }
}

int
vtkKWScubaToolEdit2DMRI::GetNewValue () {

  return mNewValue;
}

void
vtkKWScubaToolEdit2DMRI::SetNewValue ( int iNewValue ) {

  mNewValue = iNewValue;

  if ( mScaleNewValue )
    mScaleNewValue->SetValue( iNewValue );
}

void
vtkKWScubaToolEdit2DMRI::DoMouseDrag ( vtkKWScubaWindow* iWindow,
                                       vtkKWScubaView* iView,
                                       vtkKWScubaLayer* iLayer,
                                       float iRAS[3] ) {

  vtkKWScubaLayer2DMRI* layer2DMRI =
    vtkKWScubaLayer2DMRI::SafeDownCast( iLayer );
  if ( layer2DMRI ) {

    vtkFSVolumeSource* source = layer2DMRI->GetSource();

    int idx[3];
    source->ConvertRASToIndex( iRAS[0], iRAS[1], iRAS[2],
                               idx[0], idx[1], idx[2] );

    source->SetValueAtIndex( idx[0], idx[1], idx[2], mNewValue );

    iView->Render();

  }

  vtkKWScubaLayer3DMRI* layer3DMRI =
    vtkKWScubaLayer3DMRI::SafeDownCast( iLayer );
  if ( layer3DMRI ) {

    vtkFSVolumeSource* source = layer3DMRI->GetSource();

    int idx[3];
    source->ConvertRASToIndex( iRAS[0], iRAS[1], iRAS[2],
                               idx[0], idx[1], idx[2] );

    source->SetValueAtIndex( idx[0], idx[1], idx[2], mNewValue );

    iView->Render();

  }
}
