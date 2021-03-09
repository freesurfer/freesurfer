/**
 * @brief A KWWidgets progress dialog
 *
 * A simple KWWidgets dialog box with a progress bar tht listens to
 * update commands.
 */
/*
 * Original Author: Kevin Teich
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


#include "vtkKWProgressDialog.h"
#include "vtkKWApplication.h"
#include "vtkKWDialog.h"
#include "vtkKWProgressGauge.h"
#include "vtkKWLabel.h"
#include "vtkObjectFactory.h"
#include "vtkAlgorithm.h"

using namespace std;

//vtkStandardNewMacro( vtkKWProgressDialog );

vtkKWProgressDialog* vtkKWProgressDialog::New () {

  return new vtkKWProgressDialog();
}

vtkKWProgressDialog::vtkKWProgressDialog () :
    mApplication( NULL ),
    mDialog( NULL ),
    mProgressGauge ( NULL ) {

  strncpy( msWindowTitle, "", sizeof(msWindowTitle) );
}

vtkKWProgressDialog::~vtkKWProgressDialog () {}

void
vtkKWProgressDialog::SetApplication ( vtkKWApplication* iApplication ) {

  mApplication = iApplication;
}

void
vtkKWProgressDialog::SetWindowTitle ( const char* isTitle ) {

  strncpy( msWindowTitle, isTitle, sizeof(msWindowTitle) );
}

void
vtkKWProgressDialog::Execute ( vtkObject* iCaller,
			       unsigned long iEvent,
			       void* iCallData ) {

  if ( iEvent == vtkCommand::ProgressEvent ) {

    vtkAlgorithm* algo = (vtkAlgorithm*)iCaller;
    double progress = *(double*)(iCallData);

    if ( 0 == progress ) {

      if ( NULL == mDialog ) {

        mDialog = vtkKWDialog::New ();
        mDialog->SetApplication( mApplication );
        mDialog->Create();
        mDialog->SetTitle( msWindowTitle );
        mDialog->SetSize( 400, 100 );

        vtkKWLabel* label = vtkKWLabel::New();
        label->SetParent( mDialog );
        label->Create();
        label->SetText( algo->GetProgressText() );

        mProgressGauge = vtkKWProgressGauge::New ();
        mProgressGauge->SetParent( mDialog );
        mProgressGauge->Create();

        mDialog->Script( "pack %s %s -side top -expand yes",
                         label->GetWidgetName(),
                         mProgressGauge->GetWidgetName() );

        mDialog->Display();
      }

    } else if ( 1 == progress ) {

      if ( mDialog && mProgressGauge ) {

        mDialog->OK();
        mProgressGauge->Delete();
        mDialog->Delete();

        mDialog = NULL;
        mProgressGauge = NULL;
      }

    } else {

      if ( mDialog && mProgressGauge )
        mProgressGauge->SetValue( progress * 100.0 );

    }
  }
}
