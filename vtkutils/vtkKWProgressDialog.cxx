/**
 * @file  vtkKWProgressDialog.cxx
 * @brief A KWWidgets progress dialog
 *
 * A simple KWWidgets dialog box with a progress bar tht listens to
 * update commands.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/02/20 22:18:59 $
 *    $Revision: 1.2 $
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


#include "vtkKWProgressDialog.h"
#include "vtkKWApplication.h"
#include "vtkKWDialog.h"
#include "vtkKWProgressGauge.h"
#include "vtkKWLabel.h"
#include "vtkObjectFactory.h"
#include "vtkAlgorithm.h"

using namespace std;

//vtkStandardNewMacro( vtkKWProgressDialog );
vtkCxxRevisionMacro( vtkKWProgressDialog, "$Revision: 1.2 $" );

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
