/**
 * @file  test_VolumeHistogram.cpp
 * @brief test VolumeHistogram class
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2007/12/10 23:32:25 $
 *    $Revision: 1.4 $
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


#include <qapplication.h>
#include <qpushbutton.h>
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "tcl.h"
}
#include "Scuba-impl.h"
#include "QtVolumeHistogram.h"

const char* Progname = "test_VolumeHistogram";

using namespace std;

int main( int argc, char **argv ) {
  try {

    QApplication a( argc, argv );


    string fnMRI = "/Users/kteich/work/subjects/bert/mri/orig";

    char* sSubjectsDir = getenv("SUBJECTS_DIR");

    if ( NULL != sSubjectsDir ) {
      fnMRI = string(sSubjectsDir) + "/bert/mri/orig";
    }

    if ( argc == 2 ) {
      fnMRI = argv[1];
    }

    VolumeCollection vol;
    vol.SetFileName( fnMRI );
    MRI* mri = vol.GetMRI();
    if ( NULL == mri )
      exit( 1 );

    QtVolumeHistogram* histogram;
    histogram = new QtVolumeHistogram( 0, (const char*) "QtVolumeHistogram" );
    histogram->SetVolumeSource( &vol );
    histogram->SetNumberOfBins( 255 );
    histogram->SetMinIgnore( 0 );
    histogram->SetMaxIgnore( 20 );
    histogram->SetNumberOfMarkers( 4 );
    histogram->SetMarkerColor( 0, Qt::red );
    histogram->SetMarkerValue( 0, 10 );
    histogram->SetMarkerColor( 1, Qt::green );
    histogram->SetMarkerValue( 1, 30 );
    histogram->SetMarkerColor( 2, Qt::blue );
    histogram->SetMarkerValue( 2, 50 );
    histogram->SetMarkerColor( 3, Qt::yellow );
    histogram->SetMarkerValue( 3, 70 );

    histogram->resize( 600, 200 );

    a.setMainWidget( histogram );
    histogram->show();

    QApplication::setGlobalMouseTracking( true );

    return a.exec();
  } catch ( runtime_error& e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  } catch ( exception& e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  } catch (...) {
    cerr << "failed" << endl;
    exit( 1 );
  }


  exit( 0 );
}

