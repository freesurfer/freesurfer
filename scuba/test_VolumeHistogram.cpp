/**
 * @file  test_VolumeHistogram.cpp
 * @brief test VolumeHistogram class
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.5 $
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

