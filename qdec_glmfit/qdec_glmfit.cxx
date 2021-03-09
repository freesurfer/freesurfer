/**
 * @brief Wrapper for mri_glmfit: main function and option parsing
 *
 * Parses input to initialize a QdecGlmFitDesign object, runs
 * mri_glmfit, and saves the result as a .qdec file.
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

#include <getopt.h>
#include <stdlib.h> // exit

#include <string>
#include <iostream>
#include <sstream>


#include "QdecProject.h"

using namespace std;

const char* Progname = "qdec_glmfit";

void PrintUsage ();

int main ( int argc, char** argv ) {

  // The input params that we need to read in. The ones that have
  // default values are not required.
  string fnDataTable;
  string sWorkingDir = "/tmp";
  string sSubjectsDir;
  string sSubjectName = "fsaverage";
  string sAnalysisName;
  string sDiscreteFactor1 = "none";
  string sDiscreteFactor2 = "none";
  string sContinuousFactor1 = "none";
  string sContinuousFactor2 = "none";
  string sMeasurement;
  string sHemisphere;
  int smoothness = -1;
  string fnProject;

  // Try to get SUBJECTS_DIR.
  if( getenv("SUBJECTS_DIR") ) {
    sSubjectsDir = getenv("SUBJECTS_DIR");
  }

  // Our arguments to getopt_long.
  struct option aOptions[] = {
    { "data-table",        required_argument, NULL, 'd' },
    { "working-dir",       required_argument, NULL, 'w' },
    { "subjects-dir",      required_argument, NULL, 's' },
    { "average-subject",   required_argument, NULL, 'a' },
    { "analysis-name",     required_argument, NULL, 'n' },
    { "discrete-factor",   required_argument, NULL, 'f' },
    { "continuous-factor", required_argument, NULL, 'c' },
    { "measurement",       required_argument, NULL, 'm' },
    { "hemisphere",        required_argument, NULL, 'h' },
    { "smoothness",        required_argument, NULL, 't' },
    { "output",            required_argument, NULL, 'o' },
    { NULL, 0, NULL, 0 }
  };

  // Parse the command line options.
  int cDiscreteFactor = 0, cContinuousFactor = 0;
  int nOption = 0;
  int rOption = 0;
  while( -1 != (rOption = getopt_long( argc, argv, "d:w:s:a:n:f:c:m:h:t:o:",
                                       aOptions, &nOption )) ) {

    if( 'd' == rOption ) {
      fnDataTable = optarg;
    }
    else if( 'w' == rOption ) {
      sWorkingDir = optarg;
    }
    else if( 's' == rOption ) {
      sSubjectsDir = optarg;
    }
    else if( 'a' == rOption ) {
      sSubjectName = optarg;
    }
    else if( 'n' == rOption ) {
      sAnalysisName = optarg;
    }
    else if( 'f' == rOption ) {

      // If they are giving us a factor, make sure they don't give us
      // more than two.
      if( 0 == cDiscreteFactor )
        sDiscreteFactor1 = optarg;
      else if( 1 == cDiscreteFactor )
        sDiscreteFactor2 = optarg;
      else {
        cerr << "Only two discrete factors allowed; ignoring the rest."
             << endl;
        continue;
      }
      cDiscreteFactor++;
    }
    else if( 'c' == rOption ) {
      if( 0 == cContinuousFactor )
        sContinuousFactor1 = optarg;
      else if( 1 == cContinuousFactor )
        sContinuousFactor2 = optarg;
      else {
        cerr << "Only two continuous factors allowed; ignoring the rest."
             << endl;
        continue;
      }
      cContinuousFactor++;
    }
    else if( 'm' == rOption ) {
      sMeasurement = optarg;
    }
    else if( 'h' == rOption ) {
      sHemisphere = optarg;
    }
    else if( 't' == rOption ) {
      stringstream ssSmoothness;
      ssSmoothness << optarg;
      ssSmoothness >> smoothness;
    }
    else if( 'o' == rOption ) {
      fnProject = optarg;
    }
  }

  // Make sure they set our params.
  if( fnDataTable == "" ) {
    cerr << "Error: no data table specified. Use --data-table <filename>."
         << endl;
    PrintUsage();
    exit( 1 );
  }
  if( sSubjectsDir == "" ) {
    cerr << "Error: no subjects directory specified. Use --subjects-dir "
         << "<directory> or set the SUBJECTS_DIR environment variable."
         << endl;
    PrintUsage();
    exit( 1 );
  }
  if( sAnalysisName == "" ) {
    cerr << "Error: no analysis name specified. Use --analysis-name <name>."
         << endl;
    PrintUsage();
    exit( 1 );
  }
  if( sMeasurement == "" ) {
    cerr << "Error: no measurement specified. Use --measurement <filename>."
         << endl;
    PrintUsage();
    exit( 1 );
  }
  if( sHemisphere == "" ) {
    cerr << "Error: no hemisphere specified. Use --hemisphere lh|rh."
         << endl;
    PrintUsage();
    exit( 1 );
  }
  if( -1 == smoothness ) {
    cerr << "Error: no smoothness specified. Use --smoothness <smoothness>."
         << endl;
    PrintUsage();
    exit( 1 );
  }
  if( fnProject == "" ) {
    cerr << "Error: no .qdec output file specified. Use --output <filename>."
         << endl;
    PrintUsage();
    exit( 1 );
  }

  // Our working dir will be the working dir they gave us plus the
  // name of the analysis.
  sWorkingDir += "/" + sAnalysisName;

  try {

    // Load the data table.
    QdecProject project;
    if( project.LoadDataTable( fnDataTable.c_str() ) ) {
      cerr << "Error: Couldn't load data table " << fnDataTable << endl;
      exit( 1 );
    }

    // Set some values.
    project.SetSubjectsDir( sSubjectsDir.c_str() );
    project.SetAverageSubject( sSubjectName.c_str() );
    project.SetWorkingDir( sWorkingDir.c_str() );

    // Create the design based on our input params.
    if( project.CreateGlmDesign( sAnalysisName.c_str(),
                                 sDiscreteFactor1.c_str(),
                                 sDiscreteFactor2.c_str(),
                                 sContinuousFactor1.c_str(),
                                 sContinuousFactor2.c_str(),
                                 NULL, 0,
                                 sMeasurement.c_str(), sHemisphere.c_str(),
                                 smoothness, NULL ) ) {
      cerr << "Error: Couldn't create design. Make sure your parameters are "
           << "valid, including that your factors exist in the data table and "
           << "are of the right type, and that the subject exists in the "
           << "given subjects directory." << endl;

      cerr << "Input:" << endl;
      cerr << " Data table: " << fnDataTable << endl;
      cerr << " Working dir: " << sWorkingDir << endl;
      cerr << " Subjects dir: " << sSubjectsDir << endl;
      cerr << " Subject name: " << sSubjectName << endl;
      cerr << " Analysis name: " << sAnalysisName << endl;
      cerr << " Discrete factor 1: " << sDiscreteFactor1 << endl;
      cerr << " Discrete factor 2: " << sDiscreteFactor2 << endl;
      cerr << " Continuous factor 1: " << sContinuousFactor1 << endl;
      cerr << " Continuous factor 2: " << sContinuousFactor2 << endl;
      cerr << " Measurement: " << sMeasurement << endl;
      cerr << " Hemisphere: " << sHemisphere << endl;
      cerr << " Smoothness: " << smoothness << endl;
      cerr << " Output: " << fnProject << endl;
      exit( 1 );
    }

    // Run the GLM fit.
    if( project.RunGlmFit() ) {
      cerr << "Error: mri_glmfit did not return successfully. "
           << "Check the output for details." << endl;
      exit( 1 );
    }

    // Save the results.
    if( project.SaveProjectFile( fnProject.c_str() ) ) {
      cerr << "Error: Couldn't save the results to the project file "
           << fnProject << endl;
      exit( 1 );
    }

    // Delete our working directory.
    string sCommand = "rm -rf " + sWorkingDir;
    int rSystem = system( sCommand.c_str() );
    if( 0 != rSystem )
      cerr << "Warning: Couldn't remove temp working directory "
           << sWorkingDir << endl;
  }
  catch( exception& e ) {

    cerr << "Error: " << e.what() << endl;
    exit( 1 );
  }

  return 0;
}

void PrintUsage () {

  cout << "USAGE: qdec_glmfit [options...]" << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << endl;
  cout << "  --data-table, -d <filename>      Input qdec.table.dat file (reqd)"
       << endl << endl;
  cout << "  --working-dir, -w <path>         Directory in which to generate "
       << endl
       << "                                   temporary data (default /tmp)"
       << endl << endl;
  cout << "  --subjects-dir, -s <path>        Directory in which to find the "
       << endl
       << "                                   average subject (default "
       << endl
       << "                                   SUBJECTS_DIR env var)"
       << endl << endl;
  cout << "  --average-subject, -a <string>   Average subject name (reqd)"
       << endl << endl;
  cout << "  --analysis-name, -n <string>     Name for analysis (reqd)"
       << endl << endl;
  cout << "  --discrete-factor, -f <string>   Discrete factor (up to 2)"
       << endl << endl;
  cout << "  --continuous-factor, -c <string> Continuous factor (up to 2)"
       << endl << endl;
  cout << "  --measurement, -m <string>       Measurement name (reqd)"
       << endl << endl;
  cout << "  --hemisphere, -h lh|rh           Hemisphere to use (reqd)"
       << endl << endl;
  cout << "  --smoothness, -t <integer>       Smoothness to use (reqd)"
       << endl << endl;
  cout << "  --output, -o <filename>          Output .qdec filename (reqd)"
       << endl << endl;

}
