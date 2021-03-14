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

#include <iostream>
#include <sstream>

#include "QdecProject.h"

const char *Progname = "qdec_glmfit";

void PrintUsage();

int main(int argc, char **argv) {

  // The input params that we need to read in. The ones that have
  // default values are not required.
  std::string fnDataTable;
  std::string sWorkingDir = "/tmp";
  std::string sSubjectsDir;
  std::string sSubjectName = "fsaverage";
  std::string sAnalysisName;
  std::string sDiscreteFactor1   = "none";
  std::string sDiscreteFactor2   = "none";
  std::string sContinuousFactor1 = "none";
  std::string sContinuousFactor2 = "none";
  std::string sMeasurement;
  std::string sHemisphere;
  int         smoothness = -1;
  std::string fnProject;

  // Try to get SUBJECTS_DIR.
  if (getenv("SUBJECTS_DIR")) {
    sSubjectsDir = getenv("SUBJECTS_DIR");
  }

  // Our arguments to getopt_long.
  struct option aOptions[] = {
      {"data-table", required_argument, nullptr, 'd'},
      {"working-dir", required_argument, nullptr, 'w'},
      {"subjects-dir", required_argument, nullptr, 's'},
      {"average-subject", required_argument, nullptr, 'a'},
      {"analysis-name", required_argument, nullptr, 'n'},
      {"discrete-factor", required_argument, nullptr, 'f'},
      {"continuous-factor", required_argument, nullptr, 'c'},
      {"measurement", required_argument, nullptr, 'm'},
      {"hemisphere", required_argument, nullptr, 'h'},
      {"smoothness", required_argument, nullptr, 't'},
      {"output", required_argument, nullptr, 'o'},
      {nullptr, 0, nullptr, 0}};

  // Parse the command line options.
  int cDiscreteFactor = 0, cContinuousFactor = 0;
  int nOption = 0;
  int rOption = 0;
  while (-1 != (rOption = getopt_long(argc, argv, "d:w:s:a:n:f:c:m:h:t:o:",
                                      aOptions, &nOption))) {

    if ('d' == rOption) {
      fnDataTable = optarg;
    } else if ('w' == rOption) {
      sWorkingDir = optarg;
    } else if ('s' == rOption) {
      sSubjectsDir = optarg;
    } else if ('a' == rOption) {
      sSubjectName = optarg;
    } else if ('n' == rOption) {
      sAnalysisName = optarg;
    } else if ('f' == rOption) {

      // If they are giving us a factor, make sure they don't give us
      // more than two.
      if (0 == cDiscreteFactor)
        sDiscreteFactor1 = optarg;
      else if (1 == cDiscreteFactor)
        sDiscreteFactor2 = optarg;
      else {
        std::cerr << "Only two discrete factors allowed; ignoring the rest."
                  << std::endl;
        continue;
      }
      cDiscreteFactor++;
    } else if ('c' == rOption) {
      if (0 == cContinuousFactor)
        sContinuousFactor1 = optarg;
      else if (1 == cContinuousFactor)
        sContinuousFactor2 = optarg;
      else {
        std::cerr << "Only two continuous factors allowed; ignoring the rest."
                  << std::endl;
        continue;
      }
      cContinuousFactor++;
    } else if ('m' == rOption) {
      sMeasurement = optarg;
    } else if ('h' == rOption) {
      sHemisphere = optarg;
    } else if ('t' == rOption) {
      std::stringstream ssSmoothness;
      ssSmoothness << optarg;
      ssSmoothness >> smoothness;
    } else if ('o' == rOption) {
      fnProject = optarg;
    }
  }

  // Make sure they set our params.
  if (fnDataTable == "") {
    std::cerr << "Error: no data table specified. Use --data-table <filename>."
              << std::endl;
    PrintUsage();
    exit(1);
  }
  if (sSubjectsDir == "") {
    std::cerr << "Error: no subjects directory specified. Use --subjects-dir "
              << "<directory> or set the SUBJECTS_DIR environment variable."
              << std::endl;
    PrintUsage();
    exit(1);
  }
  if (sAnalysisName == "") {
    std::cerr
        << "Error: no analysis name specified. Use --analysis-name <name>."
        << std::endl;
    PrintUsage();
    exit(1);
  }
  if (sMeasurement == "") {
    std::cerr
        << "Error: no measurement specified. Use --measurement <filename>."
        << std::endl;
    PrintUsage();
    exit(1);
  }
  if (sHemisphere == "") {
    std::cerr << "Error: no hemisphere specified. Use --hemisphere lh|rh."
              << std::endl;
    PrintUsage();
    exit(1);
  }
  if (-1 == smoothness) {
    std::cerr
        << "Error: no smoothness specified. Use --smoothness <smoothness>."
        << std::endl;
    PrintUsage();
    exit(1);
  }
  if (fnProject == "") {
    std::cerr
        << "Error: no .qdec output file specified. Use --output <filename>."
        << std::endl;
    PrintUsage();
    exit(1);
  }

  // Our working dir will be the working dir they gave us plus the
  // name of the analysis.
  sWorkingDir += "/" + sAnalysisName;

  try {

    // Load the data table.
    QdecProject project;
    if (project.LoadDataTable(fnDataTable.c_str())) {
      std::cerr << "Error: Couldn't load data table " << fnDataTable
                << std::endl;
      exit(1);
    }

    // Set some values.
    project.SetSubjectsDir(sSubjectsDir.c_str());
    project.SetAverageSubject(sSubjectName.c_str());
    project.SetWorkingDir(sWorkingDir.c_str());

    // Create the design based on our input params.
    if (project.CreateGlmDesign(
            sAnalysisName.c_str(), sDiscreteFactor1.c_str(),
            sDiscreteFactor2.c_str(), sContinuousFactor1.c_str(),
            sContinuousFactor2.c_str(), nullptr, 0, sMeasurement.c_str(),
            sHemisphere.c_str(), smoothness, nullptr)) {
      std::cerr
          << "Error: Couldn't create design. Make sure your parameters are "
          << "valid, including that your factors exist in the data table and "
          << "are of the right type, and that the subject exists in the "
          << "given subjects directory." << std::endl;

      std::cerr << "Input:" << std::endl;
      std::cerr << " Data table: " << fnDataTable << std::endl;
      std::cerr << " Working dir: " << sWorkingDir << std::endl;
      std::cerr << " Subjects dir: " << sSubjectsDir << std::endl;
      std::cerr << " Subject name: " << sSubjectName << std::endl;
      std::cerr << " Analysis name: " << sAnalysisName << std::endl;
      std::cerr << " Discrete factor 1: " << sDiscreteFactor1 << std::endl;
      std::cerr << " Discrete factor 2: " << sDiscreteFactor2 << std::endl;
      std::cerr << " Continuous factor 1: " << sContinuousFactor1 << std::endl;
      std::cerr << " Continuous factor 2: " << sContinuousFactor2 << std::endl;
      std::cerr << " Measurement: " << sMeasurement << std::endl;
      std::cerr << " Hemisphere: " << sHemisphere << std::endl;
      std::cerr << " Smoothness: " << smoothness << std::endl;
      std::cerr << " Output: " << fnProject << std::endl;
      exit(1);
    }

    // Run the GLM fit.
    if (project.RunGlmFit()) {
      std::cerr << "Error: mri_glmfit did not return successfully. "
                << "Check the output for details." << std::endl;
      exit(1);
    }

    // Save the results.
    if (project.SaveProjectFile(fnProject.c_str())) {
      std::cerr << "Error: Couldn't save the results to the project file "
                << fnProject << std::endl;
      exit(1);
    }

    // Delete our working directory.
    std::string sCommand = "rm -rf " + sWorkingDir;
    int         rSystem  = system(sCommand.c_str());
    if (0 != rSystem)
      std::cerr << "Warning: Couldn't remove temp working directory "
                << sWorkingDir << std::endl;
  } catch (std::exception &e) {

    std::cerr << "Error: " << e.what() << std::endl;
    exit(1);
  }

  return 0;
}

void PrintUsage() {

  std::cout << "USAGE: qdec_glmfit [options...]" << std::endl;
  std::cout << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << std::endl;
  std::cout
      << "  --data-table, -d <filename>      Input qdec.table.dat file (reqd)"
      << std::endl
      << std::endl;
  std::cout
      << "  --working-dir, -w <path>         Directory in which to generate "
      << std::endl
      << "                                   temporary data (default /tmp)"
      << std::endl
      << std::endl;
  std::cout
      << "  --subjects-dir, -s <path>        Directory in which to find the "
      << std::endl
      << "                                   average subject (default "
      << std::endl
      << "                                   SUBJECTS_DIR env var)" << std::endl
      << std::endl;
  std::cout << "  --average-subject, -a <string>   Average subject name (reqd)"
            << std::endl
            << std::endl;
  std::cout << "  --analysis-name, -n <string>     Name for analysis (reqd)"
            << std::endl
            << std::endl;
  std::cout << "  --discrete-factor, -f <string>   Discrete factor (up to 2)"
            << std::endl
            << std::endl;
  std::cout << "  --continuous-factor, -c <string> Continuous factor (up to 2)"
            << std::endl
            << std::endl;
  std::cout << "  --measurement, -m <string>       Measurement name (reqd)"
            << std::endl
            << std::endl;
  std::cout << "  --hemisphere, -h lh|rh           Hemisphere to use (reqd)"
            << std::endl
            << std::endl;
  std::cout << "  --smoothness, -t <integer>       Smoothness to use (reqd)"
            << std::endl
            << std::endl;
  std::cout << "  --output, -o <filename>          Output .qdec filename (reqd)"
            << std::endl
            << std::endl;
}
