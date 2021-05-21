/**
 * @brief API class containing all qdec subject data and methods
 *
 * Top-level interface class containing all data associated with a users
 * subject group, and potentially mri_glmfit processed data associated with
 * that group.
 */
/*
 * Original Author: Nick Schmansky
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include <stdexcept>
#include <sstream>
#include <iostream>
#include <fstream>

#include "QdecProject.h"


// Constructors/Destructors
//

QdecProject::QdecProject ( ) :
  mDataTable( new QdecDataTable() ),
  mGlmDesign( new QdecGlmDesign( this->mDataTable ) ),
  mGlmFitter( new QdecGlmFit() ),
  msZipCommandFormat( "cd %3; zip -r %1 %2 > /dev/null" ),
  msUnzipCommandFormat( "unzip -o -d %3 %1 > /dev/null" )
{
  msStatsDataTablesDir = mGlmDesign->GetDefaultWorkingDir() + "/stats_tables/";
}

QdecProject::~QdecProject ( )
{
  delete this->mDataTable;
  delete this->mGlmDesign;
  delete this->mGlmFitter;
}

//
// Methods
//

/**
 * Load a .qdec project file (containing all necessary info to begin
 * working either on a new project, or to continue working using
 * results from a prior saved work session). isDataDir should be a
 * directory where we can expand the .qdec file (like /tmp).
 * @return int
 * @param  isFileName
 * @param  isDataDir
 */
int QdecProject::LoadProjectFile ( const char* ifnProject,
                                   const char* ifnDataDir )
{

  // If the project file name doesn't have a path, give it one.
  string fnProject( ifnProject );
  if ( fnProject.find( '/' ) == string::npos )
  {
    char sCWD[1024];
    if ( getcwd( sCWD, sizeof(sCWD) ) )
    {
      string fnProjectFull = string(sCWD) + "/" + fnProject;
      fnProject = fnProjectFull;
    }
    else
    {
      fprintf(stderr, "WARNING: QdecProject::LoadProjectFile: Can't add "
              "full path  to project file name; please specify full path." );
    }
  }

  // Find the base name of the project file.
  string fnProjectBase( fnProject );
  string::size_type nPreLastSlash = fnProject.rfind( '/' );
  if ( string::npos != nPreLastSlash )
    fnProjectBase = fnProject.substr( nPreLastSlash+1, fnProject.size() );

  // Make a target dir for the expanded file in the data dir, with a
  // directory name of the project file.
  string fnExpandedProjectBase = "qdec_project_archive";
  string fnExpandedProjectDir = string(ifnDataDir) + "/" +
    fnExpandedProjectBase;

  string sSubject;
  string sHemisphere;
  string sAnalysisName;
  string fnDataTableBase;
  string sDiscreteFactor1 = "none";
  string sDiscreteFactor2 = "none";
  string sContinuousFactor1 = "none";
  string sContinuousFactor2 = "none";
  string sMeasure;
  int smoothness = -1;

  // Check the file.
  ifstream fInput( fnProject.c_str(), std::ios::in );
  if ( !fInput || fInput.bad() )
    throw runtime_error( string("Couldn't open file " ) + fnProject );
  fInput.close();

  // Erase old working directory if present.
  string sCommand = "rm -rf " + fnExpandedProjectDir;
  int rSystem = system( sCommand.c_str() );
  if ( 0 != rSystem )
  {
    fprintf( stderr, "ERROR: QdecProject::LoadProjectFile: Couldn't "
             "remove existing temp directory (cmd=%s)\n", sCommand.c_str() );
    return -1;
  }

  // Get out command string and expand the .qdec file into the
  // destination directory.
  this->FormatCommandString( fnProject.c_str(),
                             fnExpandedProjectBase.c_str(),
                             ifnDataDir,
                             msUnzipCommandFormat.c_str(),
                             sCommand );
  rSystem = system( sCommand.c_str() );
  if ( 0 != rSystem )
  {
    fprintf( stderr, "ERROR: QdecProject::LoadProjectFile: Couldn't "
             "expand project file (cmd=%s)\n", sCommand.c_str() );
    return -1;
  }

  // Look for and check the version file.
  string fnVersion = fnExpandedProjectDir + "/Version.txt";
  ifstream fVersion( fnVersion.c_str(), ios::out );
  if ( !fVersion || fVersion.bad() )
  {
    fprintf( stderr, "ERROR: QdecProject::LoadProjectFile: Couldn't "
             "find Version file %s\n", fnVersion.c_str() );
    return -1;
  }
  int version;
  fVersion >> version;
  fVersion.close();
  if ( 1 != version )
  {
    fprintf( stderr, "ERROR: QdecProject::LoadProjectFile: Version "
             "file had wrong value (%d)\n", version );
    return -1;
  }

  // Parse the meta data file.
  string fnMetadata = fnExpandedProjectDir + "/" + this->GetMetadataFileName();
  ifstream fMetadata( fnMetadata.c_str(), ios::in );
  if ( !fMetadata || fMetadata.bad() )
  {
    fprintf( stderr, "ERROR: QdecProject::LoadProjectFile: Couldn't "
             "open metadata file %s\n", fnMetadata.c_str() );
    return -1;
  }
  // Make sure the first token is QdecProjectMetadata, and then the
  // next line is Version 1.
  string sToken;
  string asCorrectTokens[] = { "QdecProjectMetadata", "Version", "1" };
  for ( int nToken = 0; nToken < 3; nToken++ )
  {
    fMetadata >> sToken;
    if ( sToken != asCorrectTokens[nToken] )
    {
      fprintf( stderr, "ERROR: QdecProject::LoadProjectFile: Invalid "
               "metadata file %s, %s token not found\n",
               fnMetadata.c_str(), asCorrectTokens[nToken].c_str() );
      return -1;
    }
  }
  // Now we parse the file and look for names and values.
  while ( fMetadata >> sToken && !fMetadata.eof() )
  {
    if ( sToken == "Subject" ) fMetadata >> sSubject;
    else if ( sToken == "Hemisphere" ) fMetadata >> sHemisphere;
    else if ( sToken == "AnalysisName" ) fMetadata >> sAnalysisName;
    else if ( sToken == "DataTable" ) fMetadata >> fnDataTableBase;
    else if ( sToken == "Measure" ) fMetadata >> sMeasure;
    else if ( sToken == "Smoothness" ) fMetadata >> smoothness;
    else if ( sToken == "DiscreteFactor1" ) fMetadata >> sDiscreteFactor1;
    else if ( sToken == "DiscreteFactor2" ) fMetadata >> sDiscreteFactor2;
    else if ( sToken == "ContinuousFactor1" )fMetadata >> sContinuousFactor1;
    else if ( sToken == "ContinuousFactor2" )fMetadata >> sContinuousFactor2;
    else
    {
      fprintf( stderr, "WARNING: QdecProject::LoadProjectFile: Unrecognized "
               "token in QdecProjectMetadata: %s\n", sToken.c_str() );
    }
  }

  // Make sure we got some decent results.
  if ( sSubject == "" )
  {
    fprintf( stderr, "ERROR: QdecProject::LoadProjectFile: Invalid "
             "project metadata file, Subject value not found\n" );
    return -1;
  }
  if ( sHemisphere == "" )
  {
    fprintf( stderr, "ERROR: QdecProject::LoadProjectFile: Invalid "
             "project metadata file, Hemisphere value not found\n" );
    return -1;
  }
  if ( sAnalysisName == "" )
  {
    fprintf( stderr, "ERROR: QdecProject::LoadProjectFile: Invalid "
             "project metadata file, AnalysisName value not found\n" );
    return -1;
  }
  if ( fnDataTableBase == "" )
  {
    fprintf( stderr, "ERROR: QdecProject::LoadProjectFile: Invalid "
             "project metadata file, DataTable value not found\n" );
    return -1;
  }
  if ( sMeasure == "" )
  {
    fprintf( stderr, "ERROR: QdecProject::LoadProjectFile: Invalid "
             "project metadata file, Measure value not found\n" );
    return -1;
  }
  if ( -1 == smoothness )
  {
    fprintf( stderr, "ERROR: QdecProject::LoadProjectFile: Invalid "
             "project metadata file, Smoothness value not found\n" );
    return -1;
  }

  // Load our data table. Note that this might set the subjects dir,
  // but we'll set it later to our data dir.
  string fnDataTable = fnExpandedProjectDir + "/" + fnDataTableBase;
  int errorCode;
  errorCode = this->LoadDataTable( fnDataTable.c_str() );
  if ( errorCode )
    return errorCode;

  // Set the subjects dir to the data dir, so that we can find the
  // subject.
  this->SetSubjectsDir( fnExpandedProjectDir.c_str() );

  // Set the working dir to isDataDir/sAnalysisName.
  string fnWorkingDir = fnExpandedProjectDir + "/" + sAnalysisName;
  this->SetWorkingDir( fnWorkingDir.c_str() );

  // We're generating design and results here so that we can access it
  // as metadata, but we're not actually computing any new results;
  // those all exist in our data dir.

  // Create the design. This will be used in the results.
  errorCode =
    this->mGlmDesign->Create ( this->mDataTable,
                               sAnalysisName.c_str(),
                               sDiscreteFactor1.c_str(),
                               sDiscreteFactor2.c_str(),
                               sContinuousFactor1.c_str(),
                               sContinuousFactor2.c_str(),
                               NULL, 0,
                               sMeasure.c_str(),
                               sHemisphere.c_str(),
                               smoothness,
                               NULL );
  if ( errorCode )
    return errorCode;

  // Create fit results data.
  delete this->mGlmFitter;
  this->mGlmFitter = new QdecGlmFit();
  errorCode =
    mGlmFitter->CreateResultsFromCachedData ( this->mGlmDesign );
  if ( errorCode )
    return errorCode;

  return 0;
}


/**
 * Save all necessary information pertaining to this project (all subject
 * data, any results, any user preferences).
 * @return int
 * @param  isFileName
 */
int QdecProject::SaveProjectFile ( const char* ifnProject,
                                   const char* ifnDataDir )
{

  cout << "Saving project file...\n";

  // If the project file name doesn't have a path, give it one.
  string fnProject( ifnProject );
  if ( fnProject.find( '/' ) == string::npos )
  {
    char sCWD[1024];
    if ( getcwd( sCWD, sizeof(sCWD) ) )
    {
      string fnProjectFull = string(sCWD) + "/" + fnProject;
      fnProject = fnProjectFull;
    }
    else
    {
      fprintf(stderr, "WARNING: QdecProject::LoadProjectFile: Can't add "
              "full path  to project file name; please specify full path." );
    }
  }

  // If the file name doesn't end in qdec, append it now.
  if ( fnProject.find( ".qdec" ) != fnProject.size() - 5 )
  {
    fnProject += ".qdec";
  }

  /* To make our file, we create a temp directory, link in our files,
     and then compress into the destination .qdec file. This is the
     structure we want.

     $ifnWorkingDir/qdec_project_archive/
     $Subject/surf/{r,l}h.{curv,inflatd,pial,white}
     label/{r,l}h.aparc.annot
     $AnalysisName/ *
     qdec.table.dat
     QdecProjectMetadata.txt
     Version.txt
  */

  string fnSubjectsDir = this->GetSubjectsDir();
  string sSubjectName = this->GetAverageSubject();
  string fnWorkingDir = this->GetWorkingDir();
  string fnDefaultWorkingDir = this->GetDefaultWorkingDir();

  // Find the base name of the project file.
  string fnProjectBase( fnProject );
  string::size_type nPreLastSlash = fnProject.rfind( '/' );
  if ( string::npos != nPreLastSlash )
    fnProjectBase = fnProject.substr( nPreLastSlash+1, fnProject.size() );

  // Make a target dir for the expanded file in the data dir, with a
  // directory name of the project file.
  string fnExpandedProjectBase = "qdec_project_archive";
  string fnExpandedProjectDir = string(ifnDataDir) + "/" +
    fnExpandedProjectBase;

  // Erase old working directory if present.
  string sCommand = "rm -rf " + fnExpandedProjectDir;
  cout << sCommand << endl;
  int rSystem = system( sCommand.c_str() );
  if ( 0 != rSystem )
  {
    fprintf( stderr, "ERROR: QdecProject::SaveProjectFile: Couldn't "
             "remove existing temp directory (cmd=%s)\n", sCommand.c_str() );
    return -1;
  }

  // Make a temporary directory for our data.
  sCommand = "mkdir " + fnExpandedProjectDir;
  cout << sCommand << endl;
  rSystem = system( sCommand.c_str() );
  if ( 0 != rSystem )
  {
    fprintf( stderr, "ERROR: QdecProject::SaveProjectFile: Couldn't "
             "make temp directory (cmd=%s)\n", sCommand.c_str() );
    return -1;
  }

  // Write a version file to it.
  string fnVersion = fnExpandedProjectDir + "/Version.txt";
  ofstream fVersion( fnVersion.c_str(), ios::out );
  fVersion << "1" << endl;
  fVersion.close();

  // Make the average subject dir structure.
  sCommand = "mkdir -p " +
    fnExpandedProjectDir + "/" + sSubjectName + "/mri/transforms " +
    fnExpandedProjectDir + "/" + sSubjectName + "/surf " +
    fnExpandedProjectDir + "/" + sSubjectName + "/label";
  cout << sCommand << endl;
  rSystem = system( sCommand.c_str() );
  if ( 0 != rSystem )
  {
    fprintf( stderr, "ERROR: QdecProject::SaveProjectFile: Couldn't "
             "make subject dir structure (cmd=%s)\n", sCommand.c_str() );
    return -1;
  }

  // Link the necessary files. Start with mri/transforms/talairach.xfm file
  sCommand = "ln -s " +
    fnSubjectsDir + "/" + sSubjectName + 
    "/mri/transforms/talairach.xfm " +
    fnExpandedProjectDir + "/" + sSubjectName + "/mri/transforms";
  cout << sCommand << endl;
  rSystem = system( sCommand.c_str() );
  if ( 0 != rSystem )
  {
    fprintf( stderr, "ERROR: QdecProject::SaveProjectFile: Couldn't "
             "link talairach.xfm file (cmd=%s)\n", sCommand.c_str() );
    return -1;
  }

  // orig.mgz
  sCommand = "ln -s " +
    fnSubjectsDir + "/" + sSubjectName + "/mri/orig.mgz " +
    fnExpandedProjectDir + "/" + sSubjectName + "/mri";
  cout << sCommand << endl;
  rSystem = system( sCommand.c_str() );
  if ( 0 != rSystem )
  {
    fprintf( stderr, "ERROR: QdecProject::SaveProjectFile: Couldn't "
             "link orig.mgz file (cmd=%s)\n", sCommand.c_str() );
    return -1;
  }

  // Surfaces
  sCommand = "ln -s " +
    fnSubjectsDir + "/" + sSubjectName + "/surf/*.curv " +
    fnSubjectsDir + "/" + sSubjectName + "/surf/*.inflated " +
    fnSubjectsDir + "/" + sSubjectName + "/surf/*.pial " +
    fnSubjectsDir + "/" + sSubjectName + "/surf/*.white " +
    fnExpandedProjectDir + "/" + sSubjectName + "/surf";
  cout << sCommand << endl;
  rSystem = system( sCommand.c_str() );
  if ( 0 != rSystem )
  {
    fprintf( stderr, "ERROR: QdecProject::SaveProjectFile: Couldn't "
             "link surface files (cmd=%s)\n", sCommand.c_str() );
    return -1;
  }

  // Annotations.
  sCommand = "ln -s " +
    fnSubjectsDir + "/" + sSubjectName + "/label/*.aparc.annot " +
    fnExpandedProjectDir + "/" + sSubjectName + "/label";
  cout << sCommand << endl;
  rSystem = system( sCommand.c_str() );
  if ( 0 != rSystem )
  {
    fprintf( stderr, "ERROR: QdecProject::SaveProjectFile: Couldn't "
             "link annotation files (cmd=%s)\n", sCommand.c_str() );
    return -1;
  }

  // The whole working dir.
  sCommand = "ln -s " + fnWorkingDir + " " + fnExpandedProjectDir;
  cout << sCommand << endl;
  rSystem = system( sCommand.c_str() );
  if ( 0 != rSystem )
  {
    fprintf( stderr, "ERROR: QdecProject::SaveProjectFile: Couldn't "
             "link working dir (cmd=%s)\n", sCommand.c_str() );
    return -1;
  }

  // Data table and levels files.
  string fnDataTable = this->GetDataTable()->GetFileName();
  string fnDataTablePath( fnDataTable );
  string fnDataTableBase( fnDataTable );
  nPreLastSlash = fnDataTable.rfind( '/' );
  if ( string::npos != nPreLastSlash )
  {
    fnDataTableBase = fnDataTable.substr( nPreLastSlash+1, fnDataTable.size());
    fnDataTablePath = fnDataTable.substr( 0, nPreLastSlash+1);
  }

  // NOTE: Older versions of ln don't handle multiple files specified
  // in a symbolic link command properly. For example, on kani, with ln
  // version 4.5.3, this will create a dead qdec.table.dat link in /tmp:
  // ln -s $PWD/qdec.table.dat *.levels /tmp
  // so the files are linked/copied singly
  sCommand = "ln -s " + fnDefaultWorkingDir + "/" + fnDataTableBase +
    " " + fnExpandedProjectDir;
  cout << sCommand << endl;
  rSystem = system( sCommand.c_str() );
  if ( 0 != rSystem )
  {
    fprintf( stderr, "ERROR: QdecProject::SaveProjectFile: Couldn't "
             "link data table (cmd=%s)\n", sCommand.c_str() );
    return -1;
  }
  sCommand = "cp " + fnDefaultWorkingDir + "/*.levels " + fnExpandedProjectDir;
  cout << sCommand << endl;
  rSystem = system( sCommand.c_str() );
  /* .levels files may not exist, so don't check for copy error:
     if( 0 != rSystem ) {
     fprintf( stderr, "ERROR: QdecProject::SaveProjectFile: Couldn't "
     "copy *.levels (cmd=%s)\n", sCommand.c_str() );
     return -1;
     }
  */

  // Generate the meta data file.
  string fnMetadata = fnExpandedProjectDir + "/" + this->GetMetadataFileName();
  ofstream fMetadata( fnMetadata.c_str(), ios::out );
  if ( !fMetadata || fMetadata.bad() )
  {
    fprintf( stderr, "ERROR: QdecProject::SaveProjectFile: Couldn't "
             "make metadata file %s\n", fnMetadata.c_str() );
    return -1;
  }
  fMetadata << "QdecProjectMetadata" << endl;
  fMetadata << "Version 1" << endl;
  fMetadata << "Subject " << this->GetAverageSubject() << endl;
  fMetadata << "Hemisphere " << this->GetHemi() << endl;
  fMetadata << "AnalysisName " << this->GetGlmDesign()->GetName() << endl;
  fMetadata << "DataTable " << fnDataTableBase << endl;
  fMetadata << "Measure " << this->GetGlmDesign()->GetMeasure() << endl;
  fMetadata << "Smoothness " << this->GetGlmDesign()->GetSmoothness() << endl;

  // We only support the two factors of each kind, so get the vectors
  // and just write the first and second ones if they are present.
  vector<QdecFactor*> const& lDiscreteFactors =
    this->GetGlmDesign()->GetDiscreteFactors();
  if ( lDiscreteFactors.size() > 0 )
    fMetadata << "DiscreteFactor1 "
              << lDiscreteFactors[0]->GetFactorName() << endl;
  if ( lDiscreteFactors.size() > 1 )
    fMetadata << "DiscreteFactor2 "
              << lDiscreteFactors[1]->GetFactorName() << endl;

  vector<QdecFactor*> const& lContinuousFactors =
    this->GetGlmDesign()->GetContinuousFactors();
  if ( lContinuousFactors.size() > 0 )
    fMetadata << "ContinuousFactor1 "
              << lContinuousFactors[0]->GetFactorName() << endl;
  if ( lContinuousFactors.size() > 1 )
    fMetadata << "ContinuousFactor2 "
              << lContinuousFactors[1]->GetFactorName() << endl;

  fMetadata.close();

  // Get our command string and compress the directory to the
  // destination location with the .qdec filename.
  this->FormatCommandString( fnProject.c_str(),
                             fnExpandedProjectBase.c_str(),
                             ifnDataDir,
                             msZipCommandFormat.c_str(),
                             sCommand );
  cout << sCommand << endl;
  rSystem = system( sCommand.c_str() );
  if ( 0 != rSystem )
  {
    fprintf( stderr, "ERROR: QdecProject::SaveProjectFile: Couldn't "
             "compress project table (cmd=%s)\n", sCommand.c_str() );
    return -1;
  }

  // Delete the temp directory.
  sCommand = "rm -rf " + fnExpandedProjectDir;
  cout << sCommand << endl;
  rSystem = system( sCommand.c_str() );
  if ( 0 != rSystem )
  {
    fprintf( stderr, "ERROR: QdecProject::SaveProjectFile: Couldn't "
             "remove temp directory (cmd=%s)\n", sCommand.c_str() );
    return -1;
  }

  cout << "Saving project file done.\n";

  return 0;
}

int QdecProject::SetZipCommandFormat ( const char* isFormat )
{

  if ( NULL == isFormat )
    return -1;

  msZipCommandFormat = isFormat;

  return 0;
}

int QdecProject::SetUnzipCommandFormat ( const char* isFormat )
{

  if ( NULL == isFormat )
    return -1;

  msUnzipCommandFormat = isFormat;

  return 0;
}


/**
 * @return int
 * @param  isFileName
 */
int QdecProject::LoadDataTable ( const char* isFileName )
{
  char subjectsDir[3000];
  if ( this->mDataTable ) delete this->mDataTable;
  this->mDataTable = new QdecDataTable();
  int ret = 0;
  try
  {
    ret = this->mDataTable->Load ( isFileName, subjectsDir );
  }
  catch ( exception& e )
  {
    cerr << e.what() << endl;
    exit(1); // shutdown the whole shootin' match
  }
  if ( ret ) return ret;
  if ( strlen(subjectsDir) > 0 ) ret = this->SetSubjectsDir ( subjectsDir );
  if ( ret ) return ret;
  delete this->mGlmDesign;
  this->mGlmDesign = new QdecGlmDesign( this->mDataTable );
  this->VerifySubjects();
  //if ( ret )
  //{
  //  delete this->mDataTable;
  //  this->mDataTable = new QdecDataTable(); // on err, return empty table
  //}
  return ret;
}


/**
 * Check that all subjects exist in the specified subjects_dir (including the
 * specified average subject).  Print to stderr and ErrorMessage any errors
 * found (one message for each error).  Also check that thickness, sulc, curv,
 * area and jacobian_white files exist, and that their vertex numbers equal
 * their inflated surface (and that surfaces all have the same number of
 * vertices).
 * @return int
 */
int QdecProject::VerifySubjects ( )
{
  fprintf(stdout,"Verifying subject data");
  int errs=0;
  vector< QdecSubject* > theSubjects = this->mDataTable->GetSubjects();
  for (unsigned int i=0; i < theSubjects.size(); i++)
  {
    fprintf(stdout,".");
    fflush(stdout);
    string id = theSubjects[i]->GetId();
    string sCommand = "ls " +
      this->GetSubjectsDir() + "/" + id + " >& /dev/null";
    int rSystem = system( sCommand.c_str() );
    if ( 0 != rSystem )
    {
      fprintf( stderr, "ERROR: QdecProject::VerifySubjects: Couldn't "
               "find subject '%s' in SUBJECTS_DIR\n", id.c_str() );
      errs++;
      if ( errs > 30 )
      {
        fprintf( stderr, "Too many subjects not found!\n" );
        break;
      }
    }
  }
  if ( 0 == errs )
  {
    fprintf(stdout,"Subject verification complete.\n");
  }

  return errs;
}


/**
 * @return bool  true if a table has been loaded
 */
bool QdecProject::HaveDataTable ( )
{
  if ( this->mDataTable &&
       (this->mDataTable->GetNumberOfSubjects() > 0) ) return true;

  return false;
}


/**
 * @return void
 * @param  iFilePointer
 */
void QdecProject::DumpDataTable ( FILE* iFilePointer )
{
  return this->mDataTable->Dump ( iFilePointer );
}

/**
 * @return int
 * @param  isFileName
 */
int QdecProject::SaveDataTable ( const char* isFileName )
{
  return this->mDataTable->Save ( isFileName );
}


/**
 * @return QdecDataTable*
 */
QdecDataTable* QdecProject::GetDataTable ( )
{
  return this->mDataTable;
}


/**
 * @return string
 */
string QdecProject::GetSubjectsDir ( )
{
  return this->mGlmDesign->GetSubjectsDir();
}


/**
 * @param  ifnSubjectsDir
 */
int QdecProject::SetSubjectsDir ( const char* ifnSubjectsDir )
{
  return this->mGlmDesign->SetSubjectsDir( ifnSubjectsDir );
}


/**
 * @return string
 */
string QdecProject::GetAverageSubject ( )
{
  return this->mGlmDesign->GetAverageSubject();
}


/**
 * @param  isSubjectName
 */
void QdecProject::SetAverageSubject ( const char* isSubjectName )
{
  this->mGlmDesign->SetAverageSubject( isSubjectName );
}


/**
 * @return string
 */
string QdecProject::GetDefaultWorkingDir ( )
{
  return this->mGlmDesign->GetDefaultWorkingDir();
}


/**
 * @return string
 */
string QdecProject::GetWorkingDir ( )
{
  return this->mGlmDesign->GetWorkingDir();
}


/**
 * @return 0 if ok, 1 on error
 * @param  isWorkingDir
 */
int QdecProject::SetWorkingDir ( const char* isWorkingDir )
{
  return this->mGlmDesign->SetWorkingDir( isWorkingDir );
}


/**
 * @return vector< string >
 */
vector< string > QdecProject::GetSubjectIDs ( )
{
  return this->mDataTable->GetSubjectIDs();
}


/**
 * @return vector< string >
 */
vector< string > QdecProject::GetDiscreteFactorNames ( )
{
  return this->mDataTable->GetDiscreteFactorNames();
}


/**
 * @return vector< string >
 */
vector< string > QdecProject::GetContinousFactorNames ( )
{
  return this->mDataTable->GetContinuousFactorNames();
}


/**
 * @return string
 */
string QdecProject::GetHemi ( )
{
  return this->mGlmDesign->GetHemi();
}


/**
 * From the given design parameters, this creates the input data required by
 * mri_glmfit:
 *  - the 'y' data (concatenated subject volumes)
 *  - the FSGD file
 *  - the contrast vectors, as .mat files
 * and writes this data to the specified working directory.
 * @return int
 * @param  isName
 * @param  isFirstDiscreteFactor
 * @param  isSecondDiscreteFactor
 * @param  isFirstContinuousFactor
 * @param  isSecondContinuousFactor
 * @param  isNuisanceFactors
 * @parma  iNumNuisanceFactors
 * @param  isMeasure
 * @param  isHemi
 * @param  iSmoothnessLevel
 * @param  iProgressUpdateGUI
 */
int QdecProject::CreateGlmDesign ( const char* isName,
                                   const char* isFirstDiscreteFactor,
                                   const char* isSecondDiscreteFactor,
                                   const char* isFirstContinuousFactor,
                                   const char* isSecondContinuousFactor,
                                   const char** isNuisanceFactors,
                                   int iNumNuisanceFactors,
                                   const char* isMeasure,
                                   const char* isHemi,
                                   int iSmoothnessLevel,
                                   ProgressUpdateGUI* iProgressUpdateGUI )
{

  int errorCode;
  errorCode = this->mGlmDesign->Create ( this->mDataTable,
                                         isName,
                                         isFirstDiscreteFactor,
                                         isSecondDiscreteFactor,
                                         isFirstContinuousFactor,
                                         isSecondContinuousFactor,
                                         isNuisanceFactors,
                                         iNumNuisanceFactors,
                                         isMeasure,
                                         isHemi,
                                         iSmoothnessLevel,
                                         iProgressUpdateGUI );
  if ( errorCode )
    return errorCode;

  if ( iProgressUpdateGUI )
  {
    iProgressUpdateGUI->BeginActionWithProgress("Writing input files..." );
  }


  if ( mGlmDesign->WriteFsgdFile() )
  {
    fprintf( stderr, "ERROR: QdecProject::CreateGlmDesign: "
             "could not create fsgd file\n");
    return(-3);
  }

  if ( mGlmDesign->WriteContrastMatrices() )
  {
    fprintf( stderr, "ERROR: QdecProject::CreateGlmDesign: could not "
             "generate contrasts\n");
    return(-4);
  }

  if ( mGlmDesign->WriteYdataFile() )
  {
    fprintf( stderr, "ERROR: QdecProject::CreateGlmDesign: could not "
             "create y.mgh file\n");
    return(-4);
  }

  if ( iProgressUpdateGUI )
  {
    iProgressUpdateGUI->EndActionWithProgress();
  }

  return 0;
}


/**
 * From the given design parameters, this creates the input data required by
 * mri_glmfit:
 *  - the 'y' data (data points stuffed into a volume)
 *  - the FSGD file
 *  - the contrast vectors, as .mat files
 * and writes this data to the specified working directory.
 * @return int
 * @param  isName
 * @param  isFirstDiscreteFactor
 * @param  isSecondDiscreteFactor
 * @param  isFirstContinuousFactor
 * @param  isSecondContinuousFactor
 * @param  isNuisanceFactors
 * @parma  iNumNuisanceFactors
 * @param  isMeasure
 * @param  iProgressUpdateGUI
 */
int QdecProject::CreateGlmDesign ( const char* isName,
                                   const char* isFirstDiscreteFactor,
                                   const char* isSecondDiscreteFactor,
                                   const char* isFirstContinuousFactor,
                                   const char* isSecondContinuousFactor,
                                   const char** isNuisanceFactors,
                                   int iNumNuisanceFactors,
                                   const char* isMeasure,
                                   ProgressUpdateGUI* iProgressUpdateGUI )
{

  int errorCode;
  errorCode = this->mGlmDesign->Create ( this->mDataTable,
                                         isName,
                                         isFirstDiscreteFactor,
                                         isSecondDiscreteFactor,
                                         isFirstContinuousFactor,
                                         isSecondContinuousFactor,
                                         isNuisanceFactors,
                                         iNumNuisanceFactors,
                                         isMeasure,
                                         NULL,
                                         0,
                                         iProgressUpdateGUI );
  if ( errorCode )
    return errorCode;

  if ( iProgressUpdateGUI )
  {
    iProgressUpdateGUI->BeginActionWithProgress("Writing input files..." );
  }


  if ( mGlmDesign->WriteFsgdFile() )
  {
    fprintf( stderr, "ERROR: QdecProject::CreateGlmDesign: "
             "could not create fsgd file\n");
    return(-3);
  }

  if ( mGlmDesign->WriteContrastMatrices() )
  {
    fprintf( stderr, "ERROR: QdecProject::CreateGlmDesign: could not "
             "generate contrasts\n");
    return(-4);
  }

  if ( mGlmDesign->WriteYdataFile( isMeasure ) )
  {
    fprintf( stderr, "ERROR: QdecProject::CreateGlmDesign: could not "
             "create y.mgh file\n");
    return(-4);
  }

  if ( iProgressUpdateGUI )
  {
    iProgressUpdateGUI->EndActionWithProgress();
  }

  return 0;
}


/**
 * @return int
 */
int QdecProject::RunGlmFit ( )
{
  return this->mGlmFitter->Run( mGlmDesign );
}


/**
 * @return QdecGlmFitResults
 */
QdecGlmFitResults* QdecProject::GetGlmFitResults ( )
{
  return this->mGlmFitter->GetResults();
}


/**
 * Run mri_label2label on each subject, mapping the label that was drawn on
 * the average surface onto each subject. Optionally supply a GUI manager
 * to allow posting progress info.
 * @return int
 * @param  ifnLabel
 * @param  iProgressUpdateGUI
 */
int QdecProject::GenerateMappedLabelForAllSubjects
( const char* ifnLabel,
  ProgressUpdateGUI* iProgressUpdateGUI )
{
  vector< string > subjects = this->GetSubjectIDs();
  int numSubjects = this->GetSubjectIDs().size();
  float stepIncrement = 100.0 / numSubjects-1;
  int nStep = 1;

  if ( 0 == numSubjects )
    throw runtime_error( "Zero subjects! Cannot run mri_label2label\n" );

  if ( iProgressUpdateGUI )
  {
    iProgressUpdateGUI->BeginActionWithProgress
      ( "Running mri_label2label..." );
  }

  for ( int i=0; i < numSubjects; i++ )
  {
    // build a command line for this subject
    stringstream ssCommand;
    ssCommand << "mri_label2label"
              << " --srclabel " << ifnLabel
              << " --srcsubject " << this->GetAverageSubject()
              << " --trgsubject " << subjects[i]
              << " --trglabel " << ifnLabel
              << " --regmethod surface"
              << " --hemi " << this->GetHemi();

    // Now run the command.
    if ( iProgressUpdateGUI )
    {
      string status = "Running mri_label2label on subject '";
      status += subjects[i];
      status += "'...";
      iProgressUpdateGUI->UpdateProgressMessage( status.c_str() );
      iProgressUpdateGUI->UpdateProgressPercent
        ( (float)nStep++ * stepIncrement );
    }
    char* sCommand = strdup( ssCommand.str().c_str() );
    printf( "\n----------------------------------------------------------\n" );
    printf( "%s\n", sCommand );
    fflush(stdout);
    fflush(stderr);
    int rRun = system( sCommand );
    if ( -1 == rRun )
      throw runtime_error( "system call failed: " + ssCommand.str() );
    if ( rRun > 0 )
      throw runtime_error( "command failed: " + ssCommand.str() );
    free( sCommand );
  }

  if ( iProgressUpdateGUI )
  {
    iProgressUpdateGUI->UpdateProgressMessage( "Completed mri_label2label." );
    iProgressUpdateGUI->UpdateProgressPercent( 100 );
    iProgressUpdateGUI->EndActionWithProgress();
  }

  return 0;
}


/**
 * @return QdecGlmDesign
 */
QdecGlmDesign* QdecProject::GetGlmDesign ( )
{
  return this->mGlmDesign;
}

const char*
QdecProject::GetMetadataFileName () const
{

  static char fnMetadata[] = "QdecProjectMetadata.txt";
  return fnMetadata;
}

void
QdecProject::FormatCommandString ( const char* ifnProject,
                                   const char* isExpandedProjectBaseName,
                                   const char* isWorkingDir,
                                   const char* isFormat,
                                   string& iosCommand ) const
{
  assert( ifnProject );
  assert( isExpandedProjectBaseName );
  assert( isWorkingDir );
  assert( isFormat );

  // Start by copying the format string.
  iosCommand = isFormat;

  // Make our substitutions.
  string::size_type n;
  while ( string::npos != (n = iosCommand.find( "%1" )) )
    iosCommand.replace( n, 2, ifnProject );
  while ( string::npos != (n = iosCommand.find( "%2" )) )
    iosCommand.replace( n, 2, isExpandedProjectBaseName );
  while ( string::npos != (n = iosCommand.find( "%3" )) )
    iosCommand.replace( n, 2, isWorkingDir );
}



/**
 * run asegstats2table and aparcstats2table to generate fresurfer stats data
 * on each subject, for later optional inclusion into the main mDataTable.
 * creates files aseg_volume.dat, lh_aparc_thickness.dat...
 *
 * returns the names of the data that were created (aseg_volume,
 * lh_aparc_thickness...)
 *
 * @return vector< string >
 */
vector< string > QdecProject::CreateStatsDataTables ()
{
  // to be returned with names of stats data categories (files) created
  vector< string > statsDataNames;

  if ( ! this->HaveDataTable() ) return statsDataNames;

  vector< string > subjects = this->GetSubjectIDs();
  int numSubjects = this->GetSubjectIDs().size();

  if ( 0 == numSubjects )
    throw runtime_error( "Zero subjects! Cannot run asegstats2table\n" );

  // Make the sure the storage dir (/stats_tables) exists
  {
    string sCommand = "mkdir -p " + this->msStatsDataTablesDir;
    cout << sCommand << endl;
    int rRun = system( sCommand.c_str() );
    if ( -1 == rRun )
      throw runtime_error( "system call failed: " + sCommand );
    if ( rRun > 0 )
      throw runtime_error( "command failed: " + sCommand );
  }

  /*
   * start by running asegstats2table
   */
  vector< string > segs;
  segs.push_back( "aseg" );
  if (NULL == getenv("QDEC_SKIP_WMPARC_STAT"))
  {
    segs.push_back( "wmparc" );
  }

  unsigned int s;
  for (s=0; s < segs.size(); s++)
  {
    // build a command line
    stringstream name;
    name << segs[s] << ".volume";
    stringstream ssCommand;
    ssCommand << "asegstats2table --common-segs --meas volume --tablefile "
              << this->msStatsDataTablesDir
              << name.str() << ".stats.dat "
              << "--statsfile=" << segs[s] << ".stats "
              << "--subjects";
    for ( int i=0; i < numSubjects; i++ )
    {
      ssCommand << " " << subjects[i];
    }

    // and run the command...
    char* sCommand = strdup( ssCommand.str().c_str() );
    printf( "\n----------------------------------------------------------\n" );
    printf( "%s\n", sCommand );
    fflush(stdout);
    fflush(stderr);
    int rRun = system( sCommand );
    if ( -1 == rRun )
      throw runtime_error( "system call failed: " + ssCommand.str() );
    if ( rRun > 0 )
      throw runtime_error( "command failed: " + ssCommand.str() );
    free( sCommand );

    // save the name of this file
    statsDataNames.push_back( name.str() );
  }

  /*
   * now the variants of aparctats2table
   */

  vector< string > hemi;
  vector< string > parc;
  vector< string > meas;
  hemi.push_back( "lh" );
  hemi.push_back( "rh" );
  parc.push_back( "aparc" );
  parc.push_back( "aparc.a2009s" );
  meas.push_back( "area" );
  meas.push_back( "volume" );
  meas.push_back( "thickness" );
  meas.push_back( "meancurv" );
  unsigned int h,p,m;
  for ( h=0; h < hemi.size(); h++ ) // for each hemi...
  {
    for ( p=0; p < parc.size(); p++ ) // for each parcellation type
    {
      for ( m=0; m < meas.size(); m++ ) // for each measure
      {
        // construct name of output file
        stringstream ssFname;
        ssFname <<  hemi[h]
                << "." << parc[p]
                << "." << meas[m];

        // build a command line
        stringstream ssCommand;
        ssCommand << "aparcstats2table"
                  << " --hemi " << hemi[h]
                  << " --parc " << parc[p]
                  << " --meas " << meas[m]
                  << " --tablefile " << this->msStatsDataTablesDir
                  << "/" << ssFname.str() << ".stats.dat"
                  << " --subjects";
        for ( int i=0; i < numSubjects; i++ )
        {
          ssCommand << " " << subjects[i];
        }

        // and run the command...
        char* sCommand = strdup( ssCommand.str().c_str() );
        printf( "\n------------------------------------------------------\n" );
        printf( "%s\n", sCommand );
        fflush(stdout);
        fflush(stderr);
        int rRun = system( sCommand );
        free( sCommand );

        // don't exit on error, just keep trying to create tables.
        // aseg/aparcstats2table can fail if it is missing the raw data, in
        // which case the user needs to create that data.  error messages
        // will appear in the terminal showing the failures
        if ( 0 == rRun )
        {
          // save the name of this data (now that we know it was successfully
          // created
          statsDataNames.push_back( ssFname.str() );
        }
      }
    }
  }

  cout << "Completed creation of aseg and aparc stats data tables." << endl;

  return statsDataNames;
}


/**
 * Run mri_surfcluster using supplied sig.mgh (taken from contrast dir)
 * and supplied Monte Carlo threshold and sign to generate
 * cluster-wise correction for multiple comparisons results
 * making use of pre-calculated simulation data for fsaverage.
 * @return int
 * @param  isThreshold - one of: th13, th20, th23, th30, th33, th40
 * @param  isSign - one of: abs, pos, neg
 * @param  isContrast - name of contrast from which to use sig.mgh
 */
int
QdecProject::RunMonteCarloSimulation ( const char* isThreshold,
                                       const char* isSign,
                                       const char* isContrast,
                                       const char** osClusterSigFileName )
{
  stringstream ssWorkDir;
  ssWorkDir << this->GetDefaultWorkingDir()
            << "/"
            << this->GetGlmDesign()->GetName().c_str();
  char* sWorkDir = strdup( ssWorkDir.str().c_str() );

  // read fwhm.dat file, extract the value from it, and round up
  float fwhm = 0;
  stringstream fnFwhm;
  fnFwhm << sWorkDir << "/fwhm.dat";
  ifstream ifsFwhmFile( fnFwhm.str().c_str(), ios::in );
  if (ifsFwhmFile.good())
  {
    char tmpstr[1000];
    ifsFwhmFile.getline(tmpstr,1000);
    cout << "fwhm.dat: " << tmpstr;
    sscanf( tmpstr, "%f", &fwhm );
    fwhm = ceil( fwhm );
    cout << ", rounded to " << fwhm << endl;
  }
  else
  {
    cout << "ERROR: could not open " << fnFwhm.str().c_str() << endl;
    return 1;
  }
  if( (fwhm < 1) || (fwhm > 30) )
  {
    cout << "ERROR: fwhm out-of-range (< 1 or > 30)!" << endl;
    return 1;
  }
  
  // build mri_surfcluster command line
  stringstream ssCommand;
  stringstream ssClusterSigFileName;
  ssClusterSigFileName << sWorkDir << "/" << isContrast
                       << "/mc-z." << isSign << "." << isThreshold
                       << ".sig.cluster.mgh";
  ssCommand << "mri_surfcluster"
            << " --in " << sWorkDir << "/" << isContrast << "/sig.mgh"
            << " --csd " << getenv("FREESURFER_HOME")
            << "/average/mult-comp-cor/fsaverage/" << this->GetHemi()
            << "/cortex/fwhm" << (int)fwhm << "/" << isSign
            << "/" << isThreshold << "/mc-z.csd"
            << " --mask " << sWorkDir << "/mask.mgh"
            << " --cwsig " << ssClusterSigFileName.str().c_str()
            << " --vwsig " << sWorkDir << "/" << isContrast
            << "/mc-z." << isSign << "." << isThreshold << ".sig.vertex.mgh"
            << " --sum " << sWorkDir << "/" << isContrast
            << "/mc-z."<<isSign << "." << isThreshold << ".sig.cluster.summary"
            << " --ocn " << sWorkDir << "/" << isContrast
            << "/mc-z." << isSign << "." << isThreshold << ".sig.ocn.mgh"
            << " --oannot " <<sWorkDir << "/" << isContrast
            << "/mc-z." << isSign << "." << isThreshold << ".sig.ocn.annot"
            << " --csdpdf "<< sWorkDir << "/" << isContrast
            << "/mc-z." << isSign << "." << isThreshold << ".pdf.dat"
            << " --annot aparc"
            << " --cwpvalthresh 0.05 --surf white";


  // and run the command...
  char* sCommand = strdup( ssCommand.str().c_str() );
  printf( "\n------------------------------------------------------\n" );
  printf( "%s\n", sCommand );
  fflush(stdout);
  fflush(stderr);
  int rRun = system( sCommand );
  free( sCommand );

  if ( rRun )
  {
    cout << "\nERROR!\n" << endl;
  }
  else
  {
    // print cluster summary
    stringstream ssCat;
    ssCat << "cat " << sWorkDir << "/" << isContrast
          << "/mc-z." <<isSign << "." << isThreshold << ".sig.cluster.summary";
    system( ssCat.str().c_str() );

    cout << "\nSimulation complete.\n" << endl;

    // save path to results sig file
    if ( osClusterSigFileName )
    {
      *osClusterSigFileName = strdup( ssClusterSigFileName.str().c_str() );
    }
  }

  return rRun;
}
