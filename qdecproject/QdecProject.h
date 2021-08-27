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

#ifndef QDECPROJECT_H
#define QDECPROJECT_H

#include <string>
#include <vector>

#include "QdecDataTable.h"
#include "QdecGlmDesign.h"
#include "QdecGlmFit.h"
#include "QdecGlmFitResults.h"
#include "ProgressUpdateGUI.h"

using namespace std;

class QdecProject
{
public:

  // Constructors/Destructors
  //

  QdecProject ( );

  virtual ~QdecProject ( );

  /**
   * Load a .qdec project file (containing all necessary info to begin
   * working either on a new project, or to continue working using
   * results from a prior saved work session). isDataDir should be a
   * directory where we can expand the .qdec file (like /tmp).
   * @return int
   * @param  isFileName
   * @param  isDataDir
   */
  int LoadProjectFile ( const char* isFileName,
			const char* isDataDir = "/tmp" );


  /**
   * Save all necessary information pertaining to this project
   * (all subject data, any results, any user preferences).
   * @return int
   * @param  isFileName
   */
  int SaveProjectFile ( const char* isFileName,
			const char* isDataDir = "/tmp" );


  /**
   *
   * The command format strings to zip and unzip a file. Returns -1 if
   * parameter is invalid. The default is acceptable for Linux systems
   * with unzip and zip installed. The substitutions that are made
   * are:
   *
   * %1 - Full project filename
   * %2 - Expanded project base name
   * %3 - Working dir (ifnDataDir)
   *
   * Default zip format string is:
   * cd %3; zip -r %1 %2 > /dev/null
   * Default unzip format string is:
   * unzip -o -d %3 %1 > /dev/null
   *
   * @return int
   * @param isFormat
   */
  int SetZipCommandFormat ( const char* isFormat );
  int SetUnzipCommandFormat ( const char* isFormat );

  /**
   * @return int
   * @param  isFileName
   */
  int LoadDataTable ( const char* isFileName );


  /**
   * Check that all subjects exist in the specified subjects_dir (including the
   * specified average subject).  Print to stderr and ErrorMessage any errors
   * found (one message for each error).  Also check that thickness, sulc,
   * curv, area and jacobian_white files exist, and that their vertex
   * numbers equal their inflated surface (and that surfaces all have the
   * same number of vertices).
   * @return int
   */
  int VerifySubjects ( );


  /**
   * @return void
   * @param  iFilePointer
   */
  void DumpDataTable ( FILE* iFilePointer );


  /**
   * @return bool  true if a table has been loaded
   */
  bool HaveDataTable ( );


  /**
   * @return int
   * @param  isFileName
   */
  int SaveDataTable ( const char* isFileName );


  /**
   * @return QdecDataTable*
   */
  QdecDataTable* GetDataTable ( );


  /**
   * run asegstats2table and aparcstats2table to generate fresurfer stats data
   * on each subject, for later optional inclusion into the main mDataTable
   * creates file aseg_vol.dat, lh_aparc_thickness.dat...
   *
   * returns the names of the data that were created (aseg_vol, 
   * lh_aparc_thickness...)
   *
   * @return vector< string >
   */
  vector< string > CreateStatsDataTables ();


  /**
   * @return string   directory where stats data tables are stored
   */
  string GetStatsDataTablesDir () { return this->msStatsDataTablesDir; }
  

  /**
   * Merge the factor given by isFactorName from data table iDataTable
   * into our main DataTable.  Used to add selections from an asegstats
   * table or aparcstats table into the main table (which would require
   * saving the table to disk if the user wants to keep those changes)
   *
   * @return int
   * @param  isFactorName
   * @param  iDataTable
   */
  int MergeFactorIntoDataTable ( const char* isFactorName, 
                                 QdecDataTable* iDataTable )
  {
    return this->mDataTable->MergeFactor( isFactorName, iDataTable );
  }


  /**
   * Delete the factor isFactorName from data table
   *
   * @return int
   * @param  isFactorName
   */
  int RemoveFactorFromDataTable ( const char* isFactorName )
  {
    return this->mDataTable->DeleteFactor( isFactorName );
  }


  /**
   * @return string
   */
  string GetSubjectsDir ( );


  /**
   * @param  ifnSubjectsDir
   */
  int SetSubjectsDir ( const char* ifnSubjectsDir );


  /**
   * @return string
   */
  string GetAverageSubject ( );


  /**
   * @param  isSubjectName
   */
  void SetAverageSubject ( const char* isSubjectName );


  /**
   * @return string
   */
  string GetDefaultWorkingDir ( );


  /**
   * @return string
   */
  string GetWorkingDir ( );


  /**
   * @return 0 if ok, 1 on error
   * @param  isWorkingDir
   */
  int SetWorkingDir ( const char* isWorkingDir );


  /**
   * @return vector< string >
   */
  vector< string > GetSubjectIDs ( );


  /**
   * @return vector< string >
   */
  vector< string > GetDiscreteFactorNames ( );


  /**
   * @return vector< string >
   */
  vector< string > GetContinousFactorNames ( );


  /**
   * @return string
   */
  string GetHemi ( );


  /**
   *
   * For Surface-based analysis:
   * 
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
   * @param  iProgressUpdateGUI (optional)
   */
  int CreateGlmDesign ( const char* isName,
                        const char* isFirstDiscreteFactor,
                        const char* isSecondDiscreteFactor,
                        const char* isFirstContinuousFactor,
                        const char* isSecondContinuousFactor,
                        const char** isNuisanceFactors,
                        int iNumNuisanceFactors,
                        const char* isMeasure,
                        const char* isHemi,
                        int iSmoothnessLevel,
                        ProgressUpdateGUI* iProgressUpdateGUI=NULL );


  /**
   *
   * For Volume-based analysis:
   * 
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
  int CreateGlmDesign ( const char* isName,
                        const char* isFirstDiscreteFactor,
                        const char* isSecondDiscreteFactor,
                        const char* isFirstContinuousFactor,
                        const char* isSecondContinuousFactor,
                        const char** isNuisanceFactors,
                        int iNumNuisanceFactors,
                        const char* isMeasure,
                        ProgressUpdateGUI* iProgressUpdateGUI );


  /**
   * @return int
   */
  int RunGlmFit ( );


  /**
   * @return QdecGlmFitResults
   */
  QdecGlmFitResults* GetGlmFitResults ( );


  /**
   * Run mri_label2label on each subject, mapping the label that was drawn on 
   * the average surface to each subject. Optionally supply a window manager
   * to allow posting progress info
   * @return int
   * @param  ifnLabel
   * @param  iProgressUpdateGUI (optional)
   */
  int GenerateMappedLabelForAllSubjects 
    ( const char* ifnLabel,
      ProgressUpdateGUI* iProgressUpdateGUI=NULL );


  /**
   * @return QdecGlmDesign
   */
  QdecGlmDesign* GetGlmDesign ( );


  /**
   * The file name of our metadata file, for the project file archive.
   * @return const char*
   */
  const char* GetMetadataFileName () const;


  /**
   * Perform substitutions for command format strings. See
   * documentation for Set(Un)ZipCommandFormat. This will perform the
   * substitutions on isFormat and write the command to iosCommand
   * (overwriting the contents of iosCommand).
   *
   */
  void FormatCommandString ( const char* ifnProject,
			     const char* isExpandedProjectBaseName,
			     const char* isWorkingDir,
			     const char* isFormat,
			     string& iosCommand ) const;


  /**
   * Run mri_surfcluster using supplied sig.mgh (for a contrast)
   * and supplied Monte Carlo threshold and sign to generate 
   * cluster-wise correction for multiple comparisons results
   * @return int
   * @param  isThreshold - one of: th13, th20, th23, th30, th33, th40
   * @param  isSign - one of: abs, pos, neg
   * @param  isContrast - name of contrast from which to use sig.mgh
   */
  int RunMonteCarloSimulation ( const char* isThreshold,
                                const char* isSign,
                                const char* isContrast,
                                const char** osClusterSigFileName );

private:

  // private attributes
  //

  QdecDataTable* mDataTable;
  QdecGlmDesign* mGlmDesign;
  QdecGlmFit* mGlmFitter;

  // The command format to run to zip and unzip a file. 
  string msZipCommandFormat;
  string msUnzipCommandFormat;

  // directory where stats data tables are stored
  string msStatsDataTablesDir;

};

#endif // QDECPROJECT_H
