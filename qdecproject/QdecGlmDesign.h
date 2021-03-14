/**
 * @brief Contains the data and functions associated with a GLM design run
 *
 * Contains the data and functions associated with the selected GLM input
 * (design) parameters: selected discrete and continuous factors, measure,
 * hemi, smoothness level and design name.  Functions exist to create derived
 * data: .fsgd file, contrast vectors in .mat file format, and the 'y' input
 * data file (concatenated subjects file).
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

#ifndef QDECGLMDESIGN_H
#define QDECGLMDESIGN_H

#include <set>
#include <string>
#include <vector>

#include "ProgressUpdateGUI.h"
#include "QdecContrast.h"
#include "QdecDataTable.h"

class QdecGlmDesign {
public:
  // Constructors/Destructors
  //

  QdecGlmDesign(QdecDataTable *iDataTable);

  virtual ~QdecGlmDesign();

  /**
   * Returns true if this design is valid (input parameters have been set and
   * mri_glmfit input data created, ie. Create() has been called successfully)
   * @return bool
   */
  bool IsValid();

  /**
   * Initializes the design with the given design parameters.
   * @return int
   * @param  iDataTable
   * @param  isName
   * @param  isFirstDiscreteFactor
   * @param  isSecondDiscreteFactor
   * @param  isFirstContinuousFactor
   * @param  isSecondContinuousFactor
   * @param  isNuisanceFactors
   * @param  inNumNuisanceFactors
   * @param  isMeasure
   * @param  isHemi
   * @param  iSmoothnessLevel
   * @param  iProgressUpdateGUI
   */
  int Create(QdecDataTable *iDataTable, const char *isName,
             const char * isFirstDiscreteFactor,
             const char * isSecondDiscreteFactor,
             const char * isFirstContinuousFactor,
             const char * isSecondContinuousFactor,
             const char **isNuisanceFactors, int inNumNuisanceFactors,
             const char *isMeasure, const char *isHemi, int iSmoothnessLevel,
             ProgressUpdateGUI *iProgressUpdateGUI);

  /**
   *
   */
  void ClearDiscreteFactors();

  /**
   *
   */
  void AddDiscreteFactor(const char *isFactorName);

  /**
   *
   */
  void ClearContinuousFactors();

  /**
   *
   */
  void AddContinuousFactor(const char *isFactorName);

  /**
   *
   */
  void ClearNuisanceFactors();

  /**
   *
   */
  void AddNuisanceFactor(const char *isFactorName);

  /**
   * @return int
   */
  int GetDegreesOfFreedom();

  /**
   * @return string
   */
  std::string GetName();

  /**
   *
   */
  void SetName(const char *isName);

  /**
   * @return string
   */
  std::string GetHemi();

  /**
   *
   */
  void SetHemi(const char *isHemi);

  /**
   * @return string
   */
  std::string GetMeasure();

  /**
   *
   */
  void SetMeasure(const char *isMeasure);

  /**
   * @return int
   */
  int GetSmoothness();

  /**
   *
   */
  void SetSmoothness(int iVal);

  /**
   * @return string
   */
  std::string GetDesignMatrixType();

  /**
   * @param const char*
   */
  void SetDesignMatrixType(const char *isDesignMatrixType);

  /**
   * @return string
   */
  std::string GetSubjectsDir();

  /**
   * @param const char*
   */
  int SetSubjectsDir(const char *ifnSubjectsDir);

  /**
   * @return string
   */
  std::string GetAverageSubject();

  /**
   * @param const char*
   */
  void SetAverageSubject(const char *isAverageSubject);

  /**
   * returns the pathname to the fsgd file required by mri_glmfit.
   * @return string
   */
  std::string GetFsgdFileName();

  /**
   * returns the pathname to the input data, 'y', required by mri_glmfit.
   * @return string
   */
  std::string GetYdataFileName();

  /**
   * @return vector< string >
   */
  std::vector<std::string> GetContrastNames();

  /**
   * @return vector< string >
   */
  std::vector<std::string> GetContrastQuestions();

  /**
   * @return vector< string >
   */
  std::vector<std::string> GetContrastFileNames();

  /**
   * @return string
   */
  std::string GetDefaultWorkingDir();

  /**
   * @return string
   */
  std::string GetWorkingDir();

  /**
   * @return int
   * @param  isPathName
   */
  int SetWorkingDir(const char *isPathName);

  /**
   * @return ProgressUpdateGUI*
   */
  ProgressUpdateGUI *GetProgressUpdateGUI();

  /**
   * SetExcludeSubjectID ( const char* isSubjecID, bool ibExclude ) -
   * Sets a subject ID's exclusion status. If excluded, it will not be
   * included when writing the ydata file.
   * param const char* isSubjectID
   * param bool ibExclude
   */
  void SetExcludeSubjectID(const char *isSubjectID, bool ibExclude);

  /**
   * GetExcludeSubjectID ( const char* isSubjecID ) -
   * Returns a subject ID's exclusion status.
   * param const char* isSubjectID
   */
  bool GetExcludeSubjectID(const char *isSubjectID);

  /**
   * SetExcludeSubjectsFactorGT
   */
  void SetExcludeSubjectsFactorGT(const char *isFactorName, double inExcludeGT,
                                  bool ibExclude);

  /**
   * SetExcludeSubjectsFactorLT
   */
  void SetExcludeSubjectsFactorLT(const char *isFactorName, double inExcludeLT,
                                  bool ibExclude);

  /**
   * SetExcludeSubjectsFactorET
   */
  void SetExcludeSubjectsFactorET(const char *isFactorName, double inExcludeET,
                                  bool ibExclude);

  /**
   * ClearAllExcludedSubjects
   */
  void ClearAllExcludedSubjects();

  /**
   * GetNumberOfExcludedSubjects
   */
  int GetNumberOfExcludedSubjects();

  /**
   * Using the design parameters, writes FSGF file to the working directory.
   * @return int
   */
  int WriteFsgdFile();

  /**
   * Using the design parameters, writes .mat files for all our contrasts.
   * @return int
   */
  int WriteContrastMatrices();

  /**
   * Using the design parameters, creates the 'y' input data to
   * mri_glmfit, by concatenating the subject volumes, and writes it
   * to the specified filename.
   * @return int
   */
  int WriteYdataFile();

  /**
   * Creates the 'y' input data to mri_glmfit for a volume-based analysis,
   * by creating a multi-frame 'volume' of 1x1x1 size, sticking the single
   * dat point for each subject in the 'volume'.  One frame per subject.
   * @return int
   */
  int WriteYdataFile(const char *isMeasureName);

  /**
   * Access the discrete and continuous factor names.
   * Returns a const vector of QdecFactors pointers.
   */
  std::vector<QdecFactor *> const &GetDiscreteFactors() const;
  std::vector<QdecFactor *> const &GetContinuousFactors() const;
  std::vector<QdecFactor *> const &GetNuisanceFactors() const;

private:
  // private attributes
  //

  bool           mbValid;
  QdecDataTable *mDataTable;
  std::string    msName;
  // Stores seleted discrete factors.  Initially empty.
  std::vector<QdecFactor *> mDiscreteFactors;
  // Stores selected continous factors.  Initially empty.
  std::vector<QdecFactor *> mContinuousFactors;
  // Stores selected nuisance factors.  Initially empty.
  std::vector<QdecFactor *> mNuisanceFactors;
  std::string               msMeasure;
  std::string               msHemi;
  int                       mSmoothness;
  std::string               msDesignMatrixType; // dods or doss
  std::string               mfnSubjectsDir;
  std::string               msAverageSubject;
  // Stores contrasts created from an fsgdf file. Can be empty.
  std::vector<QdecContrast *> mContrasts;
  std::string                 mfnFsgdfFile;
  std::string                 mfnYdataFile;
  std::string                 mfnDefaultWorkingDir;
  std::string                 mfnWorkingDir;
  ProgressUpdateGUI *         mProgressUpdateGUI;

  // A list of excluded subjects. These will not be included when
  // writing the fsgd and ydata files. The key values are subject IDs,
  // as found in data table, and if there is a value present in the set,
  // that subject is to be excluded.
  std::set<std::string> maExcludedSubjects;

  // private methods
  //

  int GetNumberOfDiscreteFactors() { return this->mDiscreteFactors.size(); }
  int GetNumberOfContinuousFactors() { return this->mContinuousFactors.size(); }
  int GetNumberOfNuisanceFactors() { return this->mNuisanceFactors.size(); }

  /**
   * GetNumberOfClasses( ) - returns the number of classes for the design.
   * The number of classes is just all the combinations of all
   * the levels for the discrete factors.
   */
  int GetNumberOfClasses();

  /**
   * GetNumberOfRegressors() - returns the number of regressors for the
   * given design.
   */
  int GetNumberOfRegressors();

  /**
   * GetLevels2ClassName() - returns the class name given that
   * the 1st factor is at nthlevels[0],
   * the 2nd factor is at nthlevels[1], ...
   * The class name is created by appending
   *   Factor1NameLevelName-Factor2NameLevelName...
   */
  std::string GetLevels2ClassName(unsigned int *nthlevels);

  /**
   * Creates Contrast objects based on the selected factors.
   * Stores them in our 'mContrasts' QdecContrast container.
   * @return int
   */
  int GenerateContrasts();
};

#endif // QDECGLMDESIGN_H
