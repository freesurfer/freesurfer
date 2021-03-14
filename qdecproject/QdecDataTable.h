/**
 * @brief Container for the text-input file to QDEC.
 *
 * Implements loading/saving the white-space delimited data file containing
 * the list of subjects with their discrete and continuous factors.
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

#ifndef QDECDATATABLE_H
#define QDECDATATABLE_H

#include <string>
#include <vector>

#include "QdecFactor.h"
#include "QdecSubject.h"

class QdecDataTable {
public:
  // Constructors/Destructors
  //

  QdecDataTable();

  virtual ~QdecDataTable();

  /**
   * Load white-space delimited file containing subject ids and their
   * discrete and continuous factors.
   * @return int
   * @param  isFileName
   * @param  osNewSubjDir
   * @param  isFsIdColName
   */
  int Load(const char *isFileName, char *osNewSubjDir = nullptr,
           const char *isFsIdColName = nullptr);

  /**
   * @return int
   * @param  isFileName
   */
  int Save(const char *isFileName);

  /**
   * @return string
   */
  std::string GetFileName();

  /**
   * @return vector< string >
   */
  std::vector<std::string> GetSubjectIDs();

  /**
   * @return vector< QdecSubject* >
   */
  std::vector<QdecSubject *> GetSubjects();

  /**
   * @return QdecFactor*
   * @param isFactorName
   */
  QdecFactor *GetFactor(const char *isFactorName);

  /**
   * @return QdecFactor*
   * @param isSubjectName
   * @param isFactorName
   */
  QdecFactor *GetFactor(const char *isSubjectName, const char *isFactorName);

  /**
   * @return vector< string >
   */
  std::vector<std::string> GetDiscreteFactorNames();

  /**
   * @return vector< string >
   */
  std::vector<std::string> GetContinuousFactorNames();

  /**
   * GetNumberOfClasses( ) - returns the number of subjects in the table
   */
  int GetNumberOfSubjects();

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
   * dumps factors and inputs to filepointer (stdout, or file)
   * @param  iFilePointer
   */
  void Dump(FILE *iFilePointer);

  /**
   * GetMeanAndStdDev() - computes the average and stddev of continuous factor
   * @return vector< double > - first element is mean, second is the stddev
   * @param isFactorName
   */
  std::vector<double> GetMeanAndStdDev(const char *isFactorName);

  /**
   * deletes all continuous factors that have a zero mean and zero stddev.
   * those are useless factors (probably originating from a stats table).
   * @return number of factors purged
   */
  int PurgeNullFactors();

  /**
   * merge a factor from a given data table into this data table
   *
   * @return int
   * @param  isFactorName
   * @param  iDataTable
   */
  int MergeFactor(const char *isFactorName, QdecDataTable *iDataTable);

  /**
   * delete a factor from this data table
   *
   * @return int
   * @param  isFactorName
   */
  int DeleteFactor(const char *isFactorName);

private:
  // private attributes
  //

  // name of the text file from which this data is loaded
  std::string mfnFileName;

  // discrete and continuous factors as found on first line of data table
  std::vector<QdecFactor *> mFactors;

  // Stores subject data (id and factors) as read from table.dat input file.
  std::vector<QdecSubject *> mSubjects;
};

#endif // QDECDATATABLE_H
