/**
 * @brief Contains the results of a GLM fit run.
 *
 * The bulk of the result data is stored in files on disk, and this object
 * contains paths to that data, as well as some local results.
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

#ifndef QDECGLMFITRESULTS_H
#define QDECGLMFITRESULTS_H

#include <string>
#include <vector>

#include "QdecGlmDesign.h"

class QdecGlmFitResults {
public:
  // Constructors/Destructors
  //

  QdecGlmFitResults(
      QdecGlmDesign *          iGlmDesign,
      std::vector<std::string> iContrastSigFiles,      /* /<contrast>/sig.mgh */
      std::string              iConcatContrastSigFile, /* contrast.sig.mgh */
      std::string              ifnResidualErrorStdDevFile,    /* rstd.mgh */
      std::string              ifnRegressionCoefficientsFile, /* beta.mgh */
      std::string              ifnFsgdFile /* y.fsgd */);

  virtual ~QdecGlmFitResults();

  // public attribute accessor methods
  //

  /**
   * Returns the design object used as input to the GLM fitter
   * @return QdecGlmDesign*
   */
  QdecGlmDesign *GetGlmDesign();

  /**
   * Returns the names given to the contrast results produced by glmfit.
   * Example of one of the possible names: "Avg-thickness-Age-Cor"
   * @return vector< string >
   */
  std::vector<std::string> GetContrastNames();

  /**
   * Returns the human-readable questions associated with each contrast.
   * Example of one question:
   * "Does the correlation between thickness and Age differ from zero?".
   * @return vector< string >
   */
  std::vector<std::string> GetContrastQuestions();

  /**
   * Returns pathname to the concatenated contrast significance file,
   * ie sig.mgh for all contrasts.
   * @return string
   */
  std::string GetConcatContrastSigFile();

  /**
   * Returns pathnames to the contrast significance file, ie sig.mgh for that
   * contrast.
   * @return vector< string >
   */
  std::vector<std::string> GetContrastSigFiles();

  /**
   * Returns pathnames to the contrast gamma file, ie gamma.mgh for that
   * contrast.
   * @return vector< string >
   */
  std::vector<std::string> GetContrastGammaFiles();

  /**
   * Returns pathnames to the contrast F-test file, ie F.mgh for that contrast.
   * @return vector< string >
   */
  std::vector<std::string> GetContrast_F_Files();

  /**
   * Returns pathname to the beta.mgh file.
   * @return string
   */
  std::string GetRegressionCoefficientsFile();

  /**
   * Returns pathname to eres.mgh
   * @return string
   */
  std::string GetResidualErrorFile();

  /**
   * Returns pathname to rstd.mgh
   * @return string
   */
  std::string GetResidualErrorStdDevFile();

  /**
   * Returns pathname to y.fsgd
   * @return string
   */
  std::string GetFsgdFile();

private:
  // private attributes
  //

  QdecGlmDesign *          mGlmDesign;
  std::vector<std::string> mContrastNames;
  std::vector<std::string> mContrastQuestions;
  std::string              mfnConcatContrastSigFile; /* contrast.sig.mgh */
  std::vector<std::string> mfnContrastSigFiles;      /* /<contrast>/sig.mgh */
  std::vector<std::string> mfnContrastGammaFiles;
  std::vector<std::string> mfnContrast_F_Files;
  std::string              mfnRegressionCoefficientsFile; /* beta.mgh */
  std::string              mfnResidualErrorFile;          /* eres.mgh */
  std::string              mfnResidualErrorStdDevFile;    /* rstd.mgh */
  std::string              mfnFsgdFile;                   /* y.fsgd */
};

#endif // QDECGLMFITRESULTS_H
