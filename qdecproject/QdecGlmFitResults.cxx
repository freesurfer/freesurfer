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

#include "QdecGlmFitResults.h"
#include <cstring>

// Constructors/Destructors
//

QdecGlmFitResults::QdecGlmFitResults(
    QdecGlmDesign *          iGlmDesign,
    std::vector<std::string> iContrastSigFiles,      /* /<contrast>/sig.mgh */
    std::string              iConcatContrastSigFile, /* contrast.sig.mgh */
    std::string              ifnResidualErrorStdDevFile,    /* rstd.mgh */
    std::string              ifnRegressionCoefficientsFile, /* beta.mgh */
    std::string              ifnFsgdFile /* y.fsgd */) {
  assert(iGlmDesign);
  assert(iContrastSigFiles.size());

  this->mGlmDesign                    = iGlmDesign;
  this->mfnContrastSigFiles           = iContrastSigFiles;
  this->mfnConcatContrastSigFile      = iConcatContrastSigFile;
  this->mfnResidualErrorStdDevFile    = ifnResidualErrorStdDevFile;
  this->mfnRegressionCoefficientsFile = ifnRegressionCoefficientsFile;
  this->mfnFsgdFile                   = ifnFsgdFile;
}

QdecGlmFitResults::~QdecGlmFitResults() {}

//
// Methods
//

/**
 * Returns the design object used as input to the GLM fitter
 * @return QdecGlmDesign*
 */
QdecGlmDesign *QdecGlmFitResults::GetGlmDesign() { return this->mGlmDesign; }

/**
 * Returns the names given to the contrast results produced by glmfit.
 * Example of one of the possible names: "Avg-thickness-Age-Cor"
 * @return vector< string >
 */
std::vector<std::string> QdecGlmFitResults::GetContrastNames() {
  return this->mGlmDesign->GetContrastNames();
}

/**
 * Returns the human-readable questions associated with each contrast.
 * Example of one question:
 * "Does the correlation between thickness and Age differ from zero?".
 * @return vector< string >
 */
std::vector<std::string> QdecGlmFitResults::GetContrastQuestions() {
  return this->mGlmDesign->GetContrastQuestions();
}

/**
 * Returns pathname to the concatenated contrast significance file,
 * ie sig.mgh for all contrasts.
 * @return string
 */
std::string QdecGlmFitResults::GetConcatContrastSigFile() {
  return this->mfnConcatContrastSigFile;
}

/**
 * Returns pathnames to the contrast significance file, ie sig.mgh for that
 * contrast.
 * @return vector< string >
 */
std::vector<std::string> QdecGlmFitResults::GetContrastSigFiles() {
  return this->mfnContrastSigFiles;
}

/**
 * Returns pathnames to the contrast gamma file, ie gamma.mgh for
 * that contrast.
 * @return vector< string >
 */
std::vector<std::string> QdecGlmFitResults::GetContrastGammaFiles() {
  std::vector<std::string> tmp;
  return tmp; // TODO
}

/**
 * Returns pathnames to the contrast F-test file, ie F.mgh for that contrast.
 * @return vector< string >
 */
std::vector<std::string> QdecGlmFitResults::GetContrast_F_Files() {
  std::vector<std::string> tmp;
  return tmp; // TODO
}

/**
 * Returns pathname to the beta.mgh file.
 * @return string
 */
std::string QdecGlmFitResults::GetRegressionCoefficientsFile() {
  return this->mfnRegressionCoefficientsFile;
}

/**
 * Returns pathname to eres.mgh
 * @return string
 */
std::string QdecGlmFitResults::GetResidualErrorFile() {
  return this->mfnResidualErrorFile;
}

/**
 * Returns pathname to rstd.mgh
 * @return string
 */
std::string QdecGlmFitResults::GetResidualErrorStdDevFile() {
  return this->mfnResidualErrorStdDevFile;
}

/**
 * Returns pathname to y.fsgd
 * @return string
 */
std::string QdecGlmFitResults::GetFsgdFile() { return this->mfnFsgdFile; }
