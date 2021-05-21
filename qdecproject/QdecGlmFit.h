/**
 * @brief Wrapper for mri_glmfit.
 *
 * Run mri_glmfit, given its input data, and puts the output in the specified
 * working directory.
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

#ifndef QDECGLMFIT_H
#define QDECGLMFIT_H

#include <string>
#include <vector>

#include "QdecGlmDesign.h"
#include "QdecGlmFitResults.h"

using namespace std;

class QdecGlmFit
{
public:

  // Constructors/Destructors
  //

  QdecGlmFit ( );

  virtual ~QdecGlmFit ( );

  // public attribute accessor methods
  //


  /**
   * @return int
   * @param  iGlmDesign
   */
  int Run (QdecGlmDesign* iGlmDesign );

  /**
   * Creates the contrast, stddef, coefficients, and fsgdf file names
   * using the working directory from the design. Creates the results
   * object with these values. Does not generate any new data.
   * @return int
   * @param iGlmDesign
   */
  int CreateResultsFromCachedData ( QdecGlmDesign* iGlmDesign );

  /**
   * @return QdecGlmFitResults
   */
  QdecGlmFitResults* GetResults ( );

private:

  // private attributes
  //

  QdecGlmFitResults* mGlmFitResults;

};

#endif // QDECGLMFIT_H
