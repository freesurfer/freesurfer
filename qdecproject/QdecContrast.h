/**
 * @file  QdecContrast.h
 * @brief Stores a GLM contrast vector
 *
 * Stores a GLM contrast vector associated with a particular design (being
 * based on user selected factors, formulating a hypothesis to test).
 * Stores the name associated with it (to go in the .fsgd file), the
 * human-readable question (for use by GUI), and the vector itself (as an
 * int array).
 */
/*
 * Original Author: Nick Schmansky
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2007/12/13 22:41:57 $
 *    $Revision: 1.1.2.1 $
 *
 * Copyright (C) 2007,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

#ifndef QDECCONTRAST_H
#define QDECCONTRAST_H

#include <string>
#include <vector>

using namespace std;

class QdecContrast
{
public:

  // Constructors/Destructors
  //

  QdecContrast ( vector< double > iaVector,
                 string isName,
                 string isQuestion );

  virtual ~QdecContrast ( );

  /**
   * @return string
   */
  string GetName ( );


  /**
   * @return string
   */
  string GetQuestion ( );


  /**
   * @return string
   */
  string GetContrastStr ( );


  /**
   * Writes the contrast vector to a .mat file, which is readable by matlab,
   * and mri_glmfit.
   * @return int
   * @param string ifnWorkingDir
   */
  int WriteDotMatFile ( string ifnWorkingDir );


  /**
   * @return string
   */
  string GetDotMatFileName ( );

private:

  // private attributes
  //

  vector< double > maVector;
  string msName;
  string msQuestion;
  string mfnDotMatFileName;

};

#endif // QDECCONTRAST_H
