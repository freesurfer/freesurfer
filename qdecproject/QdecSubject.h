/**
 * @file  QdecSubject.h
 * @brief Stores all data associated with a subject.
 *
 * This is one row from the input data table file (qdec.table.dat).
 */
/*
 * Original Author: Nick Schmansky
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/01/21 02:56:53 $
 *    $Revision: 1.2 $
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

#ifndef QDECSUBJECT_H
#define QDECSUBJECT_H

#include <string>
#include <vector>

#include "QdecFactor.h"

using namespace std;

class QdecSubject
{
public:

  // Constructors/Destructors
  //

  QdecSubject ( string isId, vector < QdecFactor* > iFactors );

  virtual ~QdecSubject ( );

  /**
   * Get the value of msId
   * the subject identifier, as found in the 'fsid' column of the
   * table.dat input file.
   * @return the value of msId
   * @return string
   */
  string GetId ( );


  /**
   * @return string
   * @param  isFactorName
   */
  string GetDiscreteFactorValue ( const char* isFactorName );


  /**
   * @return double
   * @param  isFactorName
   */
  double GetContinuousFactorValue ( const char* isFactorName );

  /**
   * @return vector < QdecFactor* >
   */
  vector < QdecFactor* > GetContinuousFactors ( );

  /**
   * @return vector < QdecFactor* >
   */
  vector < QdecFactor* > GetFactors ( );

private:

  // private attributes
  //

  // the subject identifier, as found in the 'fsid' column
  // of the table.dat input file.
  string msId;

  // Stores factor values (either discrete or continous)
  // pertaining to this subject.
  vector < QdecFactor* > mFactors;

};

#endif // QDECSUBJECT_H
