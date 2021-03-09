/**
 * @brief Stores all data associated with a subject.
 *
 * This is one row from the input data table file (qdec.table.dat).
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
  string GetId ( ) { return this->msId; }


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
  vector < QdecFactor* > GetFactors ( ) { return this->mFactors; }


  /**
   * @return QdecFactor
   */
  QdecFactor* GetFactor ( const char* isFactorName );


  /**
   */
  void AddFactor ( QdecFactor* iFactor ) {this->mFactors.push_back(iFactor);}


  /**
   * @param  isFactorName
   */
  void DeleteFactor ( const char* isFactorName );

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
