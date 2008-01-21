/**
 * @file  QdecFactor.h
 * @brief Stores a factor, which can be either discrete or continuous
 *
 * An example of a discrete factor is gender (male or female) or
 * diagnosis (demented or nondemented).  An example continuous factor is
 * age, or volume of a subcortical structure.
 */
/*
 * Original Author: Nick Schmansky
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/01/21 02:56:53 $
 *    $Revision: 1.3 $
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

#ifndef QDECFACTOR_H
#define QDECFACTOR_H

#include <cassert>
#include <vector>
#include <string>

using namespace std;

class QdecFactor
{
public:

  static const int qdecDiscreteFactorType = 1;
  static const int qdecContinuousFactorType = 2;

  // Constructors/Destructors
  //

  QdecFactor ( const char* isName,
               int iType /* ==1 discrete or ==2 continuous */ );
  QdecFactor ( const char* isName,
               int iType, // ==1 discrete
               const char* iValue,
               vector< string > iLevelNames );
  QdecFactor ( const char* isName,
               int iType, // ==2 continuous
               double iValue );

  virtual ~QdecFactor ( );

  /**
   * @return bool
   */
  bool IsDiscrete ( );

  /**
   * 
   */
  void SetDiscrete ( ) { mType = qdecDiscreteFactorType; } ;

  /**
   * @return bool
   */
  bool IsContinuous ( );


  /**
   * @return string
   */
  string GetFactorName ( );


  /**
   * GetFactorTypeName() - returns the string name of the
   * type of the given factor: 'continuous' or 'discrete'
   * @return string
   */
  string GetFactorTypeName ( );

  /**
   * @return int
   * @param  isLevelName
   */
  void AddLevelName ( string isLevelName );

  /**
   *
   */
  bool HaveDotLevelsFile ( ) { return mHaveDotLevelsFile; };

  /**
   *
   */
  void SetHaveDotLevelsFile ( ) { mHaveDotLevelsFile = true; };
  

  /**
   * @return vector< string >
   */
  vector<string> GetLevelNames ( );

  /**
   * Returns true if the given levelName is in our list of known level names
   * @return bool
   */
  bool ValidLevelName ( const char* iLevelName );


  /**
   * Returns the value of the discrete factor stored in this instance
   * @return string
   */
  string GetDiscreteValue ( );


  /**
   * Returns the value of the continuous factor stored in this instance.
   * If this is a discrete factor, then it returns the level number.
   * @return double
   */
  double GetContinuousValue ( );

  /**
   * Returns the value as a continous even though its a  discrete factor, 
   * value is the level number of the specified discrete value
   * @return double
   */
  double GetContinuousValue ( const char* isDiscreteValue );

private:

  // private attributes
  //

  // This is the name of column in the table.dat file containing
  // this factor data for each subject.
  string msName;

  // if ==1, continuous
  // if ==2, discrete
  int mType;

  // Names of possible levels (for instance, if this factor is 'gender',
  // then the two possible names are 'Female' and 'Male').
  vector< string > mLevelNames;

  // true if user create a factor.levels file containing the valid level names
  bool mHaveDotLevelsFile;

  string msDiscreteValue;

  double mContinuousValue;

};

#endif // QDECFACTOR_H
