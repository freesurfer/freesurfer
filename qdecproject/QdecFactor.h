/**
 * @brief Stores a factor, which can be either discrete or continuous
 *
 * An example of a discrete factor is gender (male or female) or
 * diagnosis (demented or nondemented).  An example continuous factor is
 * age, or volume of a subcortical structure.
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
  static const int qdecIgnoreType = 3;

  // Constructors/Destructors
  //

  QdecFactor ( const char* isName, int iType );
  QdecFactor ( const char* isName, int iType, 
               const char* iValue,
               vector< string > iLevelNames );
  QdecFactor ( const char* isName, int iType, double iValue );  
  QdecFactor ( const char* isName, int iType, const char* iValue );
  QdecFactor ( const QdecFactor* iFactor ); //Copy constructor

  virtual ~QdecFactor ( );

  /**
   * @return bool
   */
  bool IsDiscrete ( ) 
  { if ( mType == qdecDiscreteFactorType ) return true; return false; };

  /**
   * 
   */
  void SetDiscrete ( ) { mType = qdecDiscreteFactorType; };

  /**
   * @return bool
   */
  bool IsContinuous ( ) 
  { if ( mType == qdecContinuousFactorType ) return true; return false; };

  /**
   * @return bool
   */
  bool Ignore ( ) 
  { if ( mType == qdecIgnoreType ) return true; return false; };

  /**
   * @return string
   */
  string GetFactorName ( ) { return msName; };

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
  vector<string> GetLevelNames ( ) { return mLevelNames; };

  /**
   * Returns true if the given levelName is in our list of known level names
   * @return bool
   */
  bool ValidLevelName ( const char* iLevelName );

  /**
   * @return int   Number of levels in this factor (assume its discrete)
   */
  int GetNumberOfLevels ( ) { return this->mLevelNames.size(); }

  /**
   * Returns the value of the discrete factor stored in this instance
   * @return string
   */
  string GetDiscreteValue ( ) 
  { assert(mType == qdecDiscreteFactorType); return msDiscreteValue; };

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

  /**
   * Returns the value of the 'ignore' factor stored in this instance
   * @return string
   */
  string GetIgnoreValue ( ) 
  { assert(mType == qdecIgnoreType); return msIgnoreValue; };

  /**
   *
   * Returns true if this factor is an ordinal number (not a float)
   */
  bool IsOrdinal( ) { return mOrdinal; }

  /**
   *
   * Indicate this factor is an ordinal
   */
  void SetOrdinal( ) { mOrdinal = true; };

private:

  // private attributes
  //

  // This is the name of column in the table.dat file containing
  // this factor data for each subject.
  string msName;

  // discrete, continuous or ignore
  int mType;

  string msDiscreteValue;

  double mContinuousValue;

  string msIgnoreValue;

  // Names of possible levels (for instance, if this factor is 'gender',
  // then the two possible names are 'Female' and 'Male').
  vector< string > mLevelNames;

  // true if user create a factor.levels file containing the valid level names
  bool mHaveDotLevelsFile;

  // true if user created an ordinal.factors file listing this factor,
  // in which case this factor is read as an integer (but still stored as
  // a float) and written to disk as an integer if data table is saved
  bool mOrdinal;
};

#endif // QDECFACTOR_H
