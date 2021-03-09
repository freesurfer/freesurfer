/**
 * @brief Stores a factor, which can be either discrete, continuous or ignore
 *
 * An example of a discrete factor is gender (male or female) or
 * diagnosis (demented or nondemented).  An example continuous factor is
 * age, or volume of a subcortical structure.  An example of a 'factor' to
 * ignore is a column containing a string for an alternate subject ID.
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

#include <string.h>
#include <iostream>
#include <stdexcept>
#include "QdecFactor.h"


// Constructors/Destructors
//

QdecFactor::QdecFactor ( const char* isName, int iType )
{
  msName = isName;
  mType = iType;
  assert( (mType == qdecDiscreteFactorType) || 
          (mType == qdecContinuousFactorType) || 
          (mType == qdecIgnoreType) );
  mHaveDotLevelsFile = false;
  mOrdinal = false;
}

QdecFactor::QdecFactor ( const char* isName,
                         int iType,
                         const char* iValue,
                         vector< string > iLevelNames )
{
  msName = isName;
  mType = iType;
  assert( mType == qdecDiscreteFactorType );
  msDiscreteValue = iValue;
  mLevelNames = iLevelNames;
  mHaveDotLevelsFile = false;
}


QdecFactor::QdecFactor ( const char* isName,
                         int iType,
                         double iValue )
{
  msName = isName;
  mType = iType;
  assert( mType == qdecContinuousFactorType );
  mContinuousValue = iValue;
  mHaveDotLevelsFile = false;
  mOrdinal = false;
}


QdecFactor::QdecFactor ( const char* isName,
                         int iType,
                         const char* iValue )
{
  msName = isName;
  mType = iType;
  assert( mType == qdecIgnoreType );
  msIgnoreValue = iValue;
}


//Copy constructor
QdecFactor::QdecFactor ( const QdecFactor *iFactor )
{
  msName = iFactor->msName;
  mType = iFactor->mType;
  mLevelNames = iFactor->mLevelNames;
  mContinuousValue = iFactor->mContinuousValue;
  msDiscreteValue = iFactor->msDiscreteValue;
  mHaveDotLevelsFile = iFactor->mHaveDotLevelsFile;
  mOrdinal = false;
} 

QdecFactor::~QdecFactor ( )
{ }

//
// Methods
//


/**
 * GetFactorTypeName() - returns the string name of the
 * type of the given factor: 'continuous' or 'discrete'
 * @return string
 */
string QdecFactor::GetFactorTypeName ( )
{
  if (this->IsContinuous())return("continuous");
  if (this->IsDiscrete())  return("discrete");
  if (this->Ignore())  return("ignore");
  return("type-error");
}

/**
 * @return int
 * @param  isLevelName
 */
void QdecFactor::AddLevelName ( string isLevelName )
{
  assert( mType == 1 );
  
  // check if already in our list:
  if (this->ValidLevelName( isLevelName.c_str() )) return; 

  mLevelNames.push_back( isLevelName );
}

/**
 * Returns true if the given levelName is in our list of known level names
 * @return bool
 */
bool QdecFactor::ValidLevelName ( const char* iLevelName )
{
  for ( unsigned int i=0; i < mLevelNames.size(); i++ )
  {
    if ( strcmp(iLevelName, mLevelNames[i].c_str() ) == 0 ) return true;
  }
  return false;
}

/**
 * Returns the value of the continous factor stored in this instance.
 * If this is actually a discrete factor, then return the level number.
 * @return double
 */
double QdecFactor::GetContinuousValue ( )
{
  if( mType == 2 ) return mContinuousValue;

  // else its a discrete factor
  for ( unsigned int i=0; i < mLevelNames.size(); i++ )
  {
    if ( msDiscreteValue == mLevelNames[i] ) return i+1;
  }
  throw runtime_error( "ERROR: QdecFactor::GetContinuousValue failure\n" );
  return -1.0;
}

/**
 * Returns the value as a continous even though its a  discrete factor, 
 * value is the level number of the specified discrete value
 * @return double
 */
double QdecFactor::GetContinuousValue ( const char* isDiscreteValue )
{
  assert( mType == 1 );

  for ( unsigned int i=0; i < mLevelNames.size(); i++ )
  {
    if ( isDiscreteValue == mLevelNames[i] ) return i+1;
  }
  throw runtime_error( "ERROR: QdecFactor::GetContinuousValue failure\n" );
  return -1.0;
}

