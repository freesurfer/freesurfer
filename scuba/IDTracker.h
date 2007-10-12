/**
 * @file  IDTracker.h
 * @brief Template to track classes by an ID number
 *
 * A template class that can be used by subclasses so they can be
 * automatically tracked by a unique ID number. Also provides a class
 * static function for returning a list of all instances of that
 * class. All functions are defined here in the header to avoid having
 * to explicitly declare the templated instantiations.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/12 19:56:58 $
 *    $Revision: 1.8 $
 *
 * Copyright (C) 2002-2007,
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

#ifndef IDTracker_h
#define IDTracker_h

#include <list>
#include <map>
#include <sstream>
#include <stdexcept>


template <typename ClassType> 
class IDTracker {

public:
  // Ctor adds this object to the tracked map with a unique ID
  // number. Dtor removes it from the map.
  IDTracker ( ) {
    mID = GetNextID();
    mIDMap[mID] = (ClassType*)this;
  }

  virtual ~IDTracker() {
    mIDMap.erase( mID );
  }

  // Accessor for the ID.
  int GetID() const  {
    return mID;
  }

  // Static function that will return a reference to an object given
  // its ID. Will throw an exception if not found.
  static ClassType& FindByID( int const iID ) {

    typename std::map<int,ClassType*>::const_iterator tObj;
    tObj = mIDMap.find( iID );
    if( tObj == mIDMap.end() ) {
      std::stringstream sError;
      sError << "Couldn't find object with ID " << iID;
      throw std::out_of_range( sError.str() );
    }

    return *(*tObj).second;
  }

  // Static function that will fill a list of valid IDs.
  static void GetIDList( std::list<int>& ioList ) {

    typename std::map<int,ClassType*>::const_iterator tIDObject;
    for( tIDObject = mIDMap.begin(); tIDObject != mIDMap.end(); ++tIDObject )
      ioList.push_back( (*tIDObject).first );
  }

  // Prints ID list.
  static void PrintIDList( std::ostream& os ) {

    typename std::map<int,ClassType*>::const_iterator tIDObject;
    for( tIDObject = mIDMap.begin(); tIDObject != mIDMap.end(); ++tIDObject )
      os << (*tIDObject).first << " ";
  }

protected:
  
  // Copy ctor needs to get its own ID number.
  IDTracker ( IDTracker const& i ) { 
    mID = GetNextID();
    mIDMap[mID] = (ClassType*)this;
  }

  // The map of tracked objects.
  static std::map<int,ClassType*> mIDMap;

  // ID number.
  int mID;

  // ID management.
  static int mNextID;
  static int GetNextID () {
    return mNextID++;
  }
};

template <typename T> int IDTracker<T>::mNextID = 0;
template <typename T> std::map<int,T*> IDTracker<T>::mIDMap;

#endif
