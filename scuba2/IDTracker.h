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
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.3 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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

  // Static function that will fill a list of pointers.
  static void GetPointerList( std::list<ClassType*>& iList ) {

    typename std::map<int,ClassType*>::const_iterator tIDObject;
    for ( tIDObject = mIDMap.begin(); tIDObject != mIDMap.end(); ++tIDObject ) 
      iList.push_back( (*tIDObject).second );
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
