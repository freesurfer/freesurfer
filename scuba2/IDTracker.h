/**
 * @file  IDTracker.h
 * @brief A mixin class for assigning an ID to a class
 *
 * Automatically assigns an ID to an object unique within that class,
 * and provides FindByID() and GetIDList() functions, as well as a
 * GetPointerList() function for returning all references to a tracked
 * object of a certain class.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/04/06 22:23:04 $
 *    $Revision: 1.1 $
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


// A template class that can be used by subclasses so they can be
// automatically tracked by a unique ID number. To use, declare your
// class like this:
//
//    class myClass : public IDTracker<myClass> {
//
// And in the myClass.cpp file, add these lines near the top to
// generate the template code and initialize the static variables:
//
//    template IDTracker<myClass>;
//    int IDTracker<myClass>::mNextID = 0;
//    std::map<int,myClass*> IDTracker<myClass>::mIDMap;

#define DeclareIDTracker(ClassType)
//template IDTracker<ClassType>;

template <typename ClassType> class IDTracker {

public:
  // Ctor adds this object to the tracked map with a unique ID
  // number. Dtor removes it from the map.
  IDTracker ( ) {
    mID = GetNextID();
    mIDMap[mID] = (ClassType*)this;
  }

  virtual ~IDTracker() {
    mIDMap[mID] = (ClassType*)NULL;
  }

  // Accessor for the ID.
  int GetID() const  {
    return mID;
  }

  // Static function that will return a reference to an object.
  static ClassType& FindByID( int const iID ) {

    ClassType* obj = mIDMap[iID];
    if ( NULL == obj ) {
      std::stringstream sError;
      sError << "Couldn't find object with ID " << iID;
      throw std::out_of_range( sError.str() );
    }

    return *obj;
  }

  // Static function that will fill a list of active (non-NULL) IDs.
  static void GetIDList( std::list<int>& iList ) {

    typename std::map<int,ClassType*>::iterator tIDObject;
    for ( tIDObject = mIDMap.begin(); tIDObject != mIDMap.end(); ++tIDObject ) {
      int id = (*tIDObject).first;
      ClassType* object = (*tIDObject).second;
      if ( NULL != object ) {
        iList.push_back( id );
      }
    }
  }

  // Static function that will fill a list of active (non-NULL) pointers.
  static void GetPointerList( std::list<ClassType*>& iList ) {

    typename std::map<int,ClassType*>::iterator tIDObject;
    for ( tIDObject = mIDMap.begin(); tIDObject != mIDMap.end(); ++tIDObject ) {
      ClassType* object = (*tIDObject).second;
      if ( NULL != object ) {
        iList.push_back( object );
      }
    }
  }

  // Prints ID list.
  static void PrintIDList( std::ostream& os ) {

    typename std::map<int,ClassType*>::iterator tIDObject;
    for ( tIDObject = mIDMap.begin(); tIDObject != mIDMap.end(); ++tIDObject ) {
      int id = (*tIDObject).first;
      ClassType* object = (*tIDObject).second;
      if ( NULL != object ) {
        os << id << " ";
      }
    }
  }

protected:

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
