
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
  IDTracker ( ){
    mID = GetNextID();
    mIDMap[mID] = (ClassType*)this;
  }

  virtual ~IDTracker() {
    mIDMap[mID] = (ClassType*)NULL;
  }

  // Accessor for the ID.
  int GetID() const  { return mID; }

  // Static function that will return a reference to an object.
  static ClassType& FindByID( int const iID ) {
    
    ClassType* obj = mIDMap[iID];
    if( NULL == obj ) {
      std::stringstream sError;
      sError << "Couldn't find object with ID " << iID;
      throw std::out_of_range( sError.str() );
    }
    
    return *obj;
  }

  // Static function that will fill a list of active (non-NULL) IDs.
  static void GetIDList( std::list<int>& iList ) {
    
    typename std::map<int,ClassType*>::iterator tIDObject;
    for( tIDObject = mIDMap.begin(); tIDObject != mIDMap.end(); ++tIDObject ) {
      int id = (*tIDObject).first;
      ClassType* object = (*tIDObject).second;
      if( NULL != object ) {
	iList.push_back( id );
      }
    }
  }

  // Prints ID list.
  static void PrintIDList( std::ostream& os ) {

    typename std::map<int,ClassType*>::iterator tIDObject;
    for( tIDObject = mIDMap.begin(); tIDObject != mIDMap.end(); ++tIDObject ) {
      int id = (*tIDObject).first;
      ClassType* object = (*tIDObject).second;
      if( NULL != object ) {
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
  static int GetNextID () { return mNextID++; }
};


template <typename T> int IDTracker<T>::mNextID = 0;
template <typename T> std::map<int,T*> IDTracker<T>::mIDMap; 

#endif
