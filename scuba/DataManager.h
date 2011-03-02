/**
 * @file  DataManager.h
 * @brief Reference counting data manager
 *
 * Templated object DataLoader that provides internal reference
 * counting for data that can be accessed by the client. The MRILoader
 * and MRISLoader objects, also defined here, implement that object
 * and allow access to MRI and MRIS objects. If multiple clients
 * access the same data, it is not loaded twice, and is only released
 * when all clients are not using it any more. A DataManager provides
 * singleton access to the loaders.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:36 $
 *    $Revision: 1.13 $
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

#ifndef DataManager_h
#define DataManager_h

#include <stdlib.h>
#include <list>
#include <map>
#include "string_fixed.h"
#include "DebugReporter.h"


extern "C" {
#include "mri.h"
#include "mrisurf.h"
#ifdef X
#undef X
#endif
#ifdef Y
#undef Y
#endif
#ifdef Z
#undef Z
#endif
}


// This is the generic Loader object. It provides basic algorithm for
// keeping references to loaded data of this type, and virtual functions
// for actually loading the data and freeing it.
template <typename T>
class DataLoader : public DebugReporter {
public:

  // Returns an instance to the data for the given file name, or
  // throws an error if it can't be loaded.
  T* GetData( std::string const& ifnData );

  // An object should call this function when it is done with the
  // data. It may be freed, or not if it's in use by other objects. If
  // the data i not in our list, don't do anything.
  void ReleaseData( T** ioData );

  // Returns the number of data of this type actually loaded. (Not
  // references to a specific instance of data.
  int CountLoaded() const {
    return mlData.size();
  }

  // Returns the number of references to this data.
  int CountReferences( T const* iData ) const;

  virtual ~DataLoader() {};
protected:
  DataLoader();

  // Override this function to load a specfic type of data given a
  // file name. Should throw an error if it's not loadable.
  virtual T* LoadData( std::string const& ifnData ) = 0;

  // Override this function to free a specific type of data. Should
  // set the ioData parameter to NULL if freed.
  virtual void FreeData( T** ioData ) = 0;

  // Comparison function for determining if a specific data instance
  // matches what would be loaded for this filename.
  virtual bool DoesFileNameMatchObject( T const* iData,
					std::string const& ifnData ) const = 0;

  // List of data and number of references for each data.
  static std::list<T*> mlData;
  static std::map<T*,int> maRefs;
};


// An implementation of the generic data loader using the MRI struct.
class MRILoader : public DataLoader<MRI> {
public:
  virtual ~MRILoader() {};
protected:

  // Load in an MRI structure with MRIread() from libutils.
  MRI* LoadData( std::string const& ifnData );

  // Free an MRI structure with MRIfree() from libutils.
  void FreeData( MRI** ioMRI );

  // Compare the internal file name of the MRI with ifnData.
  bool DoesFileNameMatchObject( MRI const* iData, 
				std::string const& ifnData ) const;
};


// An implementation of the generic data loader using the MRIS struct.
class MRISLoader : public DataLoader<MRIS> {
public:
  virtual ~MRISLoader() {}
protected:
  MRIS* LoadData( std::string const& ifnData );
  void FreeData( MRIS** ioMRI ) ;
  bool DoesFileNameMatchObject( MRIS const* iData,
				std::string const& ifnData ) const;
};


// A class specific for scuba to access our data loaders.
class DataManager : public DebugReporter {

public:

  static DataManager& GetManager();

  static MRILoader& GetMRILoader();
  static MRISLoader& GetMRISLoader();

protected:
  DataManager();
};

// Generate code for our structures.
template <typename T> std::list<T*> DataLoader<T>::mlData;
template <typename T> std::map<T*,int> DataLoader<T>::maRefs;


#endif

