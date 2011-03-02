/**
 * @file  ScubaTransform.cpp
 * @brief A Transform44 object wrapped for Scuba
 *
 * This is a Transform44 that behaves like a Scuba object, with
 * IDTrackeing, Tcl scripting, and Broadcasting of its changes. It
 * also manages the special case of when a transform is used as a
 * registration between two volumes, taking into account the source
 * and destination voxel transforms.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.18 $
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


#include <iomanip>
#include "ScubaTransform.h"
#include "VolumeCollection.h"

using namespace std;

ScubaTransformStaticTclListener ScubaTransform::mStaticListener;

ScubaTransform::ScubaTransform() :
  Broadcaster( "ScubaTransform" ),
  msLabel( "" ),
  mIsRegistration( false ),
  mSourceVolume( NULL ),
  mDestVolume( NULL ) {
  
  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "SetTransformLabel", 2, "transformID label",
                         "Set the label for a transform." );
  commandMgr.AddCommand( *this, "GetTransformLabel", 1, "transformID",
                         "Returns the label for a transform." );
  commandMgr.AddCommand( *this, "SetTransformValues", 2,
                         "transformID listOfValues",
                         "Sets the values of a transform. listOfValues "
                         "should be a list of 16 values in column order." );
  commandMgr.AddCommand( *this, "GetTransformValues", 1, "transformID",
                         "Returns a list of transform values in "
                         "column order." );
  commandMgr.AddCommand( *this, "LoadTransformFromLTAFile", 2,
                         "transformID LTAFileName", "Loads an LTA from a "
                         "file into an existing transform." );
  commandMgr.AddCommand( *this, "InvertTransform", 1, "transformID",
                         "Inverts a transform." );
  commandMgr.AddCommand( *this, "SetTransformToInverseOfTransform", 2,
                         "transformID inverseTransformID",
                         "Sets a transform to be the inverse of another "
                         "transform. Also sets its label appropriately." );
  commandMgr.AddCommand( *this, "TreatTransformAsRegistration", 3,
                         "transformID sourceVolumeID destVolumeID",
                         "Treats a transform as a tkregistration between "
                         "two volumes, getting the necessary geometry "
                         "information from them to register them.." );
  commandMgr.AddCommand( *this, "TreatTransformAsNative", 1, "transformID",
                         "If a transform was set to be treated as a "
                         "registration, this restores it to normal." );
  commandMgr.AddCommand( *this, "IsTransformRegistration", 1, "transformID",
                         "Returns whether or not a transform is being "
                         "treated as a registration." );
  commandMgr.AddCommand( *this, "GetTransformRegistrationSource", 1,
                         "transformID",
                         "Returns the id of the source volume if the "
                         "transform is being treated as a registration." );
  commandMgr.AddCommand( *this, "GetTransformRegistrationDest", 1,
                         "transformID",
                         "Returns the id of the dest volume if the "
                         "transform is being treated as a registration." );
}

ScubaTransform::~ScubaTransform() {}


TclCommandListener::TclCommandResult
ScubaTransform::DoListenToTclCommand ( char* isCommand,
				       int, char** iasArgv ) {

  // SetTransformLabel <transformID> <label>
  if ( 0 == strcmp( isCommand, "SetTransformLabel" ) ) {
    int transformID;
    try {
      transformID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad transform ID: ") + e.what();
      return error;
    }

    if ( mID == transformID ) {

      string sLabel = iasArgv[2];
      SetLabel( sLabel );
    }
  }

  // GetTransformLabel <transformID>
  if ( 0 == strcmp( isCommand, "GetTransformLabel" ) ) {
    int transformID;
    try {
      transformID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad transform ID: ") + e.what();
      return error;
    }

    if ( mID == transformID ) {
      sReturnFormat = "s";
      stringstream ssReturnValues;
      ssReturnValues << "\"" << GetLabel() << "\"";
      sReturnValues = ssReturnValues.str();
    }
  }

  // SetTransformValues <transformID> <listOfValues>
  if ( 0 == strcmp( isCommand, "SetTransformValues" ) ) {
    int transformID;
    try {
      transformID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad transform ID: ") + e.what();
      return error;
    }

    if ( mID == transformID ) {

      float v[4][4];
      stringstream list( iasArgv[2] );

      for ( int r = 0; r < 4; r++ ) {
        for ( int c = 0; c < 4; c++ ) {
          list >> v[c][r];
        }
      }

      SetMainTransform( v[0][0], v[1][0], v[2][0], v[3][0],
                        v[0][1], v[1][1], v[2][1], v[3][1],
                        v[0][2], v[1][2], v[2][2], v[3][2],
                        v[0][3], v[1][3], v[2][3], v[3][3] );
    }
  }

  // GetTransformValues <transformID>
  if ( 0 == strcmp( isCommand, "GetTransformValues" ) ) {
    int transformID;
    try {
      transformID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad transform ID: ") + e.what();
      return error;
    }

    if ( mID == transformID ) {

      sReturnFormat = "Lffffffffffffffffl";
      stringstream ssResult;

      for ( int r = 0; r < 4; r++ ) {
        for ( int c = 0; c < 4; c++ ) {
          ssResult << m(c,r) << " ";
        }
      }

      sReturnValues = ssResult.str();
    }
  }

  // LoadTransformFromLTAFile <transformID> <fileName>
  if ( 0 == strcmp( isCommand, "LoadTransformFromLTAFile" ) ) {
    int transformID;
    try {
      transformID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad transform ID: ") + e.what();
      return error;
    }

    if ( mID == transformID ) {

      string fnLTA( iasArgv[2] );
      try {
        LoadFromLTAFile( fnLTA );
      } catch (...) {
        sReturnFormat = "s";
        stringstream ssResult;
        ssResult << "Couldn't load transform " << fnLTA;
        sReturnValues = ssResult.str();
        return error;
      }

    }
  }

  // InvertTransform <transformID>
  if ( 0 == strcmp( isCommand, "InvertTransform" ) ) {
    int transformID;
    try {
      transformID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad transform ID: ") + e.what();
      return error;
    }

    if ( mID == transformID ) {
      Transform44 inverse = this->Inverse();
      SetMainTransform( inverse );
    }
  }

  // SetTransformToInverseOfTransform <transformID> <inverseTransformID>
  if ( 0 == strcmp( isCommand, "SetTransformToInverseOfTransform" ) ) {
    int transformID;
    try {
      transformID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad transform ID: ") + e.what();
      return error;
    }

    if ( mID == transformID ) {

      int inverseTransformID;
      try {
        inverseTransformID =
          TclCommandManager::ConvertArgumentToInt( iasArgv[2] );
      } catch ( runtime_error& e ) {
        sResult = string("bad inverseTransformID: ") + e.what();
        return error;
      }

      ScubaTransform& invertMe =
        ScubaTransform::FindByID( inverseTransformID );

      Transform44 inverse = invertMe.Inverse();
      SetMainTransform( inverse );
      SetLabel( string(invertMe.GetLabel()) + " -1" );
    }
  }

  // TreatTransformAsRegistration <transformID> <sourceVolumeID> <destVolumeID>
  if ( 0 == strcmp( isCommand, "TreatTransformAsRegistration" ) ) {

    int transformID;
    try {
      transformID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad transform ID: ") + e.what();
      return error;
    }

    if ( mID == transformID ) {

      try {
        int sourceVolumeID =
          TclCommandManager::ConvertArgumentToInt( iasArgv[2] );

        DataCollection& sourceCol =
          VolumeCollection::FindByID( sourceVolumeID );
        VolumeCollection& sourceVolume = 
	  dynamic_cast<VolumeCollection&>( sourceCol );

        int destVolumeID =
          TclCommandManager::ConvertArgumentToInt( iasArgv[3] );

        DataCollection& destCol =
          VolumeCollection::FindByID( destVolumeID );
        VolumeCollection& destVolume = 
	  dynamic_cast<VolumeCollection&>( destCol );

        TreatAsRegistration( sourceVolume, destVolume );
      } catch ( runtime_error& e ) {
        sResult = e.what();
        return error;
      }
    }
  }

  // TreatTransformAsNative <transformID>
  if ( 0 == strcmp( isCommand, "TreatTransformAsNative" ) ) {

    int transformID;
    try {
      transformID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad transform ID: ") + e.what();
      return error;
    }

    if ( mID == transformID ) {

      try {
        TreatAsNative();
      } catch ( runtime_error& e ) {
        sResult = e.what();
        return error;
      }
    }
  }

  // IsTransformRegistration <transformID>
  if ( 0 == strcmp( isCommand, "IsTransformRegistration" ) ) {

    int transformID;
    try {
      transformID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad transform ID: ") + e.what();
      return error;
    }

    if ( mID == transformID ) {

      sReturnValues =
        TclCommandManager::ConvertBooleanToReturnValue( mIsRegistration );
      sReturnFormat = "i";
    }
  }

  // GetTransformRegistrationSource <transformID>
  if ( 0 == strcmp( isCommand, "GetTransformRegistrationSource" ) ) {

    int transformID;
    try {
      transformID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad transform ID: ") + e.what();
      return error;
    }

    if ( mID == transformID ) {

      stringstream ssReturnValues;
      if ( NULL != mSourceVolume ) {
        ssReturnValues << mSourceVolume->GetID();
      } else {
        ssReturnValues << 0;
      }
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "i";
    }
  }

  // GetTransformRegistrationDest <transformID>
  if ( 0 == strcmp( isCommand, "GetTransformRegistrationDest" ) ) {

    int transformID;
    try {
      transformID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad transform ID: ") + e.what();
      return error;
    }

    if ( mID == transformID ) {

      stringstream ssReturnValues;
      if ( NULL != mDestVolume ) {
        ssReturnValues << mDestVolume->GetID();
      } else {
        ssReturnValues << 0;
      }
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "i";
    }
  }

  return ok;
}

void
ScubaTransform::SetLabel ( std::string const& isLabel ) {
  
  msLabel = isLabel;
}

string const&
ScubaTransform::GetLabel() const {
  
  return msLabel;
}

string
ScubaTransform::GetValuesAsString () const {

  stringstream ssResult;
  for ( int r = 0; r < 4; r++ ) {
    for ( int c = 0; c < 4; c++ ) {
      ssResult << m(c,r) << " ";
    }
  }
  return ssResult.str();
}

ScubaTransform& 
ScubaTransform::operator= ( ScubaTransform const& iTransform ) {

  // Check for self-assignment.
  if( this == &iTransform )
    return *this;

  // Use our superclass assignment operator.
  this->Transform44::operator=( iTransform );

  // Return a pointer to ourself.
  return *this;
}


int
ScubaTransform::GetRegistrationSource () const {

  if ( mSourceVolume ) {
    return mSourceVolume->GetID();
  } else {
    return -1;
  }
}

int
ScubaTransform::GetRegistrationDestination () const {

  if ( mDestVolume ) {
    return mDestVolume->GetID();
  } else {
    return -1;
  }
}

void
ScubaTransform::ValuesChanged () {

  // Call our superclass function.
  Transform44::ValuesChanged();

  // Broadcast our changed message so others can track us.
  SendBroadcast( "transformChanged", NULL );

}

void
ScubaTransform::TreatAsRegistration ( VolumeCollection const& iSourceVolume,
                                      VolumeCollection const& iDestVolume ) {

  // If we are already a registration, and we're using these IDs, just
  // return. Otherwise go to native to undo the current registration
  // first.
  if ( mIsRegistration ) 
  {
    if ( iSourceVolume.GetID() == mSourceVolume->GetID() &&
         iDestVolume.GetID()   == mDestVolume->GetID() ) 
      return;
    else 
      TreatAsNative();
  }

  // Get the MRIs from the volumes.
  MRI* sourceMRI = const_cast<MRI*>(iSourceVolume.GetMRI());
  if ( NULL == sourceMRI )
    throw runtime_error( "Couldn't get source MRI." );
  MRI* destMRI = const_cast<MRI*>(iDestVolume.GetMRI());
  if ( NULL == destMRI )
    throw runtime_error( "Couldn't get dest MRI." );

  MATRIX* tkReg = const_cast<MATRIX*>(GetMainMatrix().GetMatrix());
  MATRIX* native = MRItkReg2Native( sourceMRI, destMRI, tkReg );
  if ( NULL == native )
    throw runtime_error( "Couldn't get native transform." );
  native = MatrixInverse( native, NULL );
  SetMainTransform( native );
  MatrixFree( &native );

  mIsRegistration = true;
  mSourceVolume = &iSourceVolume;
  mDestVolume = &iDestVolume;
}

void
ScubaTransform::TreatAsNative () {

  // Make sure it's not already native.
  if ( !mIsRegistration ) {
    return;
  }

  // Get the MRIs from the volumes.
  MRI* sourceMRI = const_cast<MRI*>(mSourceVolume->GetMRI());
  if ( NULL == sourceMRI )
    throw runtime_error( "Couldn't get source MRI." );
  MRI* destMRI = const_cast<MRI*>(mDestVolume->GetMRI());
  if ( NULL == destMRI )
    throw runtime_error( "Couldn't get dest MRI." );

  // Convert our matrix.
  MATRIX* native = const_cast<MATRIX*>(GetMainMatrix().GetMatrix());
  native = MatrixInverse( native, NULL );
  MATRIX* tkReg = MRItkRegMtx( sourceMRI, destMRI, native );
  if ( NULL == tkReg )
    throw runtime_error( "Couldn't get tkReg transform." );
  SetMainTransform( tkReg );
  MatrixFree( &tkReg );

  mIsRegistration = false;
  mSourceVolume = NULL;
  mDestVolume = NULL;
}


ScubaTransformStaticTclListener::ScubaTransformStaticTclListener () {

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "GetTransformIDList", 0, "",
                         "Return a list of all transformIDs." );
  commandMgr.AddCommand( *this, "MakeNewTransform", 0, "",
                         "Creates a new transform and returns its id." );
}

ScubaTransformStaticTclListener::~ScubaTransformStaticTclListener () {}

TclCommandListener::TclCommandResult
ScubaTransformStaticTclListener::DoListenToTclCommand ( char* isCommand,
							int, char** ) {
  
  // GetTransformIDList
  if ( 0 == strcmp( isCommand, "GetTransformIDList" ) ) {

    // Get the ID list from the IDTracker, then make a list in the
    // output and go through and insert the IDs.
    list<int> idList;
    ScubaTransform::GetIDList( idList );
    stringstream ssFormat;
    stringstream ssResult;
    ssFormat << "L";
    list<int>::const_iterator tID;
    for ( tID = idList.begin(); tID != idList.end(); ++tID ) {
      int id = *tID;
      ssFormat << "i";
      ssResult << id << " ";
    }
    ssFormat << "l";

    sReturnFormat = ssFormat.str();
    sReturnValues = ssResult.str();
  }

  // MakeNewTransform
  if ( 0 == strcmp( isCommand, "MakeNewTransform" ) ) {

    // Make a new transform object and return the ID.
    ScubaTransform* transform = new ScubaTransform();
    sReturnFormat = "i";
    stringstream ssReturnValues;
    ssReturnValues << transform->GetID();
    sReturnValues = ssReturnValues.str();
  }

  return ok;
}

ostream&
operator <<  ( ostream& os, ScubaTransform const& iTransform ) {
  os << "Transform " << iTransform.GetLabel() << ":" << endl;
  os << setw(6) << iTransform(0,0) << " " << setw(6) << iTransform(1,0) << " "
  << setw(6) << iTransform(2,0) << " " << setw(6) << iTransform(3,0) << endl;
  os << setw(6) << iTransform(0,1) << " " << setw(6) << iTransform(1,1) << " "
  << setw(6) << iTransform(2,1) << " " << setw(6) << iTransform(3,1) << endl;
  os << setw(6) << iTransform(0,2) << " " << setw(6) << iTransform(1,2) << " "
  << setw(6) << iTransform(2,2) << " " << setw(6) << iTransform(3,2) << endl;
  os << setw(6) << iTransform(0,3) << " " << setw(6) << iTransform(1,3) << " "
  << setw(6) << iTransform(2,3) << " " << setw(6) << iTransform(3,3) << endl;
  return os;
}
