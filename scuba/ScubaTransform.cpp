#include <iomanip>
#include "ScubaTransform.h"

using namespace std;

DeclareIDTracker(ScubaTransform);

ScubaTransformStaticTclListener ScubaTransform::mStaticListener;

ScubaTransform::ScubaTransform() {
  msLabel = "";

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
}

ScubaTransform::~ScubaTransform() {

}


TclCommandListener::TclCommandResult
ScubaTransform::DoListenToTclCommand( char* isCommand, int iArgc, char** iasArgv ) {

  // SetTransformLabel <transformID> <label>
  if( 0 == strcmp( isCommand, "SetTransformLabel" ) ) {
    int transformID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad transform ID";
      return error;
    }
    
    if( mID == transformID ) {
      
      string sLabel = iasArgv[2];
      SetLabel( sLabel );
    }
  }

  // GetTransformLabel <transformID>
  if( 0 == strcmp( isCommand, "GetTransformLabel" ) ) {
    int transformID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad transform ID";
      return error;
    }
    
    if( mID == transformID ) {
      sReturnFormat = "s";
      stringstream ssReturnValues;
      ssReturnValues << "\"" << GetLabel() << "\"";
      sReturnValues = ssReturnValues.str();
    }
  }

  // SetTransformValues <transformID> <listOfValues>
  if( 0 == strcmp( isCommand, "SetTransformValues" ) ) {
    int transformID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad transform ID";
      return error;
    }
    
    if( mID == transformID ) {
      
      float v[4][4];
      stringstream list( iasArgv[2] );

      for( int r = 0; r < 4; r++ ) {
	for( int c = 0; c < 4; c++ ) {
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
  if( 0 == strcmp( isCommand, "GetTransformValues" ) ) {
    int transformID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad transform ID";
      return error;
    }
    
    if( mID == transformID ) {
      
      sReturnFormat = "Lffffffffffffffffl";
      stringstream ssResult;

      for( int r = 0; r < 4; r++ ) {
	for( int c = 0; c < 4; c++ ) {
	  ssResult << m(c,r) << " ";
	}
      }
      
      sReturnValues = ssResult.str();
    }
  }

  // LoadTransformFromLTAFile <transformID> <fileName>
  if( 0 == strcmp( isCommand, "LoadTransformFromLTAFile" ) ) {
    int transformID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad transform ID";
      return error;
    }
    
    if( mID == transformID ) {
      
      string fnLTA( iasArgv[2] );
      try {
	LoadFromLTAFile( fnLTA );
      }
      catch(...) {
	sReturnFormat = "s";
	stringstream ssResult;
	ssResult << "Couldn't load transform " << fnLTA;
	sReturnValues = ssResult.str();
	return error;
      }
      
    }
  }
  
  // InvertTransform <transformID>
  if( 0 == strcmp( isCommand, "InvertTransform" ) ) {
    int transformID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad transform ID";
      return error;
    }
    
    if( mID == transformID ) {
      SetMainTransform( this->Inverse() );
    }
  }

  return ok;
}

void
ScubaTransform::ValuesChanged () {

  SendBroadcast( "transformChanged", NULL );
  Transform44::ValuesChanged();
}

ScubaTransformStaticTclListener::ScubaTransformStaticTclListener () {

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "GetTransformIDList", 0, "", 
			 "Return a list of all transformIDs." );
  commandMgr.AddCommand( *this, "MakeNewTransform", 0, "", 
			 "Creates a new transform and returns its id." );
}

ScubaTransformStaticTclListener::~ScubaTransformStaticTclListener () {
}

TclCommandListener::TclCommandResult
ScubaTransformStaticTclListener::DoListenToTclCommand ( char* isCommand, 
					       int iArgc, char** iasArgv ) {

  // GetTransformIDList
  if( 0 == strcmp( isCommand, "GetTransformIDList" ) ) {
    list<int> idList;
    ScubaTransform::GetIDList( idList );
    stringstream ssFormat;
    stringstream ssResult;
    ssFormat << "L";
    list<int>::iterator tID;
    for( tID = idList.begin(); tID != idList.end(); ++tID ) {
      int id = *tID;
      ssFormat << "i";
      ssResult << id << " ";
    }
    ssFormat << "l";
    
    sReturnFormat = ssFormat.str();
    sReturnValues = ssResult.str();
  }

  // MakeNewTransform
  if( 0 == strcmp( isCommand, "MakeNewTransform" ) ) {

    ScubaTransform* transform = new ScubaTransform();
    sReturnFormat = "i";
    stringstream ssReturnValues;
    ssReturnValues << transform->GetID();
    sReturnValues = ssReturnValues.str();
  }

  return ok;
}

ostream& 
operator <<  ( ostream& os, ScubaTransform& iTransform ) { 
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
