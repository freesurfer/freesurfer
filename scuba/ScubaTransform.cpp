#include "ScubaTransform.h"

using namespace std;

DeclareIDTracker(ScubaTransform);

bool ScubaTransform::mbRegisteredStaticListener = false;
ScubaTransformStaticTclListener ScubaTransform::mStaticListener;

ScubaTransform::ScubaTransform() {
  msLabel = "";
  m = MatrixIdentity( 4, NULL );
  mInv = MatrixInverse( m, NULL );
  mTmpVec4Src = VectorAlloc( 4, MATRIX_REAL );
  mTmpVec4Dest = VectorAlloc( 4, MATRIX_REAL );

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
  if( !mbRegisteredStaticListener ) {
    commandMgr.AddCommand( mStaticListener, "GetTransformIDList", 0, "", 
			   "Return a list of all transformIDs." );
    commandMgr.AddCommand( mStaticListener, "MakeNewTransform", 0, "", 
			   "Creates a new transform and returns its id." );
    mbRegisteredStaticListener = true;
  }
}

ScubaTransform::~ScubaTransform() {

}

void
ScubaTransform::SetTransform ( float i0j0, float i1j0, float i2j0, float i3j0,
			  float i0j1, float i1j1, float i2j1, float i3j1,
			  float i0j2, float i1j2, float i2j2, float i3j2,
			  float i0j3, float i1j3, float i2j3, float i3j3 ) {

  SetCR(0,0,i0j0);  SetCR(1,0,i1j0);  SetCR(2,0,i2j0);  SetCR(3,0,i3j0);
  SetCR(0,1,i0j1);  SetCR(1,1,i1j1);  SetCR(2,1,i2j1);  SetCR(3,1,i3j1);
  SetCR(0,2,i0j2);  SetCR(1,2,i1j2);  SetCR(2,2,i2j2);  SetCR(3,2,i3j2);
  SetCR(0,3,i0j3);  SetCR(1,3,i1j3);  SetCR(2,3,i2j3);  SetCR(3,3,i3j3);

  ValuesChanged();
}

void 
ScubaTransform::MakeIdentity () {

  MatrixIdentity( 4, m );
}

void 
ScubaTransform::MultiplyVector3 ( float iVector[3], float oVector[3] ) {

  VECTOR_ELT( mTmpVec4Src, 1 ) = iVector[0];
  VECTOR_ELT( mTmpVec4Src, 2 ) = iVector[1];
  VECTOR_ELT( mTmpVec4Src, 3 ) = iVector[2];
  VECTOR_ELT( mTmpVec4Src, 4 ) = 1;
  MatrixMultiply( m, mTmpVec4Src, mTmpVec4Dest );
  oVector[0] = VECTOR_ELT( mTmpVec4Dest, 1 );
  oVector[1] = VECTOR_ELT( mTmpVec4Dest, 2 );
  oVector[2] = VECTOR_ELT( mTmpVec4Dest, 3 );
}

void 
ScubaTransform::InvMultiplyVector3 ( float iVector[3], float oVector[3] ) {

  VECTOR_ELT( mTmpVec4Src, 1 ) = iVector[0];
  VECTOR_ELT( mTmpVec4Src, 2 ) = iVector[1];
  VECTOR_ELT( mTmpVec4Src, 3 ) = iVector[2];
  VECTOR_ELT( mTmpVec4Src, 4 ) = 1;
  MatrixMultiply( mInv, mTmpVec4Src, mTmpVec4Dest );
  oVector[0] = VECTOR_ELT( mTmpVec4Dest, 1 );
  oVector[1] = VECTOR_ELT( mTmpVec4Dest, 2 );
  oVector[2] = VECTOR_ELT( mTmpVec4Dest, 3 );
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

      for( int r = 0; r < 4; r++ ) {
	for( int c = 0; c < 4; c++ ) {
	  SetCR( c, r, v[c][r] );
	}
      }

      ValuesChanged();
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
	  ssResult << GetCR(c,r) << " ";
	}
      }
      
      sReturnValues = ssResult.str();
    }
  }
}

void
ScubaTransform::ValuesChanged () {

  CalculateInverse();
  SendBroadcast( "transformChanged", NULL );
}

void
ScubaTransform::CalculateInverse () {

  if( NULL != m ) {
    MatrixInverse( m, mInv );
  }
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

ScubaTransformStaticTclListener::~ScubaTransformStaticTclListener () {
}

