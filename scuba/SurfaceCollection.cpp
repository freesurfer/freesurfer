#include "string_fixed.h"
#include <errno.h>
#include <stdexcept>
#include "SurfaceCollection.h"
#include "DataManager.h"

using namespace std;


SurfaceCollection::SurfaceCollection () :
  DataCollection() {

  mMRIS = NULL;

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "SetSurfaceCollectionFileName", 2, 
			 "collectionID fileName", 
			 "Sets the file name for a given surface collection.");
  commandMgr.AddCommand( *this, "GetSurfaceCollectionFileName", 1, 
			 "collectionID", 
			 "Gets the file name for a given surface collection.");
}

SurfaceCollection::~SurfaceCollection() {

  DataManager dataMgr = DataManager::GetManager();
  MRISLoader mrisLoader = dataMgr.GetMRISLoader();
  try { 
    mrisLoader.ReleaseData( &mMRIS );
  } 
  catch(...) {
    cerr << "Couldn't release data"  << endl;
  }
}


void
SurfaceCollection::SetSurfaceFileName ( string& ifnMRIS ) {

  mfnMRIS = ifnMRIS;
}

MRIS*
SurfaceCollection::GetMRIS () { 

  if( NULL == mMRIS ) {
    
    DataManager dataMgr = DataManager::GetManager();
    MRISLoader mrisLoader = dataMgr.GetMRISLoader();
    
    mMRIS = NULL;
    try { 
      mMRIS = mrisLoader.GetData( mfnMRIS );
    }
    catch( exception e ) {
      throw logic_error( "Couldn't load MRIS" );
    }


    // This is our RAS -> surfaceRAS transform.
    mDataToSurfaceTransform.SetMainTransform
      ( 1, 0, 0, -mMRIS->lta->xforms[0].src.c_r,
	0, 1, 0, -mMRIS->lta->xforms[0].src.c_a,
	0, 0, 1, -mMRIS->lta->xforms[0].src.c_s,
	0, 0, 0, 1 );


    CalcWorldToSurfaceTransform();


  }

  if( msLabel == "" ) {
    SetLabel( mfnMRIS );
  }

  return mMRIS;
}

TclCommandListener::TclCommandResult 
SurfaceCollection::DoListenToTclCommand ( char* isCommand,
					 int iArgc, char** iasArgv ) {

  // SetSurfaceCollectionFileName <collectionID> <fileName>
  if( 0 == strcmp( isCommand, "SetSurfaceCollectionFileName" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }
    
    if( mID == collectionID ) {
      
      string fnSurface = iasArgv[2];
      SetSurfaceFileName( fnSurface );
    }
  }
  
  // GetSurfaceCollectionFileName <collectionID>
  if( 0 == strcmp( isCommand, "GetSurfaceCollectionFileName" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }
    
    if( mID == collectionID ) {
      
      sReturnFormat = "s";
      sReturnValues = mfnMRIS;
    }
  }
  
  return DataCollection::DoListenToTclCommand( isCommand, iArgc, iasArgv );
}

void
SurfaceCollection::DoListenToMessage ( string isMessage, void* iData ) {
  
  if( isMessage == "transformChanged" ) {
    CalcWorldToSurfaceTransform();
  }

  DataCollection::DoListenToMessage( isMessage, iData );
}

ScubaROI*
SurfaceCollection::DoNewROI () {

  return NULL;
}

void
SurfaceCollection::RASToSurface  ( float iRAS[3], float oSurface[3] ) {

  mWorldToSurfaceTransform.MultiplyVector3( iRAS, oSurface );
}

void
SurfaceCollection::SurfaceToRAS  ( float iSurface[3], float oRAS[3] ) {

  mWorldToSurfaceTransform.InvMultiplyVector3( iSurface, oRAS );
}

int
SurfaceCollection::GetNumFaces () {

  if( NULL != mMRIS ) {
    return mMRIS->nfaces;
  }

  return 0;
}

int
SurfaceCollection::GetNumVerticesPerFace_Unsafe ( int inFace ) {

  return VERTICES_PER_FACE;
}

void
SurfaceCollection::GetNthVertexInFace_Unsafe ( int inFace, int inVertex, 
					       float oRAS[3] ) {

  VERTEX* vertex = &(mMRIS->vertices[mMRIS->faces[inFace].v[inVertex]]);
  float dataRAS[3];
  dataRAS[0] = vertex->x;
  dataRAS[1] = vertex->y;
  dataRAS[2] = vertex->z;

  SurfaceToRAS( dataRAS, oRAS );
}

void
SurfaceCollection::CalcWorldToSurfaceTransform () {

  mWorldToSurfaceTransform =
    mDataToSurfaceTransform * mDataToWorldTransform->Inverse();

  DataChanged();
}

