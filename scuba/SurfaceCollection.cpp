/**
 * @file  SurfaceCollection.cpp
 * @brief Implements DataCollection with a MRIS data set
 *
 * Very simple class that reads MRIS data and provides access to
 * geometry information.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.33 $
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


#include "string_fixed.h"
#include <errno.h>
#include <stdexcept>
#include "SurfaceCollection.h"
#include "DataManager.h"
#include "VolumeCollection.h"
#include "error.h"

using namespace std;

extern "C" {
#include "mri_circulars.h"
#include "mri.h"
}

SurfaceCollection::SurfaceCollection () :
    DataCollection(),

    mfnMRIS( "" ),
    mMRIS( NULL ),
    mbIsUsingVolumeForTransform( false ),
    mTransformVolume( NULL ),
    mHashTable( NULL ),
    mbBoundsCacheDirty( true ) {

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "SetSurfaceCollectionFileName", 2,
                         "collectionID fileName",
                         "Sets the file name for a given surface collection.");
  commandMgr.AddCommand( *this, "LoadSurfaceFromFileName", 1, "collectionID",
                         "Loads the surface from the file name.");
  commandMgr.AddCommand( *this, "GetSurfaceCollectionFileName", 1,
                         "collectionID",
                         "Gets the file name for a given surface collection.");
  commandMgr.AddCommand( *this, "SetSurfaceDataToSurfaceTransformFromVolume",
                         2, "collectionID volumeID",
                         "Gets the data to surface transform from a volume." );
  commandMgr.AddCommand( *this, "SetSurfaceDataToSurfaceTransformToDefault",
                         1, "collectionID",
                         "Sets the data to surface transform for a surface "
                         "to the default, which will be ../mri/orig or "
                         "identity." );
  commandMgr.AddCommand(*this,"IsSurfaceUsingDataToSurfaceTransformFromVolume",
                        1, "collectionID",
                        "Returns whether or not a surface collection is "
                        "using a volume to get its data to surface "
                        "transform." );
  commandMgr.AddCommand( *this, "GetSurfaceDataToSurfaceTransformVolume",
                         1, "collectionID",
                         "If a surface collection is using a volume "
                         "to get its data to surface transform, returns "
                         "the volume's collection ID." );
  commandMgr.AddCommand( *this, "LoadSurfacePatch", 2, "collectionID fileName",
                         "Loads a patch into a surface." );
  commandMgr.AddCommand( *this, "GetSurfaceUseRealRAS", 1, "collectionID",
                         "Returns whether or not a surface has its useRealRAS "
                         "flag on." );
}

SurfaceCollection::~SurfaceCollection() {

  DataManager dataMgr = DataManager::GetManager();
  MRISLoader mrisLoader = dataMgr.GetMRISLoader();
  try {
    mrisLoader.ReleaseData( &mMRIS );
  } catch (...) {
    cerr << "Couldn't release data"  << endl;
  }
}


void
SurfaceCollection::SetSurfaceFileName ( string& ifnMRIS ) {

  mfnMRIS = ifnMRIS;
}

void
SurfaceCollection::LoadSurface () {

  DataManager dataMgr = DataManager::GetManager();
  MRISLoader mrisLoader = dataMgr.GetMRISLoader();

  // If we already have data...
  if ( NULL != mMRIS ) {

    // Try to load this and see what we get. If it's the same as what
    // we already have, we're fine. If not, keep this one and release
    // the one we have.
    MRIS* newMRIS = NULL;
    try {
      newMRIS = mrisLoader.GetData( mfnMRIS );
    } catch ( exception& e ) {
      throw logic_error( "Couldn't load MRIS" );
    }

    if ( newMRIS == mMRIS ) {
      return;
    }

    /* Release old data. */
    try {
      mrisLoader.ReleaseData( &mMRIS );
    } catch (...) {
      cerr << "Couldn't release data"  << endl;
    }

    /* Save new data. */
    mMRIS = newMRIS;
    DataChanged();

    // Generate hash table.
    mHashTable = MHTfillVertexTableRes( mMRIS, NULL, CURRENT_VERTICES, 2.0 );

  } else {

    DataManager dataMgr = DataManager::GetManager();
    MRISLoader mrisLoader = dataMgr.GetMRISLoader();

    mMRIS = NULL;
    try {
      mMRIS = mrisLoader.GetData( mfnMRIS );
    } catch ( exception& e ) {
      throw logic_error( "Couldn't load MRIS" );
    }

    // Get transform.
#if 1
    if (mMRIS->vg.valid) {
      MRI *mri_cor = MRIalloc(mMRIS->vg.width, 
                              mMRIS->vg.height, 
                              mMRIS->vg.depth, 
                              MRI_UCHAR) ;
      MATRIX *m_RAS2SurfaceRAS ;
      
      MRIcopyVolGeomToMRI(mri_cor, &mMRIS->vg) ;
      MRIreInitCache(mri_cor) ;
      m_RAS2SurfaceRAS =  surfaceRASFromRAS_(mri_cor) ;
      if (getenv("FS_DEBUG_XFORMS"))
        {
          printf("RAS2SurfaceRAS:\n") ;
          MatrixPrint(stdout, m_RAS2SurfaceRAS) ;
        }
      mDataToSurfaceTransform.SetMainTransform( m_RAS2SurfaceRAS );
      MatrixFree( &m_RAS2SurfaceRAS ); 
      MRIfree(&mri_cor) ;
    }
#else
    if ( NULL != mMRIS->lta ) {

      // This is our RAS -> TkRegRAS transform.
      mDataToSurfaceTransform.SetMainTransform
      ( 1, 0, 0, -mMRIS->lta->xforms[0].src.c_r,
        0, 1, 0, -mMRIS->lta->xforms[0].src.c_a,
        0, 0, 1, -mMRIS->lta->xforms[0].src.c_s,
        0, 0, 0, 1 );
    }
#endif      

    CalcWorldToSurfaceTransform();

    // Generate hash table.
    mHashTable = MHTfillVertexTableRes( mMRIS, NULL, CURRENT_VERTICES, 2.0 );

  }

  if ( msLabel == "" ) {
    SetLabel( mfnMRIS );
  }
}

MRIS*
SurfaceCollection::GetMRIS () {

  if ( NULL == mMRIS ) {
    LoadSurface();
  }

  return mMRIS;
}

void
SurfaceCollection::LoadPatch ( string& ifnPatch ) {

  char* cfnPath;
  cfnPath = strdup( ifnPatch.c_str() );
  int rMRIS = MRISreadPatch( mMRIS, cfnPath );
  if ( rMRIS != NO_ERROR ) {
    throw runtime_error( "Error loading " + ifnPatch );
  }
  DataChanged();
  free( cfnPath );
}

void
SurfaceCollection::GetDataRASBounds ( float oRASBounds[6] ) {

  if ( mbBoundsCacheDirty ) {

    mRASBounds[0] = mRASBounds[2] = mRASBounds[4] = 999999;
    mRASBounds[1] = mRASBounds[3] = mRASBounds[5] = -999999;

    for ( int nVertex = 0; nVertex < GetNumVertices(); nVertex++ ) {

      VERTEX* vertex = &(mMRIS->vertices[nVertex]);
      float TkRegRAS[3], RAS[3];
      if (vertex->ripflag)
        continue ;
      TkRegRAS[0] = vertex->x;
      TkRegRAS[1] = vertex->y;
      TkRegRAS[2] = vertex->z;
      SurfaceToRAS( TkRegRAS, RAS );

      if ( RAS[0] < mRASBounds[0] ) mRASBounds[0] = RAS[0];
      if ( RAS[0] > mRASBounds[1] ) mRASBounds[1] = RAS[0];
      if ( RAS[1] < mRASBounds[2] ) mRASBounds[2] = RAS[1];
      if ( RAS[1] > mRASBounds[3] ) mRASBounds[3] = RAS[1];
      if ( RAS[2] < mRASBounds[4] ) mRASBounds[4] = RAS[2];
      if ( RAS[2] > mRASBounds[5] ) mRASBounds[5] = RAS[2];
    }

    mbBoundsCacheDirty = false;
  }

  oRASBounds[0] = mRASBounds[0];
  oRASBounds[1] = mRASBounds[1];
  oRASBounds[2] = mRASBounds[2];
  oRASBounds[3] = mRASBounds[3];
  oRASBounds[4] = mRASBounds[4];
  oRASBounds[5] = mRASBounds[5];
}

TclCommandListener::TclCommandResult
SurfaceCollection::DoListenToTclCommand ( char* isCommand,
    int iArgc, char** iasArgv ) {

  // SetSurfaceCollectionFileName <collectionID> <fileName>
  if ( 0 == strcmp( isCommand, "SetSurfaceCollectionFileName" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }

    if ( mID == collectionID ) {

      string fnSurface = iasArgv[2];
      SetSurfaceFileName( fnSurface );
    }
  }

  // LoadSurfaceFromFileName <collectionID>
  if ( 0 == strcmp( isCommand, "LoadSurfaceFromFileName" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }

    if ( mID == collectionID ) {
      LoadSurface();
    }
  }

  // GetSurfaceCollectionFileName <collectionID>
  if ( 0 == strcmp( isCommand, "GetSurfaceCollectionFileName" ) ) {
    int collectionID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad collection ID";
      return error;
    }

    if ( mID == collectionID ) {

      sReturnFormat = "s";
      sReturnValues = mfnMRIS;
    }
  }

  // SetSurfaceDataToSurfaceTransformFromVolume <collectionID> <volumeID>
  if ( 0 == strcmp( isCommand, "SetSurfaceDataToSurfaceTransformFromVolume") ) {

    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collection ID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      try {
        int volumeID =
          TclCommandManager::ConvertArgumentToInt( iasArgv[2] );

        DataCollection& col = DataCollection::FindByID( volumeID );
        VolumeCollection& vol = (VolumeCollection&)col;
        //   dynamic_cast<VolumeCollection&>( col );

        SetDataToSurfaceTransformFromVolume( vol );
      } catch ( runtime_error& e ) {
        sResult = e.what();
        return error;
      }
    }
  }

  // SetSurfaceDataToSurfaceTransformToDefault <collectionID>
  if ( 0 == strcmp( isCommand, "SetSurfaceDataToSurfaceTransformToDefault") ) {

    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collection ID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      SetDataToSurfaceTransformToDefault();
    }
  }

  // IsSurfaceUsingDataToSurfaceTransformFromVolume <collectionID>
  if ( 0 == strcmp( isCommand, "IsSurfaceUsingDataToSurfaceTransformFromVolume" ) ) {

    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collection ID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      sReturnValues =
        TclCommandManager::ConvertBooleanToReturnValue( mbIsUsingVolumeForTransform );
      sReturnFormat = "i";
    }
  }

  // GetSurfaceDataToSurfaceTransformVolume <collectionID>
  if ( 0 == strcmp( isCommand, "GetSurfaceDataToSurfaceTransformVolume" ) ) {

    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collection ID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      stringstream ssReturnValues;
      if ( NULL != mTransformVolume ) {
        ssReturnValues << mTransformVolume->GetID();
      } else {
        ssReturnValues << 0;
      }
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "i";
    }
  }

  // LoadSurfacePatch <collectionID> <fileName>
  if ( 0 == strcmp( isCommand, "LoadSurfacePatch" ) ) {

    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collection ID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      string fnPatch = iasArgv[2];
      LoadPatch( fnPatch );
    }
  }

  // GetSurfaceUseRealRAS <collectionID>
  if ( 0 == strcmp( isCommand, "GetSurfaceUseRealRAS" ) ) {

    int collectionID;
    try {
      collectionID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad collection ID: ") + e.what();
      return error;
    }

    if ( mID == collectionID ) {

      stringstream ssReturnValues;
      ssReturnValues << GetUseRealRAS();
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "i";
    }
  }

  return DataCollection::DoListenToTclCommand( isCommand, iArgc, iasArgv );
}

void
SurfaceCollection::DoListenToMessage ( string isMessage, void* iData ) {

  if ( isMessage == "transformChanged" ) {
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
SurfaceCollection::FindNearestVertexToRAS ( float iRAS[3], float* oDistance ) {

  float dataRAS[3];
  RASToSurface( iRAS, dataRAS );

  float minDistance = 1000;
  int nClosestVertex = -1;
  for ( int nVertex = 0; nVertex < GetNumVertices(); nVertex++ ) {

    VERTEX* vertex = &(mMRIS->vertices[nVertex]);
    float curDataRAS[3];
    if (vertex->ripflag)
      continue ;
    curDataRAS[0] = vertex->x;
    curDataRAS[1] = vertex->y;
    curDataRAS[2] = vertex->z;

    if ( !vertex->ripflag ) {

      float dx = dataRAS[0] - curDataRAS[0];
      float dy = dataRAS[1] - curDataRAS[1];
      float dz = dataRAS[2] - curDataRAS[2];
      float distance = sqrt( dx*dx + dy*dy + dz*dz );
      if ( distance < minDistance ) {
        minDistance = distance;
        nClosestVertex = nVertex;
      }
    }
  }

  if ( -1 == nClosestVertex ) {
    throw runtime_error( "No vertices found.");
  }

  if ( NULL != oDistance ) {
    *oDistance = minDistance;
  }
  return nClosestVertex;
}

int
SurfaceCollection::FindVertexAtRAS ( float iRAS[3], float* oDistance ) {

  float dataRAS[3];
  RASToSurface( iRAS, dataRAS );

  VERTEX v;
  v.x = dataRAS[0];
  v.y = dataRAS[1];
  v.z = dataRAS[2];
  float distance;
  int nClosestVertex =
    MHTfindClosestVertexNo( mHashTable, mMRIS, &v, &distance );

  if ( -1 == nClosestVertex ) {
    throw runtime_error( "No vertices found.");
  }

  if ( NULL != oDistance ) {
    *oDistance = distance;
  }

  return nClosestVertex;
}

int
SurfaceCollection::GetNumFaces () {

  if ( NULL != mMRIS ) {
    return mMRIS->nfaces;
  }

  return 0;
}

int
SurfaceCollection::GetNumVerticesPerFace_Unsafe ( int ) {

  return VERTICES_PER_FACE;
}

void
SurfaceCollection::GetNthVertexInFace_Unsafe ( int inFace, int inVertex,
    float oRAS[3], bool* oRipped ) {

  VERTEX* vertex = &(mMRIS->vertices[mMRIS->faces[inFace].v[inVertex]]);
  float dataRAS[3];
  dataRAS[0] = vertex->x;
  dataRAS[1] = vertex->y;
  dataRAS[2] = vertex->z;

  SurfaceToRAS( dataRAS, oRAS );

  if ( NULL != oRAS ) {
    *oRipped = vertex->ripflag;
  }
}

int
SurfaceCollection::GetNumVertices () {

  if ( NULL != mMRIS ) {
    return mMRIS->nvertices;
  }

  return 0;
}

void
SurfaceCollection::GetNthVertex_Unsafe ( int inVertex,
    float oRAS[3], bool* oRipped ) {

  VERTEX* vertex = &(mMRIS->vertices[inVertex]);
  float dataRAS[3];
  dataRAS[0] = vertex->x;
  dataRAS[1] = vertex->y;
  dataRAS[2] = vertex->z;

  SurfaceToRAS( dataRAS, oRAS );

  if ( NULL != oRipped ) {
    *oRipped = vertex->ripflag;
  }
}

bool
SurfaceCollection::GetUseRealRAS () {

  if ( NULL != mMRIS ) {
    return mMRIS->useRealRAS;
  }

  return false;
}

void
SurfaceCollection::CalcWorldToSurfaceTransform () {

  Transform44 worldToData = mDataToWorldTransform->Inverse();
  Transform44 tmp = mDataToSurfaceTransform * worldToData;
  mWorldToSurfaceTransform = tmp;

  DataChanged();
}

void
SurfaceCollection::SetDataToSurfaceTransformFromVolume( VolumeCollection&
    iVolume ) {

  if ( mbIsUsingVolumeForTransform ) {
    if ( NULL != mTransformVolume ) {
      if ( mTransformVolume->GetID() == iVolume.GetID() ) {
        return;
      }
    }
  }

  // We want to get the MRITkRegRASToRAS matrix from the volume.
  MRI* mri = const_cast<MRI*>(iVolume.GetMRI());
  if ( NULL == mri ) {
    throw runtime_error( "Couldn't get MRI from volume" );
  }

  // Enh... OK, this is what Tosa told me to change it to a while ago,
  // but it apparently stopped working after he left. So I changed it
  // back to what I originally had, which is just getting the matrix
  // with TkRegRASFromRAS_ and setting our dataToSurface transform
  // to it. So we'll see how that goes for now.
#if 0
  VOL_GEOM  surfaceGeometry;
  VOL_GEOM  volumeGeometry;

  memcpy( &surfaceGeometry, &(mMRIS->vg), sizeof(VOL_GEOM) );
  getVolGeom( mri, &volumeGeometry );

  if ( surfaceGeometry.valid ) {
    if ( !vg_isEqual( &volumeGeometry, &surfaceGeometry ) ) {
      printf( "Transforming surface to match volume geometry.\n" );


      MRI* transformMRI =
        MRIallocHeader( volumeGeometry.width, volumeGeometry.height,
                        volumeGeometry.depth, MRI_VOLUME_TYPE_UNKNOWN);

      useVolGeomToMRI( &volumeGeometry, transformMRI );

      int eMRIS = MRISsurf2surf( mMRIS, transformMRI, NULL );
      if (ERROR_NONE != eMRIS) { }

      MRIfree( &transformMRI );
    }
  }

#else
  MATRIX* RASToTkRegRASMatrix = surfaceRASFromRAS_( mri );
  if ( NULL == RASToTkRegRASMatrix ) {
    throw runtime_error( "Couldn't get surfaceRASFromRAS_ from MRI" );
  }

  mDataToSurfaceTransform.SetMainTransform( RASToTkRegRASMatrix );
  MatrixFree( &RASToTkRegRASMatrix );
#endif

  CalcWorldToSurfaceTransform();

  DataChanged();

  mbIsUsingVolumeForTransform = true;
  mTransformVolume = &iVolume;
}

void
SurfaceCollection::SetDataToSurfaceTransformToDefault () {

  if ( !mbIsUsingVolumeForTransform )
    return;

  // We want to get the lta from the mris or else use identity.
  if ( NULL != mMRIS->lta ) {

    // This is our RAS -> TkRegRAS transform.
    mDataToSurfaceTransform.SetMainTransform
    ( 1, 0, 0, -mMRIS->lta->xforms[0].src.c_r,
      0, 1, 0, -mMRIS->lta->xforms[0].src.c_a,
      0, 0, 1, -mMRIS->lta->xforms[0].src.c_s,
      0, 0, 0, 1 );
  } else {

    mDataToSurfaceTransform.MakeIdentity();
  }

  CalcWorldToSurfaceTransform();

  DataChanged();

  mbIsUsingVolumeForTransform = false;
  mTransformVolume = NULL;
}
