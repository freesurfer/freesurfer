#include "ScubaLayer2DMRIS.h"

using namespace std;

ScubaLayer2DMRIS::ScubaLayer2DMRIS () {
  SetOutputStreamToCerr();
  mSurface = NULL;

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "Set2DMRISLayerSurfaceCollection", 2, 
			 "layerID collectionID",
			 "Sets the surface collection for this layer." );

}

ScubaLayer2DMRIS::~ScubaLayer2DMRIS () {

}

void 
ScubaLayer2DMRIS::SetSurfaceCollection ( SurfaceCollection& iSurface ) {

  mSurface = &iSurface;
}

void
ScubaLayer2DMRIS::DrawIntoBuffer ( GLubyte* iBuffer, int iWidth, int iHeight,
				   ViewState& iViewState,
				   ScubaWindowToRASTranslator& iTranslator ) {

 if( NULL == mSurface ) {
    DebugOutput( << "No surface to draw" );
    return;
  }

  MRIS* surf = mSurface->GetMRIS();

  for( int nFace = 0; nFace < surf->nfaces; nFace++ ) {

    FACE* face = &(surf->faces[nFace]);
    for( int nVertex = 0; nVertex < VERTICES_PER_FACE; nVertex++ ) {

      int nNextVertex = nVertex + 1;
      if( nNextVertex >= VERTICES_PER_FACE ) nNextVertex = 0;

      VERTEX* vertex = &(surf->vertices[face->v[nVertex]]);
      VERTEX* nextVertex = &(surf->vertices[face->v[nNextVertex]]);

      float vertexCoord, nextVertexCoord, planeCoord;
      switch( iViewState.mInPlane ) {
      case 0: // X
	vertexCoord     = vertex->x;
	nextVertexCoord = nextVertex->x;
	planeCoord      = iViewState.mCenterRAS[0];
	break;
      case 1: // Y
	vertexCoord     = vertex->y;
	nextVertexCoord = nextVertex->y;
	planeCoord      = iViewState.mCenterRAS[1];
	break;
      case 2: // Z
	vertexCoord     = vertex->z;
	nextVertexCoord = nextVertex->z;
	planeCoord      = iViewState.mCenterRAS[2];
	break;
      }
      if( (vertexCoord - planeCoord) * (nextVertexCoord - planeCoord) <= 0.0 ) {
	float f = (planeCoord - vertexCoord) / (nextVertexCoord - vertexCoord);

	float world[3];
	switch( iViewState.mInPlane ) {
	case 0: // X
	  world[0] = iViewState.mCenterRAS[0];
	  world[1] = vertex->y + f * (nextVertex->y - vertex->y);
	  world[2] = vertex->z + f * (nextVertex->z - vertex->z);
	  break;
	case 1: // Y
	  world[0] = vertex->x + f * (nextVertex->x - vertex->x);
	  world[1] = iViewState.mCenterRAS[1];
	  world[2] = vertex->z + f * (nextVertex->z - vertex->z);
	  break;
	case 2: // Z
	  world[0] = vertex->x + f * (nextVertex->x - vertex->x);
	  world[1] = vertex->y + f * (nextVertex->y - vertex->y);
	  world[2] = iViewState.mCenterRAS[2];
	  break;
	}

	int x, y;
	iTranslator.TranslateRASToWindow( world, x, y );

	if( x >= 0 && x < iWidth && y >= 0 && y < iHeight ) { 
	  GLubyte* dest = iBuffer + (iWidth * y * 4) + (x * 4);
	  dest[0] = 0;
	  dest[1] = 255;
	  dest[2] = 0;
	  dest[3] = 255;
	}
      }
    }
  }
}
  
void
ScubaLayer2DMRIS::GetInfoAtRAS ( float inX, float inY, float inZ,
			   std::map<std::string,std::string>& iLabelValues ) {

}

TclCommandManager::TclCommandResult 
ScubaLayer2DMRIS::DoListenToTclCommand ( char* isCommand, 
					 int iArgc, char** iasArgv ) {

  // Set2DMRISLayerSurfaceCollection <layerID> <collectionID>
  if( 0 == strcmp( isCommand, "Set2DMRISLayerSurfaceCollection" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }
    
    if( mID == layerID ) {
      
      int collectionID = strtol(iasArgv[2], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad collection ID";
	return error;
      }
	
      try { 
	DataCollection& data = DataCollection::FindByID( collectionID );
	if( data.GetID() != collectionID ) {
	  cerr << "IDs didn't match" << endl;
	}
	SurfaceCollection& surface = (SurfaceCollection&)data;
	// VolumeCollection& volume = dynamic_cast<VolumeCollection&>(data);
	
	SetSurfaceCollection( surface );
      }
      catch( std::bad_cast& e ) {
	DebugOutput( << "Bad cast from DataCollection" );
	sResult = "bad collection ID, collection not a surface collection";
	return error;
      }
      catch(...) {
	sResult = "bad collection ID, collection not found";
	return error;
      }
    }
  }

  return Layer::DoListenToTclCommand( isCommand, iArgc, iasArgv );
}

