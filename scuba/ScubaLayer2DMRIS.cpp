/**
 * @file  ScubaLayer2DMRIS.cpp
 * @brief Draws a MRIS surface
 *
 * Draws the intersection of a surface and a plane into the
 * Layer. Also handles getting information about a surface at an RAS
 * point, and responding to Tcl commands requesting vertex
 * information.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:38 $
 *    $Revision: 1.40 $
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


#include <list>
#include "ScubaLayer2DMRIS.h"
#include "VectorOps.h"
#include "Utilities.h"

using namespace std;

char* const ScubaLayer2DMRIS::kaReportableInfo[ScubaLayer2DMRIS::kcReportableInfo] = 
  { (char*)"Vertex", (char*)"Distance" };

ScubaLayer2DMRIS::ScubaLayer2DMRIS () :
    mSurface( NULL ),
    mLineWidth( 1 ),
    mbDrawVertices( false ) {

  maLineColor[0] = 0;
  maLineColor[1] = 255;
  maLineColor[2] = 0;
  maVertexColor[0] = 255;
  maVertexColor[1] = 0;
  maVertexColor[2] = 255;

  SetOutputStreamToCerr();

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "Set2DMRISLayerSurfaceCollection", 2,
                         "layerID collectionID",
                         "Sets the surface collection for this layer." );
  commandMgr.AddCommand( *this, "Get2DMRISLayerSurfaceCollection", 1,
                         "layerID",
                         "Returns the surface collection for this layer." );
  commandMgr.AddCommand( *this, "Set2DMRISLayerLineColor", 4,
                         "layerID red green blue",
                         "Sets the line color for this layer. red, green, "
                         "and blue should be 0-255 integers." );
  commandMgr.AddCommand( *this, "Get2DMRISLayerLineColor", 1, "layerID",
                         "Returns the line color for this layer as a list "
                         " of red, green, and blue integers from 0-255." );
  commandMgr.AddCommand( *this, "Set2DMRISLayerVertexColor", 4,
                         "layerID red green blue",
                         "Sets the vertex color for this layer. red, green, "
                         "and blue should be 0-255 integers." );
  commandMgr.AddCommand( *this, "Get2DMRISLayerVertexColor", 1, "layerID",
                         "Returns the vertex color for this layer as a list "
                         " of red, green, and blue integers from 0-255." );
  commandMgr.AddCommand( *this, "Set2DMRISLayerLineWidth", 2, "layerID width",
                         "Sets the line width for this layer. width should "
                         "be an integer." );
  commandMgr.AddCommand( *this, "Get2DMRISLayerLineWidth", 1, "layerID",
                         "Returns the line width for this layer as an "
                         "integer." );
  commandMgr.AddCommand( *this, "Set2DMRISLayerDrawVertices", 2,
                         "layerID draw", "Sets whether or not we're drawing "
                         "vertices on the surface." );
  commandMgr.AddCommand( *this, "Get2DMRISLayerDrawVertices", 1, "layerID",
                         "Returns whether or not we're drawing vertices." );
  commandMgr.AddCommand( *this, "Get2DMRISRASCoordsFromVertexIndex", 2,
                         "layerID vertexIndex", "Returns as a list of RAS "
                         "coords the location of the vertex." );
  commandMgr.AddCommand( *this, "Get2DMRISNearestVertexIndex", 4,
                         "layerID x y z", "Returns the index of the vertex "
                         "closest to the input RAS coords and the distance "
                         "to that vertex." );

  // Add the items we can report. 
  for( int nInfo = 0; nInfo < kcReportableInfo; nInfo++ )
    AddReportableInfo( kaReportableInfo[nInfo] );

}

ScubaLayer2DMRIS::~ScubaLayer2DMRIS () {}

void
ScubaLayer2DMRIS::SetSurfaceCollection ( SurfaceCollection& iSurface ) {

  mSurface = &iSurface;

  mSurface->AddListener( *this );

  mSurface->GetMRIS();
}

void
ScubaLayer2DMRIS::DrawIntoBuffer ( GLubyte* iBuffer, int iWidth, int iHeight,
                                   ViewState& iViewState,
                                   ScubaWindowToRASTranslator& iTranslator ) {}

void
ScubaLayer2DMRIS::DrawIntoGL ( ViewState& iViewState,
                               ScubaWindowToRASTranslator& iTranslator ) {

  if ( NULL == mSurface ) {
    DebugOutput( << "No surface to draw" );
    return;
  }

  if ( !mbVisible ) {
    return;
  }

  // Get a point and a normal for our view plane.
  Point3<float> planeRAS( iViewState.GetCenterRAS() );
  Point3<float> planeN( iViewState.GetPlaneNormal() );

  // If we don't already have a cache of the draw list for this view
  // state...
  if ( !mCachedViewState.IsSameAs( iViewState ) ) {

    // Clear the cache.
    mCachedDrawList.clear();

    // Save info about the current view state.
    mCachedViewState.SetFrom( iViewState );

    int cIntersectionsInFace = 0;
    int intersectionPair[2][2];

    // We need to look for intersections of edges in a face and the
    // current plane.
    int cFaces = mSurface->GetNumFaces();
    for ( int nFace = 0; nFace < cFaces; nFace++ ) {

      cIntersectionsInFace = 0;

      // Look at each edge in this face...
      int cVerticesPerFace = mSurface->GetNumVerticesPerFace_Unsafe( nFace );
      for ( int nVertex = 0; nVertex < cVerticesPerFace; nVertex++ ) {

        int nNextVertex = nVertex + 1;
        if ( nNextVertex >= cVerticesPerFace ) nNextVertex = 0;

        // Get the vertices.
        Point3<float> vRAS, vnRAS;
        bool bRipped, bNextRipped;
        mSurface->GetNthVertexInFace_Unsafe( nFace, nVertex,
                                             vRAS.xyz(), &bRipped );
        mSurface->GetNthVertexInFace_Unsafe( nFace, nNextVertex,
                                             vnRAS.xyz(), &bNextRipped );

        // Don't draw ripped verts.
        if ( bRipped || bNextRipped ) {
          continue;
        }

        // If they cross the view's in plane coordinate...
        VectorOps::IntersectionResult rInt;
        Point3<float> intersectionRAS;
        rInt = VectorOps::SegmentIntersectsPlane
               ( vRAS, vnRAS, planeRAS, planeN, intersectionRAS );
        if ( VectorOps::intersect == rInt ) {

          // Translate intersection point to a window coord. If it's in
          // the view...
          int window[2];
          iTranslator.TranslateRASToWindow( intersectionRAS.xyz(), window );

          // Add this intersection window point to our pair. If we
          // have two intersections, they make a line of the
          // intersection of this face and the in plane, so add them
          // to the list of points to draw.
          intersectionPair[cIntersectionsInFace][0] = window[0];
          intersectionPair[cIntersectionsInFace][1] = window[1];

          cIntersectionsInFace++;

          if ( cIntersectionsInFace == 2 ) {
            cIntersectionsInFace = 0;
            mCachedDrawList.push_back( intersectionPair[0][0] );
            mCachedDrawList.push_back( intersectionPair[0][1] );
            mCachedDrawList.push_back( intersectionPair[1][0] );
            mCachedDrawList.push_back( intersectionPair[1][1] );
          }
        }
      }
    }
  }

  // Draw all the intersection points we just calced.
  bool bDraw = false;
  int window[2];
  int window1[2];
  int window2[2];
  list<int>::iterator tDrawList;
  window1[0]=0;
  window1[1]=0;
  for ( tDrawList = mCachedDrawList.begin();
        tDrawList != mCachedDrawList.end(); ++tDrawList ) {

    window[0] = *tDrawList;
    window[1] = *(++tDrawList);

    // First time around, just save a point, next time around, draw
    // the line they make. Also draw pixels for the intersection
    // points themselves.
    if ( !bDraw ) {

      window1[0] = window[0];
      window1[1] = window[1];
      bDraw = true;

    } else {

      window2[0] = window[0];
      window2[1] = window[1];
      bDraw = false;

      glLineWidth( mLineWidth );
      glColor3ub( maLineColor[0], maLineColor[1], maLineColor[2] );
      glBegin( GL_LINES );
      glVertex2d( window1[0], window1[1] );
      glVertex2d( window2[0], window2[1] );
      glEnd();

      // If we're drawing verts, draw a box around the location of the
      // vertex, the size of the line width so that it will scale
      // properly.
      if ( mbDrawVertices ) {
        glColor3ub( maVertexColor[0], maVertexColor[1], maVertexColor[2] );
        glBegin( GL_QUADS );
        glVertex2d( window1[0]-mLineWidth, window1[1]-mLineWidth );
        glVertex2d( window1[0]+mLineWidth, window1[1]-mLineWidth );
        glVertex2d( window1[0]+mLineWidth, window1[1]+mLineWidth );
        glVertex2d( window1[0]-mLineWidth, window1[1]+mLineWidth );
        glEnd();
      }
    }
  }
}

void
ScubaLayer2DMRIS::GetInfoAtRAS ( float iRAS[3],
                                 std::list<InfoAtRAS>& ioInfo ) {

  if ( !mbReportInfoAtRAS )
    return;

  if ( NULL != mSurface ) {

    InfoAtRAS info;
    try {
      float distance;
      int nVertex = mSurface->FindVertexAtRAS( iRAS, &distance );

      if( GetReportInfo( kaReportableInfo[Vertex] ) ) {

	stringstream ssVertex;
	ssVertex << nVertex;
	
	stringstream ssCallback;
	ssCallback << "SetCursorFromSurfaceVertexIndex " << mSurface->GetID();
	
	info.SetLabel( mSurface->GetLabel() + ",vertex" );
	info.SetValue( ssVertex.str() );
	info.SetInputFilter( "ui" );
	info.SetTclCallback( ssCallback.str() );
	ioInfo.push_back( info );
	info.Clear();
      }
	  
      if( GetReportInfo( kaReportableInfo[Distance] ) ) {
	
	stringstream ssDistance;
	ssDistance << distance;
	
	info.SetLabel( mSurface->GetLabel() + ",distance" );
	info.SetValue( ssDistance.str() );
	ioInfo.push_back( info );
	info.Clear();
      }

    } catch (...) {

      if( GetReportInfo( kaReportableInfo[Vertex] ) ) {
	
	stringstream ssCallback;
	ssCallback << "SetCursorFromSurfaceVertexIndex " << mSurface->GetID();
	
	info.SetLabel( mSurface->GetLabel() + ",vertex" );
	info.SetValue( "None" );
	info.SetInputFilter( "ui" );
	info.SetTclCallback( ssCallback.str() );
	ioInfo.push_back( info );
	info.Clear();
      }

      if( GetReportInfo( kaReportableInfo[Distance] ) ) {
	
	info.SetLabel( mSurface->GetLabel() + ",distance" );
	info.SetValue( "None" );
	ioInfo.push_back( info );
	info.Clear();
      }
    }
  }
}

void
ScubaLayer2DMRIS::DataChanged () {

  ClearCache();

  RequestRedisplay();
}

void
ScubaLayer2DMRIS::HandleTool ( float[3], ViewState& ,
                               ScubaWindowToRASTranslator& ,
                               ScubaToolState& , InputState&  ) {}

void
ScubaLayer2DMRIS::FindRASLocationOfVertex ( int inVertex, float oRAS[3] ) {

  if ( NULL == mSurface ) {
    throw runtime_error( "No surface loaded." );
  }

  if ( inVertex < mSurface->GetNumVertices() ) {
    mSurface->GetNthVertex_Unsafe( inVertex, oRAS, NULL );
  } else {
    throw runtime_error( "Vertex index is out of bounds." );
  }
}

TclCommandManager::TclCommandResult
ScubaLayer2DMRIS::DoListenToTclCommand ( char* isCommand,
    int iArgc, char** iasArgv ) {

  // Set2DMRISLayerSurfaceCollection <layerID> <collectionID>
  if ( 0 == strcmp( isCommand, "Set2DMRISLayerSurfaceCollection" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      int collectionID = strtol(iasArgv[2], (char**)NULL, 10);
      if ( ERANGE == errno ) {
        sResult = "bad collection ID";
        return error;
      }

      try {
        DataCollection& data = DataCollection::FindByID( collectionID );
        if ( data.GetID() != collectionID ) {
          cerr << "IDs didn't match" << endl;
        }
        SurfaceCollection& surface = (SurfaceCollection&)data;
        // VolumeCollection& volume = dynamic_cast<VolumeCollection&>(data);

        SetSurfaceCollection( surface );
      } catch (...) {
        sResult = "bad collection ID, collection not found";
        return error;
      }
    }
  }

  // Get2DMRISLayerSurfaceCollection <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRISLayerSurfaceCollection" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      stringstream ssReturnValues;
      ssReturnValues << (int) (mSurface->GetID());
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "i";
    }
  }

  // Set2DMRISLayerLineColor <layerID> <red> <green> <blue>
  if ( 0 == strcmp( isCommand, "Set2DMRISLayerLineColor" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      int red = strtol( iasArgv[2], (char**)NULL, 10);
      if ( ERANGE == errno ) {
        sResult = "bad red";
        return error;
      }

      int green = strtol( iasArgv[3], (char**)NULL, 10);
      if ( ERANGE == errno ) {
        sResult = "bad green";
        return error;
      }

      int blue = strtol( iasArgv[4], (char**)NULL, 10);
      if ( ERANGE == errno ) {
        sResult = "bad blue";
        return error;
      }

      int color[3];
      color[0] = red;
      color[1] = green;
      color[2] = blue;
      SetLineColor3d( color );
    }
  }

  // Get2DMRISLayerLineColor <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRISLayerLineColor" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {
      sReturnFormat = "Liiil";
      stringstream ssReturnValues;
      ssReturnValues << maLineColor[0] << " " << maLineColor[1] << " "
      << maLineColor[2];
      sReturnValues = ssReturnValues.str();
    }
  }

  // Set2DMRISLayerVertexColor <layerID> <red> <green> <blue>
  if ( 0 == strcmp( isCommand, "Set2DMRISLayerVertexColor" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      int red = strtol( iasArgv[2], (char**)NULL, 10);
      if ( ERANGE == errno ) {
        sResult = "bad red";
        return error;
      }

      int green = strtol( iasArgv[3], (char**)NULL, 10);
      if ( ERANGE == errno ) {
        sResult = "bad green";
        return error;
      }

      int blue = strtol( iasArgv[4], (char**)NULL, 10);
      if ( ERANGE == errno ) {
        sResult = "bad blue";
        return error;
      }

      int color[3];
      color[0] = red;
      color[1] = green;
      color[2] = blue;
      SetVertexColor3d( color );
    }
  }

  // Get2DMRISLayerVertexColor <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRISLayerVertexColor" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {
      sReturnFormat = "Liiil";
      stringstream ssReturnValues;
      ssReturnValues << maVertexColor[0] << " " << maVertexColor[1] << " "
      << maVertexColor[2];
      sReturnValues = ssReturnValues.str();
    }
  }

  // Set2DMRISLayerLineWidth <layerID> <width>
  if ( 0 == strcmp( isCommand, "Set2DMRISLayerLineWidth" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {

      int width = strtol( iasArgv[2], (char**)NULL, 10);
      if ( ERANGE == errno ) {
        sResult = "bad width";
        return error;
      }

      SetLineWidth( width );
    }
  }

  // Get2DMRISLayerLineWidth <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRISLayerLineWidth" ) ) {
    int layerID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad layer ID";
      return error;
    }

    if ( mID == layerID ) {
      sReturnFormat = "i";
      stringstream ssReturnValues;
      ssReturnValues << mLineWidth;
      sReturnValues = ssReturnValues.str();
    }
  }

  // Set2DMRISLayerDrawVertices <layerID> <draw>
  if ( 0 == strcmp( isCommand, "Set2DMRISLayerDrawVertices" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      try {
        bool bDraw =
          TclCommandManager::ConvertArgumentToBoolean( iasArgv[2] );
        SetDrawVertices( bDraw );
      } catch ( runtime_error& e ) {
        sResult = "bad draw \"" + string(iasArgv[2]) + "\"," + e.what();
        return error;
      }
    }
  }

  // Get2DMRISLayerDrawVertices <layerID>
  if ( 0 == strcmp( isCommand, "Get2DMRISLayerDrawVertices" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      bool bDraw = GetDrawVertices();
      sReturnValues =
        TclCommandManager::ConvertBooleanToReturnValue( bDraw );
      sReturnFormat = "i";
    }
  }

  // Get2DMRISRASCoordsFromVertexIndex <layerID> <vertexIndex>
  if ( 0 == strcmp( isCommand, "Get2DMRISRASCoordsFromVertexIndex" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      int nVertex;
      nVertex = TclCommandManager::ConvertArgumentToInt( iasArgv[2] );

      float ras[3];
      FindRASLocationOfVertex( nVertex, ras );

      stringstream ssReturnValues;
      ssReturnValues << ras[0] << " " << ras[1] << " " << ras[2];
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "Lfffl";

      return ok;
    }
  }

  // Get2DMRISNearestVertexIndex <layerID> <x> <y> <z>
  if ( 0 == strcmp( isCommand, "Get2DMRISNearestVertexIndex" ) ) {
    int layerID;
    try {
      layerID = TclCommandManager::ConvertArgumentToInt( iasArgv[1] );
    } catch ( runtime_error& e ) {
      sResult = string("bad layerID: ") + e.what();
      return error;
    }

    if ( mID == layerID ) {

      float ras[3];
      try {
        ras[0] = TclCommandManager::ConvertArgumentToFloat( iasArgv[2] );
        ras[1] = TclCommandManager::ConvertArgumentToFloat( iasArgv[3] );
        ras[2] = TclCommandManager::ConvertArgumentToFloat( iasArgv[4] );
      } catch ( runtime_error& e ) {
        sResult = string("bad RAS coord: ") + e.what();
        return error;
      }

      try {
        int nVertex;
        float distance;
        nVertex = mSurface->FindNearestVertexToRAS( ras, &distance );

        stringstream ssReturnValues;
        ssReturnValues << nVertex << " " << distance;
        sReturnValues = ssReturnValues.str();
        sReturnFormat = "Lifl";
      } catch ( runtime_error& e) {
        sResult = string("Couldn't find vertex: ") + e.what();
        return error;
      }

      return ok;
    }
  }

  return Layer::DoListenToTclCommand( isCommand, iArgc, iasArgv );
}

void
ScubaLayer2DMRIS::ClearCache () {

  // Just set one bogus value so it won't compare as the same next
  // time through.
  mCachedViewState.SetZoomLevel( -1 );
}

void
ScubaLayer2DMRIS::ProcessOption ( string isOption, string isValue ) {

  char sValue[1024];
  strcpy( sValue, isValue.c_str() );

  if ( 0 == isOption.compare( "linecolor" ) ||
       0 == isOption.compare( "vertexcolor" ) ) {
    // Value will be "r,g,b" so we need to separate out the three
    // values.
    vector<string> lResults;
    Utilities::SplitString( isValue, ",", lResults );
    if ( 3 != lResults.size() ) {
      throw runtime_error( "Couldn't parse three values from value string" );
    }

    // Covert each to a number 0-255.
    int color[3];
    color[0] = strtol( lResults[0].c_str(), (char**)NULL, 10 );
    if ( ERANGE == errno ) {
      throw runtime_error( "Couldn't convert first value" );
    }
    if ( color[0] < 0 || color[0] > 255 ) {
      throw runtime_error( "First value is out of range" );
    }

    color[1] = strtol( lResults[1].c_str(), (char**)NULL, 10 );
    if ( ERANGE == errno ) {
      throw runtime_error( "Couldn't convert second value" );
    }
    if ( color[1] < 0 || color[1] > 255 ) {
      throw runtime_error( "Second value is out of range" );
    }

    color[2] = strtol( lResults[2].c_str(), (char**)NULL, 10 );
    if ( ERANGE == errno ) {
      throw runtime_error( "Couldn't convert third value" );
    }
    if ( color[2] < 0 || color[2] > 255 ) {
      throw runtime_error( "Third value is out of range" );
    }

    // Set the color.
    if ( 0 == isOption.compare( "linecolor" ) ) {
      SetLineColor3d( color );
    } else if ( 0 == isOption.compare( "vertexcolor" ) ) {
      SetVertexColor3d( color );
    }

  } else if ( 0 == isOption.compare( "linewidth" ) ) {
    int width = (int) strtol( sValue, (char**)NULL, 10 );
    if ( ERANGE == errno ) {
      throw runtime_error( "Bad width value" );
    }
    SetLineWidth( width );

  } else {

    return Layer::ProcessOption( isOption, isValue );
  }
}
