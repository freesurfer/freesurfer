#include <errno.h>
#include "string_fixed.h"
#include "ScubaToolState.h"

using namespace std;

template IDTracker<ScubaToolState>;
int IDTracker<ScubaToolState>::mNextID = 0;
map<int,ScubaToolState*> IDTracker<ScubaToolState>::mIDMap;


ScubaToolState::ScubaToolState() {
  mMode = navigation;
  mBrushRadius = 0.5;
  mBrushShape = circle;
  mbBrush3D = false;
  mbFloodStopAtLines = true;
  mbFloodStopAtROIs = true;
  mFloodFuzziness = 0;
  mFloodMaxDistance = 0;
  mbFlood3D = true;
  mEdgeLineStraightBias = 0.9;
  mEdgeLineEdgeBias = 0.9;

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "SetToolMode", 2, "toolID mode",
			 "Sets the current mode of a tool." );
  commandMgr.AddCommand( *this, "GetToolMode", 1, "toolID",
			 "Gets the current mode of a tool." );
  commandMgr.AddCommand( *this, "SetToolBrushRadius", 2, "toolID radius",
			 "Sets the current brush radius of a tool." );
  commandMgr.AddCommand( *this, "GetToolBrushRadius", 1, "toolID",
			 "Gets the current brush radius of a tool." );
  commandMgr.AddCommand( *this, "SetToolBrushShape", 2, "toolID shape",
			 "Sets the current brush shape of a tool. shape "
			 "should be square or circle." );
  commandMgr.AddCommand( *this, "GetToolBrushShape", 1, "toolID",
			 "Gets the current brush shape of a tool, "
			 "square or circle." );
  commandMgr.AddCommand( *this, "SetToolBrush3D", 2, "toolID 3D",
			 "Sets the current brush 3D of a tool." );
  commandMgr.AddCommand( *this, "GetToolBrush3D", 1, "toolID",
			 "Gets the current brush 3D of a tool." );
  commandMgr.AddCommand( *this, "SetToolFloodStopAtLines", 2, "toolID stop",
			 "Specify whether a tool flood should stop at lines.");
  commandMgr.AddCommand( *this, "GetToolFloodStopAtLines", 1, "toolID",
			 "Returns whether or not a tool flood will stop "
			 "at lines." );
  commandMgr.AddCommand( *this, "SetToolFloodStopAtROIs", 2, "toolID stop",
			 "Specify whether a tool flood should stop at ROIs." );
  commandMgr.AddCommand( *this, "GetToolFloodStopAtROIs", 1, "toolID",
			 "Returns whether or not a tool flood will stop "
			 "at ROIs." );
  commandMgr.AddCommand( *this, "SetToolFloodFuzziness", 2, "toolID fuzziness",
			 "Specify a tool flood's fuzziness." );
  commandMgr.AddCommand( *this, "GetToolFloodFuzziness", 1, "toolID",
			 "Returns a tool flood's fuzziness." );
  commandMgr.AddCommand( *this, "SetToolFloodMaxDistance", 2,"toolID distance",
			 "Specify a tool flood's max distance." );
  commandMgr.AddCommand( *this, "GetToolFloodMaxDistance", 1, "toolID",
			 "Returns a tool flood's max distance." );
  commandMgr.AddCommand( *this, "SetToolFlood3D", 2, "toolID 3D",
			 "Sets the current brush 3D of a tool." );
  commandMgr.AddCommand( *this, "GetToolFlood3D", 1, "toolID",
			 "Gets the current brush 3D of a tool." );
  commandMgr.AddCommand( *this, "SetToolEdgeLineStraightBias", 2, 
			 "toolID bias", "Sets the bias (0-1) for straight "
			 "lines for the edge line tool." );
  commandMgr.AddCommand( *this, "GetToolEdgeLineStraightBias", 1, 
			 "toolID", "Returns the bias for straight "
			 "lines for the edge line tool." );
  commandMgr.AddCommand( *this, "SetToolEdgeLineEdgeBias", 2, 
			 "toolID bias", "Sets the bias (0-1) for edges "
			 "for the edge line tool." );
  commandMgr.AddCommand( *this, "GetToolEdgeLineEdgeBias", 1, 
			 "toolID", "Returns the bias for edges "
			 "for the edge line tool." );

}

ScubaToolState::~ScubaToolState() {

}

TclCommandListener::TclCommandResult
ScubaToolState::DoListenToTclCommand ( char* isCommand, 
				       int iArgc, char** iasArgv ) {

  // SetToolMode <toolID> <mode>
  if( 0 == strcmp( isCommand, "SetToolMode" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }
    
    if( GetID() == toolID ) {

      Mode newMode;
      if( 0 == strcmp( iasArgv[2], "navigation" )) {
	newMode = navigation;
      } else if ( 0 == strcmp( iasArgv[2], "plane" )) {
	newMode = plane;
      } else if ( 0 == strcmp( iasArgv[2], "marker" )) {
	newMode = marker;
      } else if ( 0 == strcmp( iasArgv[2], "voxelEditing" )) {
	newMode = voxelEditing;
      } else if ( 0 == strcmp( iasArgv[2], "roiEditing" )) {
	newMode = roiEditing;
      } else if ( 0 == strcmp( iasArgv[2], "straightLine" )) {
	newMode = straightLine;
      } else if ( 0 == strcmp( iasArgv[2], "edgeLine" )) {
	newMode = edgeLine;
      } else {
	sResult = "bad mode \"" + string(iasArgv[2]) + 
	  "\", should be navigation, plane, marker, voxelEditing, "
	  "roiEditing, straightLine, edgeLine.";
	return error;
      }
      SetMode( newMode );
    }
  }

  // GetToolMode <toolID>
  if( 0 == strcmp( isCommand, "GetToolMode" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }
    
    if( GetID() == toolID ) {

      switch( GetMode() ) {
      case navigation:
	sReturnValues = "navigation";
	break;
      case plane:
	sReturnValues = "plane";
	break;
      case marker:
	sReturnValues = "marker";
	break;
      case voxelEditing:
	sReturnValues = "voxelEditing";
	break;
      case roiEditing:
	sReturnValues = "roiEditing";
	break;
      case straightLine:
	sReturnValues = "straightLine";
	break;
      case edgeLine:
	sReturnValues = "edgeLine";
	break;
      }
      sReturnFormat = "s";
    }
  }

  // SetToolBrushRadius <toolID> <radius>
  if( 0 == strcmp( isCommand, "SetToolBrushRadius" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }
    
    if( GetID() == toolID ) {

      float radius = strtod(iasArgv[2], (char**)NULL);
      if( ERANGE == errno ) {
	sResult = "bad radius";
	return error;
      }
      SetBrushRadius( radius );
    }
  }

  // GetToolBrushRadius <toolID>
  if( 0 == strcmp( isCommand, "GetToolBrushRadius" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }
    
    if( GetID() == toolID ) {

      stringstream ssValues;
      ssValues << GetBrushRadius();
      sReturnValues = ssValues.str();
      sReturnFormat = "f";
    }
  }

  // SetToolBrushShape <toolID> <shape>
  if( 0 == strcmp( isCommand, "SetToolBrushShape" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }
    
    if( GetID() == toolID ) {

      Shape shape;
      if( 0 == strcmp( iasArgv[2], "square" )) {
	shape = square;
      } else if ( 0 == strcmp( iasArgv[2], "circle" )) {
	shape = circle;
      } else {
	sResult = "bad shape \"" + string(iasArgv[2]) + 
	  "\", should be square or circle.";
	return error;
      }
      SetBrushShape( shape );
    }
  }

  // GetToolBrushShape <toolID>
  if( 0 == strcmp( isCommand, "GetToolBrushShape" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }
    
    if( GetID() == toolID ) {

      switch( GetBrushShape() ) {
      case square:
	sReturnValues = "square";
	break;
      case circle:
	sReturnValues = "circle";
	break;
      }
      sReturnFormat = "s";
    }
  }

  // SetToolBrush3D <toolID> <3D>
  if( 0 == strcmp( isCommand, "SetToolBrush3D" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }
    
    if( GetID() == toolID ) {
      
      if( 0 == strcmp( iasArgv[2], "true" ) || 
	  0 == strcmp( iasArgv[2], "1" )) {
	SetBrush3D( true );
      } else if( 0 == strcmp( iasArgv[2], "false" ) ||
		 0 == strcmp( iasArgv[2], "0" ) ) {
	SetBrush3D( false );
      } else {
	sResult = "bad 3D\"" + string(iasArgv[2]) +
	  "\", should be true, 1, false, or 0";
	return error;	
      }
    }
  }

  // GetToolBrush3D <toolID>
  if( 0 == strcmp( isCommand, "GetToolBrush3D" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }
    
    if( GetID() == toolID ) {

      sReturnFormat = "i";
      stringstream ssReturnValues;
      ssReturnValues << (int)GetBrush3D();
      sReturnValues = ssReturnValues.str();
    }
  }

  // SetToolFloodStopAtLines <toolID> <stop>
  if( 0 == strcmp( isCommand, "SetToolFloodStopAtLines" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }
    
    if( GetID() == toolID ) {
      
      if( 0 == strcmp( iasArgv[2], "true" ) || 
	  0 == strcmp( iasArgv[2], "1" )) {
	SetFloodStopAtLines( true );
      } else if( 0 == strcmp( iasArgv[2], "false" ) ||
		 0 == strcmp( iasArgv[2], "0" ) ) {
	SetFloodStopAtLines( false );
      } else {
	sResult = "bad stop \"" + string(iasArgv[2]) +
	  "\", should be true, 1, false, or 0";
	return error;	
      }
    }
  }

  // GetToolFloodStopAtLines <toolID>
  if( 0 == strcmp( isCommand, "GetToolFloodStopAtLines" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }
    
    if( GetID() == toolID ) {

      sReturnFormat = "i";
      stringstream ssReturnValues;
      ssReturnValues << (int)GetFloodStopAtLines();
      sReturnValues = ssReturnValues.str();
    }
  }

  // SetToolFloodStopAtROIs <toolID> <stop>
  if( 0 == strcmp( isCommand, "SetToolFloodStopAtROIs" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }
    
    if( GetID() == toolID ) {
      
      if( 0 == strcmp( iasArgv[2], "true" ) || 
	  0 == strcmp( iasArgv[2], "1" )) {
	SetFloodStopAtROIs( true );
      } else if( 0 == strcmp( iasArgv[2], "false" ) ||
		 0 == strcmp( iasArgv[2], "0" ) ) {
	SetFloodStopAtROIs( false );
      } else {
	sResult = "bad stop \"" + string(iasArgv[2]) +
	  "\", should be true, 1, false, or 0";
	return error;	
      }
    }
  }

  // GetToolFloodStopAtROIs <toolID>
  if( 0 == strcmp( isCommand, "GetToolFloodStopAtROIs" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }
    
    if( GetID() == toolID ) {

      sReturnFormat = "i";
      stringstream ssReturnValues;
      ssReturnValues << (int)GetFloodStopAtROIs();
      sReturnValues = ssReturnValues.str();
    }
  }

  // SetToolFloodFuzziness <toolID> <fuzziness>
  if( 0 == strcmp( isCommand, "SetToolFloodFuzziness" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }
    
    if( GetID() == toolID ) {

      int fuzziness = strtol(iasArgv[2], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad fuzziness";
	return error;
      }
      SetFloodFuzziness( fuzziness );
    }
  }

  // GetToolFloodFuzziness <toolID>
  if( 0 == strcmp( isCommand, "GetToolFloodFuzziness" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }
    
    if( GetID() == toolID ) {

      sReturnFormat = "i";
      stringstream ssReturnValues;
      ssReturnValues << (int)GetFloodFuzziness();
      sReturnValues = ssReturnValues.str();
    }
  }

  // SetToolFloodMaxDistance <toolID> <MaxDistance>
  if( 0 == strcmp( isCommand, "SetToolFloodMaxDistance" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }
    
    if( GetID() == toolID ) {

      int maxDistance = strtol(iasArgv[2], (char**)NULL, 10);
      if( ERANGE == errno ) {
	sResult = "bad distance";
	return error;
      }
      SetFloodMaxDistance( maxDistance );
    }
  }

  // GetToolFloodMaxDistance <toolID>
  if( 0 == strcmp( isCommand, "GetToolFloodMaxDistance" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }
    
    if( GetID() == toolID ) {

      sReturnFormat = "i";
      stringstream ssReturnValues;
      ssReturnValues << (int)GetFloodMaxDistance();
      sReturnValues = ssReturnValues.str();
    }
  }

  // SetToolFlood3D <toolID> <3D>
  if( 0 == strcmp( isCommand, "SetToolFlood3D" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }
    
    if( GetID() == toolID ) {
      
      if( 0 == strcmp( iasArgv[2], "true" ) || 
	  0 == strcmp( iasArgv[2], "1" )) {
	SetFlood3D( true );
      } else if( 0 == strcmp( iasArgv[2], "false" ) ||
		 0 == strcmp( iasArgv[2], "0" ) ) {
	SetFlood3D( false );
      } else {
	sResult = "bad 3D\"" + string(iasArgv[2]) +
	  "\", should be true, 1, false, or 0";
	return error;	
      }
    }
  }

  // GetToolFlood3D <toolID>
  if( 0 == strcmp( isCommand, "GetToolFlood3D" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }
    
    if( GetID() == toolID ) {

      sReturnFormat = "i";
      stringstream ssReturnValues;
      ssReturnValues << (int)GetFlood3D();
      sReturnValues = ssReturnValues.str();
    }
  }

  // SetToolEdgeLineStraightBias <toolID> <bias>
  if( 0 == strcmp( isCommand, "SetToolEdgeLineStraightBias" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }
    
    if( GetID() == toolID ) {

      float bias = (float) strtod( iasArgv[2], (char**)NULL );
      if( ERANGE == errno ) {
	sResult = "bad bias";
	return error;
      }

      SetEdgeLineStraightBias( bias );
    }
  }

  // GetToolEdgeLineStraightBias <toolID>
  if( 0 == strcmp( isCommand, "GetToolEdgeLineStraightBias" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }
    
    if( GetID() == toolID ) {

      sReturnFormat = "f";
      stringstream ssReturnValues;
      ssReturnValues << GetEdgeLineStraightBias();
      sReturnValues = ssReturnValues.str();
    }
  }

  // SetToolEdgeLineEdgeBias <toolID> <bias>
  if( 0 == strcmp( isCommand, "SetToolEdgeLineEdgeBias" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }
    
    if( GetID() == toolID ) {

      float bias = (float) strtod( iasArgv[2], (char**)NULL );
      if( ERANGE == errno ) {
	sResult = "bad bias";
	return error;
      }

      SetEdgeLineEdgeBias( bias );
    }
  }

  // GetToolEdgeLineEdgeBias <toolID>
  if( 0 == strcmp( isCommand, "GetToolEdgeLineEdgeBias" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }
    
    if( GetID() == toolID ) {

      sReturnFormat = "f";
      stringstream ssReturnValues;
      ssReturnValues << GetEdgeLineEdgeBias();
      sReturnValues = ssReturnValues.str();
    }
  }

  return ok;
}
