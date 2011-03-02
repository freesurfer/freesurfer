/**
 * @file  ScubaToolState.cpp
 * @brief State variables for the Window-global tool
 *
 * This class simply represents the Window's tool and its current
 * settings. It is available all the way to the Layer subclasses.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:38 $
 *    $Revision: 1.31 $
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


#include <errno.h>
#include "string_fixed.h"
#include "ScubaToolState.h"
#ifdef __cplusplus
extern "C"
{
#endif
#include <string.h> // strcmp
#ifdef __cplusplus
}
#endif

using namespace std;

ScubaToolState::ScubaToolState() :
    mMode( navigation ),
    mTargetLayer( -1 ),
    mBrushRadius( 0.5 ),
    mBrushShape( voxel ),
    mbBrush3D( false ),
    mNewValue( 1 ),
    mNewValueMinThreshold( 0 ), // These are for COR volume editing
    mNewValueMaxThreshold( 5 ),
    mEraseValue( 0 ),
    mEraseValueMinThreshold( 5 ),
    mEraseValueMaxThreshold( 255 ),
    mbUseEditThreshold( false ),
    mbOnlyFillZero( false ),
    mbFloodStopAtPaths( true ),
    mbFloodStopAtROIs( true ),
    mFloodFuzziness( 0 ),
    mFloodFuzzinessType( seed ),
    mFloodMaxDistance( 0 ),
    mbFlood3D( false ),
    mFloodSourceCollection( -1 ),
    mbOnlyFloodZero( false ),
    mEdgePathEdgeBias( 0.9 ) {
  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "SetToolMode", 2, "toolID mode",
                         "Sets the current mode of a tool." );
  commandMgr.AddCommand( *this, "GetToolMode", 1, "toolID",
                         "Gets the current mode of a tool." );
  commandMgr.AddCommand( *this, "SetToolTargetLayer", 2, "toolID layerID",
                         "Sets the target layer of a tool." );
  commandMgr.AddCommand( *this, "GetToolTargetLayer", 1, "toolID",
                         "Gets the target layer of a tool." );
  commandMgr.AddCommand( *this, "SetToolNewVoxelValue", 2, "toolID value",
                         "Sets the new voxel value of a tool." );
  commandMgr.AddCommand( *this, "GetToolNewVoxelValue", 1, "toolID",
                         "Gets the new voxel value of a tool." );
  commandMgr.AddCommand( *this, "SetToolEraseVoxelValue", 2, "toolID value",
                         "Sets the erase voxel value of a tool." );
  commandMgr.AddCommand( *this, "GetToolEraseVoxelValue", 1, "toolID",
                         "Gets the erase voxel value of a tool." );
  commandMgr.AddCommand( *this, "SetToolOnlyBrushZero", 2, "toolID onlyZero",
                         "Specify whether the brush should only affect zero "
                         "values.." );
  commandMgr.AddCommand( *this, "GetToolOnlyBrushZero", 1, "toolID",
                         "Returns whether or not a brush is only affecting "
                         "zero values." );
  commandMgr.AddCommand( *this, "SetToolBrushRadius", 2, "toolID radius",
                         "Sets the current brush radius of a tool." );
  commandMgr.AddCommand( *this, "GetToolBrushRadius", 1, "toolID",
                         "Gets the current brush radius of a tool." );
  commandMgr.AddCommand( *this, "SetToolBrushShape", 2, "toolID shape",
                         "Sets the current brush shape of a tool. shape "
                         "should be voxel, square or circle." );
  commandMgr.AddCommand( *this, "GetToolBrushShape", 1, "toolID",
                         "Gets the current brush shape of a tool, "
                         "voxel, square or circle." );
  commandMgr.AddCommand( *this, "SetToolBrush3D", 2, "toolID 3D",
                         "Sets the current brush 3D of a tool." );
  commandMgr.AddCommand( *this, "GetToolBrush3D", 1, "toolID",
                         "Gets the current brush 3D of a tool." );
  commandMgr.AddCommand( *this, "SetToolFloodStopAtPaths", 2, "toolID stop",
                         "Specify whether a tool flood should stop at paths.");
  commandMgr.AddCommand( *this, "GetToolFloodStopAtPaths", 1, "toolID",
                         "Returns whether or not a tool flood will stop "
                         "at paths." );
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
                         "Sets the current flood 3D of a tool." );
  commandMgr.AddCommand( *this, "GetToolFlood3D", 1, "toolID",
                         "Gets the current flood 3D of a tool." );
  commandMgr.AddCommand( *this, "SetToolFloodSourceCollection", 2,
                         "toolID colID", "Sets the current flood "
                         "source collection of a tool." );
  commandMgr.AddCommand( *this, "SetToolFloodFuzzinessType", 2, "toolID type",
                         "Sets the tool's fuzziness type. Should be seed "
                         "or gradient." );
  commandMgr.AddCommand( *this, "GetToolFloodFuzzinessType", 1, "toolID",
                         "Returns the tool's fuzziness type: seed "
                         "or gradient." );
  commandMgr.AddCommand( *this, "GetToolFloodSourceCollection", 1,
                         "toolID", "Gets the current flood "
                         "source collection of a tool." );
  commandMgr.AddCommand( *this, "SetToolOnlyFloodZero", 2, "toolID onlyZero",
                         "Specify whether the flood should only affect zero "
                         "values.." );
  commandMgr.AddCommand( *this, "GetToolOnlyFloodZero", 1, "toolID",
                         "Returns whether or not a flood is only affecting "
                         "zero values." );
  commandMgr.AddCommand( *this, "SetToolEdgePathEdgeBias", 2,
                         "toolID bias", "Sets the bias (0-1) for edges "
                         "for the edge path tool." );
  commandMgr.AddCommand( *this, "GetToolEdgePathEdgeBias", 1,
                         "toolID", "Returns the bias for edges "
                         "for the edge path tool." );

}

ScubaToolState::~ScubaToolState() {}

TclCommandListener::TclCommandResult
ScubaToolState::DoListenToTclCommand ( char* isCommand,
                                       int, char** iasArgv ) {

  // SetToolMode <toolID> <mode>
  if ( 0 == strcmp( isCommand, "SetToolMode" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      Mode newMode;
      if ( 0 == strcmp( iasArgv[2], "navigation" )) {
        newMode = navigation;
      } else if ( 0 == strcmp( iasArgv[2], "plane" )) {
        newMode = plane;
      } else if ( 0 == strcmp( iasArgv[2], "marker" )) {
        newMode = marker;
      } else if ( 0 == strcmp( iasArgv[2], "voxelEditing" )) {
        newMode = voxelEditing;
      } else if ( 0 == strcmp( iasArgv[2], "voxelFilling" )) {
        newMode = voxelFilling;
      } else if ( 0 == strcmp( iasArgv[2], "roiEditing" )) {
        newMode = roiEditing;
      } else if ( 0 == strcmp( iasArgv[2], "roiFilling" )) {
        newMode = roiFilling;
      } else if ( 0 == strcmp( iasArgv[2], "straightPath" )) {
        newMode = straightPath;
      } else if ( 0 == strcmp( iasArgv[2], "edgePath" )) {
        newMode = edgePath;
      } else {
        sResult = "bad mode \"" + string(iasArgv[2]) +
                  "\", should be navigation, plane, marker, voxelEditing, "
                  "voxelFilling, roiEditing, roiFilling, straightPath, edgePath.";
        return error;
      }
      SetMode( newMode );
    }
  }

  // GetToolMode <toolID>
  if ( 0 == strcmp( isCommand, "GetToolMode" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      switch ( GetMode() ) {
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
      case voxelFilling:
        sReturnValues = "voxelFilling";
        break;
      case roiEditing:
        sReturnValues = "roiEditing";
        break;
      case roiFilling:
        sReturnValues = "roiFilling";
        break;
      case straightPath:
        sReturnValues = "straightPath";
        break;
      case edgePath:
        sReturnValues = "edgePath";
        break;
      }
      sReturnFormat = "s";
    }
  }

  // SetToolNewVoxelValue <toolID> <value>
  if ( 0 == strcmp( isCommand, "SetToolNewVoxelValue" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      float value = strtod(iasArgv[2], (char**)NULL);
      if ( ERANGE == errno ) {
        sResult = "bad radius";
        return error;
      }
      SetNewValue( value );
    }
  }

  // GetToolNewVoxelValue <toolID>
  if ( 0 == strcmp( isCommand, "GetToolNewVoxelValue" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      stringstream ssValues;
      ssValues << GetNewValue();
      sReturnValues = ssValues.str();
      sReturnFormat = "f";
    }
  }

  // SetToolEraseVoxelValue <toolID> <value>
  if ( 0 == strcmp( isCommand, "SetToolEraseVoxelValue" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      float value = strtod(iasArgv[2], (char**)NULL);
      if ( ERANGE == errno ) {
        sResult = "bad radius";
        return error;
      }
      SetEraseValue( value );
    }
  }

  // GetToolEraseVoxelValue <toolID>
  if ( 0 == strcmp( isCommand, "GetToolEraseVoxelValue" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      stringstream ssValues;
      ssValues << GetEraseValue();
      sReturnValues = ssValues.str();
      sReturnFormat = "f";
    }
  }

  // SetToolOnlyBrushZero <toolID> <onlyZero>
  if ( 0 == strcmp( isCommand, "SetToolOnlyBrushZero" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {
      try {
        bool bOnlyZero =
          TclCommandManager::ConvertArgumentToBoolean( iasArgv[2] );
        SetOnlyBrushZero( bOnlyZero );
      } catch ( runtime_error& e ) {
        sResult = "bad onlyZero \"" + string(iasArgv[2]) + "\"," + e.what();
        return error;
      }
    }
  }

  // GetToolOnlyBrushZero <toolID>
  if ( 0 == strcmp( isCommand, "GetToolOnlyBrushZero" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {
      sReturnValues =
        TclCommandManager::ConvertBooleanToReturnValue( GetOnlyBrushZero() );
      sReturnFormat = "i";
    }
  }

  // SetToolFloodFuzzinessType <toolID> <type>
  if ( 0 == strcmp( isCommand, "SetToolFloodFuzzinessType" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      FuzzinessType type;
      if ( 0 == strcmp( iasArgv[2], "seed" )) {
        type = seed;
      } else if ( 0 == strcmp( iasArgv[2], "gradient" )) {
        type = gradient;
      } else {
        sResult = "bad type \"" + string(iasArgv[2]) +
                  "\", should be seed or gradient.";
        return error;
      }
      SetFuzzinessType( type );
    }
  }

  // GetToolFloodFuzzinessType <toolID>
  if ( 0 == strcmp( isCommand, "GetToolFloodFuzzinessType" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      switch ( GetFuzzinessType() ) {
      case seed:
        sReturnValues = "seed";
        break;
      case gradient:
        sReturnValues = "gradient";
        break;
      }
      sReturnFormat = "s";
    }
  }

  // SetToolTargetLayer <toolID> <layerID>
  if ( 0 == strcmp( isCommand, "SetToolTargetLayer" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      int layerID = strtol(iasArgv[2], (char**)NULL, 10);
      if ( ERANGE == errno ) {
        sResult = "bad layerID";
        return error;
      }
      SetTargetLayer( layerID );
    }
  }

  // GetToolTargetLayer <toolID>
  if ( 0 == strcmp( isCommand, "GetToolTargetLayer" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      stringstream ssValues;
      ssValues << GetTargetLayer();
      sReturnValues = ssValues.str();
      sReturnFormat = "i";
    }
  }

  // SetToolBrushRadius <toolID> <radius>
  if ( 0 == strcmp( isCommand, "SetToolBrushRadius" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      float radius = strtod(iasArgv[2], (char**)NULL);
      if ( ERANGE == errno ) {
        sResult = "bad radius";
        return error;
      }
      SetBrushRadius( radius );
    }
  }

  // GetToolBrushRadius <toolID>
  if ( 0 == strcmp( isCommand, "GetToolBrushRadius" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      stringstream ssValues;
      ssValues << GetBrushRadius();
      sReturnValues = ssValues.str();
      sReturnFormat = "f";
    }
  }

  // SetToolBrushShape <toolID> <shape>
  if ( 0 == strcmp( isCommand, "SetToolBrushShape" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      Shape shape;
      if ( 0 == strcmp( iasArgv[2], "voxel" )) {
        shape = voxel;
      } else if ( 0 == strcmp( iasArgv[2], "square" )) {
        shape = square;
      } else if ( 0 == strcmp( iasArgv[2], "circle" )) {
        shape = circle;
      } else {
        sResult = "bad shape \"" + string(iasArgv[2]) +
                  "\", should be voxel, square or circle.";
        return error;
      }
      SetBrushShape( shape );
    }
  }

  // GetToolBrushShape <toolID>
  if ( 0 == strcmp( isCommand, "GetToolBrushShape" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      switch ( GetBrushShape() ) {
      case voxel:
        sReturnValues = "voxel";
        break;
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
  if ( 0 == strcmp( isCommand, "SetToolBrush3D" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      if ( 0 == strcmp( iasArgv[2], "true" ) ||
           0 == strcmp( iasArgv[2], "1" )) {
        SetBrush3D( true );
      } else if ( 0 == strcmp( iasArgv[2], "false" ) ||
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
  if ( 0 == strcmp( isCommand, "GetToolBrush3D" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      sReturnFormat = "i";
      stringstream ssReturnValues;
      ssReturnValues << (int)GetBrush3D();
      sReturnValues = ssReturnValues.str();
    }
  }

  // SetToolFloodStopAtPaths <toolID> <stop>
  if ( 0 == strcmp( isCommand, "SetToolFloodStopAtPaths" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      if ( 0 == strcmp( iasArgv[2], "true" ) ||
           0 == strcmp( iasArgv[2], "1" )) {
        SetFloodStopAtPaths( true );
      } else if ( 0 == strcmp( iasArgv[2], "false" ) ||
                  0 == strcmp( iasArgv[2], "0" ) ) {
        SetFloodStopAtPaths( false );
      } else {
        sResult = "bad stop \"" + string(iasArgv[2]) +
                  "\", should be true, 1, false, or 0";
        return error;
      }
    }
  }

  // GetToolFloodStopAtPaths <toolID>
  if ( 0 == strcmp( isCommand, "GetToolFloodStopAtPaths" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      sReturnFormat = "i";
      stringstream ssReturnValues;
      ssReturnValues << (int)GetFloodStopAtPaths();
      sReturnValues = ssReturnValues.str();
    }
  }

  // SetToolFloodStopAtROIs <toolID> <stop>
  if ( 0 == strcmp( isCommand, "SetToolFloodStopAtROIs" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      if ( 0 == strcmp( iasArgv[2], "true" ) ||
           0 == strcmp( iasArgv[2], "1" )) {
        SetFloodStopAtROIs( true );
      } else if ( 0 == strcmp( iasArgv[2], "false" ) ||
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
  if ( 0 == strcmp( isCommand, "GetToolFloodStopAtROIs" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      sReturnFormat = "i";
      stringstream ssReturnValues;
      ssReturnValues << (int)GetFloodStopAtROIs();
      sReturnValues = ssReturnValues.str();
    }
  }

  // SetToolFloodFuzziness <toolID> <fuzziness>
  if ( 0 == strcmp( isCommand, "SetToolFloodFuzziness" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      float fuzziness = strtod(iasArgv[2], (char**)NULL);
      if ( ERANGE == errno ) {
        sResult = "bad fuzziness";
        return error;
      }
      SetFloodFuzziness( fuzziness );
    }
  }

  // GetToolFloodFuzziness <toolID>
  if ( 0 == strcmp( isCommand, "GetToolFloodFuzziness" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      sReturnFormat = "i";
      stringstream ssReturnValues;
      ssReturnValues << (int)GetFloodFuzziness();
      sReturnValues = ssReturnValues.str();
    }
  }

  // SetToolFloodMaxDistance <toolID> <MaxDistance>
  if ( 0 == strcmp( isCommand, "SetToolFloodMaxDistance" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      float maxDistance = strtod(iasArgv[2], (char**)NULL);
      if ( ERANGE == errno ) {
        sResult = "bad distance";
        return error;
      }
      SetFloodMaxDistance( maxDistance );
    }
  }

  // GetToolFloodMaxDistance <toolID>
  if ( 0 == strcmp( isCommand, "GetToolFloodMaxDistance" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      sReturnFormat = "i";
      stringstream ssReturnValues;
      ssReturnValues << (int)GetFloodMaxDistance();
      sReturnValues = ssReturnValues.str();
    }
  }

  // SetToolFlood3D <toolID> <3D>
  if ( 0 == strcmp( isCommand, "SetToolFlood3D" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      if ( 0 == strcmp( iasArgv[2], "true" ) ||
           0 == strcmp( iasArgv[2], "1" )) {
        SetFlood3D( true );
      } else if ( 0 == strcmp( iasArgv[2], "false" ) ||
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
  if ( 0 == strcmp( isCommand, "GetToolFlood3D" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      sReturnFormat = "i";
      stringstream ssReturnValues;
      ssReturnValues << (int)GetFlood3D();
      sReturnValues = ssReturnValues.str();
    }
  }

  // SetToolFloodSourceCollection <toolID> <colID>
  if ( 0 == strcmp( isCommand, "SetToolFloodSourceCollection" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      int colID = strtol(iasArgv[2], (char**)NULL, 10);
      if ( ERANGE == errno ) {
        sResult = "bad colID";
        return error;
      }
      SetFloodSourceCollection( colID );
    }
  }

  // GetToolFloodSourceCollection <toolID>
  if ( 0 == strcmp( isCommand, "GetToolFloodSourceCollection" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      sReturnFormat = "i";
      stringstream ssReturnValues;
      ssReturnValues << (int)GetFloodSourceCollection();
      sReturnValues = ssReturnValues.str();
    }
  }

  // SetToolOnlyFloodZero <toolID> <onlyZero>
  if ( 0 == strcmp( isCommand, "SetToolOnlyFloodZero" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {
      try {
        bool bOnlyZero =
          TclCommandManager::ConvertArgumentToBoolean( iasArgv[2] );
        SetOnlyFloodZero( bOnlyZero );
      } catch ( runtime_error& e ) {
        sResult = "bad onlyZero \"" + string(iasArgv[2]) + "\"," + e.what();
        return error;
      }
    }
  }

  // GetToolOnlyFloodZero <toolID>
  if ( 0 == strcmp( isCommand, "GetToolOnlyFloodZero" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {
      sReturnValues =
        TclCommandManager::ConvertBooleanToReturnValue( GetOnlyFloodZero() );
      sReturnFormat = "i";
    }
  }

  // SetToolEdgePathEdgeBias <toolID> <bias>
  if ( 0 == strcmp( isCommand, "SetToolEdgePathEdgeBias" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      float bias = (float) strtod( iasArgv[2], (char**)NULL );
      if ( ERANGE == errno ) {
        sResult = "bad bias";
        return error;
      }

      SetEdgePathEdgeBias( bias );
    }
  }

  // GetToolEdgePathEdgeBias <toolID>
  if ( 0 == strcmp( isCommand, "GetToolEdgePathEdgeBias" ) ) {
    int toolID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad tool ID";
      return error;
    }

    if ( GetID() == toolID ) {

      sReturnFormat = "f";
      stringstream ssReturnValues;
      ssReturnValues << GetEdgePathEdgeBias();
      sReturnValues = ssReturnValues.str();
    }
  }

  return ok;
}

string
ScubaToolState::GetBrushShapeAsString () {

  switch ( GetBrushShape() ) {
  case voxel:
    return "voxel";
    break;
  case square:
    return "square";
    break;
  case circle:
    return "circle";
    break;
  }

  return "Unkown";
}
