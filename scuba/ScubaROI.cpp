/**
 * @file  ScubaROI.cpp
 * @brief Generic ROI with a an associated structure or color
 *
 * This is a superclass for ROIs. It can be of the type Structure or
 * Free. If Structure, it has a value associated with a ScubaColorLUT
 * and gets a color from that. If it's Free, you can give it an
 * explicit color.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:38 $
 *    $Revision: 1.11 $
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

#include "ScubaROI.h"
#include "ScubaColorLUT.h"
#ifdef __cplusplus
extern "C"
{
#endif
#include <string.h> // strcmp
#ifdef __cplusplus
}
#endif

using namespace std;

ScubaROIStaticTclListener ScubaROI::mStaticListener;


ScubaROI::ScubaROI () :
    Broadcaster( "ScubaROI" ) ,
    msLabel( "" ),
    mType( Free ),
    mLUTID( 0 ),
    mStructure( 0 ),
    mbSuspendROIChangedMessage( false ) {
  mFreeColor[0] = mFreeColor[1] = mFreeColor[2] = 0;

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "SetROILabel", 2, "roiID label",
                         "Set the label for a ROI." );
  commandMgr.AddCommand( *this, "GetROILabel", 1, "roiID",
                         "Returns the label for a ROI." );
  commandMgr.AddCommand( *this, "SetROIType", 2, "roiID type",
                         "Sets ROI type. type should be structure or free." );
  commandMgr.AddCommand( *this, "GetROIType", 1, "roiID",
                         "Returns ROI type of structure or free." );
  commandMgr.AddCommand( *this, "SetROILUTID", 2, "roiID lutID",
                         "Sets the LUT ID for an ROI." );
  commandMgr.AddCommand( *this, "GetROILUTID", 1, "roiID",
                         "Returns LUT ID for an ROI." );
  commandMgr.AddCommand( *this, "SetROIStructure", 2, "roiID structure",
                         "Sets ROI structure index." );
  commandMgr.AddCommand( *this, "GetROIStructure", 1, "roiID",
                         "Returns ROI structure index." );
  commandMgr.AddCommand( *this, "SetROIColor", 4, "roiID red green blue",
                         "Sets ROI color. red, green, and blue should be"
                         "0-255 integers." );
  commandMgr.AddCommand( *this, "GetROIColor", 1, "roiID",
                         "Returns ROI color as a list of red, green, and "
                         "blue integers from 0-255." );

}

ScubaROI::~ScubaROI () {}


TclCommandManager::TclCommandResult
ScubaROI::DoListenToTclCommand ( char* isCommand,
                                 int, char** iasArgv ) {

  // SetROILabel <transformID> <label>
  if ( 0 == strcmp( isCommand, "SetROILabel" ) ) {
    int transformID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad transform ID";
      return error;
    }

    if ( mID == transformID ) {

      string sLabel = iasArgv[2];
      SetLabel( sLabel );
    }
  }

  // GetROILabel <transformID>
  if ( 0 == strcmp( isCommand, "GetROILabel" ) ) {
    int transformID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad transform ID";
      return error;
    }

    if ( mID == transformID ) {
      sReturnFormat = "s";
      stringstream ssReturnValues;
      ssReturnValues << "\"" << GetLabel() << "\"";
      sReturnValues = ssReturnValues.str();
    }
  }

  // SetROIType <roiID> <type>
  if ( 0 == strcmp( isCommand, "SetROIType" ) ) {
    int roiID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad roi ID";
      return error;
    }

    if ( mID == roiID ) {

      Type type;
      if ( 0 == strcmp( iasArgv[2], "structure" ) ) {
        type = Structure;
      } else if ( 0 == strcmp( iasArgv[2], "free" ) ) {
        type = Free;
      } else {
        sResult = "bad type \"" + string(iasArgv[2]) +
                  "\", should be structure or free";
        return error;
      }

      SetType( type );
    }
  }

  // GetROIType <roiID>
  if ( 0 == strcmp( isCommand, "GetROIType" ) ) {
    int roiID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad roi ID";
      return error;
    }

    if ( mID == roiID ) {

      switch ( mType ) {
      case Structure:
        sReturnValues = "structure";
        break;
      case Free:
        sReturnValues = "free";
        break;
      }
      sReturnFormat = "s";
    }
  }

  // SetROILUTID <roiID> <lutID>
  if ( 0 == strcmp( isCommand, "SetROILUTID" ) ) {
    int roiID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad roi ID";
      return error;
    }

    if ( mID == roiID ) {

      int lutID = strtol(iasArgv[2], (char**)NULL, 10);
      if ( ERANGE == errno ) {
        sResult = "bad lut ID";
        return error;
      }

      try {
        ScubaColorLUT& lut = ScubaColorLUT::FindByID( lutID );
        if ( lut.GetID() != lutID ) {
          sResult = "Got wrong lut";
          return error;
        }
      } catch (...) {
        sResult = "bad lut ID, doesn't exist";
        return error;
      }

      SetColorLUT( lutID );
    }
  }

  // GetROILUTID <roiID>
  if ( 0 == strcmp( isCommand, "GetROILUTID" ) ) {
    int roiID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad roi ID";
      return error;
    }

    if ( mID == roiID ) {

      stringstream ssReturnValues;
      ssReturnValues << GetColorLUT();
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "i";
    }
  }

  // SetROIStructure <roiID> <structure>
  if ( 0 == strcmp( isCommand, "SetROIStructure" ) ) {
    int roiID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad roi ID";
      return error;
    }

    if ( mID == roiID ) {

      int structure = strtol( iasArgv[2], (char**)NULL, 10 );
      if ( ERANGE == errno ) {
        sResult = "bad structure";
        return error;
      }

      SetStructure( structure );
    }
  }

  // GetROIStructure <roiID>
  if ( 0 == strcmp( isCommand, "GetROIStructure" ) ) {
    int roiID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad roi ID";
      return error;
    }

    if ( mID == roiID ) {

      stringstream ssReturnValues;
      ssReturnValues << mStructure;
      sReturnValues = ssReturnValues.str();
      sReturnFormat = "i";
    }
  }

  // SetROIColor <roiID> <red> <green> <blue>
  if ( 0 == strcmp( isCommand, "SetROIColor" ) ) {
    int roiID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad roi ID";
      return error;
    }

    if ( mID == roiID ) {

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
      SetFreeColor( color );
    }
  }

  // GetROIColor <roiID>
  if ( 0 == strcmp( isCommand, "GetROIColor" ) ) {
    int roiID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad roi ID";
      return error;
    }

    if ( mID == roiID ) {

      sReturnFormat = "Liiil";
      stringstream ssReturnValues;
      ssReturnValues << mFreeColor[0] << " " << mFreeColor[1] << " "
      << mFreeColor[2];
      sReturnValues = ssReturnValues.str();
    }
  }

  return ok;
}

void 
ScubaROI::SetType( Type iType ) {
  mType = iType;
}

ScubaROI::Type 
ScubaROI::GetType() const {
  return mType;
}

void 
ScubaROI::SetColorLUT( int iLUTID ) {
  mLUTID = iLUTID;
}

int 
ScubaROI::GetColorLUT() const {
  return mLUTID;
}

void
ScubaROI::SetStructure( int iStructure ) {

  mStructure = iStructure;
  ROIChanged();
}

int
ScubaROI::GetStructure() const {

  return mStructure;
}

void
ScubaROI::SetFreeColor( int const iColor[3] ) {

  mFreeColor[0] = iColor[0];
  mFreeColor[1] = iColor[1];
  mFreeColor[2] = iColor[2];
  ROIChanged();
}

void 
ScubaROI::SetLabel( string const& isLabel ) {

  msLabel = isLabel;
}

string const&
ScubaROI::GetLabel() const{

  return msLabel;
}

void
ScubaROI::GetDrawColor( int oColor[3] ) const {

  // If this is a Free ROI, return the color, otherwise look up the
  // structure in the LUT and return the color.
  switch ( mType ) {

  case Free:
    oColor[0] = mFreeColor[0];
    oColor[1] = mFreeColor[1];
    oColor[2] = mFreeColor[2];
    break;

  case Structure: {
    try {
      ScubaColorLUT& lut = ScubaColorLUT::FindByID( mLUTID );
      lut.GetColorAtIndex( mStructure, oColor );
    } catch (...) {
      // Couldn't find LUT, just use red.
      oColor[0] = 255;
      oColor[1] = 0;
      oColor[2] = 0;
    }
  }
  break;
  }
}

void
ScubaROI::ROIChanged() {

  // Broadcast roiChanged message.
  if ( !mbSuspendROIChangedMessage ) {
    int id = GetID();
    SendBroadcast( "roiChanged", (void*)&id );
  }
}

void
ScubaROI::BeginBatchChanges () {

  mbSuspendROIChangedMessage = true;
}

void
ScubaROI::EndBatchChanges () {

  mbSuspendROIChangedMessage = false;
  ROIChanged();
}

ScubaROIStaticTclListener::ScubaROIStaticTclListener () {

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "GetROIIDList", 0, "",
                         "Return a list of all ROIs." );
}

ScubaROIStaticTclListener::~ScubaROIStaticTclListener () {}

TclCommandManager::TclCommandResult
ScubaROIStaticTclListener::DoListenToTclCommand ( char* isCommand,
						  int, char** ) {
  
  // GetROIIDList
  if ( 0 == strcmp( isCommand, "GetROIIDList" ) ) {
    list<int> idList;
    ScubaROI::GetIDList( idList );
    stringstream ssFormat;
    stringstream ssResult;
    ssFormat << "L";
    list<int>::iterator tID;
    for ( tID = idList.begin(); tID != idList.end(); ++tID ) {
      int id = *tID;
      ssFormat << "i";
      ssResult << id << " ";
    }
    ssFormat << "l";

    sReturnFormat = ssFormat.str();
    sReturnValues = ssResult.str();
  }

  return ok;
}
