/**
 * @file  ScubaColorLUT.cpp
 * @brief Reads a FreeSurfer LUT file and allows access to values
 *
 * Reads a FressSurfer LUT file (using the utils/colortab.c code) and
 * provides access to values, transforming from colors to indices and
 * getting label strings for indices.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:36 $
 *    $Revision: 1.19 $
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


#include <iostream>
#include <stdio.h>
#include <fstream>
#include "ScubaColorLUT.h"
#include "colortab.h"
extern "C" {
#include "error.h"
#include "colortab.h"
#include <string.h> // strcmp
}

using namespace std;

ScubaColorLUTStaticTclListener ScubaColorLUT::mStaticListener;

int const ScubaColorLUT::cDefaultEntries = 500;

ScubaColorLUT::ScubaColorLUT() {
  msLabel = "";
  mfnLUT = "";
  mHighestItemNo = 0;

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "SetColorLUTLabel", 2, "lutID label",
                         "Set the label for a color LUT." );
  commandMgr.AddCommand( *this, "GetColorLUTLabel", 1, "lutID",
                         "Returns the label for a color LUT." );
  commandMgr.AddCommand( *this, "SetColorLUTFileName", 2,
                         "lutID fileName",
                         "Set the LUT file name for a colorLUT." );
  commandMgr.AddCommand( *this, "GetColorLUTFileName", 1, "lutID",
                         "Returns the LUT file name for a transform." );
  commandMgr.AddCommand( *this, "GetColorLUTNumberOfEntries", 1, "lutID",
                         "Returns the number of entries in an LUT." );
  commandMgr.AddCommand( *this, "GetColorLUTEntryLabel", 2, "lutID entry",
                         "Returns the label for an entry in an LUT." );
  commandMgr.AddCommand( *this, "GetColorLUTEntryRGB", 2, "lutID entry",
                         "Returns the rgb values (0-255) for an entry "
                         "in an LUT." );
  commandMgr.AddCommand( *this, "IsColorLUTEntryValid", 2, "lutID entry",
                         "Returns whether or not an entry is valid. Can "
                         "be invalid if it's out of range or was missing "
                         "from the LUT." );
}

ScubaColorLUT::~ScubaColorLUT() {}

void
ScubaColorLUT::UseFile ( std::string ifnLUT ) {
  mfnLUT = ifnLUT;
  ReadFile();
}

TclCommandListener::TclCommandResult
ScubaColorLUT::DoListenToTclCommand ( char* isCommand,
                                      int, char** iasArgv ) {

  // SetColorLUTLabel <lutID> <label>
  if ( 0 == strcmp( isCommand, "SetColorLUTLabel" ) ) {
    int lutID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad lut ID";
      return error;
    }

    if ( mID == lutID ) {

      string sLabel = iasArgv[2];
      SetLabel( sLabel );
    }
  }

  // GetColorLUTLabel <lutID>
  if ( 0 == strcmp( isCommand, "GetColorLUTLabel" ) ) {
    int lutID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad lut ID";
      return error;
    }

    if ( mID == lutID ) {
      sReturnFormat = "s";
      stringstream ssReturnValues;
      ssReturnValues << "\"" << GetLabel() << "\"";
      sReturnValues = ssReturnValues.str();
    }
  }

  // SetColorLUTFileName <lutID> <fileName>
  if ( 0 == strcmp( isCommand, "SetColorLUTFileName" ) ) {
    int lutID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad lut ID";
      return error;
    }

    if ( mID == lutID ) {

      UseFile( iasArgv[2] );
    }
  }

  // GetColorLUTFileName <lutID>
  if ( 0 == strcmp( isCommand, "GetColorLUTFileName" ) ) {
    int lutID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad lut ID";
      return error;
    }

    if ( mID == lutID ) {
      sReturnFormat = "s";
      stringstream ssReturnValues;
      ssReturnValues << "\"" << mfnLUT << "\"";
      sReturnValues = ssReturnValues.str();
    }
  }

  // GetColorLUTNumberOfEntries <lutID>
  if ( 0 == strcmp( isCommand, "GetColorLUTNumberOfEntries" ) ) {
    int lutID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad lut ID";
      return error;
    }

    if ( mID == lutID ) {
      sReturnFormat = "i";
      stringstream ssReturnValues;
      ssReturnValues << mHighestItemNo + 1;
      sReturnValues = ssReturnValues.str();
    }
  }

  // GetColorLUTEntryLabel <lutID> <entry>
  if ( 0 == strcmp( isCommand, "GetColorLUTEntryLabel" ) ) {
    int lutID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad lut ID";
      return error;
    }

    if ( mID == lutID ) {

      int nEntry = strtol(iasArgv[2], (char**)NULL, 10);
      if ( ERANGE == errno ) {
        sResult = "bad entry";
        return error;
      }

      if ( nEntry < 0 || nEntry > mHighestItemNo ) {
        sResult = "No entry at this index";
        return error;
      }

      sReturnFormat = "s";
      stringstream ssReturnValues;
      ssReturnValues << "\"" << GetLabelAtIndex(nEntry) << "\"";
      sReturnValues = ssReturnValues.str();
    }
  }

  // GetColorLUTEntryRGB <lutID> <entry>
  if ( 0 == strcmp( isCommand, "GetColorLUTEntryRGB" ) ) {
    int lutID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad lut ID";
      return error;
    }

    if ( mID == lutID ) {

      int nEntry = strtol(iasArgv[2], (char**)NULL, 10);
      if ( ERANGE == errno ) {
        sResult = "bad entry";
        return error;
      }

      if ( nEntry < 0 || nEntry > mHighestItemNo ) {
        sResult = "No entry at this index";
        return error;
      }

      sReturnFormat = "Liiil";
      stringstream ssReturnValues;
      ssReturnValues << mEntries[nEntry].color[0] << " "
      << mEntries[nEntry].color[1] << " "
      << mEntries[nEntry].color[2];
      sReturnValues = ssReturnValues.str();
    }
  }

  // IsColorLUTEntryValid <lutID> <entry>
  if ( 0 == strcmp( isCommand, "IsColorLUTEntryValid" ) ) {
    int lutID = strtol(iasArgv[1], (char**)NULL, 10);
    if ( ERANGE == errno ) {
      sResult = "bad lut ID";
      return error;
    }

    if ( mID == lutID ) {

      int nEntry = strtol(iasArgv[2], (char**)NULL, 10);
      if ( ERANGE == errno ) {
        sResult = "bad entry index";
        return error;
      }

      sReturnFormat = "i";
      stringstream ssReturnValues;
      ssReturnValues << IsEntryValid( nEntry );
      sReturnValues = ssReturnValues.str();
    }
  }


  return ok;
}

void
ScubaColorLUT::GetColorAtIndex ( int iIndex, int oColor[3] ) {

  if ( iIndex >= 0 && iIndex <= mHighestItemNo &&
       mEntries[iIndex].mbValid ) {
    oColor[0] = mEntries[iIndex].color[0];
    oColor[1] = mEntries[iIndex].color[1];
    oColor[2] = mEntries[iIndex].color[2];
  } else  {
    oColor[0] = oColor[1] = oColor[2] = 0;
  }
}

void
ScubaColorLUT::GetIndexOfColor ( int[3], int& ) {
}

bool
ScubaColorLUT::IsEntryValid ( int inIndex ) {

  if ( inIndex >= 0 && inIndex <= mHighestItemNo ) {
    return mEntries[inIndex].mbValid;
  }
  return false;
}

string
ScubaColorLUT::GetLabelAtIndex ( int iIndex ) {

  if ( iIndex >= 0 && iIndex <= mHighestItemNo &&
       mEntries[iIndex].mbValid ) {
    return mEntries[iIndex].msLabel;
  } else  {
    return "No entry";
  }
}

void
ScubaColorLUT::ReadFile () {

  /* Read the file with the CTAB code. */
  char fnLUT[1024];
  strncpy( fnLUT, mfnLUT.c_str(), sizeof(fnLUT) );
  COLOR_TABLE* ctab = CTABreadASCII( fnLUT );
  if ( NULL == ctab ) {
    throw new runtime_error( string("Couldn't open color table file ") +
                             mfnLUT );
  }

  /* Go through and make entries for each valid entry we got. */
  int cEntries;
  CTABgetNumberOfTotalEntries( ctab, &cEntries );
  mHighestItemNo = -1;
  for ( int nEntry = 0; nEntry < cEntries; nEntry++ ) {

    int bValid;
    CTABisEntryValid( ctab, nEntry, &bValid );
    if ( bValid ) {

      mEntries[nEntry].mbValid = true;
      CTABrgbaAtIndexi( ctab, nEntry,
                        &mEntries[nEntry].color[0],
                        &mEntries[nEntry].color[1],
                        &mEntries[nEntry].color[2],
                        &mEntries[nEntry].alpha );
      char sName[1024];
      CTABcopyName( ctab, nEntry, sName, sizeof(sName) );
      mEntries[nEntry].msLabel = sName;

      /* Update our highest item number. */
      if ( nEntry > mHighestItemNo )
        mHighestItemNo = nEntry;

    } else {

      mEntries[nEntry].mbValid = false;
    }
  }

  CTABfree( &ctab );
}


ScubaColorLUTStaticTclListener::ScubaColorLUTStaticTclListener () {

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "GetColorLUTIDList", 0, "",
                         "Return a list of all colorLUTIDs." );
  commandMgr.AddCommand( *this, "MakeNewColorLUT", 0, "",
                         "Creates a new color LUT and returns its id." );
}

ScubaColorLUTStaticTclListener::~ScubaColorLUTStaticTclListener () {}

TclCommandListener::TclCommandResult
ScubaColorLUTStaticTclListener::DoListenToTclCommand ( char* isCommand,
    int, char** ) {

  // GetColorLUTIDList
  if ( 0 == strcmp( isCommand, "GetColorLUTIDList" ) ) {
    list<int> idList;
    ScubaColorLUT::GetIDList( idList );
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

  // MakeNewColorLUT
  if ( 0 == strcmp( isCommand, "MakeNewColorLUT" ) ) {

    ScubaColorLUT* lut = new ScubaColorLUT();
    sReturnFormat = "i";
    stringstream ssReturnValues;
    ssReturnValues << lut->GetID();
    sReturnValues = ssReturnValues.str();
  }

  return ok;
}
