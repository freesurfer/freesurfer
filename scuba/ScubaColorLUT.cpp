#include <iostream>
#include <fstream>
#include "ScubaColorLUT.h"


using namespace std;

DeclareIDTracker(ScubaColorLUT);

ScubaColorLUTStaticTclListener ScubaColorLUT::mStaticListener;

int const ScubaColorLUT::cDefaultEntries = 500;

ScubaColorLUT::ScubaColorLUT() {
  msLabel = "";
  mfnLUT = "";
  mcEntries = cDefaultEntries;

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

}

ScubaColorLUT::~ScubaColorLUT() {

}

void 
ScubaColorLUT::UseFile ( std::string ifnLUT ) {
  mfnLUT = ifnLUT;
  ReadFile();
}

TclCommandListener::TclCommandResult
ScubaColorLUT::DoListenToTclCommand ( char* isCommand, 
				      int iArgc, char** iasArgv ) {

  // SetColorLUTLabel <transformID> <label>
  if( 0 == strcmp( isCommand, "SetColorLUTLabel" ) ) {
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

  // GetColorLUTLabel <transformID>
  if( 0 == strcmp( isCommand, "GetColorLUTLabel" ) ) {
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

  // SetColorLUTFileName <transformID> <fileName>
  if( 0 == strcmp( isCommand, "SetColorLUTFileName" ) ) {
    int transformID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad transform ID";
      return error;
    }
    
    if( mID == transformID ) {
      
      UseFile( iasArgv[2] );
    }
  }

  // GetColorLUTFileName <transformID>
  if( 0 == strcmp( isCommand, "GetColorLUTFileName" ) ) {
    int transformID = strtol(iasArgv[1], (char**)NULL, 10);
    if( ERANGE == errno ) {
      sResult = "bad transform ID";
      return error;
    }
    
    if( mID == transformID ) {
      sReturnFormat = "s";
      stringstream ssReturnValues;
      ssReturnValues << "\"" << mfnLUT << "\"";
      sReturnValues = ssReturnValues.str();
    }
  }



  return ok;
}

void 
ScubaColorLUT::GetColorAtIndex ( int iIndex, int oColor[3] ) {

  if( iIndex >= 0 && iIndex < mcEntries ) {
    oColor[0] = mEntries[iIndex].color[0];
    oColor[1] = mEntries[iIndex].color[1];
    oColor[2] = mEntries[iIndex].color[2];
  } else  {
    oColor[0] = oColor[1] = oColor[2] = 0;
  }
}

void 
ScubaColorLUT::GetIndexOfColor ( int iColor[3], int& oIndex ) {

  
}

string 
ScubaColorLUT::GetLabelAtIndex ( int iIndex ) {

  if( iIndex >= 0 && iIndex < mcEntries ) {
    return mEntries[iIndex].msLabel;
  } else  {
    return "No entry";
  }
}

void 
ScubaColorLUT::ReadFile () {

  ifstream fLUT( mfnLUT.c_str() );
  if( fLUT.is_open() ) {

    int nLine = 0;
    int maxEntry = 0;
    while( !fLUT.eof() ) {

      string sLine;
      getline( fLUT, sLine );
      
      int nEntry;
      char sLabel[1024];
      int red, green, blue;
      int eRead = sscanf( sLine.c_str(), "%d %s %d %d %d %*s",
			  &nEntry, sLabel, &red, &green, &blue );
      if( 5 != eRead &&
	  -1 != eRead ) {
	DebugOutput(  "Error parsing " << mfnLUT << ": Malformed line " 
		      << nLine );
	continue;
      }

      mEntries[nEntry].color[0] = red;
      mEntries[nEntry].color[1] = green;
      mEntries[nEntry].color[2] = blue;
      mEntries[nEntry].msLabel  = sLabel;
      if( nEntry > maxEntry ) maxEntry = nEntry;
    
      nLine++;
    }

    mcEntries = maxEntry + 1;

  } else {

    for( int nEntry = 0; nEntry < mcEntries; nEntry+=1 ) {
      // Assign random colors.
      mEntries[nEntry].color[0] = rand() % 256;
      mEntries[nEntry].color[1] = rand() % 256;
      mEntries[nEntry].color[2] = rand() % 256;
      mEntries[nEntry].msLabel  = "";
    }
  }
}


ScubaColorLUTStaticTclListener::ScubaColorLUTStaticTclListener () {

  TclCommandManager& commandMgr = TclCommandManager::GetManager();
  commandMgr.AddCommand( *this, "GetColorLUTIDList", 0, "", 
			 "Return a list of all colorLUTIDs." );
  commandMgr.AddCommand( *this, "MakeNewColorLUT", 0, "", 
			 "Creates a new color LUT and returns its id." );
}

ScubaColorLUTStaticTclListener::~ScubaColorLUTStaticTclListener () {
}

TclCommandListener::TclCommandResult
ScubaColorLUTStaticTclListener::DoListenToTclCommand ( char* isCommand, 
					       int iArgc, char** iasArgv ) {

  // GetColorLUTIDList
  if( 0 == strcmp( isCommand, "GetColorLUTIDList" ) ) {
    list<int> idList;
    ScubaColorLUT::GetIDList( idList );
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

  // MakeNewColorLUT
  if( 0 == strcmp( isCommand, "MakeNewColorLUT" ) ) {

    ScubaColorLUT* lut = new ScubaColorLUT();
    sReturnFormat = "i";
    stringstream ssReturnValues;
    ssReturnValues << lut->GetID();
    sReturnValues = ssReturnValues.str();
  }

  return ok;
}
