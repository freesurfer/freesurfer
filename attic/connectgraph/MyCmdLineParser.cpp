/**
 * @file  MyCmdLineParser.h
 * @brief A simple command-line parser class.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:01 $
 *    $Revision: 1.2 $
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


#include "MyCmdLineParser.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>

#include <iostream>

#include <fstream>
#include <stdarg.h>

using namespace std;

#define CONSOLE_WIDTH 78

MyCmdLineParser::MyCmdLineParser( const char* ProgramName, CmdLineEntry* entries )
{
  SetProgramName( ProgramName );
  m_nNumberOfPureArguments = 0;
  if ( entries )
    SetValidCmdLineEntries( entries );
  
  m_bNewLineStyle = true;
}


MyCmdLineParser::~MyCmdLineParser()
{
  m_cmdLineEntries.clear();
  m_cmdLineEntriesValid.clear();
}

void MyCmdLineParser::SetValidCmdLineEntries( CmdLineEntry* entries )
{
  m_cmdLineEntriesValid.clear();
  int n = 0;
  bool bHasHelp = false;
  while ( entries[n].type != CMD_LINE_NONE )
  {
    if ( strcmp( entries[n].shortName, "h" ) == 0 )
      bHasHelp = true;
    m_cmdLineEntriesValid.push_back( entries[n] );
    n++;
  }
  if ( !bHasHelp )
  {
    CmdLineEntry e( CMD_LINE_SWITCH, "h", "help", "", "Display this help." );
    m_cmdLineEntriesValid.push_back( e );
  }
}

void MyCmdLineParser::AddValidCmdLineEntry( int nType,
    const char* shortname,
    const char* longname,
    const char* arguname,
    const char* desc,
    int min_num,
    int max_num )
{
  CmdLineEntry e( nType, shortname, longname, arguname, desc, min_num, max_num );
  m_cmdLineEntriesValid.push_back( e );
}

void MyCmdLineParser::SetProgramName( string name )
{
  m_strProgramName = name;
}

void MyCmdLineParser::SetProgramDescription( string text )
{
  m_strProgramDescription = text;
}

bool MyCmdLineParser::Found( const char* chFlag )
{
  for ( size_t i = 0; i < m_cmdLineEntries.size(); i++ )
  {
    if ( strcmp( m_cmdLineEntries[i].shortName, chFlag ) == 0 ||
         strcmp( m_cmdLineEntries[i].longName, chFlag ) == 0 )
    {
      return true;
    }
  }
  return false;
}

/*
bool MyCmdLineParser::Found( const char* ch, string* strg )
{
  if ( Found( ch ) )
  {
    *strg = GetArgument( ch, 0 );
    return true;
  }
  else
    return false;
}


bool MyCmdLineParser::Found( const char* ch, int* value )
{
  if ( Found( ch ) )
{
    string strg = GetArgument( ch, 0 );
    *value = atoi( strg.c_str() );
    return true;
}
  else
    return false;
}

bool MyCmdLineParser::Found( const char* ch, double* value )
{
  if ( Found( ch ) )
{
    string strg = GetArgument( ch, 0 );
    *value = atof( strg.c_str() );
    return true;
}
  else
    return false;
}
*/

bool MyCmdLineParser::Found( const char* ch, string_array* sa, int nIndex )
{
  if ( Found( ch ) )
  {
    *sa = GetArguments( ch, nIndex );
    return true;
  }
  else
    return false;
}

bool MyCmdLineParser::Found( const char* chFlag, CmdLineEntry* e, int nIndex )
{
  int n = 0;
  for ( size_t i = 0; i < m_cmdLineEntries.size(); i++ )
  {
    if ( strcmp( m_cmdLineEntries[i].shortName, chFlag ) == 0 ||
         strcmp( m_cmdLineEntries[i].longName, chFlag ) == 0 )
    {
      *e = m_cmdLineEntries[i];
      if ( nIndex >= 0 && n == nIndex )
        return true;
          
      n++;
    }
  }
  if ( n > 0 )
    return true;
  else
    return false;
}

int MyCmdLineParser::GetNumberOfRepeats( const char* chFlag )
{
  int n = 0; 
  for ( size_t i = 0; i < m_cmdLineEntries.size(); i++ )
  {
    if ( strcmp( m_cmdLineEntries[i].shortName, chFlag ) == 0 ||
         strcmp( m_cmdLineEntries[i].longName, chFlag ) == 0 )
    {
      n++;
    }
  }
  return n;
}

bool MyCmdLineParser::IsValid( const char* chFlag, CmdLineEntry* e )
{
  for ( size_t i = 0; i < m_cmdLineEntriesValid.size(); i++ )
  {
    if ( strcmp( m_cmdLineEntriesValid[i].shortName, chFlag ) == 0 ||
         strcmp( m_cmdLineEntriesValid[i].longName, chFlag ) == 0 )
    {
      *e = m_cmdLineEntriesValid[i];
      return true;
    }
  }
  return false;
}

int MyCmdLineParser::GetNumberOfArguments( const char* ch )
{
  CmdLineEntry e;
  if ( Found( ch, &e ) )
    return e.arguments.size();

  return 0;
}

string MyCmdLineParser::GetArgument( const char* ch, int n, const char* chDefault )
{
  CmdLineEntry e;
  if ( Found( ch, &e ) )
  {
    if ( ( int )e.arguments.size() > n )
      return e.arguments[n];
  }
  if ( chDefault )
    return chDefault;
  else
    return "";
}

string_array MyCmdLineParser::GetArguments( const char* ch, int nIndex )
{
  CmdLineEntry e;
  if ( Found( ch, &e, nIndex ) )
  {
    return e.arguments;
  }
  else
    return string_array();
}

void MyCmdLineParser::PrintHelp()
{
  cout << endl << "Usage: " << m_strProgramName.c_str() << " [OPTION <ARGUMENT:SUB-OPTION>]..." << endl;

  string desc = m_strProgramDescription;
  while ( desc.length() > CONSOLE_WIDTH )
  {
    int n = desc.rfind( " ", CONSOLE_WIDTH);
    if ( n >= 0 )
      cout << desc.substr( 0, n ).c_str() << endl;
    desc = desc.substr( n+1 );
  }
  if ( desc.length() > 0 )
    cout << desc.c_str() << endl;
  cout << endl;

  size_t nLen = 0;
  for ( size_t i = 0; i < m_cmdLineEntriesValid.size(); i++ )
  {
    CmdLineEntry e = m_cmdLineEntriesValid[i];
    string strg = string( e.shortName ) + e.longName + e.arguName;
    if ( nLen < strg.length() )
      nLen = strg.length();
  }
  nLen += 7;
  if ( m_bNewLineStyle )
    nLen = 7;
  for ( size_t i = 0; i < m_cmdLineEntriesValid.size(); i++ )
  {
    CmdLineEntry e = m_cmdLineEntriesValid[i];
    string strg( "-" );
    strg = strg + e.shortName + ", --" + e.longName + " " + e.arguName;
    cout << strg.c_str();
    if ( m_bNewLineStyle )
    {
      cout << endl;
      for ( size_t j = 0; j < nLen; j++ )
        cout << " ";
    }
    else
    {
      int nCnt = nLen - strg.length();
      for ( int j = 0; j < nCnt; j++ )
        cout << " ";
    }
    desc = e.description;
    while ( desc.length() > CONSOLE_WIDTH - nLen || desc.find( "\n" ) != string::npos )
    {
      size_t n = desc.rfind( " ", CONSOLE_WIDTH - nLen );
      size_t m = desc.substr( 0, CONSOLE_WIDTH - nLen ).find( "\n" );
      if ( m != string::npos )
        n = m;
      if ( n != string::npos )
      {
        cout << desc.substr( 0, n ).c_str() << endl;
        desc = desc.substr( n+1 );
      }
      if ( desc.size() > 0 )
      {
        for ( size_t j = 0; j < nLen; j++ )
          cout << " ";
      }
      else
        cout << endl;
    }
    if ( desc.length() > 0 )
      cout << desc.c_str() << endl;
  }
  cout << endl;
}

void MyCmdLineParser::PrintErrorMessage( string msg )
{
  PrintHelp();

  cout << msg.c_str() << endl << endl;
}

bool MyCmdLineParser::Parse( int argc, char* argv[] )
{
  // first parse the input command line into entries, don't care if they are valid or not
  vector<string_array*> entries;
  string_array pureArgs;
  string_array* sa = NULL;

  for ( int i = 1; i < argc; i++ )
  {
    if ( argv[i][0] == '-' && string( argv[i] ).length() > 1
         && !IsNumber( argv[i][1] ) && argv[i][1] != '.' )
    {
      sa = new string_array;
      sa->clear();
      if ( string( argv[i] ).length() > 2 && argv[i][1] == '-' )    // long name
        sa->push_back( argv[i]+2 );
      else
        sa->push_back( argv[i]+1 );
      entries.push_back( sa );
    }
    else if ( sa )
    {
      sa->push_back( argv[i] );
    }
    else
      pureArgs.push_back( argv[i] );
  }

  //
  m_cmdLineEntries.clear();
  CmdLineEntry e;
  bool bSucceed = true;
  string error_msg = "";
  for ( size_t i = 0; i < entries.size(); i++ )
  {
    string_array strgs = *entries[i];

    if ( !IsValid( strgs[0].c_str(), &e ) ) // && !IsValid( strgs[0].c_str() + 1, &e ) )
    {
      bSucceed = false;
      error_msg += "Option '" + strgs[0] + "' not recognized.";
      break;
    }
    if ( e.type == CMD_LINE_OPTION )
    {
      e.arguments.clear();
      for ( size_t j = 1; j < strgs.size(); j++ )
      {
        if ( j <= (size_t)e.maxArguments )
          e.arguments.push_back( strgs[j] );
        else
          pureArgs.push_back( strgs[j] );
      }
      if ( (int)e.arguments.size() < e.minArguments )
      {
        bSucceed = false;
        cout << e.arguments.size() << " " << e.minArguments << endl;
        error_msg += "Argument missing for option '" + strgs[0] + "'.";
      }
    }
    else if ( e.type == CMD_LINE_SWITCH )
    {
      for ( size_t j = 1; j < strgs.size(); j++ )
      {
        pureArgs.push_back( strgs[j] );
      }
    }
    m_cmdLineEntries.push_back( e );
  }

  // release buffers
  for ( size_t i = 0; i < entries.size(); i++ )
  {
    delete entries[i];
  }
  entries.clear();

  /*
  if ( bSucceed && (int)pureArgs.size() > m_nNumberOfPureArguments )
  {
    bSucceed = false;
    error_msg += "Option '" + pureArgs[0] + "' not recognized.";
  }
  */

  m_cmdLineFloatingArguments = pureArgs;
  
  if ( !bSucceed )
  {
    PrintErrorMessage( error_msg );
  }
  else if ( Found( "h" ) || Found( "help" ) )
  {
    PrintHelp();
    bSucceed = false;
  }

  return bSucceed;
}

string_array MyCmdLineParser::GetFloatingArguments()
{
  return m_cmdLineFloatingArguments;
}

