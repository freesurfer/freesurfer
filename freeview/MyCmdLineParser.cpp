/**
 * @file  MyCmdLineParser.h
 * @brief A simple command-line parser class.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/10/09 17:01:54 $
 *    $Revision: 1.2 $
 *
 * Copyright (C) 2002-2009,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#include "MyCmdLineParser.h"
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>

#include <iostream>

#include <fstream>
#include <stdarg.h>

using namespace std;

#define CONSOLE_WIDTH	78

MyCmdLineParser::MyCmdLineParser( const char* ProgramName, CmdLineEntry* entries )
{
	SetProgramName( ProgramName );
	m_nNumberOfPureArguments = 0;
	if ( entries )
		SetValidCmdLineEntries( entries );
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
	for ( unsigned int i = 0; i < m_cmdLineEntries.size(); i++ )
	{
		if ( strcmp( m_cmdLineEntries[i].shortName, chFlag ) == 0 || 
		     strcmp( m_cmdLineEntries[i].longName, chFlag ) == 0 )
		{
			return true;
		}   
	}
	return false;
}

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


bool MyCmdLineParser::Found( const char* ch, string_array* sa )
{
	if ( Found( ch ) )
	{
		*sa = GetArguments( ch );
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

bool MyCmdLineParser::Found( const char* chFlag, CmdLineEntry* e )
{
	for ( unsigned int i = 0; i < m_cmdLineEntries.size(); i++ )
	{
		if ( strcmp( m_cmdLineEntries[i].shortName, chFlag ) == 0 || 
				   strcmp( m_cmdLineEntries[i].longName, chFlag ) == 0 )
		{
			*e = m_cmdLineEntries[i];
			return true;
		}   
	}
	return false;
}
	
	
bool MyCmdLineParser::IsValid( const char* chFlag, CmdLineEntry* e )
{
	for ( unsigned int i = 0; i < m_cmdLineEntriesValid.size(); i++ )
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

string_array MyCmdLineParser::GetArguments( const char* ch )
{
	CmdLineEntry e;
	if ( Found( ch, &e ) )
	{
		return e.arguments;
	}
	else
		return string_array(); 
}

void MyCmdLineParser::PrintHelp()
{
	cout << endl << "Usage: " << m_strProgramName.c_str() << " [OPTION]..." << endl;
	
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
	
	unsigned int nLen = 0;
	for ( unsigned int i = 0; i < m_cmdLineEntriesValid.size(); i++ )
	{
		CmdLineEntry e = m_cmdLineEntriesValid[i];
		string strg = string( e.shortName ) + e.longName + e.arguName;
		if ( nLen < strg.length() )
			nLen = strg.length();
	}
	nLen += 7;
	for ( unsigned int i = 0; i < m_cmdLineEntriesValid.size(); i++ )
	{
		CmdLineEntry e = m_cmdLineEntriesValid[i];
		string strg( "-" );
		strg = strg + e.shortName + ", --" + e.longName + " " + e.arguName;
		int nCnt = nLen - strg.length();
		cout << strg.c_str();
		for ( int j = 0; j < nCnt; j++ )
			cout << " ";
		desc = e.description;
		while ( desc.length() > CONSOLE_WIDTH - nLen )
		{
			int n = desc.rfind( " ", CONSOLE_WIDTH - nLen );
			if ( n >= 0 )
				cout << desc.substr( 0, n ).c_str() << endl;
			desc = desc.substr( n+1 );
			for ( unsigned int j = 0; j < nLen; j++ )
				cout << " ";
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
	for ( int i = 0; i < ( int )entries.size(); i++ )
	{
		string_array strgs = *entries[i];

		if ( !IsValid( strgs[0].c_str(), &e ) && !IsValid( strgs[0].c_str() + 1, &e ) )
		{
			bSucceed = false;
			error_msg += "Option '" + strgs[0] + "' not recognized.";
			break;
		}
		if ( e.type == CMD_LINE_OPTION )
		{
			e.arguments.clear();
			for ( int j = 1; j < (int)strgs.size(); j++ )
			{
				if ( j <= e.maxArguments )
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
			for ( int j = 1; j < (int)strgs.size(); j++ )
			{
				pureArgs.push_back( strgs[j] );
			}
		}
		m_cmdLineEntries.push_back( e );
	}
	
	// release buffers
	for ( int i = 0; i < ( int )entries.size(); i++ )
	{
		delete entries[i];
	}
	entries.clear();
	
	if ( bSucceed && (int)pureArgs.size() > m_nNumberOfPureArguments )
	{
		bSucceed = false;
		error_msg += "Option '" + pureArgs[0] + "' not recognized.";
	} 
	
	if ( !bSucceed )
	{
		PrintErrorMessage( error_msg );
	}
	else if ( Found( "h" ) )
	{
		PrintHelp();
		bSucceed = false;
	}
	
	return bSucceed;
}


