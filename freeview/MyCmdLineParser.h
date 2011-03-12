/**
 * @file  MyCmdLineParser.h
 * @brief A simple command-line parser class.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/12 00:28:52 $
 *    $Revision: 1.9 $
 *
 * Copyright (C) 2008-2009,
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

#ifndef MyCmdLineParser_h
#define MyCmdLineParser_h

#include <string>
#include <vector>
#include <QStringList>

using namespace std;

typedef vector<string> string_array;

enum CmdLineEntryType
{
  CMD_LINE_SWITCH = 0,
  CMD_LINE_OPTION,
  CMD_LINE_NONE
};

struct CmdLineEntry
{
  int    type;
  const char*  shortName;
  const char*  longName;
  const char*  arguName;
  const char*  description;
  int    minArguments;
  int    maxArguments;
  string_array arguments;

  CmdLineEntry( int nType = CMD_LINE_SWITCH,
                const char* sShortName = 0,
                const char* sLongName = 0,
                const char* sArguName = 0,
                const char* desc = 0,
                int nMinArguments = 1,
                int nMaxArguments = 1 )
  {
    if ( nMaxArguments < nMinArguments )
      nMaxArguments = nMinArguments;

    type = nType;
    shortName = sShortName;
    longName = sLongName;
    arguName = sArguName;
    description = desc;
    minArguments = nMinArguments;
    maxArguments = nMaxArguments;
  }
};

class MyCmdLineParser
{
public:
  MyCmdLineParser( const char* ProgName, CmdLineEntry* entries = NULL );
  virtual ~MyCmdLineParser();

  void SetValidCmdLineEntries( CmdLineEntry* entries );
  void AddValidCmdLineEntry( int nType, const char* sShortName,
                             const char* sLongName,
                             const char* sArguName,
                             const char* desc,
                             int nMinArguments = 1,
                             int nMaxArguments = 1 );

  void SetProgramName( string name );
  void SetProgramDescription( string text );

  bool Parse( int argc, char* argv[] );
  bool Parse( const string_array& args );
  bool Parse( const QString& cmd);

  bool Found( const char* ch );
//  bool Found( const char* ch, string* strg);
//  bool Found( const char* ch, int* value );
//  bool Found( const char* ch, double* value );
  bool Found( const char* ch, string_array* sa, int nIndex = -1  ); // -1 means last one
  bool Found( const QString flag, QStringList* sa, int nIndex = -1 );
  
  int GetNumberOfRepeats( const char* ch );

  int GetNumberOfArguments( const char* ch );

  string GetArgument( const char* ch, int n, const char* chDefault = NULL );
  string_array GetArguments( const char* ch, int nIndex = -1 );

  string_array GetFloatingArguments();
  
  void PrintHelp();

  void PrintErrorMessage( string msg );

protected:
  bool Found( const char* ch, CmdLineEntry* e, int nIndex = -1 );
  bool IsValid( const char* ch, CmdLineEntry* e );

  inline bool IsNumber(const char ch)
  {
    return ch <= '9' && ch >= '0';
  }

  vector<CmdLineEntry> m_cmdLineEntries;
  vector<CmdLineEntry> m_cmdLineEntriesValid;
  string_array         m_cmdLineFloatingArguments;
  int       m_nNumberOfPureArguments;
  string    m_strProgramName;
  string    m_strProgramDescription;
  
  bool      m_bNewLineStyle;
};

#endif
