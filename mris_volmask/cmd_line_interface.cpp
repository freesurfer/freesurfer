/**
 * @brief automate the parsing of the command-line
 *
 * this class is meant to automate the parsing of the command-line
 * its main purpose is to minimize new code being written when a new
 * program is written
 *
 *  as such, the data it holds have 2 purposes:
 *    1. will process the options, as well as the input and
 *       output files the user enters
 *    2. will help automate the print_help command, which otherwise
 *       would not be centralized
 */
/*
 * Original Author: Gheorge Postelni
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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
#include <vector>
#include <string>
#include <string.h>

#include "cmd_line_interface.h"
#include "utils.h"

using namespace std;

void
CCmdLineInterface:: AddOptionInt(const char *i_cstrName,
                                 int *piVal,
                                 const char *help)
{
  CCmdLineOptionInt cmd( i_cstrName, 1, help);


  cmd.Add(piVal);

  // register that with the generic argument list
  m_lstIntOption.push_back(cmd);
}

void
CCmdLineInterface::AddOptionInt(const char *i_cstrName,
                                int *piVal1,
                                int *piVal2,
                                const char *help)
{
  CCmdLineOptionInt cmd(i_cstrName, 2, help);

  cmd.Add(piVal1);
  cmd.Add(piVal2);

  // register that with the generic argument list
  m_lstIntOption.push_back(cmd);
}

void
CCmdLineInterface::AddOptionFloat(const char *i_cstrName,
                                  float *pfVal,
                                  char *help)
{
  CCmdLineOptionFloat cmd(i_cstrName, 1, help);

  cmd.Add(pfVal);

  // register that with the generic argument list
  m_lstFloatOption.push_back(cmd);
}

void
CCmdLineInterface::AddOptionFloat(const char *i_cstrName,
                                  float *pfVal1,
                                  float *pfVal2,
                                  char *help)
{
  CCmdLineOptionFloat cmd(i_cstrName, 2, help);

  cmd.Add(pfVal1);
  cmd.Add(pfVal2);

  // register that with the generic argument list
  m_lstFloatOption.push_back(cmd);
}

void
CCmdLineInterface::AddOptionString(const char *i_cstrName,
                                   string *pstrVal,
                                   const char *help)
{
  CCmdLineOptionString cmd(i_cstrName, 1, help);

  cmd.Add(pstrVal);

  m_lstStringOption.push_back(cmd);
}

void
CCmdLineInterface::AddOptionString(const char *i_cstrName,
                                   string *pstrVal_1,
                                   string *pstrVal_2,
                                   const char *help)
{
  CCmdLineOptionString cmd(i_cstrName, 2, help);

  cmd.Add(pstrVal_1);
  cmd.Add(pstrVal_2);

  m_lstStringOption.push_back(cmd);
}


// a bool option is basically just a flag
// no parameters at all
void
CCmdLineInterface::AddOptionBool(const char *i_cstrName,
                                 bool *pbVal,
                                 const char *help)
{
  CCmdLineOptionBool cmd(i_cstrName, help);

  cmd.Add(pbVal);

  m_lstBoolOption.push_back(cmd);
}


bool
CCmdLineInterface::Parse(int argc, char *argv[])
{

  int iTotalArgs = 1;
  char *cpOption;

  argc--; // one arguments too many

  // parsing of the options
  for ( ; argc > 0 && std::string(argv[1]).substr(0,2)=="--"; argc--, argv++ )
  {
    int nargs = 0;
    cpOption = argv[1]+2;

    list<CCmdLineOptionInt>::iterator cit_int = m_lstIntOption.begin();
    list<CCmdLineOptionFloat>::iterator cit_float = m_lstFloatOption.begin();
    list<CCmdLineOptionBool>::iterator cit_bool = m_lstBoolOption.begin();
    list<CCmdLineOptionString>::iterator cit_string =
      m_lstStringOption.begin();
    bool bFound = false;

    while ( cit_int != m_lstIntOption.end() && !bFound )
    {
      if ( cit_int->Compare(cpOption) )
      {
        for (int i=0; i < cit_int->GetArgs(); i++)
        {
          int iBuf;
          iBuf = atoi(argv[2+i]);
          cit_int->Set(iBuf);
        }
        nargs += cit_int->GetArgs();
        bFound = true;
      }
      cit_int++;
    }
    while ( cit_float != m_lstFloatOption.end() && !bFound )
    {
      if ( cit_float->Compare(cpOption) )
      {
        for (int i=0; i < cit_float->GetArgs(); i++)
        {
          float fBuf;
          fBuf = atof(argv[2+i]);
          cit_float->Set(fBuf);
        }
        nargs += cit_float->GetArgs();
        bFound = true;
      }
      cit_float++;
    }
    // a bool parameter does not have options - acts like a flag

    while ( cit_bool != m_lstBoolOption.end() && !bFound )
    {
      if ( cit_bool->Compare(cpOption) )
      {
        cit_bool->Set();
        bFound = true;
      }
      cit_bool++;
    }

    while ( cit_string != m_lstStringOption.end() && !bFound )
    {
      if ( cit_string->Compare(cpOption) )
      {
        for (int i=0; i < cit_string->GetArgs(); i++ )
        {
          string strBuf = argv[2+i];
          cit_string->Set(strBuf);
        }
        nargs += cit_string->GetArgs();
        bFound = true;
      }
      cit_string++;
    }

    if ( strcmp(cpOption, "help") == 0||strcmp(cpOption, "usage") == 0 )
    {
      PrintHelp();
    }
    else if ( !bFound )
    {
      // signal the presence of a misinterpreted option
      // will almost surely generate an error in the sequel
      // most likely to occur due to a typing error for instance
      cout << " !!! unknown option " << cpOption
           << " \n\n currently available command-line arguments are \n";
      PrintHelp();
      cout << " \n ---- \n this application will now exit\n";
      exit(1);
    }
    // do some incrementing depending on the number of parameters that were parsed
    argv += nargs;
    argc -= nargs;
  }

  // parse the IO elements now
  if ( m_lstIo.empty() )
  {
    return true;
  }

  bool bOk = true;
  list<CCmdLineIo>::iterator it = m_lstIo.begin();
  while ( argc > 0 && it != m_lstIo.end() )
  {
    cpOption = argv[1];
    it->Set(cpOption);
    iTotalArgs++;
    it++;
    argv++;
    argc--;
  }
  for ( it = m_lstIo.begin(); it != m_lstIo.end(); it++ )
  {
    bOk = bOk && it->Ok();
  }

  return bOk;
}

#include "mris_volmask.help.xml.h"
void
CCmdLineInterface::PrintHelp() const
{
  outputHelpXml(mris_volmask_help_xml,
                mris_volmask_help_xml_len);
}


void CCmdLineInterface::print() const
{
  cout << "\n --- \n Parsing Status for " << m_strProgName << "\n";
  for ( std::list<CCmdLineOptionFloat>::const_iterator cit =
          m_lstFloatOption.begin(); cit != m_lstFloatOption.end(); cit++ )
  {
    cit->outStatus();
  }
  for ( std::list<CCmdLineOptionInt>::const_iterator cit
        = m_lstIntOption.begin(); cit != m_lstIntOption.end(); cit++ )
  {
    cit->outStatus();
  }
  for ( std::list<CCmdLineOptionString>::const_iterator cit
        = m_lstStringOption.begin(); cit != m_lstStringOption.end(); cit++ )
  {
    cit->outStatus();
  }
  for ( std::list<CCmdLineOptionBool>::const_iterator cit
        = m_lstBoolOption.begin(); cit != m_lstBoolOption.end(); cit++ )
  {
    cit->outStatus();
  }
  for ( std::list<CCmdLineIo>::const_iterator cit = m_lstIo.begin();
        cit != m_lstIo.end(); cit++ )
  {
    cit->outStatus();
  }
  cout << " -------------\n";
}
