/**
 * @brief Misc utilities
 *
 * A bunch of misc utilities for Qdec.
 */
/*
 * Original Author: Kevin Teich
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

#include "QdecUtilities.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stdexcept>

using namespace std;

void
QdecUtilities::AssertFileIsReadable ( string const& ifn ) {
  
  if( !QdecUtilities::IsFileReadable( ifn ) )
    throw runtime_error( string("Couldn't find file ") + ifn );

}

bool
QdecUtilities::IsFileReadable ( string const& ifn ) {

  ifstream f( ifn.c_str(), ios::in );
  if( !f || f.bad() )
    return false;

  return true;
}

const char *
QdecUtilities::FileNamePath(const char *fname, const char *pathName)
{
  char *slash ;

  strcpy((char*)pathName, fname) ;
  slash = (char*)strrchr(pathName, '/') ;
  if (slash) *slash = 0 ; /* remove file name */
  else
#ifndef Linux
    getcwd((char*)pathName, MAXPATHLEN-1) ;  /* no path at all, must be cwd */
#else
#if 0
  getcwd(pathName, MAXPATHLEN-1) ;
#else
  sprintf((char*)pathName, ".") ;
#endif
#endif

  return(pathName) ;
}


const char *
QdecUtilities::GetQdecrcResourceString(const char *key)
{
  // resource file .Qdecrc contains key = value pairs
  // the file can located in either or both of two directories:
  // $SUBJECTS_DIR/qdec or
  // ~  ($HOME)
  // the first key found is returned
  const char* value = NULL;
  if ( NULL != getenv("SUBJECTS_DIR") )
  {
    string rcFile = getenv("SUBJECTS_DIR");
    rcFile += "/qdec/.Qdecrc";
    value = QdecUtilities::GetResourceString( key, rcFile.c_str() );
  }
  if ( NULL == value ) // did not find key in $SUBJECTS_DIR/qdec/.Qdecrc ...
  {
    // ...so look for it in $HOME/.Qdecrc
    if ( NULL != getenv("HOME") )
    {
      string rcFile = getenv("HOME");
      rcFile += "/.Qdecrc";
      value = QdecUtilities::GetResourceString( key, rcFile.c_str() );
    }
  }  
  return value;
}

const char *
QdecUtilities::GetResourceString(const char *key, const char *filename)
{
  // resource file contains key = value pairs
  if ( NULL == filename ) return NULL;
  ifstream ifsFile;
  ifsFile.open(filename);
  if ( ifsFile.is_open())
  {
    size_t tmpstrMaxSize = 200000; // maximum size of one line in the file
    char *tmpstr = (char *)malloc( tmpstrMaxSize );
#define UNIX_EOL 0x0A
#define MAC_EOL  0x0D
    char eol_delim = UNIX_EOL;
    ifsFile.getline(tmpstr, tmpstrMaxSize, eol_delim);
    if (ifsFile.eof())
    {
      eol_delim =  MAC_EOL;
    }
    ifsFile.clear();
    ifsFile.seekg( 0 ); // rewind
#undef WHITESPC
#define WHITESPC " ,\"\t\n\r"
    int getLineRet = 0;
    while ( (getLineRet = 
             ifsFile.getline(tmpstr, tmpstrMaxSize, eol_delim).good()) ) 
    {
      char *token = strtok(tmpstr,WHITESPC);
      if ((token != NULL) && !strcmp(token,key))
      {
        // found our key
        token = strtok(NULL,WHITESPC); // get next token in this line
        if ((token != NULL) && !strcmp(token,"="))
        {
          // found the =, so return the value
          return strtok(NULL,WHITESPC);
        }
      }
      // else continue with next line
    }
  }
  ifsFile.close();
  return NULL;
}
