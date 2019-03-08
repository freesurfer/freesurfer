/**
 * @file  connectgraph.c
 * @brief Export connectivity data to graphviz data format
 *
 */
/*
 * Original Author: Ruopeng Wang 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:01 $
 *    $Revision: 1.3 $
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

/*
  BEGINHELP

  ENDHELP
*/

/*
  BEGINUSAGE

  ENDUSAGE
*/


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/utsname.h>
#include "MyCmdLineParser.h"
#include <string>
#include <iostream>



#include "mri.h"
#include "mrisurf.h"
#include "colortab.h"


using namespace std;

const char *Progname = "connectgraph";

CmdLineEntry cmdLineDesc[] =
{
  CmdLineEntry( CMD_LINE_OPTION, 
                "th", 
                "threshold", 
                "<THRESHOLD_VALUE>", 
                "Set threshold of the output connectivity data.", 
                1, 
                1 ),
  CmdLineEntry( CMD_LINE_NONE )
};

char progDesc[] = "Convert connectivity data to graphviz data file.";

/*---------------------------------------------------------------*/
int main( int argc, char *argv[] ) 
{
  MyCmdLineParser cmd( 
    (const char*)"connectgraph DATA_FILE ANNOTATION_FILE OUTPUT_FILE ", 
    (CmdLineEntry*)cmdLineDesc );
  cmd.SetProgramDescription( progDesc );
  if ( !cmd.Parse( argc, (char**)argv ) )
  {
    return 1;
  }

  string_array floatingArgs = cmd.GetFloatingArguments(); 
  if ( floatingArgs.size() < 3 )
  {
    cerr << "Missing argument(s)" << endl;
    cout << "Usage: connectgraph DATA_FILE ANNOTATION_FILE OUTPUT_FILE -th <THRESHOLD>" << endl << endl;
    return 1;
  }
  
  // load connectivity data
  MRI* mri = ::MRIread( floatingArgs[0].c_str() );
  if ( !mri )
  {
    cerr << "Can not load file " << floatingArgs[0].c_str() << endl;
    return 1;
  }  
  if ( mri->width != mri->height )
    cout << "Warning: input matrix is not square" << endl;
  
  // load color table
  COLOR_TABLE* ctab = NULL;
  if ( MRISreadCTABFromAnnotationIfPresent( 
         floatingArgs[1].c_str(), &ctab ) != 0 )
  {
    cerr << "Can not load color table from  " 
         << floatingArgs[1].c_str() << endl;
    return 1;
  }
    
  // find data range and also determine threshold
  double dMin = MRIgetVoxVal( mri, 0, 0, 0, 0 );
  double dMax = dMin;
  for ( int i = 0; i < mri->width; i++ )
  {
    for ( int j = i; j < mri->height; j++ )
    {
      if ( dMin > MRIgetVoxVal( mri, i, j, 0, 0 ) )
        dMin = MRIgetVoxVal( mri, i, j, 0, 0 );
      else if ( dMax < MRIgetVoxVal( mri, i, j, 0, 0 ) )
        dMax = MRIgetVoxVal( mri, i, j, 0, 0 );
    }
  }
  
  cout << "Input data value range: " << dMin << ", " << dMax << endl;
  
  double dThreshold = dMin;     
  string_array sa;
  if ( cmd.Found( "threshold", &sa ) )
  {
    dThreshold = atof( sa[0].c_str() );
  }
  if ( dThreshold > dMax )
    cout << "Warning: Threshold is out of input data range" << endl;
  
  // write data to graphviz file
  FILE* fp = fopen( floatingArgs[2].c_str(), "w" );
  if ( !fp )
  {
    cerr << "Can not write to file " << floatingArgs[2].c_str() << endl;
    return 1;
  }
     
  fprintf( fp, "graph connectivity {\n" );
  string strg_buffer = "";
  char ch[1000];
  short* nConns = new short[ mri->width ];
  memset( nConns, 0, sizeof(short)*mri->width );
  for ( int i = 0; i < mri->width; i++ )
  {
    for ( int j = i+1; j < mri->height; j++ )
    {
      if ( MRIgetVoxVal( mri, i, j, 0, 0 ) >= dThreshold )
      {
        sprintf( ch, "%d -- %d [len=1.0];\n", i, j );
        strg_buffer += ch;
        nConns[i] = nConns[j] = 1;
      }
    }
  }
  
  for ( int i = 0; i < mri->width; i++ )
  {
    if ( nConns[i] )
    {
      int r = 255, g = 255, b = 255;
      CTABrgbAtIndexi( ctab, i, &r, &g, &b );
      fprintf( fp, 
               "node [style=\"filled\", color=\"#%02x%02x%02x\"]; %d;\n", 
               r, g, b, i );
    }
  }
  fprintf( fp, "%s", strg_buffer.c_str() );
  string name = 
    floatingArgs[0].substr( floatingArgs[0].find_last_of( '/' ) + 1 );
  fprintf( fp, "label=\"%s (%f, %f)\";\n}\n", name.c_str(), dThreshold, dMax );
  
  if ( ctab )
    CTABfree( &ctab );
  
  if ( mri )
    MRIfree( &mri );
  
  delete[] nConns;
  
  return 0;
}
