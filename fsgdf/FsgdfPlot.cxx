/**
 * @file  FsgdfPlot.cxx
 * @brief C++ wrapper for fsgdfPlot Tcl/Tk functions
 *
 */
/*
 * Original Author: Nick Schmansky
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/01/18 22:02:11 $
 *    $Revision: 1.5 $
 *
 * Copyright (C) 2008-2010,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 *
 */

#include <iostream>
#include <sstream>
extern "C" {
#include <stdlib.h>
#include <errno.h>
#include "fio.h" // fio_FileExistsReadable
}
#include "FsgdfPlot.h"

using namespace std;

// Constructors/Destructors
//

FsgdfPlot::FsgdfPlot ()
{
  mGDFID = -1;
  mbGDFLoaded = false;
  mInterp = Tcl_CreateInterp();
  InitInterp();
  InitFsgdfPlot();
}

FsgdfPlot::FsgdfPlot ( Tcl_Interp *iInterp )
{
  mGDFID = -1;
  mbGDFLoaded = false;
  mInterp = iInterp;
  InitFsgdfPlot();
}

FsgdfPlot::~FsgdfPlot ()
{
  delete mInterp;
}


/**
 * Tcl/Tk/Tix/BLT initialization
 *
 * @return
 * @param
 */
int FsgdfPlot::InitInterp()
{
  if (Tcl_Init(this->mInterp) == TCL_ERROR)
  {
    cerr << "FsgdfPlot::Init:   Tcl_Init failed: " << 
      this->mInterp->result << endl;
    cerr << Tcl_GetStringResult( this->mInterp ) << endl;
  }
  if (Tk_Init(this->mInterp) == TCL_ERROR)
  {
    cerr << "FsgdfPlot::Init:  Tk_Init failed: " << 
      this->mInterp->result << endl;
    cerr << Tcl_GetStringResult( this->mInterp ) << endl;
  }

#if NEEDS_ITCL_ITK
  if (Itcl_Init(this->mInterp) == TCL_ERROR)
  {
    cerr << "FsgdfPlot::Init:  Itcl_Init failed: " << 
      this->mInterp->result << endl;
    cerr << Tcl_GetStringResult( this->mInterp ) << endl;
  }
  if (Itk_Init(this->mInterp) == TCL_ERROR)
  {
    cerr << "FsgdfPlot::Init:  Itk_Init failed: " << 
      this->mInterp->result << endl;
    cerr << Tcl_GetStringResult( this->mInterp ) << endl;
  }
#endif

  if (Tix_Init(this->mInterp) == TCL_ERROR)
  {
    cerr << "FsgdfPlot::Init:  Tix_Init failed: " << 
      this->mInterp->result << endl;
    cerr << Tcl_GetStringResult( this->mInterp ) << endl;
  }

  if (Blt_Init(this->mInterp) == TCL_ERROR)
  {
    cerr << "FsgdfPlot::Init:  Blt_Init failed: " << 
      this->mInterp->result << endl;
    cerr << Tcl_GetStringResult( this->mInterp ) << endl;
  }

  cerr << Tcl_GetStringResult( this->mInterp ) << endl;
  return( 0 );
}


/**
 * @return
 * @param
 */
int FsgdfPlot::InitFsgdfPlot()
{
  /* Initialize our Fsgdf functions. This is in fsgdf_wrap.c */
  if (Fsgdf_Init(this->mInterp) == TCL_ERROR)
  {
    cerr << "FsgdfPlot::Init:  Fsgdf_Init failed: " << 
      this->mInterp->result << endl;
    cerr << Tcl_GetStringResult( this->mInterp ) << endl;
  }

  // Look in a few places for the fsgdfPlot.tcl script.
  int code=-1;
  string fnFSGDF = "../scripts/fsgdfPlot.tcl";
  if( fio_FileExistsReadable( (char*)fnFSGDF.c_str() ) ) {
    code =  Tcl_EvalFile(this->mInterp, fnFSGDF.c_str() );
    cerr << Tcl_GetStringResult( this->mInterp ) << endl;
  }
  if( code != TCL_OK ) {
    char* pfnFreesurferDir = getenv( "FREESURFER_HOME" );
    if( NULL != pfnFreesurferDir ) {
      fnFSGDF = string(pfnFreesurferDir) + "/lib/tcl/fsgdfPlot.tcl";
      if( fio_FileExistsReadable( (char*)fnFSGDF.c_str() ) ) {
        code =  Tcl_EvalFile(this->mInterp, fnFSGDF.c_str() );
        cerr << Tcl_GetStringResult( this->mInterp ) << endl;
      }
    }
  }
  if( code != TCL_OK ) {
    char* pfnFSGDFDir = getenv( "FSGDF_DIR" );
    if( NULL != pfnFSGDFDir ) {
      fnFSGDF = string(pfnFSGDFDir) + "/fsgdfPlot.tcl";
      if( fio_FileExistsReadable( (char*)fnFSGDF.c_str() ) ) {
        code =  Tcl_EvalFile(this->mInterp, fnFSGDF.c_str() );
        cerr << Tcl_GetStringResult( this->mInterp ) << endl;
      }
    }
  }
  if( code != TCL_OK ) {
    cerr << "FsgdfPlot::InitFsgdf:  Couldn't find fsgdfPlot.tcl\n";
    return(1);
  }

  // Log what script we're using.
  string sMessage = string("Using ") + fnFSGDF;
  cout << sMessage << endl;

  // Initialize the FSGDF plotting stuff.
  code = Tcl_Eval(this->mInterp, "FsgdfPlot_Init" );
  if( code != TCL_OK ) cerr << Tcl_GetStringResult( this->mInterp ) << endl;

  // Hide the window for now
  code = Tcl_Eval(this->mInterp, "wm withdraw ." );
  if( code != TCL_OK ) cerr << Tcl_GetStringResult( this->mInterp ) << endl;

  if (code == TCL_OK) return(0); else return(1);
}


/**
 * @return
 * @param
 */
int FsgdfPlot::ReadFile( const char* ifnGDFFile )
{
  stringstream eval;
  eval << "FsgdfPlot_Read " <<  ifnGDFFile;
  int code = Tcl_Eval( this->mInterp, eval.str().c_str() );
  int id = strtol( Tcl_GetStringResult( this->mInterp ), 
                   (char**)NULL, 10 );
  if( ERANGE == errno ) {
    cerr << "FsgdfPlot::ReadFile:  FsgdfPlot_Read failure: Couldn't load GDF";
    return(1);
  }

  this->mGDFID = id;
  this->mbGDFLoaded = true;

  if (code == TCL_OK) return(0); else return(1);
}


/**
 * @return
 * @param
 */
bool FsgdfPlot::IsLoaded()
{
  return this->mbGDFLoaded;
}


/**
 * @return
 * @param
 */
int FsgdfPlot::BeginPointList()
{
  stringstream eval;
  eval << "FsgdfPlot_BeginPointList " <<  this->mGDFID;
  if (Tcl_Eval( this->mInterp, eval.str().c_str() ) == TCL_OK) return(0);
  cerr << Tcl_GetStringResult( this->mInterp ) << endl;
  return(1);
}


/**
 * @return
 * @param
 */
int FsgdfPlot::EndPointList()
{
  stringstream eval;
  eval << "FsgdfPlot_EndPointList " <<  this->mGDFID;
  if (Tcl_Eval( this->mInterp, eval.str().c_str() ) == TCL_OK) return(0);
  cerr << Tcl_GetStringResult( this->mInterp ) << endl;
  return(1);
}


/**
 * @return
 * @param
 */
int FsgdfPlot::AddPoint( int inVertex )
{
  stringstream eval;
  eval << "FsgdfPlot_AddPoint " << this->mGDFID << " " << inVertex << " 0 0";
  if (Tcl_Eval( this->mInterp, eval.str().c_str() ) == TCL_OK) return(0);
  cerr << Tcl_GetStringResult( this->mInterp ) << endl;
  return(1);
}


/**
 * @return
 * @param
 */
int FsgdfPlot::SetPoint( int inVertex )
{
  stringstream eval;
  eval << "FsgdfPlot_SetPoint " <<  this->mGDFID << " " << inVertex << " 0 0";
  int code = Tcl_Eval( this->mInterp, eval.str().c_str() );
  if (code == TCL_OK) return(0); else {
    cerr << Tcl_GetStringResult( this->mInterp ) << endl;
    return(1);
  }
}


/**
 * @return
 * @param
 */
int FsgdfPlot::SetInfo( const char* isLabel )
{
  stringstream eval;
  eval << "FsgdfPlot_SetInfo " <<  this->mGDFID << " \"" << isLabel << "\"";
  if (Tcl_Eval( this->mInterp, eval.str().c_str() ) == TCL_OK) return(0);
  cerr << Tcl_GetStringResult( this->mInterp ) << endl;
  return(1);
}


/**
 * @return
 * @param
 */
int FsgdfPlot::SaveToPostscript( const char* ifnPostscript )
{
  stringstream eval;
  eval << "FsgdfPlot_SaveToPostscript " <<  this->mGDFID << 
    " " << ifnPostscript;
  if (Tcl_Eval( this->mInterp, eval.str().c_str() ) == TCL_OK) return(0);
  cerr << Tcl_GetStringResult( this->mInterp ) << endl;
  return(1);
}


/**
 * @return
 * @param
 */
bool FsgdfPlot::IsWindowShowing()
{
  stringstream eval;
  eval << "FsgdfPlot_IsWindowShowing " <<  this->mGDFID;
  if (Tcl_Eval( this->mInterp, eval.str().c_str() ) != TCL_OK) return false;
  string showing = Tcl_GetStringResult( this->mInterp );
  if (showing == "1") return true;
  return false;
}


/**
 * @return
 * @param
 */
int FsgdfPlot::HideWindow()
{
  stringstream eval;
  eval << "FsgdfPlot_HideWindow " <<  this->mGDFID;
  if (Tcl_Eval( this->mInterp, eval.str().c_str() ) == TCL_OK) return(0);
  cerr << Tcl_GetStringResult( this->mInterp ) << endl;
  return(1);
}

