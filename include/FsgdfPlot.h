/**
 * @brief C++ wrapper for fsgdfPlot Tcl/Tk functions
 *
 */
/*
 * Original Author: Nick Schmansky
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

#ifndef FSGDFPLOT_H
#define FSGDFPLOT_H

 

#include "tix.h"
#include "fsgdf_wrap.h"

#if NEEDS_ITCL_ITK
#ifndef Itcl_Init
  int Itcl_Init(Tcl_Interp* interp);
#endif
#ifndef Itk_Init
  int Itk_Init(Tcl_Interp* interp);
#endif
#endif

#ifndef Blt_Init
  int Blt_Init ( Tcl_Interp* interp );
#endif
#ifndef Blt_SafeInit
  int Blt_SafeInit ( Tcl_Interp* interp );
#endif

#ifndef Tix_SafeInit
  int Tix_SafeInit ( Tcl_Interp* interp );
#endif


class FsgdfPlot
{
public:

  // Constructors/Destructors
  //

  FsgdfPlot ( );
  FsgdfPlot ( Tcl_Interp *iInterp );

  virtual ~FsgdfPlot ( );

  // public attribute accessor methods
  //

  /**
   * @return
   * @param
   */
  int ReadFile( const char* ifnGDFFile );

  /**
   * @return
   * @param
   */
  bool IsLoaded();

  /**
   * @return
   * @param
   */
  int BeginPointList();

  /**
   * @return
   * @param
   */
  int EndPointList();

  /**
   * @return
   * @param
   */
  int AddPoint( int inVertex );

  /**
   * @return
   * @param
   */
  int SetPoint( int inVertex );

  /**
   * @return
   * @param
   */
  int SetInfo( const char* isLabel );

  /**
   * @return
   * @param
   */
  int SaveToPostscript( const char* ifnPostscript );

  /**
   * @return
   * @param
   */
  bool IsWindowShowing();

  /**
   * @return
   * @param
   */
  int HideWindow();

private:

  // private methods
  //
  int InitInterp();
  int InitFsgdfPlot();

  // private attributes
  //
  int mGDFID;
  int mbGDFLoaded;
  Tcl_Interp *mInterp;

};

#endif // FSGDFPLOT_H
