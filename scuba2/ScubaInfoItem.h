/**
 * @file  ScubaInfoItem.h
 * @brief Interface function for a displayable info item
 *
 * Provides interface for an item that can be displayed in the main
 * window's 'info item area'.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
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


#ifndef ScubaInfoItem_h
#define ScubaInfoItem_h

#include <string>

class ScubaInfoItem {

public:
  ScubaInfoItem ();

  // The label portion of the display. 
  void SetLabel ( char const* is );
  char const* GetLabel () const;

  // The value portion of the display.
  void SetValue ( char const* is );
  char const* GetValue () const;

  // Tcl function to be called when the value is changed (Or null if
  // it can't be changed).
  void SetCallback ( char const* is );
  char const* GetCallback () const;

  // A sprintf-style filter for validating the input.
  void SetInputFilter ( char const* is );
  char const* GetInputFilter () const;

  // Whether or not the value should be shortened.
  void SetShortenHint ( bool ib );
  bool GetShortenHint () const;

  void Clear();

protected:
  std::string msLabel, msValue;
  std::string msCallback;
  std::string msInputFilter;
  bool mbShortenHint;
};

#endif
