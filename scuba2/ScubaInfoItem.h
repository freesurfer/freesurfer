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
 *    $Author: kteich $
 *    $Date: 2007/04/06 22:23:04 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2002-2007,
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
