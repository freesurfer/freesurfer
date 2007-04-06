/**
 * @file  ScubaInfoItem.cxx
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


#include "ScubaInfoItem.h"

ScubaInfoItem::ScubaInfoItem () :
    msLabel(""),
    msValue(""),
    msCallback(""),
    msInputFilter(""),
    mbShortenHint(false) {}

void
ScubaInfoItem::SetLabel ( char const* is ) {
  msLabel = is;
}

char const*
ScubaInfoItem::GetLabel () const {
  return msLabel.c_str();
}

void
ScubaInfoItem::SetValue ( char const* is ) {
  msValue = is;
}

char const*
ScubaInfoItem::GetValue () const {
  return msValue.c_str();
}


void
ScubaInfoItem::SetCallback ( char const* is ) {
  msCallback = is;
}

char const*
ScubaInfoItem::GetCallback () const {
  return msCallback.c_str();
}


void
ScubaInfoItem::SetInputFilter ( char const* is ) {
  msInputFilter = is;
}

char const*
ScubaInfoItem::GetInputFilter () const {
  return msInputFilter.c_str();
}


void
ScubaInfoItem::SetShortenHint ( bool ib ) {
  mbShortenHint = ib;
}

bool
ScubaInfoItem::GetShortenHint () const {
  return mbShortenHint;
}

void
ScubaInfoItem::Clear() {
  msLabel = "";
  msValue = "";
  msCallback = "";
  msInputFilter = "";
  mbShortenHint = false;
}

