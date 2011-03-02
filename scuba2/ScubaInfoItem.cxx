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

