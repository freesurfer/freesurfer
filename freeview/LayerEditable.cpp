/**
 * @file  LayerEditable.cpp
 * @brief Base Layer class for editable volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:02 $
 *    $Revision: 1.10 $
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

#include "LayerEditable.h"

LayerEditable::LayerEditable() : Layer(),
    m_nMaxUndoSteps( 100 ),
    m_bModified( false ),
    m_bEditable( true )
{
  m_strTypeNames.push_back( "Editable" );
}

LayerEditable::~LayerEditable()
{}

void LayerEditable::SetModified()
{
  m_bModified = true;
  this->SendBroadcast( "LayerModified", this );
}
