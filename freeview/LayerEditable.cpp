/**
 * @brief Base Layer class for editable volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
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
 *
 */

#include "LayerEditable.h"
#include <QDateTime>
#include <QVariant>

LayerEditable::LayerEditable( QObject* parent ) : Layer( parent ),
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
  setProperty("last_modified", QDateTime::currentMSecsSinceEpoch());
  emit Modified();
}
