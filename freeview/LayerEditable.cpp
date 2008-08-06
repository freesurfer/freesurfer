/**
 * @file  LayerEditable.cpp
 * @brief Layer data object for MRI surface.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/08/06 21:07:44 $
 *    $Revision: 1.6 $
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

#include "LayerEditable.h"

LayerEditable::LayerEditable() : Layer(),
		m_nMaxUndoSteps( 100 ),
		m_bModified( false ),
		m_bEditable( true )
{
	m_strTypeNames.push_back( "Editable" );
}

LayerEditable::~LayerEditable()
{	
}
 		
void LayerEditable::SetModified()
{
	m_bModified = true;
	this->SendBroadcast( "LayerModified", this );
}
