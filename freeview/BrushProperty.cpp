/**
 * @file  BrushProperty.cpp
 * @brief class to hold brush properties
 *
 * Simple mix-in class for use with the Listener class so text
 * messages with a pointer data can be sent to a list of listeners.
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/06/04 20:43:23 $
 *    $Revision: 1.1.2.1 $
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


#include "BrushProperty.h"
#include "LayerEditable.h"

using namespace std;

BrushProperty::BrushProperty () :
		m_nBrushSize( 1 ), m_nBrushTolerance( 0 ), m_layerRef( NULL )
{
}
	
int BrushProperty::GetBrushSize()
{
	return m_nBrushSize;
}

void BrushProperty::SetBrushSize( int nSize )
{
	m_nBrushSize = nSize;
}
	
int BrushProperty::GetBrushTolerance()
{
	return m_nBrushTolerance;
}

void BrushProperty::SetBrushTolerance( int nTolerance )
{
	m_nBrushTolerance = nTolerance;
}

LayerEditable* BrushProperty::GetReferenceLayer()
{
	return m_layerRef;
}
	
void BrushProperty::SetReferenceLayer( LayerEditable* layer )
{
	m_layerRef = layer;
}

