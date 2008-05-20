/**
 * @file  BrushProperty.h
 * @brief class to hold brush properties
 *
 * Simpleclass for use with the Listener class so text
 * messages with a pointer data can be sent to a list of listeners.
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/05/20 16:50:17 $
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


#ifndef BrushProperty_h
#define BrushProperty_h

class LayerEditable;

class BrushProperty 
{
public:

  	BrushProperty ();
  	virtual ~BrushProperty () {}

	int GetBrushSize();
	void SetBrushSize( int nSize );
	
	int GetBrushTolerance();
	void SetBrushTolerance( int nTolerance );
	
	LayerEditable* GetReferenceLayer();
	void SetReferenceLayer( LayerEditable* layer );
	
protected:
  	int	m_nBrushSize;
  	int	m_nBrushTolerance;
	LayerEditable*	m_layerRef;
};


#endif
