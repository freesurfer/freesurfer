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
 *    $Date: 2008/10/08 19:14:34 $
 *    $Revision: 1.4 $
 *
 * Copyright (C) 2002-2008,
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

class LayerVolumeBase;

class BrushProperty 
{
public:

  	BrushProperty ();
	virtual ~BrushProperty (); 

	int 	GetBrushSize();
	void	SetBrushSize( int nSize );
	
	int 	GetBrushTolerance();
	void 	SetBrushTolerance( int nTolerance );
	
	LayerVolumeBase* 	GetReferenceLayer();
	void 				SetReferenceLayer( LayerVolumeBase* layer );
	
	double* GetDrawRange();
	void 	SetDrawRange( double* range );
	void 	SetDrawRange( double low, double high );
	
	bool	GetDrawRangeEnabled();
	void	SetDrawRangeEnabled( bool bEnable ); 
	 
	double* GetExcludeRange();
	void	SetExcludeRange( double* range );
	void 	SetExcludeRange( double low, double high );
	
	bool	GetExcludeRangeEnabled();
	void	SetExcludeRangeEnabled( bool bEnable );
	 
	bool	GetDrawConnectedOnly();
	void	SetDrawConnectedOnly( bool bEnable );
	
protected:
  	int		m_nBrushSize;
  	int		m_nBrushTolerance;
	double	m_dDrawRange[2];
	bool	m_bEnableDrawRange;
	double	m_dExcludeRange[2];
	bool	m_bEnableExcludeRange;
	bool	m_bDrawConnectedOnly;
	
	LayerVolumeBase*	m_layerRef;
};


#endif
