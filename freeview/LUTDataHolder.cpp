/**
 * @file  LUTDataHolder.h
 * @brief LUT data object.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/03/27 20:38:59 $
 *    $Revision: 1.2 $
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
 
#include <wx/wx.h>
#include <wx/filename.h>
#include "LUTDataHolder.h"

LUTDataHolder::LUTDataHolder()
{
	ColorTableData ctd;
	wxFileName fn = wxFileName::FileName( "$FREESURFER_HOME/FreeSurferColorLUT.txt" );
	fn.Normalize();
	ctd.table = CTABreadASCII( fn.GetFullPath().char_str() );
	ctd.name = "FreeSurferColorLUT";
	if ( ctd.table )
		m_tables.push_back( ctd );
	
	fn = wxFileName::FileName( "$FREESURFER_HOME/tkmeditParcColorsCMA" );
	fn.Normalize();
	ctd.table = CTABreadASCII( fn.GetFullPath().char_str() );
	ctd.name = "tkmeditParcColorsCMA";
	if ( ctd.table )
		m_tables.push_back( ctd );
	
	fn = wxFileName::FileName( "$FREESURFER_HOME/surface_labels.txt" );
	fn.Normalize();
	ctd.table = CTABreadASCII( fn.GetFullPath().char_str() );
	ctd.name = "surface_labels";
	if ( ctd.table )
		m_tables.push_back( ctd );
}
	
LUTDataHolder::~LUTDataHolder()
{
	for ( int i = 0; i < GetCount(); i++ )
	{
		CTABfree( &m_tables[i].table );
	}
}

const char* LUTDataHolder::GetName( int i )
{
	if ( i < GetCount() )
		return m_tables[i].name.c_str();
	else
		return NULL;
}
		
COLOR_TABLE* LUTDataHolder::GetColorTable( int i )
{
	if ( i < GetCount() )
		return m_tables[i].table;
	else
		return NULL;	
}
		
int LUTDataHolder::GetCount()
{
	return m_tables.size();
}

int LUTDataHolder::GetIndex( COLOR_TABLE* ct )
{
	for ( int i = 0; i < GetCount(); i++ )
	{
		if ( ct == m_tables[i].table )
			return i;
	}
	
	return -1;
}
