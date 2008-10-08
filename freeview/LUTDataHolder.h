/**
 * @file  LUTDataHolder.h
 * @brief LUT data object.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/10/08 19:14:35 $
 *    $Revision: 1.3 $
 *
 * Copyright (C) 2002-2009,
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
 
#ifndef LUTDataHolder_h
#define LUTDataHolder_h

#include <string>
#include <vector>

extern "C" {
#include "colortab.h"
}

class LUTDataHolder
{
	public:
		LUTDataHolder();
		virtual ~LUTDataHolder();	
		
		const char* GetName( int i );
		
		COLOR_TABLE* GetColorTable( int i );
		
		int GetIndex( COLOR_TABLE* ct );
		
		int GetCount();
		
	protected:
		struct ColorTableData
		{
			COLOR_TABLE*	table;
			std::string		name;
		};
		
		std::vector<ColorTableData>	m_tables;
};

#endif 


