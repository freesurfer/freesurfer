/**
 * @file  CommonDataStruct.h
 * @brief Mix-in class for message receiving
 *
 * Simple mix-in class for use with the Broadcaster class so text
 * messages with a pointer data can be sent to a list of listeners.
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2008/10/08 19:14:34 $
 *    $Revision: 1.2 $
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


#ifndef CommonDataStruct_h
#define CommonDataStruct_h

#include <vector>


struct Settings2D
{
	bool	SyncZoomFactor;
};

struct SettingsScreenshot
{
	bool	HideCursor;
	bool	HideCoords;
	int		Magnification;
};

struct RotationElement
{
	int 	Plane;
	double 	Angle;
	double	Point[3];
};

#endif
